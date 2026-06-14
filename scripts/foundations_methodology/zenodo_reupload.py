#!/usr/bin/env python3
"""zenodo_reupload.py — recompile + re-upload the STALE Zenodo PDFs flagged by zenodo_parity_audit.py.

The 2026-05-31 pass shipped PDFs that predated the audit note (it uploaded stale build artifacts
without recompiling). This fixes that, per record, in the only order that cannot re-ship a stale PDF:

    1. recompile  .tex -> fresh PDF (latexmk)
    2. VERIFY the fresh PDF actually contains the audit note (pdftotext) -- ABORT this record if not
    3. new-version on the existing Zenodo record (same concept-DOI; new version-DOI)
    4. delete the old file, upload the fresh PDF under the SAME filename (stable file URL)
    5. update metadata.description to carry the audit note (fixes the abstract gap)
    6. publish

Worklist = the STALE-PDF records in zenodo_parity_report.json. Re-run zenodo_parity_audit.py after.

SAFETY
------
* DRY-RUN BY DEFAULT. Steps 1-2 run (compile + verify); steps 3-6 are only PRINTED. Nothing is
  written to Zenodo without --publish.
* Published Zenodo versions CANNOT be deleted, only superseded. So: dry-run first, then do ONE real
  record (--only <stem> --publish), confirm it on the Zenodo web UI, then batch the rest.
* --sandbox targets sandbox.zenodo.org (use a sandbox token) to rehearse the whole flow harmlessly.
* The Zenodo deposit-API flow below is the standard one; it has NOT been run end-to-end from here
  (no token), so the single-record test is not optional.

Requires: latexmk + pdftotext on PATH; network; $ZENODO_TOKEN (or $ZENODO_SANDBOX_TOKEN for --sandbox).
Token scopes: deposit:write, deposit:actions. Pure stdlib otherwise.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path
from urllib.parse import quote

NOTE_SIGS = re.compile(r"audit note|predates the framework|methodology audit", re.I)
PROD = "https://zenodo.org"
SANDBOX = "https://sandbox.zenodo.org"
MIN_READABLE = 200


# ---- local helpers (kept self-contained so this runs standalone on the box) -------------------

def has_note(text: str) -> bool:
    return bool(NOTE_SIGS.search(text or ""))


def local_tex(stem: str, root: Path) -> Path | None:
    direct = root / f"{stem}.tex"
    if direct.is_file():
        return direct
    d = root / Path(stem).parent
    if d.is_dir():
        texs = sorted(d.glob("*.tex"))
        if texs:
            return texs[0]
    return None


def pdftotext(pdf: Path) -> tuple[str, bool]:
    out = subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True, timeout=120)
    t = out.stdout.decode("utf-8", "replace")
    return t, len(t.strip()) >= MIN_READABLE


def recompile(tex: Path) -> tuple[Path | None, str]:
    """latexmk -pdf in the .tex's directory. Returns (pdf_path|None, message)."""
    pdf = tex.with_suffix(".pdf")
    r = subprocess.run(
        ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", tex.name],
        cwd=tex.parent, capture_output=True, timeout=600,
    )
    if not pdf.is_file():
        tail = r.stdout.decode("utf-8", "replace")[-400:]
        return None, f"latexmk produced no PDF (rc={r.returncode}). tail: {tail}"
    return pdf, f"compiled (latexmk rc={r.returncode})"


def extract_note_html(tex_text: str) -> str | None:
    """Best-effort plain-text of the \\paragraph*{Audit note ...} block, for the Zenodo abstract."""
    m = re.search(r"\\paragraph\*?\{[^}]*[Aa]udit note[^}]*\}(.+?)(?:\\section|\\paragraph|\n\s*\n)",
                  tex_text, re.S)
    if not m:
        return None
    s = m.group(1)
    s = re.sub(r"\\(textbf|emph|textit|texttt)\{([^}]*)\}", r"\2", s)
    s = re.sub(r"\\cite\{[^}]*\}", "", s)
    s = s.replace(r"\S", "§").replace(r"\,", " ").replace(r"\\", " ").replace("~", " ")
    s = re.sub(r"\$([^$]*)\$", r"\1", s)
    s = re.sub(r"\s+", " ", s).strip()
    return ("<p><strong>Audit note.</strong> " + s[:900] + "</p>") if s else None


# ---- Zenodo deposit API ------------------------------------------------------------------------

class Zenodo:
    def __init__(self, base: str, token: str):
        self.base, self.token = base, token

    def _req(self, method: str, url: str, json_body=None, raw: bytes | None = None):
        headers = {"Authorization": f"Bearer {self.token}"}
        data = None
        if json_body is not None:
            data, headers["Content-Type"] = json.dumps(json_body).encode(), "application/json"
        elif raw is not None:
            data, headers["Content-Type"] = raw, "application/octet-stream"
        req = urllib.request.Request(url, data=data, method=method, headers=headers)
        with urllib.request.urlopen(req, timeout=300) as r:
            body = r.read()
            return json.loads(body) if body and r.headers.get("Content-Type", "").startswith("application/json") else None

    def new_version_draft(self, rec_id: str) -> dict:
        r = self._req("POST", f"{self.base}/api/deposit/depositions/{rec_id}/actions/newversion")
        draft_url = (r.get("links") or {}).get("latest_draft")
        if not draft_url:
            raise RuntimeError("newversion returned no latest_draft link")
        return self._req("GET", draft_url)

    def replace_pdf(self, draft: dict, fresh_pdf: Path):
        newid = draft["id"]
        for f in draft.get("files", []):                     # drop inherited (stale) files
            fid = f.get("id")
            if fid:
                self._req("DELETE", f"{self.base}/api/deposit/depositions/{newid}/files/{fid}")
        bucket = (draft.get("links") or {}).get("bucket")
        if not bucket:
            raise RuntimeError("draft has no bucket link for upload")
        self._req("PUT", f"{bucket}/{quote(fresh_pdf.name)}", raw=fresh_pdf.read_bytes())

    def set_description(self, draft: dict, note_html: str):
        meta = dict(draft.get("metadata") or {})
        desc = meta.get("description", "") or ""
        if not has_note(desc) and note_html:                 # don't double-prepend
            meta["description"] = note_html + "\n" + desc
            self._req("PUT", f"{self.base}/api/deposit/depositions/{draft['id']}",
                      json_body={"metadata": meta})

    def publish(self, newid) -> dict:
        return self._req("POST", f"{self.base}/api/deposit/depositions/{newid}/actions/publish")


# ---- driver ------------------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--report", default="zenodo_parity_report.json", help="output of zenodo_parity_audit.py")
    ap.add_argument("--only", default="", help="comma-separated stem substrings (TEST ONE FIRST)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--publish", action="store_true", help="actually write to Zenodo (default: dry-run)")
    ap.add_argument("--sandbox", action="store_true", help="target sandbox.zenodo.org")
    ap.add_argument("--no-fix-abstract", action="store_true", help="re-upload PDF only; leave metadata")
    ap.add_argument("--keep-going", action="store_true", help="continue after a per-record error")
    args = ap.parse_args()

    for tool in ("latexmk", "pdftotext"):
        if not shutil.which(tool):
            print(f"ERROR: {tool} not found on PATH.", file=sys.stderr); return 2
    root = Path(args.root).resolve()
    report = root / args.report
    if not report.is_file():
        print(f"ERROR: report not found: {report}. Run zenodo_parity_audit.py --json first.", file=sys.stderr)
        return 2

    base = SANDBOX if args.sandbox else PROD
    token = os.environ.get("ZENODO_SANDBOX_TOKEN" if args.sandbox else "ZENODO_TOKEN", "")
    if args.publish and not token:
        var = "ZENODO_SANDBOX_TOKEN" if args.sandbox else "ZENODO_TOKEN"
        print(f"ERROR: --publish needs ${var} set (scopes: deposit:write, deposit:actions). "
              f"export {var}=... first.", file=sys.stderr)
        return 2

    stale = [r for r in json.loads(report.read_text()) if r.get("verdict") == "STALE-PDF"]
    if args.only:
        subs = [s.strip() for s in args.only.split(",") if s.strip()]
        stale = [r for r in stale if any(s in r["stem"] for s in subs)]
    if args.limit:
        stale = stale[: args.limit]
    if not stale:
        print("No STALE-PDF records match.", file=sys.stderr); return 2

    mode = f"PUBLISH -> {base}" if args.publish else "DRY-RUN (no writes)"
    print(f"# zenodo_reupload — {len(stale)} stale record(s) | {mode}\n")
    zen = Zenodo(base, token) if args.publish else None
    done, failed = [], []

    for r in stale:
        stem, rec_id = r["stem"], r["record_id"]
        print(f"--- {stem}  (record {rec_id})")
        try:
            tex = local_tex(stem, root)
            if tex is None:
                raise RuntimeError("no local .tex found")
            pdf, msg = recompile(tex)
            print(f"    recompile: {msg}")
            if pdf is None:
                raise RuntimeError("recompile failed")
            text, readable = pdftotext(pdf)
            if not (readable and has_note(text)):
                raise RuntimeError("fresh PDF still lacks the audit note — NOT uploading (check the .tex note renders)")
            print(f"    verify: fresh PDF contains the audit note  [{pdf.name}]")
            note_html = None if args.no_fix_abstract else extract_note_html(tex.read_text(errors="replace"))

            if not args.publish:
                print(f"    DRY-RUN: would new-version {rec_id}, upload {pdf.name}, "
                      f"{'set description' if note_html else 'leave description'}, publish.")
                done.append(stem); continue

            draft = zen.new_version_draft(rec_id)
            print(f"    new draft: {draft['id']}")
            zen.replace_pdf(draft, pdf)
            print(f"    uploaded fresh PDF")
            if note_html:
                zen.set_description(draft, note_html)
                print(f"    updated description")
            pub = zen.publish(draft["id"])
            doi = (pub.get("metadata") or {}).get("doi") or (pub.get("links") or {}).get("doi", "?")
            print(f"    PUBLISHED new version: {doi}")
            done.append(stem)
        except (urllib.error.URLError, urllib.error.HTTPError, RuntimeError,
                subprocess.SubprocessError, OSError) as e:
            print(f"    ERROR: {e}")
            failed.append((stem, str(e)))
            if not args.keep_going and args.publish:
                print("\nStopping (use --keep-going to continue past errors).")
                break

    print(f"\n# Summary: {len(done)} {'processed' if args.publish else 'ready (dry-run)'}, {len(failed)} failed")
    for stem, e in failed:
        print(f"  FAILED {stem}: {e}")
    if not args.publish:
        print("\nNext: `--only <one-stem> --publish` (with $ZENODO_TOKEN) to test ONE record, verify it "
              "on zenodo.org, then re-run without --only to batch the rest. Re-run zenodo_parity_audit.py to confirm.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
