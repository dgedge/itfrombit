#!/usr/bin/env python3
"""zenodo_parity_audit.py — does each paper's PUBLISHED Zenodo version carry the audit note
its local .tex already has?

Background (discovered 2026-06-07): the 2026-05-31 scripted re-upload ("audit-note pass") shipped
STALE PDFs for the entire numbered "Part NN" series — the audit/retraction notes are present in the
local .tex but ABSENT from the published Zenodo PDFs — while standalone papers shipped fresh. This
re-checks the whole corpus so a re-upload targets only what is actually broken, instead of blind.

For each paper (parsed from zenodo_reconciliation.md so the paper->DOI map stays single-sourced):
    local_note : the local .tex contains an audit-note signature
    pdf_note   : the PUBLISHED Zenodo PDF (downloaded + pdftotext'd) contains it
    abs_note   : the Zenodo record's abstract/metadata.description contains it
    verdict    : one of
        STALE-PDF      local note present, PDF readable, note MISSING from PDF   <-- re-upload target
        FRESH          local note present, note present in PDF
        NO-LOCAL-NOTE  .tex has no audit-note signature (uninformative — may simply not need one)
        NO-TEX | NO-RECORD | NO-PDF | PDF-UNREADABLE | ERROR

The abstract gap is reported separately (it was universal in the spot-check, even on fresh PDFs),
because it is orthogonal: a fresh PDF can still have a stale abstract.

Requires: `pdftotext` (poppler-utils) on PATH, and network access to zenodo.org. Pure stdlib otherwise.
Exit 1 if any STALE-PDF is found (usable as a re-upload gate); exit 0 if none; exit 2 on setup error.
No aggregate score — itemised verdicts only.

Usage:
    python3 zenodo_parity_audit.py                       # audit every paper in the reconciliation
    python3 zenodo_parity_audit.py --only part_05,narrow # spot-check a known-stale + known-fresh pair
    python3 zenodo_parity_audit.py --json report.json    # also dump full machine-readable results
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

API = "https://zenodo.org/api/records/{}"
FILE_URL = "https://zenodo.org/records/{id}/files/{key}"            # verified-working public download URL
UA = {"User-Agent": "zenodo-parity-audit/1.0 (corpus self-audit)"}

# An audit note is identified by any of these (case-insensitive). Tune if the note wording changes.
NOTE_SIGS = [r"audit note", r"predates the framework", r"methodology audit"]
_SIG_RE = re.compile("|".join(NOTE_SIGS), re.I)

_ROW_PATH = re.compile(r"`([^`]+)`")            # first backtick group on a table row = the paper path
_ZENODO_ID = re.compile(r"zenodo\.(\d+)")       # all DOIs on the row; the LAST is the current/new one
MIN_READABLE_CHARS = 200                          # below this, treat extraction as failed (guard vs false STALE)


def has_note(text: str) -> bool:
    return bool(_SIG_RE.search(text or ""))


def parse_reconciliation(path: Path) -> list[dict]:
    """Pull (paper_stem, record_id) from each table row that has a backtick path + a zenodo DOI.

    Section A rows carry prev+new DOIs (take the LAST = current); section B rows carry one. Rows with
    no DOI (section C / deferred) and the headerless Strong-CP line are skipped."""
    papers, seen = [], set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.lstrip().startswith("|"):
            continue
        mp = _ROW_PATH.search(line)
        ids = _ZENODO_ID.findall(line)
        if not mp or not ids:
            continue
        stem = mp.group(1).strip()
        rid = ids[-1]                            # current DOI for this paper
        if rid in seen:
            continue
        seen.add(rid)
        papers.append({"stem": stem, "record_id": rid})
    return papers


def local_tex(stem: str, root: Path) -> Path | None:
    """Resolve the paper's .tex. The reconciliation stem is `dir/filestem` (no extension); fall back
    to globbing the directory (handles odd entries like `.../.tex`)."""
    direct = root / f"{stem}.tex"
    if direct.is_file():
        return direct
    d = root / Path(stem).parent
    if d.is_dir():
        texs = sorted(d.glob("*.tex"))
        if texs:
            return texs[0]
    return None


def http_get(url: str, timeout: float = 60.0) -> bytes:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:       # follows redirects for GET
        return r.read()


def fetch_record(record_id: str) -> dict:
    """Return {'description': str, 'pdf_key': str|None} from the Zenodo REST record."""
    data = json.loads(http_get(API.format(record_id)).decode("utf-8", "replace"))
    desc = (data.get("metadata") or {}).get("description", "") or ""
    pdf_key = next(
        (f.get("key") for f in (data.get("files") or [])
         if str(f.get("key", "")).lower().endswith(".pdf")),
        None,
    )
    return {"description": desc, "pdf_key": pdf_key}


def pdf_text(pdf_bytes: bytes) -> tuple[str, bool]:
    """Extract text via pdftotext. Returns (text, readable). readable=False guards against treating an
    extraction failure as a missing note."""
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=True) as tf:
        tf.write(pdf_bytes)
        tf.flush()
        out = subprocess.run(["pdftotext", tf.name, "-"], capture_output=True, timeout=120)
    text = out.stdout.decode("utf-8", "replace")
    return text, len(text.strip()) >= MIN_READABLE_CHARS


def audit_one(paper: dict, root: Path) -> dict:
    stem, rid = paper["stem"], paper["record_id"]
    res = {"stem": stem, "record_id": rid, "local_note": None,
           "pdf_note": None, "abs_note": None, "verdict": "ERROR", "detail": ""}
    try:
        tex = local_tex(stem, root)
        if tex is None:
            res["verdict"] = "NO-TEX"; res["detail"] = "no .tex found locally"; return res
        res["local_note"] = has_note(tex.read_text(encoding="utf-8", errors="replace"))

        try:
            rec = fetch_record(rid)
        except (urllib.error.URLError, urllib.error.HTTPError, ValueError) as e:
            res["verdict"] = "NO-RECORD"; res["detail"] = f"record fetch: {e}"; return res
        res["abs_note"] = has_note(rec["description"])

        if not rec["pdf_key"]:
            res["verdict"] = "NO-PDF"; res["detail"] = "no PDF attached to record"; return res
        try:
            blob = http_get(FILE_URL.format(id=rid, key=rec["pdf_key"]))
            text, readable = pdf_text(blob)
        except (urllib.error.URLError, urllib.error.HTTPError, subprocess.SubprocessError) as e:
            res["verdict"] = "PDF-UNREADABLE"; res["detail"] = f"download/extract: {e}"; return res
        if not readable:
            res["verdict"] = "PDF-UNREADABLE"; res["detail"] = "pdftotext returned <200 chars"; return res
        res["pdf_note"] = has_note(text)

        if not res["local_note"]:
            res["verdict"] = "NO-LOCAL-NOTE"
        elif res["pdf_note"]:
            res["verdict"] = "FRESH"
        else:
            res["verdict"] = "STALE-PDF"
    except Exception as e:                                         # never let one paper crash the run
        res["detail"] = f"{type(e).__name__}: {e}"
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="repo root (default: cwd)")
    ap.add_argument("--reconciliation", default="zenodo_reconciliation.md")
    ap.add_argument("--only", default="", help="comma-separated substrings; audit only matching paths")
    ap.add_argument("--limit", type=int, default=0, help="audit at most N papers (0 = all)")
    ap.add_argument("--delay", type=float, default=1.0, help="seconds between records (be polite to Zenodo)")
    ap.add_argument("--json", default="", help="write full results to this JSON file")
    args = ap.parse_args()

    if not shutil.which("pdftotext"):
        print("ERROR: pdftotext not found. Install poppler-utils "
              "(macOS: brew install poppler; Debian/Ubuntu: apt install poppler-utils).", file=sys.stderr)
        return 2
    root = Path(args.root).resolve()
    recon = root / args.reconciliation
    if not recon.is_file():
        print(f"ERROR: reconciliation file not found: {recon}", file=sys.stderr)
        return 2

    papers = parse_reconciliation(recon)
    if args.only:
        subs = [s.strip() for s in args.only.split(",") if s.strip()]
        papers = [p for p in papers if any(s in p["stem"] for s in subs)]
    if args.limit:
        papers = papers[: args.limit]
    if not papers:
        print("No papers matched.", file=sys.stderr); return 2

    print(f"# Zenodo parity audit — {len(papers)} papers from {recon.name}\n")
    print(f"{'paper':40} {'rec':>9}  loc  pdf  abs  verdict")
    print("-" * 78)
    results = []
    for i, p in enumerate(papers):
        r = audit_one(p, root)
        results.append(r)
        b = lambda v: "—" if v is None else ("Y" if v else "·")
        short = (Path(p["stem"]).parent.name or p["stem"])[:40]
        print(f"{short:40} {p['record_id']:>9}  {b(r['local_note']):^3}  "
              f"{b(r['pdf_note']):^3}  {b(r['abs_note']):^3}  {r['verdict']}"
              + (f"   ({r['detail']})" if r["detail"] else ""))
        if args.delay and i < len(papers) - 1:
            time.sleep(args.delay)

    stale = [r for r in results if r["verdict"] == "STALE-PDF"]
    counts: dict[str, int] = {}
    for r in results:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    abs_missing = sum(1 for r in results if r["abs_note"] is False)

    print("\n# Summary")
    for k in sorted(counts):
        print(f"  {k:16} {counts[k]}")
    print(f"  {'abstract-note missing':16} {abs_missing}  (orthogonal: applies even to FRESH PDFs)")
    if stale:
        print("\n# STALE-PDF — re-upload these record IDs (recompile .tex -> fresh PDF):")
        for r in stale:
            print(f"  {r['record_id']}  {r['stem']}")
    else:
        print("\n# No STALE-PDF papers — every local note is present in its published PDF.")

    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2), encoding="utf-8")
        print(f"\nWrote {args.json}")

    return 1 if stale else 0


if __name__ == "__main__":
    sys.exit(main())
