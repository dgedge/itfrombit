#!/usr/bin/env python3
"""Re-version the two dark-sector records after the 2026-06-15 canon split.

Default mode stages new-version drafts only.  Use ``--publish`` to publish the
new versions and update the local Zenodo state files.

Records:
  * papers/dark_sector          -> consolidated current dark-sector paper
  * dark_sector_cosmology       -> older corpus paper, now historical/superseded
"""
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BASE = os.environ.get("ZENODO_URL", "https://zenodo.org").rstrip("/")
TOKEN = os.environ.get("ZENODO_TOKEN", "")
STAGE_FILE = HERE / "zenodo_reversion_dark_sector_drafts.json"
STATE_MAIN = HERE / "zenodo_state.json"
STATE_CORPUS = HERE / "zenodo_corpus_state.json"

RECORDS = {
    "dark_sector": {
        "state_file": STATE_MAIN,
        "state_key": "dark_sector",
        "id": 20671862,
        "pdf": HERE / "dark_sector/dark_sector.pdf",
        "version": "2.0",
        "erratum": (
            "<p><b>Version 2.0 (2026-06-15).</b> Updated the dark-sector "
            "taxonomy after the K04 mobility/equation-of-state and CMB "
            "third-peak audits. K04 debris is now treated as a pinned, "
            "substrate-static fossil or upper-bound component, not the mobile "
            "pressureless halo medium and not the component that completes the "
            "CMB third acoustic peak. Mobile halo and recombination-era cold "
            "clustering are charged to the R4/MOND line-current sector, the "
            "sterile-&nu;<sub>R</sub> branch, or an AeST-class pressureless "
            "completion. The island-floor surface is retained as a fossil "
            "bookkeeping diagnostic rather than a full dark-matter budget.</p>"
        ),
    },
    "dark_sector_cosmo": {
        "state_file": STATE_CORPUS,
        "state_key": "dark_sector_cosmo",
        "id": 20672999,
        "pdf": ROOT / "dark_sector_cosmology/dark_sector_cosmology.pdf",
        "version": "3.0",
        "erratum": (
            "<p><b>Version 3.0 (2026-06-15).</b> Added a canon supersession "
            "notice. This older dark-sector cosmology paper is retained as a "
            "historical audit trail: its fluid-coupled equation of state, "
            "horizon-input dark-energy magnitude, graph ratio "
            "&Omega;<sub>DE</sub>/&Omega;<sub>DM</sub>=12/5, and fixed 80/20 "
            "composition are no longer current canon. Current claims are "
            "carried by the consolidated cosmology, dark-sector, "
            "defect-network, and falsification papers. The 17.7 keV "
            "sterile-&nu;<sub>R</sub> branch remains a named testable "
            "component, but not a complete dark-sector budget by itself; the "
            "CMB third acoustic peak remains an open no-go without an "
            "AeST-class or equivalent pressureless completion.</p>"
        ),
    },
}


def api(method, path, body=None, content_type="application/json", raw=False):
    url = path if path.startswith("http") else BASE + path
    data = body if raw else (json.dumps(body).encode() if body is not None else None)
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {TOKEN}")
    if body is not None:
        req.add_header("Content-Type", content_type)
    try:
        with urllib.request.urlopen(req) as resp:
            text = resp.read().decode()
            return json.loads(text) if text else {}
    except urllib.error.HTTPError as err:
        raise SystemExit(
            f"FATAL {method} {url} -> HTTP {err.code}\n{err.read().decode()}"
        ) from None


def load_json(path):
    return json.loads(path.read_text()) if path.exists() else {}


def save_json(path, data):
    path.write_text(json.dumps(data, indent=4, sort_keys=False) + "\n")


def upload_pdf(dep, pdf):
    bucket = dep["links"]["bucket"]
    data = pdf.read_bytes()
    info = api(
        "PUT",
        f"{bucket}/{urllib.parse.quote(pdf.name)}",
        body=data,
        content_type="application/octet-stream",
        raw=True,
    )
    size = info.get("size", 0)
    if size != len(data):
        raise SystemExit(f"FATAL upload size mismatch for {pdf}: {size} != {len(data)}")
    return len(data)


def find_unsubmitted_draft_by_title(title):
    """Zenodo sometimes leaves latest_draft pointing at the parent record.
    Search the user's depositions for an unsubmitted draft with the same title.
    """
    for dep in api("GET", "/api/deposit/depositions?status=draft&size=100"):
        if dep.get("submitted") or dep.get("state") == "done":
            continue
        if dep.get("metadata", {}).get("title") == title:
            return dep
    return None


def stage_one(key, rec):
    if not rec["pdf"].exists():
        raise SystemExit(f"FATAL missing PDF: {rec['pdf']}")
    try:
        dep = api("POST", f"/api/deposit/depositions/{rec['id']}/actions/newversion")
    except SystemExit as exc:
        # Zenodo keeps an unpublished latest_draft when a previous staging run
        # was interrupted after cloning the record.
        msg = str(exc).lower()
        if (
            "already" not in msg
            and "draft" not in msg
            and "remove all files first" not in msg
        ):
            raise
        parent = api("GET", f"/api/deposit/depositions/{rec['id']}")
        latest = parent.get("links", {}).get("latest_draft")
        if latest and not latest.rstrip("/").endswith(f"/{rec['id']}"):
            dep = api("GET", latest)
        else:
            dep = find_unsubmitted_draft_by_title(parent.get("metadata", {}).get("title", ""))
            if not dep:
                raise
    if "latest_draft" in dep.get("links", {}):
        draft_id = int(dep["links"]["latest_draft"].rstrip("/").split("/")[-1])
    else:
        draft_id = int(dep["id"])
    draft = api("GET", f"/api/deposit/depositions/{draft_id}")
    metadata = draft["metadata"]
    metadata["description"] = metadata.get("description", "") + rec["erratum"]
    metadata["version"] = rec["version"]
    metadata.pop("prereserve_doi", None)

    for file_info in api("GET", f"/api/deposit/depositions/{draft_id}/files"):
        try:
            api("DELETE", f"/api/deposit/depositions/{draft_id}/files/{file_info['id']}")
        except SystemExit as exc:
            print(
                f"[{key}] warning: could not delete inherited file "
                f"{file_info.get('filename', file_info.get('id'))}; "
                "continuing with bucket overwrite"
            )

    dep = api("PUT", f"/api/deposit/depositions/{draft_id}", body={"metadata": metadata})
    nbytes = upload_pdf(dep, rec["pdf"])
    doi = (dep["metadata"].get("prereserve_doi") or {}).get("doi", "")
    print(f"[{key}] staged draft {draft_id}; uploaded {rec['pdf'].name} ({nbytes} bytes); DOI {doi}")
    return {
        "draft_id": draft_id,
        "parent": rec["id"],
        "doi": doi,
        "bytes": nbytes,
        "pdf": rec["pdf"].name,
        "version": rec["version"],
    }


def publish_one(key, rec, staged):
    draft_id = staged["draft_id"]
    pub = api("POST", f"/api/deposit/depositions/{draft_id}/actions/publish")
    entry = {
        "id": draft_id,
        "published": True,
        "file": rec["pdf"].name,
        "bytes": staged["bytes"],
        "version": rec["version"],
        "previous": load_json(rec["state_file"]).get(rec["state_key"], {"id": rec["id"]}),
        "doi": pub["doi"],
        "record": pub["links"].get("record_html", ""),
        "conceptdoi": pub.get("conceptdoi", ""),
    }
    state = load_json(rec["state_file"])
    state[rec["state_key"]] = entry
    save_json(rec["state_file"], state)
    print(f"[{key}] published {pub['doi']} ({entry['record']})")
    return entry


def main():
    if not TOKEN:
        raise SystemExit("FATAL: ZENODO_TOKEN is not set")
    keys = list(RECORDS)
    if "--records" in sys.argv:
        raw = sys.argv[sys.argv.index("--records") + 1]
        keys = [k.strip() for k in raw.split(",") if k.strip()]
    unknown = [k for k in keys if k not in RECORDS]
    if unknown:
        raise SystemExit(f"FATAL unknown records: {unknown}")

    publish = "--publish" in sys.argv
    staged = {}
    if publish and STAGE_FILE.exists():
        staged = load_json(STAGE_FILE)

    for key in keys:
        rec = RECORDS[key]
        if key not in staged:
            staged[key] = stage_one(key, rec)
            save_json(STAGE_FILE, staged)
        if publish:
            publish_one(key, rec, staged[key])
            staged.pop(key, None)
            save_json(STAGE_FILE, staged)

    if not publish:
        print(f"Drafts staged in {STAGE_FILE}; rerun with --publish after inspection.")


if __name__ == "__main__":
    main()
