#!/usr/bin/env python3
"""Batch-submit the 21 revised corpus papers to Zenodo as v2.0.

Per record: newversion -> remove inherited file -> upload revised PDF ->
metadata preserved verbatim EXCEPT version=2.0, publication_date=2026-06-12,
and one additive description line explaining the v2 -> publish.
Special case: a record that is still an unpublished/file-less draft gets its
file uploaded and is published as-is (first publication, no newversion).

Idempotent: state in zenodo_corpus_state.json; re-runs skip published entries.
Usage: ZENODO_TOKEN=... python3 zenodo_corpus_submit.py [--publish]
Without --publish it stops after preparing each draft (prints pre-reserved DOI).
"""
import json
import os
import sys
import urllib.parse
import urllib.request
from pathlib import Path

BASE = os.environ.get("ZENODO_URL", "https://zenodo.org").rstrip("/")
TOKEN = os.environ.get("ZENODO_TOKEN", "")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STATE_FILE = HERE / "zenodo_corpus_state.json"

V2_NOTE = ("<p><strong>v2.0 (2026-06-12):</strong> an in-PDF dated status/"
           "erratum note has been added reflecting the June 2026 canon audit "
           "(DRIFT/ANCHOR ledger); see the paper's status note for the "
           "specific corrections, supersessions, or upgrades.</p>")

PAPERS = [
    ("forPhysicists",       20356640, "forPhysicists/forPhysicists.pdf"),
    ("velocity_note",       20323169, "technical_notes/velocity_unification_note.pdf"),
    ("velocity_rg",         20394235, "velocity_unification_rg/velocity_unification_rg.pdf"),
    ("part_12",             20374338, "part_12_fine _structure/part_12_fine_structure.pdf"),
    ("part_20",             20580699, "part_20_planck_mass/part_20_planck_mass.pdf"),
    ("holographic_dilution",20580672, "gravity/holographic_dilution.pdf"),
    ("strong_gravity",      20478725, "strong_gravity/gravity_paper.pdf"),
    ("emergent_entropic",   20478678, "emergent_entropic_gravity/emergent_entropic_gravity.pdf"),
    ("part_16",             20580696, "part_16_dark_energy/part_16_dark_energy.pdf"),
    ("part_17",             20478717, "part_17_energy_trajectory/part_17_energy_trajectory.pdf"),
    ("cosmological_qec",    20478670, "cosmological_qec_engine/cosmological_qec_engine.pdf"),
    ("dark_sector_cosmo",   20478672, "dark_sector_cosmology/dark_sector_cosmology.pdf"),
    ("ckm_hierarchy",       20478668, "CKM-hierarchy_as_quantum_error_correction/The_CKM_HierarchyAsQuantumError-Correction.pdf"),
    ("origin_of_mass",      20478694, "origin_of_mass/origin_of_mass.pdf"),
    ("part_05",             20580678, "part_05_neutrino_masses/part_05_neutrino_masses.pdf"),
    ("glueball_capstone",   20478680, "glueball_capstone/glueball_capstone.pdf"),
    ("crystallisation",     20478671, "crystallize/crystallisation.pdf"),
    ("entropy_q3",          20478679, "entropy/entropy_q3_substrate.pdf"),
    ("lindbladian",         20478690, "lindbladian_closure/lindbladian_closure.pdf"),
    ("baryogenesis_eta",    20478664, "baryogenesis_eta/baryogenesis_eta.pdf"),
    ("summary_synthesis",   20082547, "summary_May_26/summary.pdf"),
    ("eight_easy_pieces",   19834762, 'eight_easy_pieces/An accessible account of an "It for bit" quantum theory.pdf'),
]


def api(method, path, body=None, content_type="application/json", raw=False):
    url = path if path.startswith("http") else BASE + path
    data = None
    if body is not None:
        data = body if raw else json.dumps(body).encode()
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {TOKEN}")
    if body is not None:
        req.add_header("Content-Type", content_type)
    try:
        with urllib.request.urlopen(req) as r:
            txt = r.read().decode()
            return json.loads(txt) if txt else {}
    except urllib.error.HTTPError as e:
        raise SystemExit(f"FATAL {method} {url} -> HTTP {e.code}\n{e.read().decode()[:800]}")


def load_state():
    return json.loads(STATE_FILE.read_text()) if STATE_FILE.exists() else {}


def save_state(st):
    STATE_FILE.write_text(json.dumps(st, indent=1) + "\n")


def swap_file_and_meta(dep_id, pdf, version_label):
    dep = api("GET", f"/api/deposit/depositions/{dep_id}")
    for f in api("GET", f"/api/deposit/depositions/{dep_id}/files"):
        api("DELETE", f"/api/deposit/depositions/{dep_id}/files/{f['id']}")
    meta = dep["metadata"]
    meta.pop("prereserve_doi", None)
    meta["version"] = version_label
    meta["publication_date"] = "2026-06-12"
    desc = meta.get("description", "") or ""
    if "v2.0 (2026-06-12)" not in desc:
        meta["description"] = desc + "\n" + V2_NOTE
    dep = api("PUT", f"/api/deposit/depositions/{dep_id}", body={"metadata": meta})
    bucket = dep["links"]["bucket"]
    pdf_bytes = pdf.read_bytes()
    info = api("PUT", f"{bucket}/{urllib.parse.quote(pdf.name)}",
               body=pdf_bytes, content_type="application/octet-stream", raw=True)
    assert info.get("size", 0) == len(pdf_bytes), f"upload size mismatch ({pdf.name})"
    doi = (dep["metadata"].get("prereserve_doi") or {}).get("doi", "")
    return dep, doi, len(pdf_bytes)


def process(key, rec_id, rel_pdf, st, do_publish):
    entry = st.get(key, {})
    if entry.get("published"):
        print(f"  [{key}] already published v2: {entry['doi']} — skipping")
        return entry
    pdf = ROOT / rel_pdf
    if not pdf.exists():
        raise SystemExit(f"FATAL: missing PDF {pdf}")

    if not entry.get("draft_id"):
        dep = api("GET", f"/api/deposit/depositions/{rec_id}")
        if dep.get("submitted"):
            nv = api("POST", f"/api/deposit/depositions/{rec_id}/actions/newversion")
            draft_id = int(nv["links"]["latest_draft"].rstrip("/").split("/")[-1])
            mode = "newversion"
        else:
            draft_id = rec_id          # file-less draft: first publication
            mode = "first-publication"
        entry = {"old_id": rec_id, "draft_id": draft_id, "mode": mode}
        print(f"  [{key}] {mode}: draft {draft_id}")
    else:
        print(f"  [{key}] resuming draft {entry['draft_id']}")

    vlabel = "2.0" if entry["mode"] == "newversion" else "1.0"
    dep, predoi, nbytes = swap_file_and_meta(entry["draft_id"], pdf, vlabel)
    entry.update(predoi=predoi, bytes=nbytes)
    print(f"  [{key}] uploaded {pdf.name} ({nbytes} b), pre-reserved DOI {predoi}")

    if do_publish:
        pub = api("POST", f"/api/deposit/depositions/{entry['draft_id']}/actions/publish")
        entry.update(doi=pub["doi"], published=True,
                     conceptdoi=pub.get("conceptdoi", ""),
                     record=pub["links"].get("record_html", ""))
        print(f"  [{key}] PUBLISHED: https://doi.org/{pub['doi']}")
    return entry


def main():
    if not TOKEN:
        raise SystemExit("FATAL: ZENODO_TOKEN not set")
    do_publish = "--publish" in sys.argv
    st = load_state()
    for key, rec_id, rel_pdf in PAPERS:
        st[key] = process(key, rec_id, rel_pdf, st, do_publish)
        save_state(st)
    print("\n== Final table ==")
    for key, rec_id, _ in PAPERS:
        e = st[key]
        print(f"  {key:<22s} {e.get('doi', e.get('predoi','?')):<32s} "
              f"{'published' if e.get('published') else 'draft'}")


if __name__ == "__main__":
    main()
