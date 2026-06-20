#!/usr/bin/env python3
"""Re-version the records touched by the 2026-06-15 CMB/GW canon shift, and redirect
the accidental cosmology duplicate. Stdlib only; mirrors zenodo_submit.py's api().

SAFETY MODEL
  * Default = STAGE: create new-version DRAFTS (private, deletable), preserve the
    record's EXISTING metadata, append only the erratum paragraph, swap in the new PDF.
    Nothing public changes. Prints each draft's pre-reserved DOI for inspection.
  * --publish: publish the staged drafts (PERMANENT) + apply the duplicate redirect.
  * --sandbox: run a throwaway FRESH-upload smoke test on sandbox.zenodo.org (validates
    metadata schema + PDF upload + description HTML end-to-end). new-version-from-existing
    cannot be sandboxed (the production records do not exist there), so the real rehearsal
    is the production STAGE step.
  * --records a,b : restrict to these keys (default: the HIGH three).

Usage:
  export ZENODO_TOKEN=... ; export ZENODO_SANDBOX_TOKEN=...
  python3 zenodo_reversion_cmb_2026_06_15.py --sandbox --records defect_network   # smoke test
  python3 zenodo_reversion_cmb_2026_06_15.py                                       # stage HIGH drafts
  python3 zenodo_reversion_cmb_2026_06_15.py --publish                            # publish + redirect
"""
import json, os, sys, urllib.error, urllib.parse, urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
SANDBOX = "--sandbox" in sys.argv
BASE = ("https://sandbox.zenodo.org" if SANDBOX else
        os.environ.get("ZENODO_URL", "https://zenodo.org")).rstrip("/")
TOKEN = os.environ.get("ZENODO_SANDBOX_TOKEN" if SANDBOX else "ZENODO_TOKEN", "") \
        or os.environ.get("ZENODO_TOKEN", "")
STAGE_FILE = HERE / "zenodo_reversion_cmb_drafts.json"

COPY = "Copyright (C) 2026 D. G. Elliman."

# Each record: published parent id, the rebuilt PDF, and the HTML erratum to APPEND.
RECORDS = {
    "defect_network": dict(
        id=20688584, pdf=HERE / "defect_network/defect_network.pdf",
        erratum=("<p><b>Revision (2026-06-15).</b> Added an equation-of-state note: a frozen "
                 "p-dimensional defect network has w = &minus;p/3, so the dominant frustrated "
                 "walls are w = &minus;2/3 (strings &minus;1/3) &mdash; a pinned, negative-pressure "
                 "relic, NOT the w = 0 clustering cold dark matter; the cold / CMB-anchoring role "
                 "belongs to the sterile-neutrino relic. The CMB third peak is a parameter-free "
                 "open no-go (matter&ndash;radiation equality z_eq &asymp; z_rec), closeable only "
                 "by an external AeST-class field.</p>")),
    "cosmology": dict(
        id=20671860, pdf=HERE / "cosmology/cosmology.pdf",
        erratum=("<p><b>Version 2.0 (2026-06-15).</b> Added the CMB acoustic third-peak open "
                 "no-go: with only baryons + the sterile-neutrino relic clustering "
                 "(&Omega;<sub>m</sub>h&sup2; &asymp; 0.046), matter&ndash;radiation equality is "
                 "z_eq &asymp; z_rec instead of &asymp;3400; both substrate-native completions "
                 "(a conserved boundary-QEC charge; a glassy line-current freeze) fail, so the "
                 "standard CMB requires an external relativistic-MOND (AeST-class) pressureless "
                 "field or stands as a falsification. Also folds in the 2026-06-14 tensor-ledger "
                 "update.</p>")),
    "falsification_sheet": dict(
        id=20693693, pdf=HERE / "falsification_sheet/falsification_sheet.pdf",
        erratum=("<p><b>Revision (2026-06-15).</b> Added the CMB third-peak parameter-free no-go "
                 "(z_eq &asymp; z_rec) and two new gravitational-wave channels: a pulsar-timing-band "
                 "stochastic background (sub-threshold) and a ringdown-echo template "
                 "(&Delta;t &asymp; 56 ms at 30 M&#8857;, reflectivity R &asymp; 1, spacing set by "
                 "the lattice cutoff a&#8320; = 1/&Lambda;<sub>QCD</sub>).</p>")),
    "gravity_blackholes": dict(
        id=20671864, pdf=HERE / "gravity_blackholes/gravity_blackholes.pdf",
        erratum=("<p><b>Version 2.0 (2026-06-15).</b> Added a gravitational-wave channels section: "
                 "the Einstein field <i>form</i> G<sub>&mu;&nu;</sub>=8&pi;G T<sub>&mu;&nu;</sub> is "
                 "intrinsic (entanglement thermodynamics, Jacobson &delta;Q=T&delta;S; emergent "
                 "M<sub>Pl</sub>&asymp;1.7 GeV at the lattice scale); a pulsar-timing-band stochastic "
                 "background (sub-threshold, T<sub>*</sub>&sim;&Lambda;<sub>QCD</sub>); and ringdown "
                 "echoes (&Delta;t&asymp;56 ms at 30 M&#8857;, reflectivity R&asymp;1, spacing set by "
                 "a&#8320;=1/&Lambda;<sub>QCD</sub>) as a LIGO/LISA test.</p>")),
    "blackhole_deep": dict(
        id=20683867, pdf=HERE / "blackhole_deep/blackhole_deep.pdf",
        erratum=("<p><b>Version 2.0 (2026-06-15).</b> Added a ringdown-echo prediction: the "
                 "firewall-over-rigid-core horizon reflects the graviton, so the remnant rings with "
                 "post-merger echoes (general relativity has none). With the reflecting surface at the "
                 "lattice cutoff a&#8320;=1/&Lambda;<sub>QCD</sub>, the spacing is &Delta;t&asymp;56 ms "
                 "at 30 M&#8857; (about half a generic Planck-scale echo) with R&asymp;1 &mdash; a sharp "
                 "LIGO/LISA test.</p>")),
    "boundary_printing": dict(
        id=20693231, pdf=HERE / "boundary_printing/boundary_printing.pdf",
        erratum=("<p><b>Version 2.0 (2026-06-15).</b> Added a dark-matter caveat: the printer fixes a "
                 "dark-matter abundance but does not complete the CMB third peak &mdash; with only "
                 "baryons + the sterile-neutrino relic clustering (&Omega;<sub>m</sub>h&sup2;&asymp;0.046), "
                 "matter&ndash;radiation equality is z_eq&asymp;z_rec rather than &asymp;3400; the R4 "
                 "line-current (w&rarr;&minus;1) and the frozen defect network (w=&minus;p/3) are both "
                 "w&lt;0, so closing the peak needs an external AeST-class pressureless mode.</p>")),
}

DUPLICATE = dict(
    id=20693186, canonical_doi="10.5281/zenodo.20671859",
    note=("<p><b>DUPLICATE RECORD &mdash; please cite 10.5281/zenodo.20671859.</b> This is an "
          "accidental 2026-06-14 re-upload of &lsquo;Cosmology, Dark Energy, and Inflation in the "
          "Finite-QEC Substrate&rsquo;, which already exists as concept DOI 10.5281/zenodo.20671859 "
          "(created 2026-06-12). Please use that record; this one is superseded.</p>"))

VERSION = "2.0"

# --echo-v21: re-version the 3 echo-affected records to v2.1 (parents = their v2.0 ids),
# correcting the now-withdrawn R~1 ringdown-echo claim to the canonical large-echo null.
ECHO_V21 = {
    "gravity_blackholes": dict(
        id=20701782, pdf=HERE / "gravity_blackholes/gravity_blackholes.pdf",
        erratum=("<p><b>Version 2.1 (2026-06-15).</b> Echo correction: the v2.0 near-unit "
                 "ringdown-echo (R~1) reading is WITHDRAWN. The QEC horizon is a one-way "
                 "record-writing channel, not a coherent mirror, so the canonical prediction is a "
                 "LARGE-ECHO NULL; LIGO/LISA echo searches are upper-bound tests on extra horizon "
                 "structure, not a positive prediction. The PTA stochastic-background and "
                 "intrinsic-Einstein-form additions are unchanged.</p>")),
    "blackhole_deep": dict(
        id=20701787, pdf=HERE / "blackhole_deep/blackhole_deep.pdf",
        erratum=("<p><b>Version 2.1 (2026-06-15).</b> Echo correction: the v2.0 near-unit "
                 "ringdown-echo (R~1) reading is WITHDRAWN -- the QEC horizon channel is one-way "
                 "record writing, not a coherent mirror; the canonical prediction is a large-echo "
                 "null (echo searches retained as upper-bound tests on extra horizon structure). "
                 "The earlier 0.31-Eddington horizon-bandwidth reading is likewise withdrawn.</p>")),
    "falsification_sheet": dict(
        id=20701528, pdf=HERE / "falsification_sheet/falsification_sheet.pdf",
        erratum=("<p><b>Version 2.1 (2026-06-15).</b> Echo correction: the v2.0 near-unit "
                 "ringdown-echo (R~1) channel is WITHDRAWN; the canonical prediction is a large-echo "
                 "null (one-way record channel, not a coherent mirror), echo searches retained as "
                 "upper-bound tests. The CMB third-peak no-go and the PTA stochastic-background "
                 "entry are unchanged.</p>")),
}


def api(method, path, body=None, ctype="application/json", raw=False):
    url = path if path.startswith("http") else BASE + path
    data = body if raw else (json.dumps(body).encode() if body is not None else None)
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {TOKEN}")
    if body is not None:
        req.add_header("Content-Type", ctype)
    try:
        with urllib.request.urlopen(req) as r:
            txt = r.read().decode()
            return json.loads(txt) if txt else {}
    except urllib.error.HTTPError as e:
        raise SystemExit(f"FATAL {method} {url} -> HTTP {e.code}\n{e.read().decode()}") from None


def upload_pdf(dep, pdf: Path):
    bucket = dep["links"]["bucket"]
    data = pdf.read_bytes()
    info = api("PUT", f"{bucket}/{urllib.parse.quote(pdf.name)}", body=data,
               ctype="application/octet-stream", raw=True)
    assert info.get("size", 0) == len(data), f"upload size mismatch {pdf.name}"
    return len(data)


def stage(keys):
    state = {}
    for k in keys:
        r = RECORDS[k]
        if not r["pdf"].exists():
            raise SystemExit(f"FATAL missing PDF {r['pdf']}")
        dep = api("POST", f"/api/deposit/depositions/{r['id']}/actions/newversion")
        new_id = int(dep["links"]["latest_draft"].rstrip("/").split("/")[-1])
        draft = api("GET", f"/api/deposit/depositions/{new_id}")
        meta = draft["metadata"]
        meta["description"] = meta.get("description", "") + r["erratum"]
        meta["version"] = VERSION
        meta.pop("prereserve_doi", None)
        for f in api("GET", f"/api/deposit/depositions/{new_id}/files"):
            api("DELETE", f"/api/deposit/depositions/{new_id}/files/{f['id']}")
        dep = api("PUT", f"/api/deposit/depositions/{new_id}", body={"metadata": meta})
        nbytes = upload_pdf(dep, r["pdf"])
        doi = (dep["metadata"].get("prereserve_doi") or {}).get("doi", "")
        state[k] = dict(draft_id=new_id, parent=r["id"], doi=doi, bytes=nbytes,
                        pdf=r["pdf"].name)
        print(f"  [{k}] staged draft {new_id} <- parent {r['id']}; PDF {r['pdf'].name} "
              f"({nbytes} B); pre-reserved DOI {doi}")
        # verify the erratum is in the staged description
        chk = api("GET", f"/api/deposit/depositions/{new_id}")
        assert "2026-06-15" in chk["metadata"]["description"], f"{k}: erratum not in draft"
        print(f"  [{k}] verified: erratum present, file attached, version {VERSION}")
    STAGE_FILE.write_text(json.dumps(state, indent=2) + "\n")
    print(f"\nStaged {len(state)} draft(s) -> {STAGE_FILE.name}. Drafts are PRIVATE and deletable.")
    return state


def publish_staged():
    if not STAGE_FILE.exists():
        raise SystemExit("FATAL: no staged drafts; run stage first.")
    state = json.loads(STAGE_FILE.read_text())
    for k, s in state.items():
        if s.get("published"):
            print(f"  [{k}] already published {s['doi']}"); continue
        pub = api("POST", f"/api/deposit/depositions/{s['draft_id']}/actions/publish")
        s.update(published=True, doi=pub["doi"], conceptdoi=pub.get("conceptdoi", ""),
                 record=pub["links"].get("record_html", ""))
        print(f"  [{k}] PUBLISHED v{VERSION}: https://doi.org/{pub['doi']}  "
              f"(concept {s['conceptdoi']})")
        STAGE_FILE.write_text(json.dumps(state, indent=2) + "\n")
    return state


def redirect_duplicate():
    d = DUPLICATE
    api("POST", f"/api/deposit/depositions/{d['id']}/actions/edit")
    dep = api("GET", f"/api/deposit/depositions/{d['id']}")
    meta = dep["metadata"]
    if "DUPLICATE RECORD" not in meta.get("description", ""):
        meta["description"] = d["note"] + meta.get("description", "")
    rels = [r for r in meta.get("related_identifiers", [])
            if not (r.get("identifier") == d["canonical_doi"]
                    and r.get("relation") == "isIdenticalTo")]
    rels.append({"identifier": d["canonical_doi"], "relation": "isIdenticalTo"})
    meta["related_identifiers"] = rels
    meta.pop("prereserve_doi", None)
    api("PUT", f"/api/deposit/depositions/{d['id']}", body={"metadata": meta})
    api("POST", f"/api/deposit/depositions/{d['id']}/actions/publish")
    print(f"  [duplicate] redirected {d['id']} -> isIdenticalTo {d['canonical_doi']} "
          f"(metadata edit, same DOI)")


def sandbox_smoke(keys):
    """Throwaway fresh-upload of each PDF + erratum metadata to sandbox; validates mechanics."""
    for k in keys:
        r = RECORDS[k]
        dep = api("POST", "/api/deposit/depositions", body={})
        meta = {"upload_type": "publication", "publication_type": "preprint",
                "title": f"[SANDBOX TEST] {k}", "description": "<p>smoke test.</p>" + r["erratum"],
                "creators": [{"name": "Elliman, David"}], "version": VERSION,
                "publication_date": "2026-06-15", "license": "cc-by-4.0", "language": "eng"}
        dep = api("PUT", f"/api/deposit/depositions/{dep['id']}", body={"metadata": meta})
        nbytes = upload_pdf(dep, r["pdf"])
        pub = api("POST", f"/api/deposit/depositions/{dep['id']}/actions/publish")
        print(f"  [sandbox {k}] OK: published throwaway {pub['doi']} ({nbytes} B) "
              f"-- metadata+PDF+description accepted")


def main():
    global RECORDS, VERSION
    if "--echo-v21" in sys.argv:
        RECORDS = ECHO_V21
        VERSION = "2.1"
    if not TOKEN:
        raise SystemExit("FATAL: token not set (ZENODO_TOKEN / ZENODO_SANDBOX_TOKEN).")
    keys = list(RECORDS)
    if "--records" in sys.argv:
        keys = sys.argv[sys.argv.index("--records") + 1].split(",")
    print(f"== target {BASE} | records: {keys} ==")
    if SANDBOX:
        sandbox_smoke(keys); return
    if "--publish" in sys.argv:
        publish_staged()
        if "--redirect" in sys.argv:
            redirect_duplicate()
        return
    stage(keys)
    print("\nReview the pre-reserved DOIs above, then re-run with --publish to publish "
          "(PERMANENT) and apply the duplicate redirect.")


if __name__ == "__main__":
    main()
