#!/usr/bin/env python3
"""Submit the five finite-QEC companion papers to Zenodo as preprints.

Usage:
    export ZENODO_TOKEN=...          # personal token with deposit:write + deposit:actions
    python3 zenodo_submit.py         # phase 1: create/refresh DRAFT depositions,
                                     #          upload PDFs, set metadata,
                                     #          print the pre-reserved DOIs
    python3 zenodo_submit.py --publish   # phase 2: publish ALL five (irreversible)

The overview paper is deliberately NOT in the list: it will cite these five
DOIs in its references and is submitted afterwards.

Design notes:
  * stdlib only (urllib) — no pip dependencies.
  * Two-phase: --publish refuses to run unless every paper has a complete
    draft (file uploaded + metadata accepted).  Publishing is permanent on
    Zenodo (files frozen; metadata stays editable), so all validation
    happens at draft stage.
  * State is recorded in zenodo_state.json next to this script; reruns are
    idempotent (existing depositions are updated, not duplicated).
  * The DOI printed at draft stage is pre-reserved: it is the same DOI that
    becomes active at publish, so it can be quoted in the overview's
    references immediately.
  * Set ZENODO_URL=https://sandbox.zenodo.org (with a sandbox token) to
    rehearse against the sandbox.
"""

import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

BASE = os.environ.get("ZENODO_URL", "https://zenodo.org").rstrip("/")
TOKEN = os.environ.get("ZENODO_TOKEN", "")
HERE = Path(__file__).resolve().parent
STATE_FILE = HERE / "zenodo_state.json"

CREATOR = {
    "name": "Elliman, David",
    "affiliation": "Neuro-symbolic Ltd",
    "orcid": "0009-0001-9792-4036",
}
COPYRIGHT = "Copyright (C) 2026 D. G. Elliman."
REPO = "https://github.com/dgedge/itfrombit"
SERIES_NOTE = (
    COPYRIGHT
    + " Part of the finite-QEC substrate paper series; the overview paper "
    + "cites the DOIs of the five companions. Code and reproducibility "
    + f"repository: {REPO} — project page: https://neusym.ai"
)

COMMON_META = {
    "upload_type": "publication",
    "publication_type": "preprint",
    "publication_date": "2026-06-12",
    "version": "1.0",
    "language": "eng",
    "license": "cc-by-4.0",
    "creators": [CREATOR],
    "notes": SERIES_NOTE,
    "related_identifiers": [
        {"identifier": REPO, "relation": "isSupplementedBy"},
    ],
    "prereserve_doi": True,
}

PAPERS = [
    {
        "key": "foundations",
        "pdf": HERE / "foundations/foundations.pdf",
        "title": (
            "Foundations and Methodology for a Finite-QEC Substrate: "
            "Code, Crystallisation, Ledgers, and Audit Protocol"
        ),
        "description": (
            "<p>States the foundational objects and audit methodology of the "
            "finite-QEC (quantum-error-correction) substrate program: the "
            "eight-bit local code, the bi-cubic crystallisation picture, the "
            "distinction between syndrome and strain readouts, the 28 = 2×14 "
            "service-clock construction, and the methodological rules that "
            "classify every claim as locked, computed, conditional, retired, "
            "or open. The goal is to make the companion physics papers "
            "auditable: every coefficient must name its carrier, event unit, "
            "scheduler, observable map, and reproducibility check.</p>"
            "<p><b>Version 3.0.</b> Added a radiation-sector domain-of-validity "
            "caveat (the emergent gauge sector is an effective theory below the "
            "substrate's Brillouin scale; trans-cutoff single photons are an "
            "open foundational question) and the relativity companion "
            "citation.</p>"
            f"<p>{COPYRIGHT}</p>"
        ),
        "keywords": [
            "quantum error correction", "finite substrate", "discrete spacetime",
            "it from bit", "extended Hamming code", "truncated cubic honeycomb",
            "crystallisation", "audit methodology", "reproducible research",
            "service clock", "syndrome decoding", "strain ledger",
            "emergent physics", "digital physics",
        ],
    },
    {
        "key": "matter_gauge",
        "pdf": HERE / "matter_gauge/matter_gauge.pdf",
        "title": (
            "Matter, Gauge Structure, and Spectroscopy in the Finite-QEC "
            "Substrate"
        ),
        "description": (
            "<p>States the matter and gauge canon of the finite-QEC substrate "
            "program: the recovery of Standard-Model-like charge, colour, and "
            "three-generation structure from an eight-bit local code on the "
            "truncated cubic honeycomb; the gauge-web construction; hadron "
            "spectroscopy including the nucleon mass at the 10<sup>-3</sup> "
            "level; and the named open problems, including the continuum "
            "chiral/mirror-gapping question.</p>"
            f"<p>{COPYRIGHT}</p>"
        ),
        "keywords": [
            "Standard Model", "gauge theory", "quantum error correction",
            "discrete spacetime", "lattice gauge theory", "fermion doubling",
            "chiral fermions", "hadron spectroscopy", "proton mass",
            "generation structure", "fine-structure constant",
            "emergent particle physics", "finite substrate", "it from bit",
        ],
    },
    {
        "key": "cosmology",
        "pdf": HERE / "cosmology/cosmology.pdf",
        "title": "Cosmology, Dark Energy, and Inflation in the Finite-QEC Substrate",
        "description": (
            "<p>States the cosmology canon of the finite-QEC substrate "
            "program: a Landauer/QEC accounting of the cosmological constant "
            "landing on the observed dark-energy density; the dynamic "
            "dark-energy law w(a) = -1 + a/28 with its discriminating "
            "w<sub>0</sub> = -27/28; the holographic-boundary-crystallisation "
            "route to the spectral index n<sub>s</sub> = 27/28 = 0.96429 "
            "(-0.15 sigma against Planck 2018); and named exposures, "
            "including the tensor-amplitude question for the candidate "
            "printing scale.</p>"
            f"<p>{COPYRIGHT}</p>"
        ),
        "keywords": [
            "cosmology", "dark energy", "cosmological constant",
            "equation of state", "w0-wa parametrisation", "spectral index",
            "inflation alternative", "Landauer principle",
            "quantum error correction", "Hubble constant",
            "primordial power spectrum", "discrete spacetime",
            "finite substrate", "emergent cosmology",
        ],
    },
    {
        "key": "dark_sector",
        "pdf": HERE / "dark_sector/dark_sector.pdf",
        "title": "Dark Matter, MOND, and K04 Debris in the Finite-QEC Substrate",
        "description": (
            "<p>States the dark-sector canon of the finite-QEC substrate "
            "program: a Kibble-Zurek debris mechanism in which rapid boot "
            "crystallisation traps durable domain-wall fragments whose "
            "gravitating ledger is the recorded boundary shadow; an R4/MOND "
            "line-current mechanism producing the cubic AQUAL action and the "
            "baryonic Tully-Fisher relation under a Poisson line-current "
            "theorem; and the per-branch kill conditions — the island-floor "
            "abundance surface, the chi<sub>R</sub> = 1 normalisation, and "
            "cluster-scale lensing discrimination of pinned versus "
            "collisionless wall debris.</p>"
            f"<p>{COPYRIGHT}</p>"
        ),
        "keywords": [
            "dark matter", "MOND", "modified Newtonian dynamics",
            "baryonic Tully-Fisher relation", "Kibble-Zurek mechanism",
            "domain walls", "topological defects", "AQUAL",
            "galaxy rotation curves", "Bullet Cluster",
            "pseudo-isothermal halo", "quantum error correction",
            "discrete spacetime", "dark sector",
        ],
    },
    {
        "key": "gravity_blackholes",
        "pdf": HERE / "gravity_blackholes/gravity_blackholes.pdf",
        "title": "Gravity, Horizons, and Black Holes in the Finite-QEC Substrate",
        "description": (
            "<p>States the gravity and horizon canon of the finite-QEC "
            "substrate program: the proton-primary route in which the "
            "measured proton mass fixes the chiral scale and the accounting "
            "chain predicts G, the Planck mass, and H<sub>0</sub> = 67.27 "
            "km/s/Mpc; the Bekenstein severing-channel count C = 55/8 and "
            "its alpha-free image Omega<sub>Lambda</sub> = 12&pi;/55; a "
            "finite-cell horizon isometry, a Schwarzschild shell channel, "
            "and a finite Hawking ladder on the 208-state invalid subspace; "
            "and explicit falsification surfaces — a 22-ppm CODATA test of G "
            "once the alpha-convention gate resolves, and death of the "
            "H<sub>0</sub> prediction under a SH0ES-side resolution of the "
            "Hubble tension.</p>"
            f"<p>{COPYRIGHT}</p>"
        ),
        "keywords": [
            "emergent gravity", "Newton's constant", "black hole entropy",
            "Bekenstein-Hawking entropy", "Hawking radiation",
            "black hole information", "holography", "KMS condition",
            "Hubble tension", "proton mass", "quantum error correction",
            "discrete spacetime", "horizon thermodynamics", "finite substrate",
        ],
    },
    {
        "key": "blackhole_deep",
        "pdf": HERE / "blackhole_deep/blackhole_deep.pdf",
        "publication_date": "2026-06-13",
        "title": (
            "Going Deeper into a Black Hole: The Horizon as a "
            "Quantum-Error-Correcting Record"
        ),
        "description": (
            "<p>A focused deep-dive companion to the gravity and black-hole "
            "paper (DOI 10.5281/zenodo.20671864) in the finite-QEC "
            "(quantum-error-correction) substrate series. It develops one "
            "thread: the black-hole horizon as a finite quantum-error-correcting "
            "record governed by a single object, the three-cube coboundary "
            "&delta;.</p>"
            "<p>Three results. (1) One &delta; controls the boundary strain "
            "record, the firewall isometry V<sub>cell</sub>, and the single "
            "blind degree of freedom; the Hawking degeneracies and the area "
            "coefficient additionally require a register-validity/monogamy "
            "ingredient that is provably independent of &delta;. (2) The "
            "Bekenstein area law, in substrate units (node area "
            "a<sub>0</sub><sup>2</sup>/4, proton-primary G), is exactly "
            "equivalent to a microscopic records rate: each horizon node is "
            "assigned about 10<sup>38</sup> nats while a cell holds at most "
            "8&nbsp;ln&nbsp;2 &approx; 5.5, so the area entropy cannot be "
            "stored (standing storage fails by 37&ndash;45 orders of magnitude) "
            "and must instead flow, at rate H<sub>0</sub>M<sub>P</sub><sup>2</sup>"
            "/(16&Lambda;<sub>QCD</sub><sup>3</sup>) = C&alpha;<sub>0</sub>"
            "<sup>2</sup>. (3) The coefficient is C = 55/8, with the "
            "value-level-versus-address-level direction-tag fork discharged by "
            "two independent arguments: the record channel is the syndrome "
            "itself, and AGL(3,2) covariance forbids an address-level "
            "orientation.</p>"
            "<p>Predictions: &Omega;<sub>&Lambda;</sub> = 12&pi;/55 = 0.6854 "
            "(&alpha;-free, +0.1&sigma;, Planck branch of the Hubble tension); a "
            "universal, mass-independent horizon information bandwidth of "
            "&approx;0.31 Eddington at the base service rate (&approx;43 "
            "Eddington at once-per-tick), with super-Eddington accretors as the "
            "falsification surface; and a discrete Hawking line spectrum with "
            "fixed integer-strain degeneracies and gap F<sub>min</sub> = 3. The "
            "absolute Planck-mass scale is not settled here: it rests on an "
            "un-derived &alpha;<sup>2</sup>, deferred to the gravity "
            "companion.</p>"
            f"<p>{COPYRIGHT}</p>"
        ),
        "keywords": [
            "black hole entropy", "Bekenstein-Hawking entropy",
            "Hawking radiation", "black hole information", "holography",
            "horizon thermodynamics", "quantum error correction", "dark energy",
            "Hubble tension", "Eddington limit", "KMS condition",
            "discrete spacetime", "finite substrate", "emergent gravity",
        ],
        "notes": (
            COPYRIGHT
            + " A deep-dive companion to the gravity/black-holes paper "
            + "(DOI 10.5281/zenodo.20671864) in the finite-QEC substrate "
            + f"series. Code and reproducibility repository: {REPO} — "
            + "project page: https://neusym.ai"
        ),
        "related": [
            {"identifier": "10.5281/zenodo.20671864",
             "relation": "isSupplementTo"},
        ],
    },
    {
        "key": "relativity",
        "pdf": HERE / "relativity/relativity.pdf",
        "title": (
            "Special and General Relativity from the Finite-QEC Substrate: "
            "The Propagation Clock, the Equivalence Principle, and the "
            "Horizon Ledger"
        ),
        "description": (
            "<p>States how relativity, in both forms, is recovered from the "
            "finite-QEC substrate. Special relativity is the universal "
            "reversible propagation clock: the invariant speed c is one "
            "lattice step per service tick, and observers built from the same "
            "walk/QEC clock infer Lorentz symmetry despite a preferred update "
            "order (1D walk theorem; an exactly isotropic massless Weyl cone "
            "H<sup>2</sup> = (&Sigma; sin<sup>2</sup>k<sub>d</sub>)I from an "
            "anticommuting chirality Clifford triple; no bare Dirac mass). "
            "General relativity is the coarse-grained covariance of that clock "
            "under strain: straining the cone gives the metric directly, "
            "g<sup>ij</sup> = &delta;<sup>ij</sup> + 2&epsilon;<sup>ij</sup>, "
            "so the metric perturbation IS the strain field; the matter "
            "strain-response is a conserved symmetric stress tensor; the "
            "linearised Bianchi identity and the substrate's strain-ledger "
            "conservation are one statement, making G = 8&pi;G T a consistent "
            "field equation. The equivalence principle follows once gravity "
            "couples to the full energy (kinetic plus confinement/Yukawa "
            "mass), not the bare hopping (which would be anti-mass); and the "
            "four horizon objects (area law, Hawking ladder, firewall "
            "isometry, gravitational source) are one register syndrome read "
            "four ways. Tiers are explicit: the kinematic spine is computed "
            "and machine-verified; the Sakharov and "
            "Weinberg steps are flagged imports.</p>"
            "<p>The physical photon is the Gauss-projected transverse "
            "SC-Maxwell mode on the dual simple-cubic gauge web, with an "
            "irrelevant O((a<sub>0</sub>k)<sup>2</sup>) velocity anisotropy, so "
            "Lorentz invariance holds at accessible energies; the open residual "
            "is the high-energy substrate cutoff.</p>"
            "<p><b>Version 3.0.</b> Photon section reconciled: the physical "
            "photon is the Gauss-projected SC-Maxwell mode (not the K6 "
            "T<sub>1u</sub>/E<sub>g</sub> band); this retracts the v2 "
            "Collins&ndash;Perez&ndash;Sudarsky 'open marginal anisotropy' "
            "framing, and the open residual is the high-energy substrate "
            "cutoff.</p>"
            f"<p>{COPYRIGHT}</p>"
        ),
        "keywords": [
            "special relativity", "general relativity", "Lorentz invariance",
            "equivalence principle", "emergent gravity", "Weyl fermion",
            "Einstein field equations", "stress-energy tensor",
            "induced gravity", "Nielsen-Ninomiya", "quantum error correction",
            "discrete spacetime", "graviton", "finite substrate",
        ],
    },
    {
        "key": "overview",
        "pdf": HERE / "overview/overview.pdf",
        "title": (
            "A Finite-QEC Substrate Program for Particle Physics and "
            "Cosmology: Current Canon and Audit Methodology"
        ),
        "description": (
            "<p>Entry-point overview of a six-paper series in which particle "
            "physics and cosmology are modelled as consequences of a "
            "quantum-error-correcting register geometry rather than as "
            "fields on a pre-assumed continuum. Summarises the current "
            "canon under an explicit audit methodology — every claim "
            "classified as locked, computed, conditional, retired, or open, "
            "with a standing retraction ledger. Headline themes: recovery "
            "of Standard-Model-like matter structure from an eight-bit "
            "local code; the 28 = 2&times;14 service clock; a "
            "proton-primary route predicting G, the Planck mass, and "
            "H<sub>0</sub>; the dynamic dark-energy law w(a) = -1 + a/28; "
            "a spectral-index route to n<sub>s</sub> = 27/28; and a "
            "crystallisation/debris dark sector. The paper states its "
            "standing falsification targets, and cites the six companion "
            "preprints (DOIs 10.5281/zenodo.20672042, .20671858, "
            ".20671860, .20671862, .20671864, .20677639) that carry the "
            "sector-by-sector detail.</p>"
            "<p><b>Version 5.0.</b> Relativity-section photon paragraph "
            "reconciled to the companion's v3: the photon is the "
            "Gauss-projected SC-Maxwell mode, Lorentz-invariant at accessible "
            "energies, with the high-energy substrate cutoff as the open "
            "residual.</p>"
            f"<p>{COPYRIGHT}</p>"
        ),
        "keywords": [
            "quantum error correction", "discrete spacetime", "it from bit",
            "finite substrate", "Standard Model", "emergent gravity",
            "dark energy", "cosmological constant", "dark matter", "MOND",
            "black hole entropy", "Hubble constant", "spectral index",
            "proton mass", "audit methodology", "reproducible research",
        ],
        "notes": (
            COPYRIGHT
            + " Entry-point overview of the finite-QEC substrate paper "
            + "series; the six companion preprints are DOIs "
            + "10.5281/zenodo.20672042, 10.5281/zenodo.20671858, "
            + "10.5281/zenodo.20671860, 10.5281/zenodo.20671862, "
            + "10.5281/zenodo.20671864, 10.5281/zenodo.20677639. Code and "
            + f"reproducibility repository: {REPO} — project page: https://neusym.ai"
        ),
        "related": [
            {"identifier": "10.5281/zenodo.20672042", "relation": "cites"},
            {"identifier": "10.5281/zenodo.20671858", "relation": "cites"},
            {"identifier": "10.5281/zenodo.20671860", "relation": "cites"},
            {"identifier": "10.5281/zenodo.20671862", "relation": "cites"},
            {"identifier": "10.5281/zenodo.20671864", "relation": "cites"},
            {"identifier": "10.5281/zenodo.20677639", "relation": "cites"},
        ],
    },
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
        detail = e.read().decode()
        raise SystemExit(
            f"FATAL {method} {url} -> HTTP {e.code}\n{detail}"
        ) from None


def load_state():
    if STATE_FILE.exists():
        return json.loads(STATE_FILE.read_text())
    return {}


def save_state(state):
    STATE_FILE.write_text(json.dumps(state, indent=2) + "\n")


def build_meta(paper, version=None):
    meta = dict(COMMON_META)
    meta.update(
        title=paper["title"],
        description=paper["description"],
        keywords=paper["keywords"],
    )
    if paper.get("notes"):
        meta["notes"] = paper["notes"]
    if paper.get("related"):
        meta["related_identifiers"] = (
            COMMON_META["related_identifiers"] + paper["related"]
        )
    if paper.get("publication_date"):
        meta["publication_date"] = paper["publication_date"]
    if version is not None:
        meta["version"] = version
    return meta


def ensure_draft(paper, state):
    """Create or refresh one draft deposition; return its state entry."""
    key = paper["key"]
    entry = state.get(key, {})

    if not paper["pdf"].exists():
        raise SystemExit(f"FATAL: missing PDF {paper['pdf']}")

    if entry.get("id"):
        dep = api("GET", f"/api/deposit/depositions/{entry['id']}")
        if dep.get("submitted"):
            print(f"  [{key}] already PUBLISHED as {dep['doi']} — skipping")
            entry.update(doi=dep["doi"], published=True,
                         record=dep["links"].get("record_html", ""))
            return entry
        print(f"  [{key}] reusing existing draft {entry['id']}")
    else:
        dep = api("POST", "/api/deposit/depositions", body={})
        entry["id"] = dep["id"]
        print(f"  [{key}] created draft deposition {dep['id']}")

    dep = api("PUT", f"/api/deposit/depositions/{entry['id']}",
              body={"metadata": build_meta(paper)})

    doi = (dep["metadata"].get("prereserve_doi") or {}).get("doi", "")
    entry["doi"] = doi

    # upload (or replace) the PDF via the bucket API
    bucket = dep["links"]["bucket"]
    fname = paper["pdf"].name
    pdf_bytes = paper["pdf"].read_bytes()
    info = api("PUT", f"{bucket}/{urllib.parse.quote(fname)}",
               body=pdf_bytes, content_type="application/octet-stream",
               raw=True)
    size = info.get("size", 0)
    assert size == len(pdf_bytes), f"upload size mismatch for {key}"
    print(f"  [{key}] uploaded {fname} ({size} bytes)")
    print(f"  [{key}] pre-reserved DOI: {doi}")

    entry.update(published=False, file=fname, bytes=size,
                 draft_url=dep["links"].get("html", ""))
    return entry


def publish(paper, state):
    key = paper["key"]
    entry = state[key]
    if entry.get("published"):
        print(f"  [{key}] already published: {entry['doi']}")
        return entry
    dep = api("POST",
              f"/api/deposit/depositions/{entry['id']}/actions/publish")
    entry.update(doi=dep["doi"], published=True,
                 record=dep["links"].get("record_html", ""))
    print(f"  [{key}] PUBLISHED: https://doi.org/{dep['doi']}")
    return entry


def new_version(paper, state, version_string, do_publish):
    """Mint a new version of an already-published record: clone via the
    newversion action, swap the file, bump version metadata, optionally
    publish. The concept DOI is unchanged and always resolves to latest."""
    key = paper["key"]
    entry = state.get(key, {})
    assert entry.get("published"), f"{key} has no published version to bump"
    if not paper["pdf"].exists():
        raise SystemExit(f"FATAL: missing PDF {paper['pdf']}")

    dep = api("POST",
              f"/api/deposit/depositions/{entry['id']}/actions/newversion")
    new_id = int(dep["links"]["latest_draft"].rstrip("/").split("/")[-1])
    print(f"  [{key}] new-version draft {new_id} (from {entry['id']})")

    for f in api("GET", f"/api/deposit/depositions/{new_id}/files"):
        api("DELETE", f"/api/deposit/depositions/{new_id}/files/{f['id']}")
        print(f"  [{key}] removed inherited file {f.get('filename', '?')}")

    dep = api("PUT", f"/api/deposit/depositions/{new_id}",
              body={"metadata": build_meta(paper, version=version_string)})

    bucket = dep["links"]["bucket"]
    fname = paper["pdf"].name
    pdf_bytes = paper["pdf"].read_bytes()
    info = api("PUT", f"{bucket}/{urllib.parse.quote(fname)}",
               body=pdf_bytes, content_type="application/octet-stream",
               raw=True)
    assert info.get("size", 0) == len(pdf_bytes), "upload size mismatch"
    print(f"  [{key}] uploaded {fname} ({len(pdf_bytes)} bytes)")

    new_entry = {
        "id": new_id, "published": False, "file": fname,
        "bytes": len(pdf_bytes), "version": version_string,
        "previous": entry,
        "doi": (dep["metadata"].get("prereserve_doi") or {}).get("doi", ""),
    }
    print(f"  [{key}] pre-reserved DOI: {new_entry['doi']}")

    if do_publish:
        pub = api("POST",
                  f"/api/deposit/depositions/{new_id}/actions/publish")
        new_entry.update(doi=pub["doi"], published=True,
                         record=pub["links"].get("record_html", ""),
                         conceptdoi=pub.get("conceptdoi", ""))
        print(f"  [{key}] PUBLISHED v{version_string}: "
              f"https://doi.org/{pub['doi']}")
        if new_entry["conceptdoi"]:
            print(f"  [{key}] concept DOI (always latest): "
                  f"{new_entry['conceptdoi']}")
    return new_entry


def link_cited_by(state):
    """Add an isCitedBy -> overview-concept-DOI relation to each companion
    record. Pure metadata edit (edit -> PUT -> publish): no new version is
    minted and the files stay frozen."""
    ov = state.get("overview", {})
    concept = ov.get("conceptdoi", "")
    if not concept:
        dep = api("GET", f"/api/deposit/depositions/{ov['id']}")
        concept = dep.get("conceptdoi", "")
    assert concept, "overview concept DOI unavailable"
    print(f"  overview concept DOI: {concept}")
    for key in ["foundations", "matter_gauge", "cosmology", "dark_sector",
                "gravity_blackholes", "relativity"]:
        entry = state[key]
        dep = api("GET", f"/api/deposit/depositions/{entry['id']}")
        meta = dep["metadata"]
        rels = meta.get("related_identifiers", [])
        if any(r.get("identifier") == concept
               and r.get("relation") == "isCitedBy" for r in rels):
            print(f"  [{key}] isCitedBy already present — skipping")
            continue
        api("POST", f"/api/deposit/depositions/{entry['id']}/actions/edit")
        meta.pop("prereserve_doi", None)
        meta["related_identifiers"] = rels + [
            {"identifier": concept, "relation": "isCitedBy"}]
        api("PUT", f"/api/deposit/depositions/{entry['id']}",
            body={"metadata": meta})
        api("POST", f"/api/deposit/depositions/{entry['id']}/actions/publish")
        print(f"  [{key}] isCitedBy {concept} added (metadata edit, "
              f"same version)")


def retarget_repo(state, old, new):
    """Swap a repository-URL substring in the metadata of every PUBLISHED record
    (edit -> PUT -> publish; files stay frozen, no new version, DOI unchanged).
    Used to point the isSupplementedBy link and notes at the renamed code repo."""
    for key, entry in state.items():
        if not entry.get("published") or not entry.get("id"):
            continue
        try:
            dep = api("GET", f"/api/deposit/depositions/{entry['id']}")
            meta = dep["metadata"]
            swapped = json.loads(json.dumps(meta).replace(old, new))
            if swapped == meta:
                print(f"  [{key}] no '{old}' in metadata -- skipping")
                continue
            swapped.pop("prereserve_doi", None)
            api("POST", f"/api/deposit/depositions/{entry['id']}/actions/edit")
            api("PUT", f"/api/deposit/depositions/{entry['id']}",
                body={"metadata": swapped})
            api("POST", f"/api/deposit/depositions/{entry['id']}/actions/publish")
            print(f"  [{key}] retargeted -> {new} (metadata edit, DOI {entry.get('doi')})")
        except SystemExit as e:
            print(f"  [{key}] FAILED: {e}")


def main():
    if not TOKEN:
        raise SystemExit("FATAL: ZENODO_TOKEN is not set in the environment.")
    do_publish = "--publish" in sys.argv

    if "--retarget-repo" in sys.argv:
        i = sys.argv.index("--retarget-repo")
        retarget_repo(load_state(), sys.argv[i + 1], sys.argv[i + 2])
        return

    if "--link-cited-by" in sys.argv:
        link_cited_by(load_state())
        return

    if "--new-version" in sys.argv:
        key = sys.argv[sys.argv.index("--new-version") + 1]
        vstr = "2.0"
        if "--version-string" in sys.argv:
            vstr = sys.argv[sys.argv.index("--version-string") + 1]
        paper = next((p for p in PAPERS if p["key"] == key), None)
        if paper is None:
            raise SystemExit(f"FATAL: unknown paper key {key!r}")
        state = load_state()
        state[key] = new_version(paper, state, vstr, do_publish)
        save_state(state)
        return

    state = load_state()

    print(f"== Phase 1: drafts on {BASE} ==")
    for paper in PAPERS:
        state[paper["key"]] = ensure_draft(paper, state)
        save_state(state)

    if not do_publish:
        print("\nAll five drafts are complete. Pre-reserved DOIs above are "
              "final (they activate on publish).\nRun again with --publish "
              "to publish all five (PERMANENT: files freeze; metadata stays "
              "editable).")
        return

    print("\n== Phase 2: publishing all five ==")
    for paper in PAPERS:
        state[paper["key"]] = publish(paper, state)
        save_state(state)

    print("\n== Final DOI table ==")
    for paper in PAPERS:
        e = state[paper["key"]]
        print(f"  {paper['key']:<20s} {e['doi']}   {e.get('record','')}")


if __name__ == "__main__":
    main()
