#!/usr/bin/env python3
"""Generate publication-refresh inventories from local TeX sources.

This script is intentionally mechanical. It identifies candidate paper roots from
main TeX files under papers/ and book chapter sources under
holographic_circlette_book/. It also inventories legacy standalone TeX papers
from the repository root. The legacy status labels are routing triage, not final
canon judgments.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "publication_refresh" / "inventory"

TOPIC_BY_FOLDER = {
    "overview": "methodology_falsification",
    "foundations": "foundations",
    "born_rule": "foundations",
    "matter_gauge": "matter_gauge_spectroscopy",
    "baryogenesis": "matter_gauge_spectroscopy",
    "constant_ledger": "cosmology_dark_energy",
    "cosmology": "cosmology_dark_energy",
    "boundary_printing": "cosmology_dark_energy",
    "holographic_bridge": "cosmology_dark_energy",
    "dark_sector": "dark_matter_mond_k04",
    "defect_network": "dark_matter_mond_k04",
    "gravity_blackholes": "gravity_horizons_black_holes",
    "blackhole_deep": "gravity_horizons_black_holes",
    "relativity": "relativity_trans_lambda",
    "past_hypothesis": "methodology_falsification",
    "falsification_sheet": "methodology_falsification",
    "canon_update_recent_derivations": "methodology_falsification",
}

BOOK_TOPIC_BY_PART = {
    "part_1_foundations": "foundations",
    "part_2_substrate": "foundations",
    "part_3_standard_model": "matter_gauge_spectroscopy",
    "part_4_emergence": "gravity_horizons_black_holes",
    "part_5_unification": "methodology_falsification",
    "appendices": "book",
    "frontmatter": "book",
    "backmatter": "book",
}

SCAN_EXCLUDE_ROOTS = {
    ".git",
    "papers",
    "holographic_circlette_book",
    "publication_refresh",
}

PDF_ONLY_IGNORE_PARTS = {
    "figures",
    "images",
    "python_code",
}


LEGACY_TRIAGE = {
    "CKM-hierarchy_as_quantum_error_correction": (
        "matter_gauge_spectroscopy", "rewrite", "high",
        "CKM/CP status changed; recovery holonomy remains open"
    ),
    "TCH_Coulomb_law": (
        "matter_gauge_spectroscopy", "rewrite", "medium",
        "Older electrodynamics/gauge-web claim; compare to current photon and LSZ status"
    ),
    "accessible_infor_first": (
        "methodology_falsification", "rewrite", "high",
        "Accessible overview should be replaced by record-reconstruction opening"
    ),
    "ai_methodology": (
        "methodology_falsification", "brief errata", "low",
        "Methodology mostly stable; check PTMS/R1 wording"
    ),
    "alpha_decay": (
        "matter_gauge_spectroscopy", "rewrite", "medium",
        "Nuclear/second-scale sector remains a current open wall"
    ),
    "atomic_shells": (
        "matter_gauge_spectroscopy", "rewrite", "medium",
        "Atomic/nuclear-scale claims need second-scale audit"
    ),
    "baryogenesis_eta": (
        "matter_gauge_spectroscopy", "rewrite", "high",
        "Baryogenesis now depends on R14 billing and R15 CP-holonomy status"
    ),
    "baryon_resonance_ladder": (
        "matter_gauge_spectroscopy", "brief errata", "medium",
        "Matter/spectroscopy mostly held but needs canon status check"
    ),
    "birefringence_skeleton": (
        "relativity_trans_lambda", "rewrite", "medium",
        "Photon/relativity status changed via Gauss-projected and null-chain layers"
    ),
    "black_holes": (
        "gravity_horizons_black_holes", "rewrite", "high",
        "Superseded by current black-hole horizon-QEC papers"
    ),
    "cosmological_qec_engine": (
        "cosmology_dark_energy", "rewrite", "high",
        "Dark-energy/CC route changed; Landauer-exhaust framing needs audit"
    ),
    "crystallize": (
        "dark_matter_mond_k04", "brief errata", "medium",
        "Crystallisation remains important but toy/embedded and K33 artefact status changed"
    ),
    "dark_sector_cosmology": (
        "dark_matter_mond_k04", "rewrite", "high",
        "Superseded by current dark-sector and CMB/zero-mode statuses"
    ),
    "double_split_arXiv": (
        "relativity_trans_lambda", "brief errata", "low",
        "Numerical diffraction/decoherence note; likely mostly independent"
    ),
    "electrodynamics": (
        "relativity_trans_lambda", "rewrite", "high",
        "Photon identity and trans-Lambda representation substantially changed"
    ),
    "emergent_entropic_gravity": (
        "gravity_horizons_black_holes", "withdraw", "high",
        "Old entropic/induced-gravity quantitative route superseded/retracted"
    ),
    "entropy": (
        "foundations", "rewrite", "high",
        "Thermodynamic claims now governed by R13/R14 record-action discipline"
    ),
    "forPhysicists": (
        "methodology_falsification", "rewrite", "high",
        "Broad summary pre-dates record reconstruction and current frontier map"
    ),
    "full_explanation_in_3D": (
        "foundations", "rewrite", "high",
        "Broad model introduction pre-dates forced-byte/reconstruction framing"
    ),
    "gauge_coupling": (
        "matter_gauge_spectroscopy", "rewrite", "high",
        "Gauge-coupling/alpha story changed after R14 billing and dressed-alpha split"
    ),
    "glueball_capstone": (
        "matter_gauge_spectroscopy", "brief errata", "medium",
        "Gauge/spectroscopy claim; check against SMG/TCH continuum status"
    ),
    "glueball_decay_widths": (
        "matter_gauge_spectroscopy", "brief errata", "medium",
        "Gauge/spectroscopy claim; check against SMG/TCH continuum status"
    ),
    "glueball_ladder": (
        "matter_gauge_spectroscopy", "brief errata", "medium",
        "Gauge/spectroscopy programme note; likely targeted update"
    ),
    "gravity": (
        "gravity_horizons_black_holes", "withdraw", "high",
        "Holographic-dilution Planck route superseded by proton-primary/selector framing"
    ),
    "hartman": (
        "methodology_falsification", "rewrite", "medium",
        "Older part-series note; needs canon alignment"
    ),
    "helium": (
        "matter_gauge_spectroscopy", "rewrite", "medium",
        "Nuclear/second-scale sector remains a current open wall"
    ),
    "higgs": (
        "matter_gauge_spectroscopy", "rewrite", "high",
        "EW/Higgs second scale remains open"
    ),
    "information_basis": (
        "foundations", "brief errata", "medium",
        "QEC-code basis likely useful but should reference forced-byte theorem"
    ),
    "information_to_geometry": (
        "foundations", "rewrite", "high",
        "Foundational opening should be recast around record reconstruction"
    ),
    "lindbladian_closure": (
        "foundations", "brief errata", "medium",
        "Open-system/QEC closure likely still relevant; check R13/R14 status"
    ),
    "narrow_higgs": (
        "matter_gauge_spectroscopy", "rewrite", "high",
        "EW/Higgs second scale remains open"
    ),
    "nyquist_hqet_strange": (
        "matter_gauge_spectroscopy", "brief errata", "medium",
        "Matter/spectroscopy mostly held but needs canon status check"
    ),
    "octagonal_honeycomb": (
        "matter_gauge_spectroscopy", "rewrite", "medium",
        "Gauge-geometry claim needs comparison to current TCH/SMG status"
    ),
    "origin_of_mass": (
        "matter_gauge_spectroscopy", "rewrite", "high",
        "Mass/EW second-scale sector remains open"
    ),
    "pauli_antisymmetry": (
        "matter_gauge_spectroscopy", "brief errata", "low",
        "Structural matter result likely stable; verify current wording"
    ),
    "photons": (
        "relativity_trans_lambda", "rewrite", "high",
        "Photon identity and trans-Lambda status changed substantially"
    ),
    "pivot": (
        "foundations", "rewrite", "medium",
        "Bridge note predates forced-byte/current Q3 framing"
    ),
    "sm_matter_from_substrate": (
        "matter_gauge_spectroscopy", "brief errata", "medium",
        "Matter content mostly held; update record-reconstruction references"
    ),
    "spinor_algebra": (
        "relativity_trans_lambda", "brief errata", "medium",
        "Dirac/spinor structure likely useful; check causal-set spin frontier"
    ),
    "strong_cp": (
        "matter_gauge_spectroscopy", "rewrite", "high",
        "CP/holonomy sector remains open; old strong-CP framing needs audit"
    ),
    "strong_gravity": (
        "gravity_horizons_black_holes", "withdraw", "high",
        "Old gravity/mass-hierarchy route superseded by current gravity status"
    ),
    "summary_May_26": (
        "methodology_falsification", "rewrite", "high",
        "Broad summary is stale after record reconstruction and canon updates"
    ),
    "technical_notes": (
        "methodology_falsification", "brief errata", "medium",
        "Technical notes need per-note audit; not merged with current paper set"
    ),
    "velocity_unification_rg": (
        "relativity_trans_lambda", "brief errata", "medium",
        "Velocity-unification result partly recovered; update vacuum-polarisation status"
    ),
    "willow": (
        "relativity_trans_lambda", "rewrite", "medium",
        "Older birefringence/communication note; check photon status"
    ),
    "wilson_strings_markov": (
        "matter_gauge_spectroscopy", "rewrite", "medium",
        "Gauge/string Markov framing needs TCH/SMG current-status audit"
    ),
    "yang_mills_mass_gap": (
        "matter_gauge_spectroscopy", "rewrite", "high",
        "SMG/TCH continuum frontier changed; likely stale"
    ),
}


PART_TRIAGE = {
    "part_01_encoding_dynamics": ("foundations", "rewrite", "high", "Early part-series foundation; superseded by forced-byte/R14 framing"),
    "part_02_composites": ("matter_gauge_spectroscopy", "rewrite", "medium", "Early part-series matter note; needs canon alignment"),
    "part_04_ckm_mixing": ("matter_gauge_spectroscopy", "rewrite", "high", "CKM/CP holonomy status changed"),
    "part_05_neutrino_masses": ("matter_gauge_spectroscopy", "rewrite", "high", "Neutrino/sterile/Majorana status changed"),
    "part_06_pmns_matrix": ("matter_gauge_spectroscopy", "rewrite", "high", "PMNS delta/holonomy status changed"),
    "part_07_anomaly_cancellation": ("matter_gauge_spectroscopy", "brief errata", "medium", "Anomaly arithmetic mostly held; update universal-2/9 split"),
    "part_08_hadron_topology": ("matter_gauge_spectroscopy", "brief errata", "medium", "Hadron topology mostly held; update spectroscopy audit"),
    "part_09_pmns_dactorisation": ("matter_gauge_spectroscopy", "rewrite", "high", "PMNS factorisation status changed"),
    "part_10_algorithmic_inertia": ("foundations", "rewrite", "medium", "Older mechanism note; compare to record-action measure"),
    "part_11_baryon_octet": ("matter_gauge_spectroscopy", "brief errata", "medium", "Baryon octet/spectroscopy mostly held; update status"),
    "part_12_fine _structure": ("matter_gauge_spectroscopy", "rewrite", "high", "Bare alpha0 now derived differently; dressed alpha remains open"),
    "part_13_chiral_pion": ("matter_gauge_spectroscopy", "brief errata", "medium", "Pion derivation mostly held; update status"),
    "part_14_meson_splitting": ("matter_gauge_spectroscopy", "brief errata", "medium", "Meson splitting needs current spectroscopy audit"),
    "part_15_emergent_gravity": ("gravity_horizons_black_holes", "rewrite", "high", "Emergent-geometry/gravity framing superseded"),
    "part_16_dark_energy": ("cosmology_dark_energy", "rewrite", "high", "Dark-energy equation-of-state status changed"),
    "part_17_energy_trajectory": ("cosmology_dark_energy", "rewrite", "high", "Energy trajectory/w(0) route superseded by current selector/CC status"),
    "part_18_feshback_mechanism": ("matter_gauge_spectroscopy", "withdraw", "high", "Feshbach trace route non-constructible/superseded in canon"),
    "part_19_pati_salem": ("matter_gauge_spectroscopy", "rewrite", "high", "Pati-Salam/broken-load status changed"),
    "part_20_planck_mass": ("gravity_horizons_black_holes", "rewrite", "high", "Planck-mass route superseded by proton-primary/selector analysis"),
    "part_21_vector_mesons": ("matter_gauge_spectroscopy", "brief errata", "medium", "Vector-meson spectroscopy needs current audit"),
    "part_22_summary": ("methodology_falsification", "rewrite", "high", "Summary is stale after record reconstruction"),
    "part_23_force_coupling_constants": ("matter_gauge_spectroscopy", "rewrite", "high", "Force-coupling constants and dressed alpha status changed"),
}


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return path.read_text(encoding="latin-1")


def tex_title(text: str) -> str:
    m = re.search(r"\\title(?:\[[^\]]*\])?\s*\{(.+?)\}", text, re.S)
    if m:
        return clean_tex(m.group(1))
    m = re.search(r"\\chapter(?:\[[^\]]*\])?\s*\{(.+?)\}", text, re.S)
    if m:
        return clean_tex(m.group(1))
    m = re.search(r"\\section\*?\s*\{(.+?)\}", text, re.S)
    if m:
        return clean_tex(m.group(1))
    return ""


def clean_tex(value: str) -> str:
    value = re.sub(r"%.*", "", value)
    value = value.replace("\n", " ")
    value = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^\]]*\])?\{([^{}]*)\}", r"\1", value)
    value = re.sub(r"[{}]", "", value)
    value = re.sub(r"\s+", " ", value)
    return value.strip()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def find_pdf(folder: Path, stem: str) -> str:
    candidates = [folder / f"{stem}.pdf"] + sorted(folder.glob("*.pdf"))
    for candidate in candidates:
        if candidate.exists():
            return rel(candidate)
    return ""


def find_bib(folder: Path) -> str:
    candidates = sorted(folder.glob("*.bib"))
    if candidates:
        return ";".join(rel(p) for p in candidates)
    common = ROOT / "papers" / "common" / "references.bib"
    if common.exists():
        return rel(common)
    return ""


def is_main_tex(path: Path, text: str) -> bool:
    if "\\documentclass" not in text:
        return False
    name = path.name.lower()
    if name in {"macros.tex"}:
        return False
    if "/common/" in rel(path):
        return False
    return True


def infer_paper_topic(folder: Path) -> str:
    return TOPIC_BY_FOLDER.get(folder.name, "methodology_falsification")


def legacy_triage_for(path: Path) -> tuple[str, str, str, str]:
    top = path.relative_to(ROOT).parts[0]
    if top in PART_TRIAGE:
        return PART_TRIAGE[top]
    if top.startswith("part_"):
        return (
            "methodology_falsification",
            "rewrite",
            "high",
            "Older part-series paper; needs explicit canon comparison",
        )
    if top in LEGACY_TRIAGE:
        return LEGACY_TRIAGE[top]
    return (
        "methodology_falsification",
        "unreviewed",
        "unrated",
        "Legacy standalone source found by deep scan; no triage rule yet",
    )


def keep_status(value: str, default: str) -> str:
    return value if value and value not in {"unreviewed", "unrated"} else default


def generate_paper_status() -> None:
    existing = read_existing_csv(OUT / "paper_status.csv", key="folder")
    zenodo = read_existing_csv(OUT / "zenodo_records.csv", key="local_folder")
    rows = []
    for tex in sorted((ROOT / "papers").glob("**/*.tex")):
        text = read_text(tex)
        if not is_main_tex(tex, text):
            continue
        folder = tex.parent
        old = existing.get(rel(folder), {})
        z = zenodo.get(rel(folder), {})
        rows.append({
            "topic": infer_paper_topic(folder),
            "folder": rel(folder),
            "title": tex_title(text) or folder.name.replace("_", " ").title(),
            "main_tex": rel(tex),
            "pdf": find_pdf(folder, tex.stem),
            "bib": find_bib(folder),
            "doi": z.get("doi", ""),
            "conceptdoi": z.get("conceptdoi", ""),
            "zenodo_record": z.get("record_id", ""),
            "publication_date": z.get("publication_date", ""),
            "status": keep_status(old.get("status", ""), "unreviewed"),
            "risk": keep_status(old.get("risk", ""), "unrated"),
            "last_reviewed": old.get("last_reviewed", ""),
            "canon_baseline": old.get("canon_baseline", ""),
            "notes": old.get("notes", ""),
        })
    write_csv(OUT / "paper_status.csv", rows)


def generate_legacy_paper_status() -> None:
    existing = read_existing_csv(OUT / "legacy_paper_status.csv", key="main_tex")
    zenodo = read_existing_csv(OUT / "zenodo_records.csv", key="local_folder")
    rows = []
    for tex in sorted(ROOT.glob("**/*.tex")):
        rel_parts = tex.relative_to(ROOT).parts
        if not rel_parts or rel_parts[0] in SCAN_EXCLUDE_ROOTS:
            continue
        text = read_text(tex)
        if not is_main_tex(tex, text):
            continue
        folder = tex.parent
        old = existing.get(rel(tex), {})
        z = zenodo.get(rel(tex), {}) or zenodo.get(rel(folder), {})
        topic, status, risk, default_note = legacy_triage_for(tex)
        rows.append({
            "topic": old.get("topic") or topic,
            "folder": rel(folder),
            "title": tex_title(text) or folder.name.replace("_", " ").title(),
            "main_tex": rel(tex),
            "pdf": find_pdf(folder, tex.stem),
            "bib": find_bib(folder),
            "doi": z.get("doi", ""),
            "conceptdoi": z.get("conceptdoi", ""),
            "zenodo_record": z.get("record_id", ""),
            "publication_date": z.get("publication_date", ""),
            "status": keep_status(old.get("status", ""), status),
            "risk": keep_status(old.get("risk", ""), risk),
            "last_reviewed": old.get("last_reviewed") or "2026-06-20",
            "canon_baseline": old.get("canon_baseline") or "rolling-2026-06-20",
            "notes": old.get("notes") or default_note,
        })
    write_csv(OUT / "legacy_paper_status.csv", rows)


def generate_book_chapter_status() -> None:
    existing = read_existing_csv(OUT / "book_chapter_status.csv", key="source")
    rows = []
    book_root = ROOT / "holographic_circlette_book"
    for tex in sorted(book_root.glob("**/*.tex")):
        if tex.name in {"book.tex", "bookkindle.tex"}:
            continue
        text = read_text(tex)
        part = tex.relative_to(book_root).parts[0]
        old = existing.get(rel(tex), {})
        rows.append({
            "topic": BOOK_TOPIC_BY_PART.get(part, "book"),
            "part": part,
            "source": rel(tex),
            "title": tex_title(text) or tex.stem.replace("_", " ").title(),
            "status": keep_status(old.get("status", ""), "unreviewed"),
            "risk": keep_status(old.get("risk", ""), "unrated"),
            "last_reviewed": old.get("last_reviewed", ""),
            "canon_baseline": old.get("canon_baseline", ""),
            "notes": old.get("notes", ""),
        })
    write_csv(OUT / "book_chapter_status.csv", rows)


def generate_pdf_only_status() -> None:
    existing = read_existing_csv(OUT / "pdf_only_status.csv", key="pdf")
    rows = []
    tex_pdfs = {row["pdf"] for row in read_existing_csv(OUT / "paper_status.csv", "pdf").values() if row.get("pdf")}
    legacy_pdfs = {row["pdf"] for row in read_existing_csv(OUT / "legacy_paper_status.csv", "pdf").values() if row.get("pdf")}
    book_pdfs = {"holographic_circlette_book/book.pdf", "holographic_circlette_book/bookkindle.pdf"}
    known = tex_pdfs | legacy_pdfs | book_pdfs
    for pdf in sorted(ROOT.glob("**/*.pdf")):
        rel_pdf = rel(pdf)
        parts = pdf.relative_to(ROOT).parts
        if rel_pdf in known:
            continue
        if any(part in PDF_ONLY_IGNORE_PARTS for part in parts):
            continue
        if parts[0] == "papers":
            continue
        old = existing.get(rel_pdf, {})
        z = read_existing_csv(OUT / "zenodo_records.csv", key="local_folder").get(rel_pdf, {})
        rows.append({
            "topic": old.get("topic") or "pdf_only",
            "pdf": rel_pdf,
            "doi": z.get("doi", ""),
            "conceptdoi": z.get("conceptdoi", ""),
            "zenodo_record": z.get("record_id", ""),
            "publication_date": z.get("publication_date", ""),
            "status": keep_status(old.get("status", ""), "unreviewed"),
            "risk": keep_status(old.get("risk", ""), "unrated"),
            "last_reviewed": old.get("last_reviewed") or "",
            "canon_baseline": old.get("canon_baseline") or "",
            "notes": old.get("notes") or "PDF-only candidate found by deep scan; source not identified in this pass",
        })
    write_csv(OUT / "pdf_only_status.csv", rows)


def read_existing_csv(path: Path, key: str) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    return {row.get(key, ""): row for row in rows if row.get(key)}


def generate_zenodo_placeholder() -> None:
    path = OUT / "zenodo_records.csv"
    if path.exists() and path.stat().st_size > 0:
        return
    rows = []
    write_csv(path, rows, fieldnames=[
        "title", "doi", "conceptdoi", "record_id", "conceptrecid",
        "version", "publication_date", "local_folder", "notes"
    ])


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    generate_paper_status()
    generate_legacy_paper_status()
    generate_book_chapter_status()
    generate_pdf_only_status()
    generate_zenodo_placeholder()
    print(f"wrote {OUT / 'paper_status.csv'}")
    print(f"wrote {OUT / 'legacy_paper_status.csv'}")
    print(f"wrote {OUT / 'book_chapter_status.csv'}")
    print(f"wrote {OUT / 'pdf_only_status.csv'}")
    print(f"wrote {OUT / 'zenodo_records.csv'}")


if __name__ == "__main__":
    main()
