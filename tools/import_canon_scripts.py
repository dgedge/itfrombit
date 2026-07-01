#!/usr/bin/env python3
"""
Import canon source scripts from the working canon tree into this source-code repo.

The importer is deliberately conservative:
- copies only source-like files, not run outputs or numerical artifacts;
- preserves original filenames and subpaths where useful;
- classifies scripts into broad paper/topic folders;
- writes per-topic README files and a manifest with purpose/run commands.

Run from the repo root:
    python tools/import_canon_scripts.py

Optionally pass a different canon source root.  The public `scripts/` mirror is
normally refreshed from the canonical script directory itself:
    python tools/import_canon_scripts.py /path/to/octahedrons/python_code
"""

from __future__ import annotations

import ast
import csv
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


DEFAULT_SOURCE = Path(
    "/Users/davidelliman/Library/Mobile Documents/com~apple~CloudDocs/octahedrons"
)

SOURCE_SUFFIXES = {".py", ".sh", ".ipynb", ".jl", ".m", ".jsx", ".ts", ".tsx"}

EXCLUDED_DIR_PARTS = {
    ".git",
    ".agents",
    ".claude",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".venv",
    "venv",
    "node_modules",
    "dist",
    "build",
    "smg_dmrg_runs",
    "results",
    "submission",
}

EXCLUDED_FILE_PATTERNS = (
    re.compile(r".*\.min\.js$"),
)

TOPIC_TITLES = {
    "foundations_methodology": "Foundations and Methodology",
    "matter_gauge_spectroscopy": "Matter, Gauge Structure, and Spectroscopy",
    "smg_tch_gauge_cell": "SMG, TCH Gauge Cell, and Lattice Certification",
    "cosmology_dark_energy_inflation": "Cosmology, Dark Energy, and Inflation",
    "dark_matter_mond_k04": "Dark Matter, MOND, and K04 Debris",
    "gravity_horizons_blackholes": "Gravity, Horizons, and Black Holes",
    "relativity_photon_qed": "Relativity, Photon, and Causal-Set QED",
    "legacy_misc": "Legacy and Miscellaneous Source",
}

TOPIC_DESCRIPTIONS = {
    "foundations_methodology": (
        "Finite-record foundations, code/ledger audits, claim-status machinery, "
        "and reproducibility infrastructure."
    ),
    "matter_gauge_spectroscopy": (
        "Standard Model register checks, anomaly arithmetic, hadron/lepton "
        "spectroscopy, and gauge/spectroscopy support scripts."
    ),
    "smg_tch_gauge_cell": (
        "TCH gauge-cell, Peter-Weyl/CG, mirror-Fock, SMG gap, and sparse-lattice "
        "certification scripts."
    ),
    "cosmology_dark_energy_inflation": (
        "Cosmological-constant queue route, boot residuals, HBC inflation, "
        "printer saturation, and large-scale audits."
    ),
    "dark_matter_mond_k04": (
        "Embedded K04 crystallisation, debris relics, orphan/island policy, "
        "R4/MOND line-current checks, and halo audits."
    ),
    "gravity_horizons_blackholes": (
        "Gravity input/output ledgers, horizon accounting, black-hole isometry, "
        "Hawking ladder, and Bekenstein-channel audits."
    ),
    "relativity_photon_qed": (
        "Photon identity, Lorentz/clock limits, causal-set QED, trans-QCD "
        "service-packet checks, and Maxwell-loop scaffolds."
    ),
    "legacy_misc": (
        "Older exploratory scripts and material not yet classified into the "
        "paper-aligned topics."
    ),
}


@dataclass(frozen=True)
class ScriptRecord:
    source: Path
    destination: Path
    topic: str
    purpose: str
    run: str
    has_docstring: bool


def should_copy(path: Path) -> bool:
    if path.suffix not in SOURCE_SUFFIXES:
        return False
    if any(part in EXCLUDED_DIR_PARTS for part in path.parts):
        return False
    return not any(pattern.fullmatch(path.name) for pattern in EXCLUDED_FILE_PATTERNS)


def extract_docstring(path: Path) -> str | None:
    if path.suffix != ".py":
        return None
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    doc = ast.get_docstring(tree)
    if not doc:
        return None
    first_para = doc.strip().split("\n\n", 1)[0]
    first_para = " ".join(line.strip() for line in first_para.splitlines())
    return first_para[:240]


def words_from_name(path: Path) -> str:
    stem = path.stem
    stem = re.sub(r"^(item|part)_?0*", r"item ", stem)
    stem = stem.replace("_", " ").replace("-", " ")
    stem = re.sub(r"\s+", " ", stem).strip()
    return stem[:1].upper() + stem[1:] if stem else path.name


def classify(rel: Path) -> str:
    parts = [part.lower() for part in rel.parts]
    name = rel.name.lower()
    rel_s = str(rel).lower()

    if parts[0] == "ai_methodology" or name.startswith("zenodo_") or "ptms" in parts:
        return "foundations_methodology"

    if rel_s.startswith("python_code/smg_dmrg/") or rel_s.startswith("python_code/smg_validation/"):
        return "smg_tch_gauge_cell"

    if any(token in rel_s for token in (
        "css_", "smg", "tch", "peter_weyl", "mirror_fock", "gauss_proxy",
        "wilson_twist", "openlink", "qutrit", "lanczos", "krylov",
        "string_tension", "luscher", "gauge_dressing", "gap_under_twist",
        "matter_coupled_proxy", "mirror_invalid_subspace", "bulk_domainwall",
        "hgp_gauss", "gw_", "ginsparg", "overlap_gw",
        "record_grammar_cubic_su3", "record_grammar_dual_cubic_su3",
        "record_grammar_su3", "record_grammar_tch", "polyakov", "creutz",
        "steiner", "colour_singlet", "hadron_taxonomy", "face_action",
        "z3_", "su2_gauss", "sc_gauge_web", "tqo_tee",
    )):
        return "smg_tch_gauge_cell"

    if any(token in rel_s for token in (
        "causet", "photon", "trans_lambda", "null_causal", "bn_service",
        "bundle_irreducibility", "recoil_vertex", "hodge", "maxwell",
        "item102", "item115", "velocity", "graphene",
        "relativity_", "birefringence", "wavefront", "dressed_alpha",
        "vacpol", "alpha_graph_and_zeta", "zeta_lattice",
        "walk_band", "walk_kernel", "spin_coin", "spin_statistics",
    )):
        return "relativity_photon_qed"

    if any(token in rel_s for token in (
        "bh_", "black_hole", "hawking", "schwarzschild", "isometry",
        "bekenstein", "horizon", "gravity_", "g_route", "planck", "alpha_power",
        "delta_unification", "proton_anchor_g", "mp_", "item7_", "keff_",
        "r4_activation_identification", "emergent_einstein", "modular_hamiltonian",
        "qec_echo", "substrate_rt_wedge", "item_geo_a_curvature",
    )):
        return "gravity_horizons_blackholes"

    if any(token in rel_s for token in (
        "item132", "debris", "dark_matter", "dark_sector", "mond", "halo",
        "depinning", "k04", "orphan", "island", "bridge2_energy",
        "d_to_p", "foundation_annealing", "instrumentation_onset",
        "pairing_from_dynamics", "onset_alignment", "defect_", "item118",
        "sterile_release", "willow_q4", "ulx_", "advection_nonlocal",
        "wall_syndrome",
    )):
        return "dark_matter_mond_k04"

    if any(token in rel_s for token in (
        "cc_", "cosmology", "lambda", "item123", "item131", "inflation",
        "hbc", "printer", "saturation", "gamma_tee", "boot_residual",
        "structure_wz", "w_to_28", "register_handoff", "operator_map",
        "large_scale_dirac", "ledger_sky", "lroute_pth", "noise_annealing",
        "a1_", "areq", "boundary_printing", "boltzmann_brain", "cmb_",
        "cosmological_selector", "correlated_nulls", "gamma_mon",
        "handoff_action", "item057", "item144", "m5_57", "w0_ns",
    )):
        return "cosmology_dark_energy_inflation"

    if any(token in rel_s for token in (
        "anomaly", "baryon", "ckm", "decuplet", "glueball", "koide",
        "hadron", "spectral", "amu", "g_factor", "alpha_ds", "item086",
        "item126", "item48", "item75", "boson_mass", "charm", "bottom",
        "feshbach", "pseudospin", "su6", "charge", "pion", "atomic",
        "omega_", "node_trace", "item79", "l1_eigenbasis", "ew_", "item55",
        "item55f", "item87", "item97", "item99", "item113", "fermion_",
        "heavy_quark", "neutrino_sector", "so10", "sm_singlet", "v_phase",
        "v_program", "verify_neutrino", "verify_rk_lfu", "finalboss",
        "item_geo_a_gaugeweb", "bridge_web",
    )):
        return "matter_gauge_spectroscopy"

    if any(token in rel_s for token in (
        "foundation", "foundations", "record", "audit", "methodology",
        "crackpot", "input_count", "axis_booking", "engine_readout",
        "e1_", "e2_", "d2_consumption", "equipartition", "numerology_baseline",
        "alpha0_", "constant_ledger", "continuum_lift", "entropy_arrow",
        "measurement_service", "qca_", "r12_", "r13_", "r15_", "r4_complex",
        "r4_auxiliary", "substrate_", "target_a_", "two_sector",
        "unified_prediction", "precision_phenomenology", "scheduler_alpha",
    )):
        return "foundations_methodology"

    return "legacy_misc"


def destination_for(repo_root: Path, source_root: Path, path: Path) -> tuple[str, Path]:
    rel = path.relative_to(source_root)
    topic = classify(rel)

    if rel.parts[0] in {"python_code", "ai_methodology", "papers", "willow"}:
        rel_tail = Path(*rel.parts[1:]) if len(rel.parts) > 1 else Path(rel.name)
    else:
        rel_tail = rel

    # Keep package/subdirectory structure for multi-file modules, but avoid
    # nesting everything under a redundant python_code directory.
    destination = repo_root / "scripts" / topic / rel_tail
    return topic, destination


def run_command(destination: Path, repo_root: Path, topic: str) -> str:
    rel = destination.relative_to(repo_root)
    if destination.suffix == ".py":
        return f"PYTHONPATH=scripts/{topic}:$PYTHONPATH python {rel.as_posix()}"
    if destination.suffix == ".sh":
        return f"bash {rel.as_posix()}"
    if destination.suffix == ".ipynb":
        return f"jupyter notebook {rel.as_posix()}"
    return f"# inspect or run with the appropriate tool: {rel.as_posix()}"


def collect(source_root: Path, repo_root: Path) -> list[ScriptRecord]:
    records: list[ScriptRecord] = []
    for path in sorted(source_root.rglob("*")):
        if not path.is_file() or not should_copy(path):
            continue
        topic, destination = destination_for(repo_root, source_root, path)
        doc = extract_docstring(path)
        purpose = doc or f"{words_from_name(path)}. Purpose inferred from filename; needs a fuller script docstring."
        records.append(
            ScriptRecord(
                source=path,
                destination=destination,
                topic=topic,
                purpose=purpose,
                run=run_command(destination, repo_root, topic),
                has_docstring=bool(doc),
            )
        )
    return records


def write_topic_readmes(repo_root: Path, records: list[ScriptRecord]) -> None:
    by_topic: dict[str, list[ScriptRecord]] = {topic: [] for topic in TOPIC_TITLES}
    for record in records:
        by_topic.setdefault(record.topic, []).append(record)

    for topic, topic_records in by_topic.items():
        topic_dir = repo_root / "scripts" / topic
        topic_dir.mkdir(parents=True, exist_ok=True)
        readme = topic_dir / "README.md"
        rows = [
            f"# {TOPIC_TITLES[topic]}",
            "",
            TOPIC_DESCRIPTIONS[topic],
            "",
            "Each row gives the migrated script, its current description, and a direct run command.",
            "Descriptions marked as inferred should be refined as the source code is curated.",
            "",
            "| Script | Description | Run |",
            "|---|---|---|",
        ]
        for record in sorted(topic_records, key=lambda r: r.destination.as_posix()):
            rel = record.destination.relative_to(topic_dir)
            purpose = record.purpose.replace("|", "\\|")
            run = record.run.replace("|", "\\|")
            rows.append(f"| `{rel.as_posix()}` | {purpose} | `{run}` |")
        readme.write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_manifest(repo_root: Path, source_root: Path, records: list[ScriptRecord]) -> None:
    docs_dir = repo_root / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)
    manifest = docs_dir / "script_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "topic",
                "script",
                "source_path",
                "description",
                "run_command",
                "has_docstring",
            ]
        )
        for record in sorted(records, key=lambda r: (r.topic, r.destination.as_posix())):
            writer.writerow(
                [
                    record.topic,
                    record.destination.relative_to(repo_root).as_posix(),
                    record.source.relative_to(source_root).as_posix(),
                    record.purpose,
                    record.run,
                    "yes" if record.has_docstring else "no",
                ]
            )


def write_scripts_readme(repo_root: Path, records: list[ScriptRecord]) -> None:
    counts: dict[str, int] = {}
    missing_docstrings: dict[str, int] = {}
    for record in records:
        counts[record.topic] = counts.get(record.topic, 0) + 1
        if not record.has_docstring:
            missing_docstrings[record.topic] = missing_docstrings.get(record.topic, 0) + 1

    lines = [
        "# Canon Source Scripts",
        "",
        "This directory contains a curated migration of source scripts from the working canon tree.",
        "It excludes bulky numerical outputs, cached files, and run directories.",
        "",
        "The source-of-truth manifest is [`../docs/script_manifest.csv`](../docs/script_manifest.csv).",
        "",
        "| Topic | Scripts | Needs docstring pass |",
        "|---|---:|---:|",
    ]
    for topic in TOPIC_TITLES:
        lines.append(
            f"| [{TOPIC_TITLES[topic]}]({topic}/) | {counts.get(topic, 0)} | {missing_docstrings.get(topic, 0)} |"
        )
    lines.extend(
        [
            "",
            "Run commands in the topic READMEs assume the repository root as the working directory.",
            "Many scripts are audit/reproduction scripts rather than reusable libraries; read the topic README and script docstring before interpreting output.",
            "",
            "Large run outputs should not be committed here. Put regenerated artifacts in `results/` only when they are small, named, and documented.",
        ]
    )
    (repo_root / "scripts" / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_records(records: list[ScriptRecord]) -> None:
    for record in records:
        record.destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(record.source, record.destination)


def clean_managed_script_dirs(repo_root: Path) -> None:
    scripts_root = repo_root / "scripts"
    scripts_root.mkdir(exist_ok=True)
    for topic in TOPIC_TITLES:
        topic_dir = scripts_root / topic
        if topic_dir.exists():
            shutil.rmtree(topic_dir)
    readme = scripts_root / "README.md"
    if readme.exists():
        readme.unlink()


def main(argv: list[str]) -> int:
    repo_root = Path.cwd()
    source_root = Path(argv[1]).expanduser() if len(argv) > 1 else DEFAULT_SOURCE
    source_root = source_root.resolve()
    if not source_root.exists():
        raise SystemExit(f"source root does not exist: {source_root}")
    if not (repo_root / ".git").exists():
        raise SystemExit("run from the itfrombit repository root")

    records = collect(source_root, repo_root)
    clean_managed_script_dirs(repo_root)
    copy_records(records)
    write_topic_readmes(repo_root, records)
    write_scripts_readme(repo_root, records)
    write_manifest(repo_root, source_root, records)

    print(f"Imported {len(records)} source files from {source_root}")
    for topic in TOPIC_TITLES:
        count = sum(1 for record in records if record.topic == topic)
        missing = sum(1 for record in records if record.topic == topic and not record.has_docstring)
        print(f"{topic:38s} {count:4d} scripts, {missing:4d} inferred descriptions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
