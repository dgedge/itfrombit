#!/usr/bin/env python3
r"""Build the unified prediction/falsification ledger.

The existing falsification sheet is a prose attack surface.  This script turns
the same material into a ranked triage sheet: every row names the canonical
prediction, the current experimental status, the observation that would kill or
demote it, and the damage scope.

This is not a derivation script.  It is public-facing bookkeeping, intended to
keep outreach claims honest and to make the most decisive tests visible.

Run from the repository root:

    ~/bin/py13_7/bin/python python_code/unified_prediction_falsification_ledger.py
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from neutrino_sector_falsification_ledger import (
    active_mass_readouts,
    sterile_readouts,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "papers" / "unified_prediction_falsification_ledger.md"
DATE = "2026-06-15"


@dataclass(frozen=True)
class Entry:
    sector: str
    observable: str
    prediction: str
    current_status: str
    falsifier: str
    scope: str
    horizon: str
    decisiveness: int
    near_term: int
    pressure: int
    confidence: int
    source: str
    canon_hook: str
    note: str = ""

    @property
    def score(self) -> int:
        return (
            3 * self.decisiveness
            + 2 * self.near_term
            + 2 * self.pressure
            + self.confidence
        )

    @property
    def band(self) -> str:
        if self.score >= 34:
            return "kill switch"
        if self.score >= 30:
            return "high priority"
        if self.score >= 25:
            return "watch closely"
        if self.score >= 20:
            return "longer range"
        return "background"


def md_escape(text: str) -> str:
    return text.replace("|", r"\|").replace("\n", "<br>")


def entry_table(entries: list[Entry]) -> str:
    header = (
        "| rank | score | band | sector | observable | canonical prediction | "
        "current status | falsifier | scope |\n"
        "|---:|---:|---|---|---|---|---|---|---|\n"
    )
    rows: list[str] = []
    for rank, entry in enumerate(entries, start=1):
        rows.append(
            "| {rank} | {score} | {band} | {sector} | {observable} | {prediction} | "
            "{status} | {falsifier} | {scope} |".format(
                rank=rank,
                score=entry.score,
                band=md_escape(entry.band),
                sector=md_escape(entry.sector),
                observable=md_escape(entry.observable),
                prediction=md_escape(entry.prediction),
                status=md_escape(entry.current_status),
                falsifier=md_escape(entry.falsifier),
                scope=md_escape(entry.scope),
            )
        )
    return header + "\n".join(rows) + "\n"


def compact_table(entries: list[Entry], predicate) -> str:
    rows = [entry for entry in entries if predicate(entry)]
    if not rows:
        return "_No rows in this category._\n"
    header = "| observable | prediction | current status | falsifier |\n|---|---|---|---|\n"
    body = "\n".join(
        "| {obs} | {pred} | {status} | {falsifier} |".format(
            obs=md_escape(entry.observable),
            pred=md_escape(entry.prediction),
            status=md_escape(entry.current_status),
            falsifier=md_escape(entry.falsifier),
        )
        for entry in rows
    )
    return header + body + "\n"


def sources_table(entries: list[Entry]) -> str:
    header = "| observable | source / status hook | canon hook |\n|---|---|---|\n"
    body = "\n".join(
        "| {obs} | {source} | {hook} |".format(
            obs=md_escape(entry.observable),
            source=md_escape(entry.source),
            hook=md_escape(entry.canon_hook),
        )
        for entry in entries
    )
    return header + body + "\n"


def build_entries() -> list[Entry]:
    active = active_mass_readouts()
    sterile = sterile_readouts()

    sum_m_mev = active["sum_m"] * 1.0e3
    mbb_min = active["mbb_min"] * 1.0e3
    mbb_max = active["mbb_max"] * 1.0e3
    m_beta_mev = active["m_beta"] * 1.0e3
    line_kev = sterile["line_kev"]
    sin2_2theta = sterile["sin2_2theta"]

    return [
        Entry(
            sector="cosmology / dark energy",
            observable="late dark-energy equation of state",
            prediction="w(a) = -1 + a/28; w0 = -27/28; wa = -1/28",
            current_status=(
                "DESI DR2 reports robust low-z preference for evolving dark energy "
                "in combined BAO+CMB+SN analyses; current pressure is high but not settled."
            ),
            falsifier=(
                "robust reconstruction excluding the line, especially a confirmed phantom "
                "epoch without a derived negative-rate channel"
            ),
            scope="kills or rewrites R4 late dark-energy branch",
            horizon="now / DESI DR2-DR3",
            decisiveness=5,
            near_term=5,
            pressure=5,
            confidence=4,
            source="DESI DR2 guide: https://www.desi.lbl.gov/2025/03/19/desi-dr2-results-march-19-guide/",
            canon_hook="papers/falsification_sheet.md; ANCHOR item 123/131 cluster",
        ),
        Entry(
            sector="CMB / R4 completion",
            observable="early R4 equation of state and CMB fit",
            prediction=(
                "R4 support resolves to d=1, w=-1/3 for the line branch; it is not early dust "
                "unless a distinct derived mode is added"
            ),
            current_status=(
                "internal consistency gate is sharp after the R4 EoS audit; LCDM CMB remains "
                "the external standard that any completion must reproduce"
            ),
            falsifier=(
                "failure to derive an AeST-like/dust-like early component while keeping the "
                "late R4 line branch intact"
            ),
            scope="major cosmology-sector incompleteness, not a matter-sector failure",
            horizon="internal theorem before public CMB claim",
            decisiveness=5,
            near_term=4,
            pressure=5,
            confidence=5,
            source="internal scripts: python_code/r4_eos_cmb_resolution.py; python_code/item123_cmb_completion.py",
            canon_hook="DRIFT M15; ANCHOR item 123/132 consistency notes",
        ),
        Entry(
            sector="dark matter / halo phenomenology",
            observable="mobile halo mass versus pinned K04 debris",
            prediction=(
                "K04 debris is pinned, gauge-blind, substrate-static; mobile halo phenomenology "
                "is charged to the R4 zero-mode reservoir, active R4/MOND response, and/or nu_R"
            ),
            current_status=(
                "cluster mergers and galaxy halos require mobile gravitational structure; "
                "K04 is therefore an upper-bound/null component, not the main moving halo"
            ),
            falsifier=(
                "dominant claimed K04 component required to follow galaxies through mergers, "
                "or no R4 zero-mode/nu_R component can carry the observed mobile mass"
            ),
            scope="kills K04-as-halo branch; may leave R4 zero-mode/nu_R dark sector alive",
            horizon="now, with cluster and strong-lensing reanalyses",
            decisiveness=4,
            near_term=5,
            pressure=5,
            confidence=4,
            source="dark-defect guide: papers/defect_network/dark_defect_phenomenology.md",
            canon_hook="ANCHOR K04/debris and item 132 MOND notes",
        ),
        Entry(
            sector="gravity / proton-primary route",
            observable="G and H0 from the proton anchor",
            prediction=(
                "m_p fixes Lambda_p; route outputs G and H0 = 67.27 km/s/Mpc "
                "rather than consuming H0"
            ),
            current_status=(
                "Planck base-LCDM gives H0 = 67.36 +/- 0.54, consistent; "
                "SH0ES gives 73.04 +/- 1.04, strongly discrepant with the route"
            ),
            falsifier=(
                "resolved alpha-convention plus precision G/H0/proton audit outside the predicted "
                "window; a settled SH0ES-side H0 without a derived late-universe correction is fatal"
            ),
            scope="major gravity/cosmology route failure; matter code may survive",
            horizon="H0 tension resolution and G precision",
            decisiveness=5,
            near_term=4,
            pressure=4,
            confidence=4,
            source=(
                "Planck H0: https://www.aanda.org/articles/aa/pdf/2020/09/aa33910-18.pdf; "
                "SH0ES H0: https://arxiv.org/abs/2112.04510"
            ),
            canon_hook="python_code/proton_anchor_g_prediction.py; python_code/g_route_input_ledger.py",
        ),
        Entry(
            sector="inflation / boundary printing",
            observable="primordial tensor-to-scalar ratio",
            prediction=(
                "r_linear = 0; scalar-induced floor r_induced ~ 2e-9; "
                "no observable primordial B-modes"
            ),
            current_status=(
                "BICEP/Keck gives r_0.05 < 0.036; LiteBIRD targets delta r = 0.001 "
                "and CMB-S4 aims at r < 0.001 if no detection"
            ),
            falsifier=(
                "robust primordial B-mode detection at r >= 1e-3 after dust, lensing, "
                "and systematics are removed"
            ),
            scope="kills scalar boundary-printer / no-squeezing branch",
            horizon="LiteBIRD / CMB-S4 era",
            decisiveness=5,
            near_term=5,
            pressure=1,
            confidence=5,
            source=(
                "BICEP/Keck r<0.036: https://bicepkeck.org/figures/bk-xiii.html; "
                "LiteBIRD delta r=0.001: https://arxiv.org/abs/2406.02724; "
                "CMB-S4 r reach: https://arxiv.org/abs/2008.12619"
            ),
            canon_hook=(
                "recent_papers/boundary_printing/boundary_printing.tex; "
                "python_code/boundary_printing_tensor_prediction_audit.py"
            ),
        ),
        Entry(
            sector="gravitational waves / crystallisation",
            observable="PTA-band stochastic background from the boot transition",
            prediction=(
                "Lambda_QCD-scale crystallisation peaks in the PTA band (~20-50 nHz), "
                "but pinned/vacuum-grade K04 walls make the amplitude sub-threshold"
            ),
            current_status=(
                "NANOGrav 15-year data show evidence for an all-sky gravitational-wave background; "
                "current interpretation is compatible with an astrophysical SMBH-binary background"
            ),
            falsifier=(
                "a strong cosmological tens-of-nHz background at NANOGrav amplitude with the wrong "
                "suppression/signature, or a confirmed cosmological peak only in LISA/LIGO bands"
            ),
            scope="kills the K04 pinning/GW-amplitude reading; frequency scale is separate",
            horizon="PTA spectral-shape discrimination",
            decisiveness=3,
            near_term=4,
            pressure=3,
            confidence=4,
            source="NANOGrav 15-year summary: https://nanograv.org/15yr/Summary/Background",
            canon_hook="python_code/gw_stochastic_background_nanograv.py; dark-defect canon",
        ),
        Entry(
            sector="foundations / Lorentz bridge",
            observable="trans-Lambda_QCD high-energy quanta",
            prediction=(
                "TeV photons are not lattice Bloch modes; they require framed causal-set/null-chain "
                "service events without a UV oscillator vacuum tower"
            ),
            current_status=(
                "ordinary high-energy QED phenomenology is observed; the framework now has a "
                "normalized null-chain endpoint/interface audit, while precision process and "
                "astrophysical-transfer calculations remain standard-EFT work"
            ),
            falsifier=(
                "a required high-energy QED amplitude or propagation effect that cannot be "
                "represented with the normalized null-chain external leg without introducing a "
                "fine UV oscillator tower and reopening the CC problem"
            ),
            scope="global single-scale / cosmological-constant compatibility bridge",
            horizon="theory gate, then gamma-ray/QED phenomenology",
            decisiveness=5,
            near_term=2,
            pressure=5,
            confidence=3,
            source="canon notes around trans-Lambda_QCD null-chain QED and normalized-leg audits",
            canon_hook="ANCHOR item 150; DRIFT K12; python_code/foundations_trans_lambda_qed_phenomenology.py",
        ),
        Entry(
            sector="MOND / R4 line current",
            observable="RAR, BTFR, cored profiles, Jeans consistency",
            prediction=(
                "R4 one-dimensional line support gives BTFR/RAR-like behaviour only if the "
                "matched-rate, nonexclusive line-occupancy lemma holds"
            ),
            current_status=(
                "galaxy phenomenology is promising but structurally gated; CMB and cluster "
                "completion charged away from K04"
            ),
            falsifier=(
                "matched-rate/nonexclusive R4 Kraus theorem fails, or high-quality halo data "
                "reject the predicted line-density law in its claimed regime"
            ),
            scope="kills MOND/R4 halo branch, not automatically gravity or matter",
            horizon="near-term theory plus galaxy-data tests",
            decisiveness=4,
            near_term=4,
            pressure=4,
            confidence=4,
            source="papers/defect_network/dark_defect_phenomenology.md; MOND/R4 canon scripts",
            canon_hook="ANCHOR item 132",
        ),
        Entry(
            sector="inflation / HBC",
            observable="scalar spectral index",
            prediction="n_s = 27/28 = 0.9642857",
            current_status="Planck 2018 gives n_s = 0.9649 +/- 0.0042, so the value is currently consistent.",
            falsifier="future CMB/LSS result excluding 0.9643 after running, foreground, and model-extension checks",
            scope="kills HBC mode-local scalar-clock tilt branch",
            horizon="CMB-S4 / LSS precision era",
            decisiveness=4,
            near_term=4,
            pressure=1,
            confidence=5,
            source="Planck inflation constraints: https://inspirehep.net/literature/1682899",
            canon_hook="papers/boundary_printing/boundary_printing.tex; ANCHOR item 131",
        ),
        Entry(
            sector="strong CP / correlated nulls",
            observable="neutron electric dipole moment",
            prediction="d_n around 1e-31 e cm",
            current_status="current best limit is |d_n| < 1.8e-26 e cm; next searches target lower orders.",
            falsifier="confirmed d_n >> 1e-30 e cm without an added CP source",
            scope="kills bare strong-CP phase-null closure",
            horizon="next-generation nEDM searches",
            decisiveness=4,
            near_term=4,
            pressure=1,
            confidence=5,
            source="nEDM status: https://exnuc.wordpress.ncsu.edu/neutron-edm/",
            canon_hook="python_code/correlated_nulls_ledger.py; strong-CP canon notes",
        ),
        Entry(
            sector="constant ledger",
            observable="secular alpha drift",
            prediction="alpha_dot/alpha = 0 for code-level alpha; ordinary energy running only",
            current_status=(
                "clock benchmark gives alpha_dot/alpha = (-1.6 +/- 2.3)e-17 / yr; "
                "no confirmed secular drift"
            ),
            falsifier="clock or astrophysical drift requiring a changing code coupling",
            scope="global code-ledger failure",
            horizon="ongoing clock networks",
            decisiveness=5,
            near_term=3,
            pressure=1,
            confidence=5,
            source="Rosenband et al. clock ratio: https://pubmed.ncbi.nlm.nih.gov/18323415/",
            canon_hook="python_code/correlated_nulls_ledger.py; constant-ledger canon",
        ),
        Entry(
            sector="constant ledger",
            observable="secular G drift",
            prediction="G_dot/G = 0 after lock-in",
            current_status=(
                "LLR-scale analyses report no confirmed drift; published limits are model-dependent "
                "and at roughly 1e-13 to 1e-14 / yr order"
            ),
            falsifier="reproducible secular G_dot/G after ephemeris, solar-mass-loss, and local modelling checks",
            scope="global frozen-gravity ledger failure",
            horizon="ongoing LLR/ephemeris improvements",
            decisiveness=5,
            near_term=3,
            pressure=1,
            confidence=4,
            source="LLR review/status: https://pmc.ncbi.nlm.nih.gov/articles/PMC5253913/",
            canon_hook="python_code/correlated_nulls_ledger.py; constant-ledger canon",
        ),
        Entry(
            sector="neutrino / sterile branch",
            observable="17.7 keV sterile-neutrino X-ray line",
            prediction=(
                f"m_nuR = {2.0 * line_kev:.2f} keV; E_gamma = {line_kev:.2f} keV; "
                f"sin^2(2theta) = {sin2_2theta:.1e}"
            ),
            current_status=(
                "XRISM cluster stacking now constrains unidentified 2.5-15 keV lines; "
                "no required bright line is established"
            ),
            falsifier=(
                "X-ray bounds excluding any required 17.7 keV abundance/mixing, or a secure "
                "incompatible dark-sector line"
            ),
            scope="kills sterile-nu_R dark component, not whole dark sector",
            horizon="XRISM / Athena line searches",
            decisiveness=3,
            near_term=5,
            pressure=3,
            confidence=4,
            source="XRISM line constraints: https://inspirehep.net/literature/3074298",
            canon_hook="python_code/neutrino_sector_falsification_ledger.py",
        ),
        Entry(
            sector="inflation / amplitude",
            observable="scalar amplitude",
            prediction="A_s = (3/4) alpha_0^4 under the saturation stop rule",
            current_status=(
                "observed A_s is close, but external measurement is not the hard gate; "
                "the live risk is the internal stop-rule/current-volume theorem"
            ),
            falsifier="failure of N_shell alpha_0^4 = 4/3 or the scalar-current event-unit map",
            scope="kills HBC amplitude derivation; tilt may survive separately",
            horizon="internal derivation gate",
            decisiveness=4,
            near_term=3,
            pressure=3,
            confidence=4,
            source="boundary-printing/HBC scripts and papers; Planck A_s benchmark",
            canon_hook="ANCHOR item 131",
        ),
        Entry(
            sector="neutrino / active branch",
            observable="neutrinoless double beta decay",
            prediction=f"active Majorana allowed; normal-ordering m_bb = {mbb_min:.1f}-{mbb_max:.1f} meV",
            current_status=(
                "KamLAND-Zen full data gives T_1/2 > 3.8e26 yr and m_bb < 28-122 meV; "
                "still above the canonical meV band"
            ),
            falsifier=(
                "positive 0nuBB implying an incompatible mass/phase pattern, or null after experiments "
                "reach the meV normal-ordering band"
            ),
            scope="neutrino identity / mass-operator branch",
            horizon="LEGEND/nEXO/CUPID/KamLAND2-Zen generation",
            decisiveness=3,
            near_term=3,
            pressure=1,
            confidence=5,
            source="KamLAND-Zen full dataset: https://arxiv.org/abs/2406.11438",
            canon_hook="python_code/neutrino_sector_falsification_ledger.py; ANCHOR neutrino sector",
        ),
        Entry(
            sector="neutrino / direct mass",
            observable="beta-endpoint neutrino mass",
            prediction=f"m_beta = {m_beta_mev:.1f} meV; sum m_nu = {sum_m_mev:.1f} meV",
            current_status="KATRIN direct upper limit is m_beta < 0.45 eV, far above the prediction.",
            falsifier="direct mass or cosmological mass sum requiring a scale incompatible with the normal-ordering branch",
            scope="neutrino absolute-scale branch",
            horizon="KATRIN final / cosmology",
            decisiveness=2,
            near_term=3,
            pressure=1,
            confidence=5,
            source="KATRIN latest: https://www.science.org/doi/10.1126/science.adq9592",
            canon_hook="python_code/neutrino_sector_falsification_ledger.py",
        ),
        Entry(
            sector="black holes / horizon QEC",
            observable="horizon echoes and Hawking ladder",
            prediction=(
                "no large coherent echoes in the minimal V_cell channel; Hawking thermality should "
                "come from localized QEC steady state with finite strain-ladder ratios"
            ),
            current_status=(
                "no decisive observational echo signal; Hawking spectroscopy is inaccessible except "
                "for hypothetical primordial/micro black holes"
            ),
            falsifier="large coherent echoes with no derived reflection channel, or failure of localized-mass KMS steady state",
            scope="black-hole QEC channel, not the full matter/gauge sector",
            horizon="GW catalogues plus theory gate",
            decisiveness=3,
            near_term=3,
            pressure=1,
            confidence=3,
            source="papers/blackhole_deep/black_hole_qec_observational_guide.md",
            canon_hook="gravity_blackholes paper and BH scripts",
        ),
        Entry(
            sector="matter / spectroscopy",
            observable="Koide tau and hadronic structural checks",
            prediction="precision mass relations remain consistency anchors, not the main kill switches",
            current_status="current tau average is close to the Koide readout; spectroscopy side is comparatively stable.",
            falsifier="future precision mass/lattice results break the claimed locked identities beyond uncertainties",
            scope="matter/spectroscopy sector demotion, not cosmology",
            horizon="Belle II / lattice-QCD precision",
            decisiveness=2,
            near_term=2,
            pressure=1,
            confidence=4,
            source="ANCHOR matter/spectroscopy notes; PDG/Belle II status in canon",
            canon_hook="ANCHOR matter/spectroscopy section",
        ),
        Entry(
            sector="quantum information / code identity",
            observable="single-fermion entanglement ceiling",
            prediction="single-fermion reduced density matrix entropy saturates at 8 bits",
            current_status=(
                "multi-qubit entanglement experiments can probe the regime in principle; "
                "no dedicated substrate-code analogue test is recorded"
            ),
            falsifier=(
                "properly encoded single-fermion analogue shows unbounded growth beyond the "
                "8-bit code ceiling rather than saturation"
            ),
            scope="damages finite-code identity picture, not directly cosmology",
            horizon="analogue quantum-information experiment",
            decisiveness=3,
            near_term=2,
            pressure=1,
            confidence=2,
            source="ANCHOR section 2 falsifiable entanglement-saturation note",
            canon_hook="ANCHOR around the 8-bit single-fermion entropy prediction",
        ),
    ]


def build_markdown(entries: list[Entry]) -> str:
    top_five = entries[:5]
    source_count = len({entry.source for entry in entries})
    top_lines = "\n".join(
        f"{idx}. **{entry.observable}** ({entry.sector}) - {entry.band}; "
        f"score {entry.score}.  Falsifier: {entry.falsifier}"
        for idx, entry in enumerate(top_five, start=1)
    )
    null_table = compact_table(
        entries,
        lambda e: "drift" in e.observable
        or "tensor" in e.observable
        or "electric dipole" in e.observable,
    )
    return "\n".join(
        [
            "# Unified Prediction/Falsification Ledger",
            "",
            f"Generated by `python_code/unified_prediction_falsification_ledger.py` on {DATE}.",
            "",
            "This sheet consolidates the framework's forward-facing falsifiers into one",
            "ranked triage ledger.  It is deliberately attack-surface first: agreements and",
            "retrodictions are secondary to observations that would force a retraction,",
            "demotion, or major rewrite.",
            "",
            "Scope rule: this is not every claim in ANCHOR.  It is the unified list of",
            "current, externally or internally falsifiable prediction classes.  Precision",
            "retrodictions and already-retired numerology are left out unless they remain",
            "live tests.",
            "",
            "## Scoring Rule",
            "",
            "Score = `3*decisiveness + 2*near_term + 2*current_pressure + status_confidence`.",
            "",
            "- `decisiveness`: how much canon would be damaged by a failure.",
            "- `near_term`: how likely the test is to sharpen with existing or planned work.",
            "- `current_pressure`: how much present data or internal consistency already",
            "  stresses the claim.",
            "- `status_confidence`: how cleanly the current external status is known.",
            "",
            "Bands:",
            "",
            "- `kill switch`: score >= 34",
            "- `high priority`: 30-33",
            "- `watch closely`: 25-29",
            "- `longer range`: 20-24",
            "- `background`: below 20",
            "",
            "## Executive Triage",
            "",
            "The top five are the places to spend attention first:",
            "",
            top_lines,
            "",
            "## Ranked Ledger",
            "",
            entry_table(entries),
            "## Near-Term Kill Switches",
            "",
            compact_table(entries, lambda e: e.band == "kill switch"),
            "## At-Risk But Not Yet Killed",
            "",
            "These rows carry real pressure now.  They should be treated as active",
            "research targets rather than public victories.",
            "",
            compact_table(entries, lambda e: e.pressure >= 4 and e.band != "kill switch"),
            "## Correlated Null Signature",
            "",
            "The framework's most recognizable nulls are not independent absences.",
            "They are tied to frozen ledgers, absence of a high-scale tensor oscillator,",
            "and the strong-CP phase-null.",
            "",
            null_table,
            "## Neutrino Cluster",
            "",
            "The neutrino sector is compact and experimentally legible: active Majorana",
            "normal ordering at the meV 0nuBB scale, a separate 17.7 keV sterile branch,",
            "and a direct-mass scale below current KATRIN reach.",
            "",
            compact_table(entries, lambda e: e.sector.startswith("neutrino")),
            "## Source And Canon Hooks",
            "",
            f"External status was spot-checked against {source_count} source/status hooks.",
            "Refresh this section before public use if new experimental releases appear.",
            "",
            sources_table(entries),
            "## Public One-Paragraph Version",
            "",
            "The framework can be attacked cleanly.  Its most decisive live tests are",
            "a specific non-phantom dark-energy line `w(a)=-1+a/28`, a CMB completion",
            "that must not contradict the resolved R4 equation of state, dark matter",
            "morphology that separates pinned K04 debris from mobile R4-zero/nu_R halos,",
            "the proton-primary `G/H0` route, an exact linear primordial-tensor",
            "null `r_linear = 0` with induced floor `r_induced ~ 2e-9`, and the",
            "global trans-`Lambda_QCD` Lorentz bridge.",
            "The strongest correlated nulls are no secular drift of `G` or code-level",
            "`alpha`, a neutron EDM near `1e-31 e cm`, and no observable primordial",
            "B-modes.  The compact neutrino tests are the 17.7 keV sterile X-ray",
            "branch, meV-scale neutrinoless double beta decay, and sub-0.1 eV",
            "normal-ordering mass scale.",
            "",
        ]
    )


def main() -> None:
    entries = sorted(build_entries(), key=lambda e: (-e.score, e.observable))

    assert len(entries) >= 15
    assert entries[0].score >= entries[-1].score
    assert all(1 <= entry.decisiveness <= 5 for entry in entries)
    assert all(1 <= entry.near_term <= 5 for entry in entries)
    assert all(1 <= entry.pressure <= 5 for entry in entries)
    assert all(1 <= entry.confidence <= 5 for entry in entries)
    assert all(entry.prediction and entry.current_status and entry.falsifier for entry in entries)
    assert any("m_nuR" in entry.prediction and "keV" in entry.prediction for entry in entries)
    assert any("r_linear = 0" in entry.prediction for entry in entries)
    assert any("w(a) = -1 + a/28" in entry.prediction for entry in entries)

    OUT.write_text(build_markdown(entries), encoding="utf-8")

    print("UNIFIED PREDICTION/FALSIFICATION LEDGER")
    print("=" * 80)
    print(f"rows: {len(entries)}")
    print(f"output: {OUT.relative_to(ROOT)}")
    print("\nTop triage rows:")
    for idx, entry in enumerate(entries[:8], start=1):
        print(f"{idx:2d}. {entry.score:2d} {entry.band:<13} {entry.observable}")
    print("\nexit 0 -- unified ledger generated.")


if __name__ == "__main__":
    main()
