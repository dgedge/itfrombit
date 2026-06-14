# Publication Plan: Canon Snapshot Series, June 2026

This plan records the publication strategy for the current canon state. It is meant to guide drafting, not to serve as a finished abstract.

## Immediate Goal

Prepare one readable overview paper for serious external review, followed by five focused technical companion papers. The overview should be suitable for a numerate graduate-level reader and for an expert reader who wants a compact map of the framework without hunting through the whole canon.

The immediate motivating reader is John Baez, who has offered to read a summary paper. The first document therefore needs to be concise, prose-led, honest about status, and explicit about which results are derived, computed, conditional, retired, or still open.

## Author And Public Pointers

Author block for all papers:

```tex
\author{David Elliman\\
Neuro-symbolic Ltd\\
\texttt{dave@neusym.ai}}
```

Public information pointer:

- https://neusym.ai

Code and reproducibility repository:

- https://github.com/dgedge/itforbit_model.git

Each paper should cite the exact canon/code commit used for the draft. New derivations should be claimed by canon item, script, expected output, and commit hash.

## LaTeX And Bibliography Standard

Use a simple universal LaTeX format unless a target venue later requires something else.

Recommended default:

```tex
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{hyperref}
\usepackage[numbers,sort&compress]{natbib}
\bibliographystyle{plainnat}
```

Compile with LuaLaTeX:

```bash
lualatex paper.tex
bibtex paper
lualatex paper.tex
lualatex paper.tex
```

Shared bibliography:

- `papers/common/references.bib`

Optional shared macros:

- `papers/common/macros.tex`

Do not use a journal-specific class such as RevTeX until a target venue is chosen. The first priority is portability and readability.

## Common Status Vocabulary

Every headline result should carry one of these labels:

- **Locked**: finite derivation or exact theorem with no named load-bearing gap.
- **Computed**: reproducible numerical or finite-algebra result, with script and expected output.
- **Conditional**: derivation reduces to named premises or lemmas that remain open.
- **Horizon / Dirac bridge**: relation uses the QCD-to-cosmological large-number bridge. The paper must state whether the horizon quantity is an empirical input, a scaffolding variable that cancels, or an output of the route. The current proton-primary gravity route is in the last category: it predicts `G` and `H_0` from the proton anchor plus dimensionless microphysics.
- **Retired / retracted**: retained for provenance, no longer a live canonical result.
- **Open frontier**: named unsolved target with a concrete closure/falsification condition.

This vocabulary should appear in the overview and be reused in each companion paper.

## Series Structure

### Paper 0: Overview

Working title:

> A Finite-QEC Substrate Program for Particle Physics and Cosmology: Current Canon and Audit Methodology

Purpose:

- Give the shortest faithful map of the framework.
- Explain the substrate, the code structure, the crystallising geometry, and the methodology.
- Present the strongest current results and the remaining open fronts.
- Make clear that the framework has an active retraction/audit ledger, not a promotional one-way ratchet.

Target length:

- 15-25 pages main text.
- Appendices only for headline equations, status tables, and reproduction pointers.

Required sections:

1. Motivation and scope.
2. The finite substrate picture.
3. The right bi-cubic crystallisation picture.
4. The `[8,4,4]` code and the `28 = 2 x 14` service clock.
5. Matter and gauge structure: what is locked, what is conditional.
6. Cosmology and dark sector: dark energy, HBC/inflation, debris dark matter, MOND.
7. Gravity and horizons: proton-primary prediction of `G`/`H_0`, and the historical Dirac-class inversions.
8. Black holes and horizon QEC.
9. Methodology: DRIFT, PTMS, claim tiers, and reproducibility.
10. Open problems and falsification targets.

### Required Overview Explanation: Bi-Cubic Crystallisation

The overview must give a clear, non-technical but precise explanation of how the right bi-cubic columns/volumes crystallise into being.

The prose goal is that a graduate reader can picture the mechanism before seeing the appendices:

- The substrate is not initially assumed to be a smooth continuum.
- Local register-bearing cells form an embedded cubic/crystalline phase through degree-constrained bonding.
- Vertices carry the local register/event content.
- Branches/bonds carry adjacency, line-current support, and the finite move structure.
- Faces/plaquettes carry the loop and cell constraints that make cube tilings, gauge-cell structure, and QEC boundary checks meaningful.
- Crystallisation selects the ordered visible substrate; mis-ordered embedded domain walls become the debris sector.
- The physical claims should always distinguish:
  - ordered crystal,
  - wall/debris defects,
  - instrumented cells,
  - uninstrumented or shadow-read wall regions.

The overview should explain this qualitatively. The foundations paper should supply the full lattice definitions, cycle counts, move set, sector/parity constraints, Kibble-Zurek ramp protocol, and island/orphan-policy bookkeeping.

### Paper 1: Foundations, Canon, And Methodology

Focus:

- Finite substrate and `[8,4,4]` code.
- Decoder theorem: strain ledger vs syndrome ledger.
- `28 = 2 x 14` service-clock structure.
- Claim-tier protocol, DRIFT, PTMS, reproducibility gates.
- The thermodynamic-claim protocol and why many older entropy claims were demoted.
- Detailed crystallisation mechanics: vertices, faces, branches, cube tilings, embedded move sets, sectors, and defect taxonomy.

This paper establishes the rules of evidence.

### Paper 2: Matter, Gauge Structure, And Spectroscopy

Focus:

- Standard Model matter content from the code structure.
- Exact anomaly arithmetic, including recovered chiral `-2/9` identity.
- Hadronic and spectroscopic results that survived audit.
- Koide tau status.
- Strong CP: bare-substrate closure and continuum caveats.
- SMG/mirror-gap status, including the extended-rep and neutral-hopping computations.

This paper should distinguish exact finite-code results from continuum-lift problems.

### Paper 3: Cosmology, Dark Energy, And Inflation

Focus:

- Dynamic dark energy `w(a) = -1 + a/28`.
- Activation-identification result and non-phantom branch.
- HBC scalar-clock route to `n_s = 27/28`.
- Inflation amplitude status: candidate event-unit route, still conditional unless the service-current correlation volume closes.
- Cosmological constant route: active-address demux, post-service readout, generation-vertex loop, and remaining next-order residual.

This paper must be especially strict about conditional labels because the dark-sector/thermodynamic cluster had the largest historical over-claim risk.

### Paper 4: Dark Matter, MOND, And Debris

Focus:

- Embedded K04 crystallisation.
- Toy configuration-model failure and `K_{3,3}` artifact.
- Kibble-Zurek trapping, durability, percolation, aging, and domain-wall spectra.
- Boundary-local rescue / orphan policy and island-floor surface.
- R4 exhaust and MOND/BTFR: Poisson line-current theorem, matched-rate/nonexclusive-line premise, remaining Jeans/cored-profile targets.

This paper should keep the debris mechanism and R4/MOND mechanism distinct until a quantitative reconciliation is derived.

### Paper 5: Gravity, Horizons, And Black Holes

Focus:

- Proton-primary gravity route: `m_p` fixes `Lambda_p = m_p/(2 sqrt(2))`; the route predicts `G`, `M_P`, and `H_0` rather than consuming a Hubble measurement.
- Input/output ledger for the gravity route, including the distinction between `H_0` appearing in derivational scaffolding and cancelling from the final input list (`python_code/g_route_input_ledger.py`).
- Historical `M_P`, `a0`, and `rho_Lambda` horizon/Dirac consolidation: retain as the audit trail and inverse checks, not as the preferred presentation.
- What is derived, what remains conditional, and which older formulations consumed the horizon as an empirical input.
- Bekenstein severing-channel count status.
- `V_cell` and Schwarzschild/radial shell-channel construction.
- Hawking ladder wavepacket verification as an open dynamical target.

This paper should be conservative but updated: the honest headline is no longer "gravity consumes `H_0`." It is "proton-first predicts `G` and `H_0` inside a still-conditional QEC/horizon-accounting chain." The conditions and historical inversions must be shown explicitly.

## Reproducibility Appendix Standard

Each paper should include an appendix table with this schema:

| Claim | Canon item | Status | Script | Invocation | Expected output | Commit |
|---|---:|---|---|---|---|---|

Each referenced script should have:

- one-line purpose,
- input assumptions,
- invocation,
- expected pass/fail or numerical output,
- dependency/environment note,
- whether it is a proof-enumeration, numerical scan, toy model, or physical-ensemble run.

Preferred Python invocation:

```bash
~/bin/py13_7/bin/python path/to/script.py
```

For PTMS checks:

```bash
PYTHONPATH=ai_methodology ~/bin/py13_7/bin/python -m ptms extract --write
PYTHONPATH=ai_methodology ~/bin/py13_7/bin/python -m ptms check
```

## Drafting Rules

- Prose first, mathematics second.
- Put headline mathematics in the main text.
- Put detailed derivations in appendices.
- Never hide a premise in prose if it is load-bearing.
- Mark every conditional result with the exact missing lemma.
- Do not present horizon-consuming Dirac-class relations as intrinsic predictions. Conversely, when a later route has an input/output ledger showing the horizon variable cancels, present it as such and cite the ledger.
- Keep historical/retracted claims only when they explain why the current statement is sharper.
- Every new result should be accompanied by a canon item, script, and commit hash.

## Priority Strategy

1. Freeze a canon snapshot commit.
2. Draft Paper 0 first.
3. Build the shared bibliography while drafting Paper 0.
4. Convert existing script/canon references into the reproducibility table.
5. Draft companion papers in the order:
   - Foundations/methodology,
   - Cosmology/dark energy/inflation,
   - Dark matter/MOND/debris,
   - Matter/gauge/spectroscopy,
   - Gravity/black holes.
6. After the overview is externally reviewed, update existing older papers only by explicit erratum or supersession notes, not silent rewrites.

## Immediate Next Files To Create

- `papers/common/references.bib`
- `papers/common/macros.tex`
- `papers/overview/overview.tex`
- `papers/overview/script_index.md`
- `papers/foundations/foundations.tex`
- `papers/cosmology/cosmology.tex`
- `papers/dark_sector/dark_sector.tex`
- `papers/matter_gauge/matter_gauge.tex`
- `papers/gravity_blackholes/gravity_blackholes.tex`

The current file is the planning checkpoint before creating those drafts.
