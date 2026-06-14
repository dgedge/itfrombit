# Derivation Metro Map Generation Prompt

Use this prompt to regenerate or revise the canon overview metro map.

```text
Create a clean, publication-quality "derivation metro map" for the It-for-Bit canon.

Purpose:
- Show the major derivation spines of the framework, not every citation.
- Findings are stations.
- Proofs, scripts, reductions, and assumptions are edge labels.
- Open blockers are visible as terminal diamond/dashed stations.
- Avoid a hairball: use a curated transitive reduction with no more than 45 nodes.

Visual style:
- Wide landscape figure suitable for a graduate-level overview paper.
- Use metro-map grammar: colored routes, transfer stations, small script labels.
- Background white, restrained palette, readable at A4 width.
- Left-to-right flow from foundations to physical sectors.
- Use consistent status encoding:
  - blue filled station: derived / script-backed
  - amber station: conditional closure
  - purple dashed station: open frontier
  - grey station: imported standard theorem
  - red station: retired/refuted, only if needed
- Use small monospaced labels for scripts and proof artifacts.
- Do not draw every ANCHOR citation. Draw load-bearing derivation edges only.

Routes:
1. Substrate/code spine:
   pre-geometric record soup -> F_2^3 cube/XOR geometry -> [8,4,4] code/AGL(3,2)
   -> bi-cubic/TCH cell geometry -> QEC service scheduler -> E_{1/2} service frame/tetrad.

2. Matter/spectroscopy:
   [8,4,4] code -> SM register/anomaly identities -> hadron/lepton spectroscopy
   -> Lambda_QCD/proton anchor -> service scheduler.

3. Gauge/photon/causal-set QED:
   TCH cell geometry -> Wilson/Gauss photon identity -> causal-set QED Ward/recoil
   -> framed loop Maxwell F^2 -> E_{1/2} frame.

4. Cosmology/dark energy:
   QEC scheduler -> cosmological-constant queue route -> G,H0 outputs/proton-first route
   -> HBC inflation tilt/amplitude -> open a=1 completion selector.

5. Dark matter/MOND/defects:
   QEC scheduler -> K04 crystallisation/Kibble-Zurek trapping -> dark defects/shadow energy
   -> open orphan policy/island floor -> open R4 line MOND proof.

6. Gravity/horizons/black holes:
   E_{1/2} frame -> metric strain/clock covariance -> horizon delta ledger
   -> Hawking ladder/KMS shell channel -> open freeze shell r_F(F;M).

Required cross-sector arrows:
- proton anchor -> G,H0 outputs, label "input/output ledger".
- horizon delta ledger -> CC route, label "55/8 severing".
- service frame -> loop Maxwell F^2, label "Spin(3,1) Hodge sign".
- K04 crystallisation -> CC route, label "boot cooling".
- dark defects -> horizon delta ledger, label "wall shadow".
- HBC inflation -> scheduler, label "saturation premise".

Important script/proof labels:
- record_alphabet_derivation.py
- spin_coin_2O.py
- foundations_causet_loop_f2_lorentzian_no_go.py
- foundations_causet_hodge_sign_from_frame.py
- foundations_causet_mesoscopic_continuum_action.py
- relativity_TR6_horizon_ledger_identity.py
- item131_saturation_closure.py
- k04_embedded_sweep.py / k04_orphan_policy_audit.py

Caption idea:
"Curated derivation metro map of the canon. Stations are major findings; route colors are sectors; dashed stations are live blockers. Script labels identify reproducibility artifacts. The map is a transitive reduction, not a complete citation graph."
```

## Render Commands

From the repo root:

```bash
cd papers/common/figures
lualatex -interaction=nonstopmode derivation_metro_map.tex
pdftocairo -svg derivation_metro_map.pdf derivation_metro_map.svg
```

Use `derivation_metro_map.pdf` in LaTeX papers and `derivation_metro_map.svg` for the repository or web.
