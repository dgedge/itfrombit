# Seven-Page Derivation Metro Atlas Prompt

Use this prompt to regenerate or revise the multi-page canon derivation atlas.

```text
Create a seven-page publication-quality derivation metro atlas for the
It-for-Bit canon. The atlas should match the overview paper plus the six
technical companion papers.

General design:
- Use a clean metro-map visual grammar.
- Each page is landscape, white background, readable at A4 or letter width.
- Findings are stations; proofs, scripts, reductions, or bridge assumptions are
  edge labels.
- Use colored routes for sectors and grey arrows for cross-sector dependencies.
- Keep each page to roughly 12--18 stations so it remains readable on paper.
- Use a consistent status grammar:
  - blue station: derived or script-backed
  - amber station: conditional closure
  - purple dashed station: open frontier
  - grey station: imported standard theorem
  - red dashed station: retired or superseded route
- Do not draw every ANCHOR citation. Draw the curated transitive reduction: the
  load-bearing dependencies that a reader needs to understand the proof flow.

Pages:
1. Overview:
   Show the shared spine:
   F_2^3 cube/XOR records -> [8,4,4] code -> bi-cubic/TCH cell
   -> QEC service scheduler -> E_{1/2} service frame -> causal-set null records.
   Branch to matter/spectroscopy, CC/HBC cosmology, gravity/horizon ledger,
   photon/Maxwell/QED, dark defects/MOND, and black-hole shell channel.

2. Foundations and Methodology:
   record soup -> XOR closure -> right bi-cubic volume -> [8,4,4] stabilizer
   -> AGL(3,2)/28-clock -> TCH gauge cell -> QEC service ledger -> T1--T9
   claim protocol and PTMS reproducibility.
   Show the fine-UV representation as a named frontier.

3. Matter, Gauge Structure, and Spectroscopy:
   [8,4,4] code -> SM register -> anomaly identities -> fermion families.
   In parallel: Lambda_QCD anchor -> hadron sweep -> Koide/tau conditional
   status. Then TCH cell -> Wilson/CG extended reps -> SMG certified domain
   -> weak-beta frontier. Include Gauss-projected photon identity as a transfer.

4. Cosmology, Dark Energy, and Inflation:
   QEC scheduler -> billed-ledger generation -> active-address demux
   -> generation-vertex loop dressing -> rho_Lambda queue slot.
   In parallel: proton-primary Lambda_QCD route -> G,H0 outputs -> historical
   horizon audit. Then saturated printer premise -> n_s=27/28 -> A_nu amplitude,
   with a=1 completion selector as open.

5. Dark Matter, MOND, and K04 Debris:
   K04 embedded lattice -> cube crystal -> Kibble-Zurek trapping -> aging and
   wall spectra -> boundary shadow readout -> orphan-rescue policy -> island
   floor surface. Keep the toy K_{3,3} phase visibly retired. Add separate R4
   line occupancy -> Kraus matched rates -> cored profiles/Jeans -> BTFR/MOND.

6. Gravity, Horizons, and Black Holes:
   E_{1/2} service frame -> metric as clock covariance -> alpha^1 erasure
   prefactor -> M_P/G/H0 proton-primary route.
   In parallel: finite-cell isometry V -> Schwarzschild shell channel
   -> horizon delta ledger -> localized QEC KMS steady state -> Hawking ladder,
   ending at wavepacket dynamics and freeze shell as open frontiers.

7. Relativity, Photon, and Causal-Set QED:
   frictionless clock ticks -> special relativity IR limit and general
   relativity strain covariance. Then Gauss-projected Maxwell cohomology
   -> QED Ward/recoil -> B_N null service packet -> one LSZ leg with total P^mu.
   In parallel: mesoscopic Alexandrov loops -> framed Wick/Hodge sign
   -> gauge-invariant Maxwell F^2 scaffold -> causal-set Dirac spin and
   finite-density variance control as open frontiers.

Caption:
"Seven-page derivation metro atlas of the canon. Stations are major findings;
route colors are sectors; dashed stations are open blockers; red dashed stations
are retired routes. The atlas is a curated transitive reduction of the PTMS/ANCHOR
claim graph, not a complete citation graph."
```

## Render Commands

From the repo root:

```bash
cd papers/common/figures
lualatex -interaction=nonstopmode derivation_metro_atlas.tex
mkdir -p atlas_svg
for i in 1 2 3 4 5 6 7; do
  pdftocairo -f "$i" -l "$i" -svg derivation_metro_atlas.pdf "atlas_svg/derivation_metro_atlas-${i}.svg"
done
```

`dvisvgm` can also export the pages if its PDF backend is available. On this
machine it reports that PDF input requires `mutool` or Ghostscript older than
10.01, so `pdftocairo` is the working converter.
