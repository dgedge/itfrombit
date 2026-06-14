# Derivation Tree Atlas Prompt

Use this prompt when regenerating or revising the tree-style canon atlas.

```text
Create a seven-page derivation tree atlas for the It-for-Bit canon.

The graph is not a pure tree, but it should be drawn as a tree-first structure:
approximately 95% of the visual flow should be parent -> child proof
dependencies. Only add dashed cross-arrows for genuinely load-bearing
non-tree dependencies.

Purpose:
- Make the proof structure readable on paper.
- Avoid the visual density of a metro map.
- Use large labels; do not let boxes overlap.
- Prefer 10--16 stations per page.
- Use captions and surrounding prose for nuance; keep the diagram itself sparse.

Visual grammar:
- Landscape pages.
- One root node on the left.
- Two or three horizontal proof branches per page.
- Heavy solid arrows: main derivation tree.
- Dashed grey arrows: cross-sector dependency.
- Status encoding:
  - blue: derived / script-backed
  - amber: conditional closure
  - purple dashed: open frontier
  - red dashed: retired or superseded
- Put the status legend on the overview page only.

Pages:
1. Overview:
   finite record substrate and QEC ledger as root; branches to foundations,
   matter, cosmology, dark sector, relativity/photon, gravity, and black holes.

2. Foundations and Methodology:
   finite record soup -> XOR/F_2^3 cube -> [8,4,4] code/AGL(3,2) -> TCH cell
   -> QEC service scheduler. Parallel branch: K04 crystallisation -> T1--T9
   protocol -> PTMS/script reproducibility -> finer-UV frontier.

3. Matter, Gauge Structure, and Spectroscopy:
   [8,4,4]/TCH root. Branches: SM register/anomaly identities/families;
   Lambda_QCD/hadron sweep/Koide status; Wilson-CG basis/SMG certified
   domain/weak-beta frontier/Gauss photon identity.

4. Cosmology, Dark Energy, and Inflation:
   QEC scheduler root. Branches: billed-ledger generation -> active demux
   -> generation vertex -> rho_Lambda; proton-primary Lambda_QCD -> G,H0 outputs;
   saturated printer -> n_s=27/28 -> A_nu threshold -> a=1 completion selector.

5. Dark Matter, MOND, and K04 Debris:
   embedded K04 substrate root. Branches: cube crystal -> Kibble-Zurek trapping
   -> durability -> shadow readout -> orphan policy -> island-floor surface;
   retired K33 toy phase as negative control; R4 support -> Kraus matched rates
   -> cored profile/Jeans -> BTFR/MOND.

6. Gravity, Horizons, and Black Holes:
   E_{1/2} frame root. Branches: metric covariance -> alpha^1 erasure -> proton
   primary M_P/G/H0 route -> historical horizon routes; finite-cell V -> V_Sch(M)
   -> horizon delta ledger -> KMS shell channel -> Hawking ladder -> wavepacket
   and freeze-shell frontiers.

7. Relativity, Photon, and Causal-Set QED:
   service clock/null record root. Branches: SR IR propagation -> Gauss Maxwell
   cohomology -> QED Ward/recoil -> B_N service packet -> one LSZ leg; GR strain
   covariance -> Alexandrov loops -> Hodge sign -> Maxwell F^2 scaffold -> Dirac
   spin and finite-density variance frontiers.
```

## Render Commands

From the repo root:

```bash
cd papers/common/figures
lualatex -interaction=nonstopmode derivation_tree_atlas.tex
mkdir -p tree_atlas_svg
for i in 1 2 3 4 5 6 7; do
  pdftocairo -f "$i" -l "$i" -svg derivation_tree_atlas.pdf "tree_atlas_svg/derivation_tree_atlas-${i}.svg"
done
```

Preview PNGs used during layout checks:

```bash
mkdir -p tree_atlas_preview
pdftoppm -png -r 150 -f 1 -l 7 derivation_tree_atlas.pdf tree_atlas_preview/derivation_tree_atlas
```
