# Canon Overview Material

This folder mirrors selected explanatory material from the working canon tree:

`/Users/davidelliman/Library/Mobile Documents/com~apple~CloudDocs/octahedrons`

It is included here so the public source-code repo has enough orientation
material for readers to understand the code organization and the derivation
graph.

## Included

| Path | Purpose |
|---|---|
| `publication_plan_2026_06.md` | Seven-paper publication plan and topic split. |
| `model_for_teaching.md` | Teaching-oriented conceptual model. |
| `figures/derivation_tree_atlas.*` | Preferred seven-page tree-style derivation atlas. |
| `figures/tree_atlas_svg/` | Per-page SVG exports for the tree atlas. |
| `figures/derivation_metro_map.*` | Earlier one-page metro-map overview. |
| `figures/derivation_metro_atlas.*` | Earlier seven-page metro atlas, retained as an alternative visual grammar. |
| `figures/metro_atlas_svg/` | Per-page SVG exports for the metro atlas. |

## Not Included

`ANCHOR.md` and `DRIFT.md` are not copied here. They remain canonical in the
working canon tree, and copying them into this repo would create a drift risk.

## Regenerating Diagrams

From this directory:

```bash
cd figures
lualatex -interaction=nonstopmode derivation_tree_atlas.tex
mkdir -p tree_atlas_svg
for i in 1 2 3 4 5 6 7; do
  pdftocairo -f "$i" -l "$i" -svg derivation_tree_atlas.pdf "tree_atlas_svg/derivation_tree_atlas-${i}.svg"
done
```

The tree atlas is the preferred version for paper use because it is more
readable than the denser metro-map version.
