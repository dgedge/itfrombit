# Canon Source Scripts

This directory contains a curated migration of source scripts from the working canon tree.
It excludes bulky numerical outputs, cached files, and run directories.

The source-of-truth manifest is [`../docs/script_manifest.csv`](../docs/script_manifest.csv).

| Topic | Scripts | Needs docstring pass |
|---|---:|---:|
| [Foundations and Methodology](foundations_methodology/) | 57 | 5 |
| [Matter, Gauge Structure, and Spectroscopy](matter_gauge_spectroscopy/) | 56 | 0 |
| [SMG, TCH Gauge Cell, and Lattice Certification](smg_tch_gauge_cell/) | 106 | 14 |
| [Cosmology, Dark Energy, and Inflation](cosmology_dark_energy_inflation/) | 74 | 0 |
| [Dark Matter, MOND, and K04 Debris](dark_matter_mond_k04/) | 50 | 7 |
| [Gravity, Horizons, and Black Holes](gravity_horizons_blackholes/) | 28 | 0 |
| [Relativity, Photon, and Causal-Set QED](relativity_photon_qed/) | 58 | 7 |
| [Legacy and Miscellaneous Source](legacy_misc/) | 43 | 5 |

Run commands in the topic READMEs assume the repository root as the working directory.
Many scripts are audit/reproduction scripts rather than reusable libraries; read the topic README and script docstring before interpreting output.

Large run outputs should not be committed here. Put regenerated artifacts in `results/` only when they are small, named, and documented.
