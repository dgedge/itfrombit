#!/usr/bin/env python3
r"""Compatibility entry point for the Item 118/123 CMB completion audit.

The original version of this script recorded the live 2026-06-15 issue:
R4 appeared in canon as both w=+1/3 and w=-1/3, and neither value could
anchor the CMB third peak.  The current audit in r4_eos_cmb_resolution.py
settles the status split:

* w=+1/3 is the active microscopic/unbound Landauer-exhaust reading;
* w_eff>0 is the active bound halo line-current/topological-fluid reading;
* w=-1/3 was the retired string-gas/negative-pressure branch.

The CMB no-go remains: only the 20% sterile-nu_R component supplies early
w=0 clustering matter, giving a 5x third-peak-anchor shortfall unless a new
R4-as-early-dust or AeST-like theorem is derived.
"""

from r4_eos_cmb_resolution import main


if __name__ == "__main__":
    main()
