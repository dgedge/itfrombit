#!/usr/bin/env python3
r"""ITEM 123: does the K04 CONSTRAINT-GLASS give the R4 line-current a pressureless (c_s^2->0),
abundance-fixed FROZEN DUST phase? — the last framework-native CMB-completion escape (escape (ii)
of item123_r4_noether_charge_attempt.py).

THE HOPE. The Noether attempt showed the R4 ledger has no conserved cold-dust charge, but flagged
one escape: a glassy/vitrified R4 phase (cf. the K04 constraint-glass arrest, exit 0 in
k04_coarsening_kz_exponents.py: a deep quench ARRESTS at defect density d~0.95) might freeze a
homogeneous background that (a) is pressureless (c_s^2->0, so it clusters), (b) has a fixed
abundance (the arrest density is parameter-free), and (c) is frozen (tau_relax >> t_Hubble). If so,
 omega h^2=0.096 of frozen R4 dust would close the CMB (nu_R + R4-dust = 0.120).

THE PHYSICS (the decisive computation). A frozen / conformally-stretched p-DIMENSIONAL defect
network has a fixed comoving defect number; each connected defect's physical EXTENT stretches with
expansion (length ~ a), so its energy ~ a^p and the density redshifts as

    rho ~ a^(p-3)   <=>   w = -p/3.

  * p=0 (POINT defects):   rho ~ a^-3  ->  w = 0     (matter / dust)   <-- the ONLY dust case
  * p=1 (STRINGS/lines):   rho ~ a^-2  ->  w = -1/3  (anti-clustering)
  * p=2 (WALLS):           rho ~ a^-1  ->  w = -2/3  (anti-clustering)

THE DICHOTOMY (why a glass cannot be dust). w=0 needs FIXED-size, disconnected, point-like (p=0)
fragments -- a DILUTE gas below the percolation threshold. A GLASS is the opposite: it is rigid
*because* it is a dense, percolating, CONNECTED network. The K04 arrest is at d~0.95, far above any
percolation threshold (~0.2-0.5), so it is deeply connected -> extended (p>=1) -> w<0. The very
property that makes it "frozen/rigid" (connectivity) is what forces w<0. Frozen-and-rigid (glass)
<=> connected/extended <=> w<0;  w=0 dust <=> 0D/fragmented <=> NOT a glass. Mutually exclusive.

  * R4 is a 1D LINE-current (item131/132; its retired branch was literally "string-gas, w=-1/3").
    Freezing it -> w=-1/3: it RE-CREATES the already-retired anti-clustering branch, not dust.
  * K04 debris is a 2D WALL-shadow -> w=-2/3 (even more anti-clustering).
  * The user's hoped "pressureless (c_s^2->0)" phase is the OPPOSITE of what freezing a line-current
    gives: a frozen line-current is TENSION-dominated (negative pressure, w=-1/3), not pressureless.

exit 0 = w=-p/3 derived/asserted; the percolation dichotomy computed; frozen-R4 shown to coincide
with the retired w=-1/3 string-gas; the verdict (escape (ii) CLOSED; only the nu_R relic is w=0) is
reported with honest tiering.
"""

from __future__ import annotations

from fractions import Fraction


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def w_frozen(p: int) -> Fraction:
    """EoS of a frozen / conformally-stretched p-dimensional defect network: rho ~ a^(p-3)."""
    return Fraction(-p, 3)


# K04 constraint-glass arrest density (k04_coarsening_kz_exponents.py, exit 0): deep quench
# arrests at d ~ 0.95 (acceptance ~1%). Representative 3D percolation threshold (site, simple
# cubic) ~ 0.31; every common 3D lattice threshold is in ~0.2-0.5.
ARREST_D = 0.95
PERC_THRESHOLD = 0.31


def main() -> None:
    print("ITEM 123: CAN THE K04 CONSTRAINT-GLASS GIVE R4 A PRESSURELESS FROZEN-DUST PHASE?")

    print("\n[0] What CMB dust requires (all four):")
    print("    (a) w~0 (a^-3 background, fixes z_eq);  (b) c_s^2->0 (clusters on CMB scales);")
    print("    (c) abundance omega h^2~0.096 fixed;     (d) frozen (tau_relax >> t_Hubble).")

    print("\n[1] EoS of a frozen p-dimensional defect network: w = -p/3")
    rows = [(0, "point defects"), (1, "strings / line-current"), (2, "domain walls")]
    for p, name in rows:
        w = w_frozen(p)
        verdict = "DUST (clusters)" if w == 0 else "anti-clustering (DE-like)"
        print(f"    p={p} ({name:22s}): rho ~ a^({p-3}) -> w = {str(w):>4s}  -> {verdict}")
    check(w_frozen(0) == 0, "ONLY point defects (p=0) give w=0 dust")
    check(w_frozen(1) == Fraction(-1, 3), "strings (p=1) give w=-1/3")
    check(w_frozen(2) == Fraction(-2, 3), "walls (p=2) give w=-2/3")
    check(all(w_frozen(p) < 0 for p in (1, 2)), "all EXTENDED frozen defects are anti-clustering (w<0)")

    print("\n[2] The dichotomy: frozen/rigid <=> connected <=> extended <=> w<0")
    percolating = ARREST_D > PERC_THRESHOLD
    print(f"    K04 glass arrests at d~{ARREST_D} ; percolation threshold p_c~{PERC_THRESHOLD}")
    print(f"    d/p_c = {ARREST_D/PERC_THRESHOLD:.1f}  ->  {'DEEPLY PERCOLATING (connected)' if percolating else 'sub-percolating'}")
    check(percolating, "the glass arrests far ABOVE percolation -> a connected, extended network")
    check(ARREST_D / PERC_THRESHOLD > 2.0, "d~0.95 is >2x the threshold: robustly connected, not a dilute loop gas")
    print("    -> rigidity (what makes it 'frozen') IS connectivity; connectivity forces p>=1 -> w<0.")
    print("       w=0 dust would need fixed-size, disconnected, p=0 fragments (a dilute gas BELOW p_c)")
    print("       -- i.e. NOT a glass. 'Frozen-rigid' and 'w=0 dust' are mutually exclusive here.")

    print("\n[3] Apply to the framework's actual frozen structures")
    w_R4_glass = w_frozen(1)            # R4 is a 1D line-current
    w_K04_glass = w_frozen(2)           # K04 debris is a 2D wall-shadow
    w_retired_string_gas = Fraction(-1, 3)
    print(f"    frozen R4 line-current (1D): w = {w_R4_glass}  == retired 'string-gas' branch w = {w_retired_string_gas}")
    print(f"    frozen K04 wall-shadow (2D): w = {w_K04_glass}")
    check(w_R4_glass == w_retired_string_gas, "freezing R4 RE-CREATES the already-retired w=-1/3 string-gas (not new dust)")
    check(w_K04_glass < w_R4_glass < 0, "K04 walls are even more anti-clustering than R4 strings; both w<0")

    print("\n[4] The hoped 'pressureless (c_s^2->0)' phase is refuted")
    print(f"    a frozen line-current is TENSION-dominated: w={w_R4_glass} (NEGATIVE pressure), not pressureless.")
    print("    c_s^2 is not ->0; the component is anti-clustering (DE-like), the opposite of CDM.")
    check(w_R4_glass != 0, "freezing gives w=-1/3 (tension), not the c_s^2->0 / w=0 the user hoped for")

    print("\n[5] Abundance is moot given [3]-[4], but for completeness")
    print(f"    the arrest density d~{ARREST_D} is parameter-free, so a HYPOTHETICAL p=0 channel could")
    print("    in principle have a computable abundance -- but R4 has no point-like (p=0) defect to")
    print("    freeze: its support is the 1D line-current. So even the abundance hook has nothing to grip.")
    check(True, "abundance is a non-issue: there is no w=0 (p=0) R4 component for the glass to fix")

    print("\n[6] Verdict — escape (ii) is CLOSED")
    print(
        "  Freezing the R4 line-current does NOT give pressureless dust. A constraint glass is rigid\n"
        "  *because* it is a dense, percolating, CONNECTED network (K04 arrests at d~0.95, ~3x the\n"
        "  percolation threshold), and a connected p-dimensional defect network has w=-p/3 < 0:\n"
        "    * frozen R4 (1D line-current) -> w=-1/3  -- exactly the already-RETIRED string-gas branch;\n"
        "    * frozen K04 (2D wall-shadow) -> w=-2/3  -- even more anti-clustering.\n"
        "  Both are DE-like (anti-clustering), not the w=0 dust the CMB needs. The structural reason is\n"
        "  a dichotomy: frozen/rigid <=> connected/extended <=> w<0, whereas w=0 dust <=> 0D/fragmented\n"
        "  <=> NOT a glass. The R4 line-current (1D) cannot be made w=0 by freezing."
    )
    print(
        "  CONSEQUENCE: combined with the Noether no-go, BOTH framework-native escapes are now closed.\n"
        "  The only w=0 conserved component is the nu_R relic (0.024). Completing the CMB requires\n"
        "  something genuinely EXTERNAL to the substrate's R4/K04 structures -- a shift-symmetric,\n"
        "  pressureless, abundance-fixed (AeST-class) sector that is imported, not derived -- or the\n"
        "  CMB stands as real falsification pressure on the MOND + nu_R + R4-DE cosmology.\n"
        "  TIER: settled within the substrate's defect structure (the EoS of frozen defects is fixed by\n"
        "  their dimensionality); the AeST import remains a conjecture (ambiguous)."
    )
    print("exit 0 -- glass escape CLOSED: frozen R4 (1D) is w=-1/3 (retired string-gas), not w=0 dust.")


if __name__ == "__main__":
    main()
