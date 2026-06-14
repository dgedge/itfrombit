#!/usr/bin/env python3
"""
E1 "staged follow-up" (§16.3 search-space accounting): competitor counts for the
per-section simple-ratio matches whose ANCHOR headers still read "derivation /
parameter-free / zero-parameter prediction":
   - string tension  sqrt(sigma)/Lambda = 4/3   (§7.17)
   - glueball        m(2++)/m(0++)      = 7/5   (glueball section)
   - dark sector     Omega_DE/Omega_DM  = 12/5  (§16 dark-sector)
   - rho mass        m_rho/Lambda       = sqrt2*phi  (§9.1)

§16.3 (ANCHOR) already classifies 4/3, 7/5, 12/5, sqrt2*phi as class-1
"consistency checks, not predictions." This supplies the VERIFIED competitor
count at each ratio's ACTUAL match tolerance, reusing the canonical alphabet in
numerology_baseline.py. A ratio with >=~2 equally-simple competitors inside its
own match tolerance carries near-zero standalone evidential weight (class 1).

Reports counts; asserts the script reproduces the canonical density (the
clearest class-1 case, the 7%-tolerance dark-sector ratio, has many competitors)
and that numbers are computed. Nothing fitted.
"""
import numerology_baseline as nb
L = 332.0  # Lambda_QCD (MeV), framework canonical chiral scale (§1.4)
fam = nb.build_family(12)   # canonical alphabet, integer ceiling 12

phi = (1 + 5**0.5)/2
cases = [
    # label, framework expr value, empirical target value
    ("string tension  4/3   (§7.17)",  4/3,            440.0/L),     # sqrt(sigma)_QCD ~440 MeV
    ("glueball 2++/0++ 7/5  (glueball)", 7/5,           1.397),       # LQCD 1.397
    ("dark sector  12/5     (§16)",     12/5,           2.58),        # observed Omega_DE/Omega_DM
    ("rho  sqrt2*phi        (§9.1)",    2**0.5 * phi,   775.26/L),    # m_rho/Lambda, exp 775.26 MeV
]

print(f"{'ratio':<34}{'expr':>8}{'target':>9}{'match%':>8}{'#comp':>7}   sample competitors")
results = []
for label, expr, target in cases:
    tol = abs(expr - target)/target            # the ACTUAL match tolerance
    comps = nb.competitors(fam, target, tol)   # equally/less-simple exprs within that band
    n = len(comps)
    results.append((label, expr, target, tol, n))
    sample = ", ".join(f"{v:.4f}" for v, t in comps[:6])
    print(f"{label:<34}{expr:>8.4f}{target:>9.4f}{tol*100:>7.2f}%{n:>7}   {sample}")

print()
print("Interpretation per §16.3 (class-1 = many equally-simple competitors in the")
print("match band => near-zero standalone weight; the evidential weight, if any,")
print("rests on whether the *derivation* of the integers is free of formula-choice):")
for label, expr, target, tol, n in results:
    cls = "class-1 (consistency check)" if n >= 2 else "tight (0-1 competitor; class-2, weight in derivation)"
    print(f"  {label:<34} {n:>3} competitors @ {tol*100:.2f}%  -> {cls}")

# robust assertions: the script computed counts, and the loosest-tolerance ratio
# (dark sector, ~7%) is firmly class-1 with many competitors.
dark = [r for r in results if "dark sector" in r[0]][0]
assert dark[4] >= 5, f"dark-sector should have many competitors at 7%, got {dark[4]}"
assert all(isinstance(r[4], int) for r in results), "counts not computed"
print("\nexit 0 == competitor counts computed; dark-sector (7% band) firmly class-1.")
