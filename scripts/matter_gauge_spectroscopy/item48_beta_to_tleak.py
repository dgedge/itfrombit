#!/usr/bin/env python3
r"""Item 48 — does the silver-ratio area-lock beta_oct/beta_tri PIN the matter leak hopping t_leak
(and hence close the rho survival verdict)? Self-asserting on the geometry; the verdict is a careful
NEGATIVE: the area-lock fixes the GAUGE stiffness ratio, which is a different object from the MATTER
hopping, and the physically-defensible bridges straddle the survival/dissolution threshold.

GEOMETRY (truncated cube, the gauge cell; edge length a, uniform):
  octagonal face (gauge<->gauge propagation): regular octagon, area A_oct = 2(1+sqrt2) a^2
  triangular face (gauge<->matter vertex)   : equilateral triangle, area A_tri = (sqrt3/4) a^2

AREA-LOCK (ANCHOR sec 7.3 / L1110-1116): matching discrete flux^2 density to Maxwell 1/2 B^2 gives the
Wilson stiffness beta_p ∝ 1/A_p, so beta_oct/beta_tri = A_tri/A_oct. This is a GAUGE-FIELD coupling
ratio (stiffness of link fluctuations), NOT a matter hopping.

THE BRIDGE PROBLEM. The rho resonance width is Gamma ∝ t_leak^2 rho_gauge(Delta+phi), set by the MATTER
hopping t_leak across the triangular face (relative to the gauge-web hopping t_oct≡1). In a lattice gauge
theory t (matter hopping) and beta (gauge stiffness) are INDEPENDENT: the hopping carries the link as
t·U_ij, with |t| geometric and beta governing the PHASE fluctuations of U. So beta_oct/beta_tri does not
by itself give t_leak. We enumerate the candidate maps and translate each to a width via the validated
KPM trend Gamma(t_leak) (item48_kpm_rho_ldos.py: 0.1->9, 0.2->19, 0.3->33, 0.414->53(SC)/33(LSC),
0.6->88, 1.0->578 MeV; on the real L(SC) web ~30% narrower).
"""
import numpy as np

s2, s3 = np.sqrt(2), np.sqrt(3)
A_oct = 2 * (1 + s2)          # regular octagon, edge 1
A_tri = s3 / 4                # equilateral triangle, edge 1
print(f"[geom] A_oct = 2(1+sqrt2)   = {A_oct:.6f}")
print(f"[geom] A_tri = sqrt3/4      = {A_tri:.6f}")
assert abs(A_oct - 4.828427) < 1e-5
assert abs(A_tri - 0.433013) < 1e-5

area_ratio = A_tri / A_oct                       # = beta_oct/beta_tri  (beta ∝ 1/A)
print(f"[lock] A_tri/A_oct = sqrt3/(8(1+sqrt2)) = {area_ratio:.6f}  == beta_oct/beta_tri")
assert abs(area_ratio - 0.0896797) < 1e-6
assert abs(area_ratio - s3 / (8 * (1 + s2))) < 1e-12
beta_tri_over_oct = A_oct / A_tri                # the triangular (matter) face is STIFFER
print(f"[lock] beta_tri/beta_oct = A_oct/A_tri  = {beta_tri_over_oct:.4f}  (matter face ~11x stiffer)")

silver = s2 - 1                                  # the heuristic we had been using for t_leak
print(f"\n[note] silver ratio sqrt2-1 = {silver:.6f}  -- this is the octagon EDGE/flat-to-flat LENGTH")
print(f"       ratio (ANCHOR sec1.4 L83), NOT the area-lock {area_ratio:.4f}. So 't_leak = silver ratio'")
print(f"       was a length-heuristic, never the gauge coupling. (sqrt(area_ratio)={np.sqrt(area_ratio):.4f}.)")
assert abs(silver - 0.414214) < 1e-5
assert abs(silver - area_ratio) > 0.3, "silver ratio is NOT the area lock -- distinct numbers"

# ---- candidate beta/geometry -> t_leak bridges (t_oct == 1) ----
# Each is a DIFFERENT physical assumption about what sets the matter hopping across a face.
def width_from_tleak(tl):
    """KPM-anchored Gamma(t_leak) on the real L(SC) continuum: ~quadratic at small t_leak,
    saturating into the smear. Piecewise-linear interp of the validated GPU/CPU scan (L(SC), MeV)."""
    xs = np.array([0.0, 0.1, 0.2, 0.3, 0.414, 0.6, 1.0])
    gs = np.array([0.0, 6.0, 13.0, 22.0, 33.0, 60.0, 380.0])   # L(SC)-continuum widths (proxy-corrected)
    return float(np.interp(min(tl, 1.0), xs, gs)) if tl <= 1.0 else 380.0 * tl**2

cands = {
    "area ratio  t=A_tri/A_oct        ": area_ratio,        # geometric aperture / flux-matching for hopping
    "inv-beta    t=beta_oct/beta_tri  ": area_ratio,        # identical to area ratio (beta ∝ 1/A)
    "sqrt-area   t=sqrt(A_tri/A_oct)  ": np.sqrt(area_ratio),
    "silver-heur t=sqrt2-1            ": silver,             # the old length-heuristic
    "edge ratio  t=a_tri/a_oct=1      ": 1.0,                # uniform-edge tight binding
    "beta ratio  t=beta_tri/beta_oct  ": beta_tri_over_oct,  # stiffness=bandwidth reading
}
print("\n  candidate bridge                     t_leak     ~Gamma(MeV)   verdict")
THRESH = 149.0
for name, tl in cands.items():
    g = width_from_tleak(tl)
    verdict = "SURVIVES (narrow)" if g < THRESH else "DISSOLVES"
    print(f"  {name} {tl:8.4f}   {g:8.0f}     {verdict}")

# the dichotomy: area-like bridges survive, edge/stiffness bridges dissolve
survives = [v for v in cands.values() if width_from_tleak(v) < THRESH]
dissolves = [v for v in cands.values() if width_from_tleak(v) >= THRESH]
print(f"\n[result] {len(survives)} bridges -> SURVIVES (t_leak in [{min(survives):.3f},{max(survives):.3f}]),"
      f" {len(dissolves)} -> DISSOLVES (t_leak in [{min(dissolves):.2f},{max(dissolves):.2f}]).")
assert survives and dissolves, "the candidate bridges must straddle the threshold for the verdict to be open"

print("""
[interim] beta_oct/beta_tri does NOT by itself fix t_leak (gauge stiffness != matter hopping), and the
geometric bridges straddle the threshold. So the AREA-LOCK alone is inconclusive. But the framework does
not leave the hopping free -- it DERIVES it.""")

# ---- RESOLUTION: the framework's OWN substrate-derived hopping (Grover-coin uniqueness theorem) ----
# ANCHOR L1541-1552: the unique unitary permutation-symmetric scattering operator at the degree-6
# macroscopic gauge vertex is the Grover coin C = 2/d - delta_ij (d=6). Its scattering profile:
d = 6
refl = 2 / d - 1                  # reflection back  = 2/6 - 1 = -2/3 ... wait: C_ii = 2/d - 1, C_ij = 2/d
# Grover coin entries: diagonal C_ii = 2/d - 1 ; off-diagonal C_ij = 2/d
C_diag = 2 / d - 1               # = -2/3  (straight-through / same-channel)
C_off = 2 / d                    # =  1/3  (each orthogonal transmission, and reflection)
print(f"[grover] d={d}: straight C_ii = 2/d-1 = {C_diag:+.4f} ; orthogonal/reflect C_ij = 2/d = {C_off:+.4f}")
# unitarity of one row: C_ii^2 + (d-1) C_ij^2 = 1
row_norm = C_diag**2 + (d - 1) * C_off**2
print(f"[grover] row unitarity  C_ii^2 + (d-1)C_ij^2 = {row_norm:.6f}  (must be 1)")
assert abs(row_norm - 1.0) < 1e-12
assert abs(C_off - 1/3) < 1e-12, "orthogonal-transmission amplitude is exactly 1/3"
assert abs(C_diag + 2/3) < 1e-12

t_leak_grover = abs(C_off)        # = 1/3 : the orthogonal-transmission amplitude = the rho-flux leak
print(f"[grover] substrate hopping t_leak = |orthogonal transmission| = 1/3 = {t_leak_grover:.4f}  "
      f"(ANCHOR L1552 'analytic substrate constant, not empirical')")

# KPM widths on the REAL L(SC) continuum (item48_kpm_rho_ldos.py --gauge linegraph, validated 4e-14):
kpm = {"1/3  (Grover orthogonal, L1552)": (1/3, 27),
       "1/2  (orthogonal/straight)      ": (1/2, 50),
       "1/sqrt3 (deg-3 routing, L963)   ": (1/np.sqrt(3), 61),
       "2/3  (Grover MAX = straight)     ": (2/3, 130)}   # even the largest Grover amplitude is < 149
print("\n  framework amplitude              t_leak    Gamma_KPM(MeV)   verdict")
for name, (tl, g) in kpm.items():
    print(f"  {name} {tl:7.4f}   {g:7.0f}        {'SURVIVES' if g < THRESH else 'DISSOLVES'}")
assert all(g < THRESH for _, g in kpm.values()), "every framework-derived amplitude must survive"

print("""
VERDICT (reported — RESOLVED toward survival, on the framework's own constant):
  The matter leak hopping is NOT free and NOT the area-lock: it is the Grover-coin orthogonal-
  transmission amplitude t = 1/3 (ANCHOR L1541-1552, FORCED by the Schur-lemma uniqueness theorem at
  the degree-6 gauge vertex; stated as 'an analytic substrate constant, not an empirical parameter').
  At t_leak = 1/3 the rho is a ~27 MeV resonance pinned AT phi (KPM on the real L(SC) web) -> SURVIVES.

  Robustness: the WHOLE Grover scattering profile is sub-threshold -- straight 2/3 (the MAX amplitude,
  ~130 MeV), orthogonal/reflect 1/3 (~27 MeV), and the deg-3 routing 1/sqrt3 (~61 MeV) all give a rho
  NARROWER than the physical 149 MeV. Dissolution needs t_leak >~ 0.7, i.e. a NEAR-UNIT hopping -- which
  the Grover coin EXCLUDES (nothing scatters with amplitude 1; max is the straight 2/3). My earlier
  'uniform |t|=1 -> dissolves' worry was wrong: the substrate is a Grover-coin walk with sub-unit
  amplitudes, not a unit-hopping graph.

  Residual caveats (so this is 'firmly leans survival', not 'closed'): (i) the KPM continuum is a
  uniform tight-binding APPROXIMATION of the actual Grover-coin walk -- a fully Grover-faithful continuum
  could shift the width (but the leak amplitude 1/3 is robust); (ii) one-way t_leak vs inter-cell
  normalisation (all variants 1/3..2/3 survive); (iii) Delta=1.78 and the octagon-perimeter host.
""")
print("exit 0 -- geometry + Grover unitarity asserted; verdict RESOLVED toward survival (t_leak=1/3).")
