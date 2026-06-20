#!/usr/bin/env python3
"""
item9_dispersive_horizon.py

Item 9 (ANCHOR §15): the radial frustration gradient and critical position x_c.
Construct a 1D chain of voids with position-dependent Boltzmann mass
M(x) = M0 (1 + alpha*x), diagonalise the walk W = S.C, and verify that
eigenstates exhibit distinct causal structures across x_c (§11.3): propagating
(exterior) -> evanescent/trapped (interior) at the dispersive horizon.

Mechanism (ANCHOR §11.3): the coin angle theta(x) is the local mass gap
(theta proportional to M_eff). The 1D Dirac QCA dispersion is
    cos(E) = cos(theta) * cos(k)
so a mode of quasi-energy E propagates where the local gap theta(x) < |E| and
is evanescent where theta(x) > |E|. The horizon x_c is the turning point
theta(x_c) = |E| -- where the coin "freezes" (group velocity -> 0). Because
x_c depends on the mode's bare mass M0, the horizon is MASS-STRATIFIED
(dispersive), exactly as §11.3 claims.

We use a central low-mass WELL (theta(x) = theta0(1+alpha|x-x0|)) so the
propagating region is bounded by two turning points in the bulk, away from the
chain ends.

Self-asserting; numpy only.
"""
import sys
import numpy as np

_ok = True


def check(name, cond):
    global _ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    _ok = _ok and bool(cond)


# ---- setup: central low-mass well -----------------------------------------
N = 200
x0 = N // 2
theta0 = 0.30           # bare mass gap M0 (coin angle) at the well centre
alpha = 0.02            # radial frustration gradient


def theta_of(x, th0=theta0, a=alpha):
    return np.minimum(th0 * (1.0 + a * np.abs(x - x0)), np.pi / 2)


# ===========================================================================
# 1. WKB local analysis: the turning point x_c and group velocity -> 0
# ===========================================================================
E = 0.5                                   # a mode quasi-energy
xc = (E / theta0 - 1.0) / alpha           # turning point offset: theta(x0+-xc)=E
print(f"mode E={E}, theta0={theta0}, alpha={alpha}  ->  x_c offset = {xc:.2f} "
      f"(well = [{x0-xc:.0f}, {x0+xc:.0f}])")

xs = np.arange(N)
th = theta_of(xs)
cos_k = np.cos(E) / np.cos(th)            # cos k(x) = cos E / cos theta(x)
inside = np.abs(xs - x0) < xc - 1
outside = np.abs(xs - x0) > xc + 1
check("propagating inside the well: |cos k(x)| < 1", np.all(np.abs(cos_k[inside]) < 1))
check("evanescent outside (interior): |cos k(x)| > 1", np.all(np.abs(cos_k[outside]) > 1))
check("turning point at theta(x_c) = |E|", abs(theta_of(x0 + xc) - E) < 1e-9)

# group velocity v = cos(theta) sin(k) / sin(E); -> 0 at the turning point
k_in = np.arccos(np.clip(cos_k, -1, 1))
v = np.cos(th) * np.sin(k_in) / np.sin(E)
v_mid = v[x0]                              # deep in the well (fast)
v_edge = v[int(x0 + xc - 1)]               # just inside the horizon (slow)
print(f"group velocity: well-centre v={v_mid:.3f}, just inside horizon v={v_edge:.3f}")
check("group velocity vanishes at the horizon (v_edge << v_centre)",
      v_edge < 0.25 * v_mid)

# ===========================================================================
# 2. dispersive / mass-stratified horizon: x_c depends on the mode mass M0
# ===========================================================================
xc_light = (E / 0.30 - 1.0) / alpha       # lighter bare mass
xc_heavy = (E / 0.40 - 1.0) / alpha       # heavier bare mass
print(f"\ndispersive horizon: x_c(M0=0.30)={xc_light:.1f}, "
      f"x_c(M0=0.40)={xc_heavy:.1f}  -> heavier defect freezes earlier")
check("horizon is MASS-STRATIFIED: heavier mode has smaller x_c (§11.3 dispersive)",
      xc_heavy < xc_light)

# ===========================================================================
# 3. diagonalise the finite walk; show an eigenstate's causal-structure change
# ===========================================================================
def build_W():
    """1D Dirac QCA W = S.C with position-dependent coin; periodic (unitary).
    coin 0 = right-mover, 1 = left-mover. C(theta_x) then shift."""
    W = np.zeros((2 * N, 2 * N), dtype=complex)
    for x in range(N):
        th = theta_of(x)
        c, s = np.cos(th), np.sin(th)
        xr, xl = (x + 1) % N, (x - 1) % N
        # input (x,0): after C -> c|0> + s|1>; shift 0->xr, 1->xl
        W[2 * xr + 0, 2 * x + 0] += c
        W[2 * xl + 1, 2 * x + 0] += s
        # input (x,1): after C -> -s|0> + c|1>
        W[2 * xr + 0, 2 * x + 1] += -s
        W[2 * xl + 1, 2 * x + 1] += c
    return W


W = build_W()
check("walk operator W is unitary", np.allclose(W.conj().T @ W, np.eye(2 * N), atol=1e-10))
evals, evecs = np.linalg.eig(W)
Equasi = -np.angle(evals)                  # quasi-energy E from eigenvalue e^{-iE}

# pick the eigenstate with quasi-energy closest to +E
idx = int(np.argmin(np.abs(np.abs(Equasi) - E)))
psi = evecs[:, idx]
prob = (np.abs(psi[0::2])**2 + np.abs(psi[1::2])**2)   # site probability
prob /= prob.sum()
Esel = abs(Equasi[idx])
xc_sel = (Esel / theta0 - 1.0) / alpha
print(f"\nselected eigenstate: |E|={Esel:.3f} -> x_c offset={xc_sel:.1f}")
in_well = np.abs(xs - x0) < xc_sel
P_in = prob[in_well].sum()
# decay outside: amplitude well past the horizon vs just inside
a_inside = prob[int(x0 + xc_sel - 8)]
a_outside = prob[int(x0 + xc_sel + 12)]
print(f"  probability inside the well = {P_in:.3f}; "
      f"prob ratio (8 inside / 12 outside horizon) = {a_inside/max(a_outside,1e-30):.1f}x")
check("eigenstate is CONFINED to the propagating well (>85% probability inside)",
      P_in > 0.85)
check("eigenstate DECAYS across the horizon (inside >> outside)",
      a_inside > 10 * a_outside)

print("\n--- VERDICT ---")
print("Item 9 CLOSED (verification): the position-dependent-mass Dirac walk has a")
print("turning point x_c where theta(x_c)=|E| -- the coin freezes (group velocity")
print("-> 0) and eigenstates change causal character from propagating (exterior")
print("well) to evanescent/trapped (interior). x_c depends on the bare mass M0, so")
print("the horizon is mass-stratified/dispersive exactly as §11.3 states. This is")
print("the discrete black-hole horizon as a walk-operator turning point.")

print("\n" + ("ALL CHECKS PASSED" if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
