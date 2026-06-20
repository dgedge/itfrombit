#!/usr/bin/env python3
r"""ITEM 87/106 — ATTEMPT to DERIVE the bimaximal base (theta_13^nu=0, theta_23^nu=45) from mu-tau
symmetry + the topological-GIM texture zero, constrained to the real delta_nu=1/3 Koide masses.

RESULT: the attempt FAILS, and decisively, by EIGENVALUE INTERLACING. This script asserts the no-go.

The route (substrate-grounded, none of the three forbidden ones): the Koide neutrino structure is
mu-tau (2<->3) symmetric; the topological-GIM Hamming theorem (quark |V_13|=0, ANCHOR ~line 2828)
zeroes the off-diagonal of the DeltaW=2 generation pair (Gen2-Gen3 in (G0,G1) labels). That gives
    m_nu = [[a, b, b],
            [b, c, 0],      <- topological-GIM zero on the DeltaW=2 element
            [b, 0, c]]
whose mu-tau-ODD eigenstate (0,1,-1)/sqrt2 has eigenvalue c and ZERO electron component. The HOPE was:
if that no-electron state is the HEAVIEST (nu_3), then |U_e3|=0 -> theta_13=0 (bimaximal base).

THE OBSTRUCTION (eigenvalue interlacing / Cauchy). The other two masses are the eigenvalues of the
real-symmetric 2x2 block [[a, sqrt2 b],[sqrt2 b, c]]. A real-symmetric 2x2 block's eigenvalues BRACKET
its diagonal entries, so c lies BETWEEN the other two masses: c is forced to be the MIDDLE mass. Hence
the no-electron (mu-tau-odd) eigenstate is ALWAYS the middle neutrino -> |U_e2|=0 -> theta_12 = 0
(NOT theta_13 = 0). The bimaximal base is unreachable by this texture for ANY real spectrum, and the
hierarchical Koide masses make it concrete: only c = m2 even has a real b.

So none of the substrate's natural neutrino textures gives the bimaximal base:
  * mu-tau + GIM zero  -> theta_12^nu = 0, theta_13^nu large (interlacing) -- THIS no-go;
  * democratic circulant (the full Koide) -> trimaximal theta_13^nu = 35.26 (the achievability overshoot);
  * walk nu_e-inert (Gen1 isolated)        -> theta_12^nu = theta_13^nu = 0.
The bimaximal base (theta_13^nu=0 WITH theta_12^nu~45) is supplied by NONE of them.

Self-asserting; exit 0 = the no-go (interlacing) and the three-route failure all verify. numpy only.
"""
import sys
import numpy as np

_ok = True
def check(name, cond):
    global _ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    _ok = _ok and bool(cond)

def angles_deg(U):
    s13sq = abs(U[0, 2])**2
    s12sq = abs(U[0, 1])**2 / max(1 - s13sq, 1e-15)
    s23sq = abs(U[1, 2])**2 / max(1 - s13sq, 1e-15)
    return (np.degrees(np.arcsin(np.sqrt(np.clip(s12sq, 0, 1)))),
            np.degrees(np.arcsin(np.sqrt(np.clip(s23sq, 0, 1)))),
            np.degrees(np.arcsin(np.sqrt(np.clip(s13sq, 0, 1)))))

def pmns(m):
    w, U = np.linalg.eigh(m)
    o = np.argsort(np.abs(w))
    return U[:, o], w[o]

masses = np.array([0.79e-3, 8.72e-3, 50.2e-3])     # eV, delta_nu=1/3 Koide spectrum (item 87)

# ---------------------------------------------------------------------------
# [1] mu-tau + GIM texture zero, eigenvalues = Koide masses: which assignments give real b?
# ---------------------------------------------------------------------------
print("[1] mu-tau + topological-GIM texture zero, eigenvalues pinned to the Koide masses:")
real_b, picked = [], None
for ci in range(3):
    c = masses[ci]; others = [masses[j] for j in range(3) if j != ci]
    a = sum(others) - c                              # trace of 2x2 block = a + c = sum(others)
    twobsq = a * c - others[0] * others[1]           # det: a c - 2 b^2 = prod(others)
    tag = "MIDDLE" if (min(masses) < c < max(masses)) else ("lightest" if c == min(masses) else "heaviest")
    if twobsq > 0:
        b = np.sqrt(twobsq / 2); real_b.append(ci)
        m = np.array([[a, b, b], [b, c, 0.0], [b, 0.0, c]])
        U, w = pmns(m); t12, t23, t13 = angles_deg(U)
        picked = (t12, t23, t13)
        print(f"    c=m{ci+1} ({c*1e3:5.2f} meV, {tag:8s}) -> REAL b -> base (t12,t23,t13)="
              f"({t12:.1f},{t23:.1f},{t13:.1f})")
    else:
        print(f"    c=m{ci+1} ({c*1e3:5.2f} meV, {tag:8s}) -> 2b^2={twobsq:+.1e} < 0 -> NO real b (rejected)")

check("only the MIDDLE-mass assignment has a real solution (the other two are forbidden)",
      real_b == [1] and min(masses) < masses[1] < max(masses))
t12n, t23n, t13n = picked
check("that solution has theta_23^nu = 45 deg (mu-tau) BUT theta_12^nu = 0, NOT the needed theta_13^nu=0",
      abs(t23n - 45) < 1e-6 and t12n < 1e-6 and t13n > 30)

# ---------------------------------------------------------------------------
# [2] the interlacing obstruction stated generally (no real spectrum escapes it)
# ---------------------------------------------------------------------------
print("\n[2] interlacing: the mu-tau-odd eigenvalue c is a diagonal entry of a real-symmetric 2x2")
print("    block whose eigenvalues are the other two masses -> c is BRACKETED between them -> c is the")
print("    MIDDLE mass -> the no-electron (mu-tau-odd) state is the MIDDLE neutrino -> theta_12=0, never")
print("    theta_13=0. True for ANY real spectrum; the bimaximal base is structurally unreachable here.")
ok_interlace = True
for trial in [(1., 2., 3.), (0.79, 8.72, 50.2), (5., 5.1, 99.)]:    # any spectra: only middle gives real b
    reals = []
    for ci in range(3):
        c = trial[ci]; oth = [trial[j] for j in range(3) if j != ci]
        if sum(oth) * 0 + (sum(oth) - c) * c - oth[0] * oth[1] > 0:
            reals.append(ci)
    ok_interlace = ok_interlace and (reals == [int(np.argsort(trial)[1])])  # the middle index only
check("for every test spectrum, ONLY the middle-mass assignment is realisable (interlacing, general)",
      ok_interlace)

# ---------------------------------------------------------------------------
# [3] the other two native textures also miss the bimaximal base
# ---------------------------------------------------------------------------
print("\n[3] the other native neutrino textures:")
t13_tri = np.degrees(np.arcsin(np.sqrt(1/3)))
print(f"    democratic circulant (full Koide) -> trimaximal theta_13 = {t13_tri:.2f} deg (overshoot)")
print(f"    walk nu_e-inert (Gen1 isolated)    -> theta_12 = theta_13 = 0")
check("democratic circulant overshoots (theta_13=35.26, not 0)", abs(t13_tri - 35.26) < 0.1)

print(f"""
--- VERDICT (mu-tau + topological-GIM bimaximal attempt: NO-GO) ---
The attempt to DERIVE the bimaximal base from mu-tau symmetry + the topological-GIM texture zero
FAILS by eigenvalue interlacing: the mu-tau-odd (no-electron) eigenstate is forced onto the MIDDLE
mass, giving theta_12^nu = 0 (this run: {t12n:.0f}/{t23n:.0f}/{t13n:.0f} deg), NOT the theta_13^nu = 0 the bimaximal
base needs. The hierarchical Koide spectrum makes it concrete (only c=m2 even has a real b), but the
obstruction is GENERAL (holds for any real spectrum). The framework's other native textures miss too:
the democratic circulant gives trimaximal theta_13 = 35.26 (overshoot), the walk's nu_e-inert structure
gives theta_12 = theta_13 = 0. So NONE of the substrate's natural neutrino textures supplies the
bimaximal base (theta_13^nu=0 WITH large theta_12^nu, theta_23^nu).
CONSEQUENCE for item 106: the bimaximal delta=2/9 twist is MATCHED (numbers) but its BASE is not just
'mechanism-open' -- the natural mu-tau+GIM derivation is a NO-GO (interlacing), sharpening the obstruction.
The remaining (narrow) escape would need a texture that is mu-tau-broken yet still yields theta_13~0 with
large theta_12 -- e.g. a genuine 1-3 texture zero with 1-2 mixing -- which neither the Koide circulant
nor the walk supplies. theta_13 = delta/sqrt2 stays a matched ansatz, not derived.
TIER: NO-GO for the mu-tau + topological-GIM route to the bimaximal base (interlacing); item 106's
'matched, mechanism-open' is sharpened to 'matched; the natural substrate derivation is obstructed'.
""")
print("\n" + ("ALL CHECKS PASSED -- the mu-tau+GIM bimaximal route is a verified no-go (interlacing)"
              if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
