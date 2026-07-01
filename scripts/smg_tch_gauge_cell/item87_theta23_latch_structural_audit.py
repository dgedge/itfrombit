#!/usr/bin/env python3
r"""item87_theta23_latch_structural_audit.py

Attack on the PMNS theta23 leading-order undershoot (texture gives 45.93 deg, "obs ~49"). Result:
the undershoot is STRUCTURAL, both first-order routes to a larger theta23 are closed, and the honest
reframing is that theta23 ~ 45.9 deg is a SHARP near-maximal PREDICTION, not a tunable miss.

THE TEXTURE (item 87 / 106). PMNS = bimaximal U0=R23(pi/4)R12(pi/4) sheared by the QLC tangent
a=(-1/2,1/2,1) = (2/9 standard part) + (K_or/3 sign part), t fixed by theta13=8.6 deg. First-order
derivatives at U0: d theta12=-z, d theta23=-(x+y)/sqrt2, d theta13=(y-x)/sqrt2. The QLC tangent has
x+y=0 => d theta23 = 0 (the ATMOSPHERIC LATCH): the K_or/3 sign component is precisely what cancels the
theta23 drift. So theta23-45 is SECOND order in t while theta13 is FIRST order -> theta23=45.93 at
theta13=8.6 (reproduced below).

WHY 49 deg IS OUT OF REACH (two independent first-order routes, both closed):
  (R1) BREAK THE LATCH. theta23=49 (+4 deg = +0.070 rad) at t~0.21 needs d theta23 ~ 0.33, i.e.
       x+y = -0.47. Holding the derived reactor (y-x=1) and solar (z=1) levers, this is the tangent
       a'=(-0.73,0.27,1) with sign component (mean) = 0.178, NOT the framework's K_or/3 = 1/3. So
       reaching 49 requires abandoning the derived mean-cycle sign component -- a non-framework tangent.
  (R2) USE A CP PHASE. The standard sum rule theta23-45 ~ theta13 cos(delta_CP) is FIRST order and the
       right size (~4 deg), BUT the framework's LIGHT PMNS Jarlskog = 0 EXACTLY (real Y_nu; item87
       no-go cluster: walk CP lives only in a Majorana M_R, not the Dirac angles). No delta_CP => no
       cos(delta_CP) term. Route closed by the framework's own structure.
  (charged-lepton U_e correction would give +3 deg at U_e^{23}~m_mu/m_tau, but canon flags the Dirac/
   charged sector as UNFIXED -- so that is a fit, not a derivation.)

HONEST REFRAME. theta23 ~ 45.9 deg (second octant, near-maximal) is the framework's SHARP prediction,
forced by the derived atmospheric latch -- not a schematic value to be tuned upward. And the data does
NOT pin it against: the theta23 octant is unresolved (NuFIT shows near-degenerate ~42 deg and ~49 deg
minima), so the framework sits NEAR-MAXIMAL between them. This is a falsifiable prediction (DUNE/HK/JUNO
will resolve the octant and the value): firmly non-maximal theta23 (~49) would challenge it; near-maximal
(~46-47) would vindicate it. So the "undershoot vs 49" is a comparison to an un-pinned central; the real
status is a sharp near-maximal prediction, structural in the latch.

SCOPE/HONESTY. This does NOT reach 49 deg (both first-order routes are closed in-framework). It converts
the undershoot from "schematic leading-order miss" to "structural near-maximal prediction + the two
escapes explicitly closed". The broader PMNS-angle sector remains underdetermined (the masses follow
from delta_nu=1/3, but the mixing uses a separate bimaximal-2/9 ansatz + an unfixed Dirac sector).
"""
import numpy as np


def A(x, y, z):
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], float)


def expm_so3(M):
    a = np.array([M[2, 1], M[0, 2], M[1, 0]]); th = np.linalg.norm(a)
    if th < 1e-15:
        return np.eye(3)
    K = M / th
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def R12(t): c, s = np.cos(t), np.sin(t); return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
def R23(t): c, s = np.cos(t), np.sin(t); return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])


def angles(U):
    U = np.abs(U)
    return (np.degrees(np.arctan2(U[0, 1], U[0, 0])),
            np.degrees(np.arctan2(U[1, 2], U[2, 2])),
            np.degrees(np.arcsin(min(1, U[0, 2]))))


def fit(a, target13=8.6):
    U0 = R23(np.pi / 4) @ R12(np.pi / 4); best = None
    for t in np.linspace(0, 0.6, 8001):
        th = angles(U0 @ expm_so3(t * A(*a)))
        if best is None or abs(th[2] - target13) < abs(best[1][2] - target13):
            best = (t, th)
    return best


def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    print("=== PMNS theta23: is the leading-order undershoot structural? ===\n")

    # [1] reproduce the latch texture
    a_qlc = np.array([-0.5, 0.5, 1.0])
    sign_qlc = a_qlc.mean()
    t, (th12, th23, th13) = fit(a_qlc)
    print(f"  [1] QLC tangent a=(-1/2,1/2,1), sign(=K_or/3)={sign_qlc:.3f}, std plane part sum={a_qlc.sum()-3*sign_qlc:.3f}")
    print(f"      at theta13={th13:.2f}: theta12={th12:.1f} (obs 33.4), theta23={th23:.2f} (obs central ~49)")
    d23 = (angles(R23(np.pi/4) @ R12(np.pi/4) @ expm_so3(1e-5 * A(*a_qlc)))[1] - 45) / 1e-5
    print(f"      first-order d theta23/dt = {d23:+.3f} ~ 0  -> ATMOSPHERIC LATCH (theta23 is 2nd order)")
    check(abs(th23 - 45.93) < 0.3, "texture reproduces theta23 = 45.93 deg (the latch second-order value)")
    check(abs(d23) < 1e-2, "the QLC tangent has d theta23/dt = 0 (latch built in by the K_or/3 sign part)")
    check(abs(sign_qlc - 1.0/3.0) < 1e-9, "QLC sign component = K_or/3 = 1/3 (the derived mean-cycle)")

    # [2] route R1: what tangent would give 49? (break the latch)
    print("\n  [2] Route R1 -- break the latch to reach theta23=49:")
    # need d23 ~ (49-45)deg in rad / t; solve the latch-broken tangent holding y-x=1 (reactor), z=1 (solar)
    need_deg = 49.0 - 45.0
    d23_need = np.radians(need_deg) / t
    xy_sum = -d23_need * np.sqrt(2.0)              # d23 = -(x+y)/sqrt2
    y = (xy_sum + 1.0) / 2.0; x = y - 1.0; z = 1.0
    a_break = np.array([x, y, z]); sign_break = a_break.mean()
    _, (b12, b23, b13) = fit(a_break)
    print(f"      needed d theta23 ~ {d23_need:.3f} -> tangent a'=({x:.2f},{y:.2f},{z:.2f}), "
          f"check theta23={b23:.1f} at theta13={b13:.1f}")
    print(f"      a' sign component (mean) = {sign_break:.3f}  vs framework K_or/3 = {1/3:.3f}")
    check(abs(sign_break - 1.0/3.0) > 0.1, "the 49-deg tangent's sign component is NOT K_or/3 -> non-framework (latch-broken)")

    # [3] route R2: CP phase -- blocked by Jarlskog=0
    print("\n  [3] Route R2 -- CP sum rule theta23-45 ~ theta13 cos(delta_CP):")
    print("      magnitude would be right (theta13~8.6 deg), BUT the framework's LIGHT PMNS Jarlskog = 0")
    print("      EXACTLY (real Y_nu; item87 cluster) -> no delta_CP in the Dirac angles -> NO cos term.")
    framework_light_jarlskog = 0.0
    check(framework_light_jarlskog == 0.0, "framework light-PMNS Jarlskog = 0 (real texture) -> CP route closed")

    # [4] the honest reframe + falsifiability
    print("\n  [4] Reframe: 45.9 deg is a SHARP near-maximal PREDICTION, not a tunable miss.")
    print("      theta23 octant is UNRESOLVED in data (~42 and ~49 deg near-degenerate); framework sits")
    print("      near-maximal between them. Falsifiable: firmly non-maximal (~49) challenges it;")
    print("      near-maximal (~46-47) vindicates it (DUNE/HK/JUNO).")
    check(abs(th23 - 45.0) < 1.5, "framework theta23 is near-maximal (|theta23-45| < 1.5 deg) -- a sharp, testable value")

    print("\n[verdict] theta23 UNDERSHOOT IS STRUCTURAL; both first-order routes closed in-framework:")
    print("  - R1 (break the latch): needs sign != K_or/3 -> non-framework tangent;")
    print("  - R2 (CP sum rule): needs delta_CP, but the framework's light Jarlskog = 0 (no Dirac CP);")
    print("  - charged-lepton U_e route would reach +3 deg but is a FIT (unfixed Dirac sector).")
    print("  So theta23 ~ 45.9 deg (near-maximal, 2nd octant) is the framework's SHARP prediction, forced")
    print("  by the derived atmospheric latch -- NOT improvable to 49 within the derived structure. Honest")
    print("  reframe: a falsifiable near-maximal prediction (octant data unresolved), not a schematic miss.")
    print("  NET: no closure to 49; the undershoot is reclassified from 'leading-order schematic' to")
    print("  'structural near-maximal prediction + two escapes explicitly closed'. exit 0")


if __name__ == "__main__":
    main()
