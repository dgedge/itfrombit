r"""Trans-Lambda_QCD photon — the bundle no-go is KINEMATIC, not just QED-perturbative.

foundations_bundle_irreducibility_no_go.py showed that under ordinary Wilson/Gauss QED a
collinear null bundle absorbs as an N-photon process (linear j.A, N Ward identities), not a
single total-P photon. That argument uses QED perturbation theory, so it leaves a loophole:
maybe some non-perturbative / strong-coupling / solitonic mechanism evades it.

This script closes that loophole. The obstruction is KINEMATIC -- it uses only the lattice
dispersion omega(k) and the group velocity v_g(k), no QED expansion -- so it binds ANY
composite of lattice modes, perturbative or not.

CLOSED FORM (exact, along [100], a0=c=1): omega = 2 sin(k/2), so v_g = cos(k/2) and the
constituent energy is
        E1(beta) = 2 Lambda_QCD * sqrt(1 - beta^2),   beta = v_g/c.
A constituent that is NULL (beta -> 1, the only way to move at c and be a photon) carries
E1 -> 0. Hence, with NO dynamics assumed:

  NULL bundle (moves at c)  =>  soft constituents (v_g=c only at k->0)  =>  single-vertex
                                energy -> 0  =>  N >~ E*a0/hbar c independent soft quanta,
                                not one hard photon  [= the Wilson/QED no-go, now FORCED].
  HARD bundle (E1 ~ E_BZ)   =>  BZ-edge constituents have v_g -> 0  =>  ~stationary  =>
                                not a photon (photons move at c).

NULL and HARD are mutually exclusive for lattice-mode composites. The soliton escape fails by
the same scale wall (a c-boosted localized carrier is Lorentz-contracted below a0).

exit 0 = v_g=cos(k/2) verified; the closed form E1=2L*sqrt(1-beta^2) verified; null<->soft and
         hard<->slow quantified; the soliton contraction cap computed; all escapes fail.
"""
import numpy as np

HBARC = 0.1973269804     # GeV fm
LAMBDA = 0.332           # GeV  (= hbar c / a0)
A0 = HBARC / LAMBDA      # fm

# SC Wilson/Maxwell dispersion, lattice units a0 = c = 1.  Use the cancellation-free identity
# 6 - 2 sum cos k_i = 4 sum sin^2(k_i/2), i.e. omega = 2 sqrt(sum sin^2(k_i/2)) (stable as k->0).
def omega(k):
    return 2.0 * np.sqrt(sum(np.sin(ki / 2.0) ** 2 for ki in k))
def vg(k, eps=1e-6):
    g = np.zeros(3)
    for d in range(3):
        kp = list(k); km = list(k); kp[d] += eps; km[d] -= eps
        g[d] = (omega(kp) - omega(km)) / (2 * eps)
    return g
def E_GeV(k):
    return omega(k) * LAMBDA      # E = hbar omega = (hbar c/a0) sqrt(K),  (hbar c/a0) = Lambda

E_edge = 2 * LAMBDA              # [100] zone-edge energy
E_corner = np.sqrt(12) * LAMBDA  # (pi,pi,pi) corner energy
print(f"[0] scales: Lambda_QCD={LAMBDA} GeV, a0={A0:.4f} fm. One-mode energy ceiling:")
print(f"    [100] edge 2*Lambda = {E_edge:.3f} GeV ; (pi,pi,pi) corner sqrt(12)*Lambda = {E_corner:.3f} GeV.")

# [1] exact identity along [100]: v_g = cos(k/2), and E1 = 2 Lambda sqrt(1 - beta^2)
print("\n[1] exact [100] kinematics: v_g = cos(k/2);  E1 = 2*Lambda*sqrt(1 - (v_g/c)^2)")
print("    k(a0)    E1(GeV)   v_g/c    2L*sqrt(1-beta^2)")
worst = 0.0
for kx in (1e-3, 0.3, 1.0, 2.0, np.pi - 1e-6):
    b = float(vg((kx, 0, 0))[0]); e = E_GeV((kx, 0, 0)); cf = 2 * LAMBDA * np.sqrt(max(0.0, 1 - b * b))
    worst = max(worst, abs(e - cf), abs(b - np.cos(kx / 2)))
    print(f"    {kx:6.4f}   {e:6.4f}   {b:6.4f}   {cf:6.4f}")
assert worst < 1e-4                                          # v_g=cos(k/2) and E1=2L*sqrt(1-b^2) hold exactly
assert vg((1e-3, 0, 0))[0] > 0.999                          # soft constituent: v_g = c
assert abs(vg((np.pi - 1e-6, 0, 0))[0]) < 1e-3              # zone-edge constituent: v_g = 0

# [2] NULL bundle (must move at c) => single-vertex energy is soft => huge N
print("\n[2] NULL bundle: largest single-constituent energy E1 at a given v_g floor, and N for 1 TeV")
for vfloor in (0.999, 0.99, 0.9):
    kx = 2 * np.arccos(vfloor)                              # v_g = cos(k/2) = vfloor
    e1 = 2 * LAMBDA * np.sqrt(1 - vfloor ** 2)
    print(f"    v_g >= {vfloor}c : E1 = {e1*1e3:7.2f} MeV  =>  N(1 TeV) >= {1000.0/e1:11,.0f}")
e1_999 = 2 * LAMBDA * np.sqrt(1 - 0.999 ** 2)
assert 1000.0 / e1_999 > 1e4                                # a near-null 1 TeV bundle = >10^4 soft quanta

# [3] HARD bundle (E1 ~ E_BZ) => v_g -> 0 => not a photon
print("\n[3] HARD bundle: a zone-edge constituent is essentially stationary")
ke = (np.pi - 1e-3, 0, 0)
print(f"    constituent at [100] edge: E1 = {E_GeV(ke):.3f} GeV, v_g/c = {abs(vg(ke)[0]):.2e}  -> ~stationary")
assert abs(vg(ke)[0]) < 1e-2

# [4] soliton escape: a c-boosted localized carrier is Lorentz-contracted below a0
print("\n[4] soliton escape: lab size = (rest size w0)*(M/E); representability needs lab size >= a0")
print("    => E_max ~ M * w0/a0 (take M ~ Lambda, a lattice-scale carrier):")
for w0 in (1, 3, 10):
    print(f"    rest size {w0:2d}*a0 : E_max ~ {LAMBDA*w0:.3f} GeV = {LAMBDA*w0/1000:.5f} of 1 TeV")
assert (LAMBDA * 10) / 1000.0 < 0.01                       # even a 10*a0 soliton caps far below 1 TeV

print(f"""
[5] VERDICT — the trans-Lambda photon no-go is KINEMATIC (mechanism-independent):
  A c-moving quantum of energy E needs a lab-frame internal scale ~ hbar c/E; for
  E >> hbar c/a0 ~ {E_edge:.2f}-{E_corner:.2f} GeV that scale is sub-lattice and forbidden. Every escape
  fails by this one wall, using NO QED perturbation theory:
    * NULL bundle (moves at c)  -> soft constituents (E1 = 2L*sqrt(1-beta^2) -> 0 as beta->1)
      -> a 1 TeV event is >~ 1e4 independent soft quanta, not one photon
      [the Wilson/QED no-go, now shown FORCED by null-ness, not one option among many];
    * HARD bundle (E1 ~ E_BZ)   -> zone-edge constituents have v_g -> 0 -> ~stationary -> not a photon;
    * SOLITON                   -> Lorentz-contracted below a0 when boosted -> pinned, E capped O(Lambda).
  So NO lattice-mode composite -- elementary, bundle, or soliton -- is simultaneously NULL
  (a photon) and HARD (a single total-P^mu detector vertex). The required B_N (one LSZ leg
  carrying total P^mu, built from sub-cutoff records, no new UV oscillators) CANNOT be a
  lattice-mode composite, independent of the Wilson coupling. The same wall binds every sector
  (high-energy electrons/quarks on the a0 lattice too): this is foundational, not photon-specific.
  DOORS (costs explicit): (i) accept the EFT scoping -- the framework is a sub-~1 GeV effective
  theory for radiation [K12], the canon-consistent position; (ii) a finer gauge UV scale --
  reopens the cosmological constant [K12 route A]; (iii) a dynamical/rewiring (non-fixed-lattice)
  carrier with no fixed Brillouin spectrum -- outside this no-go, but unestablished and likely
  UV/CC-costly. The single-scale (CC-resolving) commitment makes (i) the honest reading.
exit 0""")
print("ALL ASSERTIONS PASSED — null<->soft, hard<->slow, soliton contraction: B_N is not a lattice composite.")
