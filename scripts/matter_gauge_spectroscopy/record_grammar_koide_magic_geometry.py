#!/usr/bin/env python3
r"""Rung 3, the attack: is the Koide amplitude sqrt(2) FORCED by the record cell's
magic-state (octahedron) geometry, or is it free?

The companion guard-rail (record_grammar_koide_walk_probe.py) established the
negative: a forced symmetric walk coin gives a degenerate/massless spectrum, and
reproducing Koide Q=2/3 with three free eigenphases is a fit.  The geometric
content of the target isolates *one* number:

    For a Z3-symmetric spectrum  sqrt(m_i) = M * (1 + b*cos(phi + 2*pi*i/3)),
        Q = (sum m_i)/(sum sqrt(m_i))^2 = (1 + b^2/2)/3,
    which is independent of the scale M and the phase phi.  Hence
        Q = 2/3   <=>   b^2 = 2   <=>   b = sqrt(2).
The whole of Koide reduces to the single amplitude ratio b = sqrt(2)
(Foot 1994: sqrt(m) at 45 deg to the democratic axis).  So: does the framework
FORCE b = sqrt(2)?

The framework's matter substrate is the single-qubit record cell, whose state
space is the Bloch ball with the stabilizer OCTAHEDRON (the L1 ball |x|+|y|+|z|<=1)
inscribed.  Magic = distance outside the octahedron; the L1 norm of a pure state
measures it.  Three Koide-BLIND geometric extrema (computed below, not recalled):

    on-axis stabilizer extent           = 1      (octahedron vertex, e.g. |0>)
    max magic L1 in a 2-plane           = sqrt(2)  (equatorial T-state, face diag)
    max magic L1 in full 3-space        = sqrt(3)  ((1,1,1)/sqrt3 T-state, body diag)

The Z3 generation symmetry is a 120-degree ROTATION; a rotation fixes one axis
(the democratic baseline) and acts in the orthogonal 2-PLANE.  So the generation-
distinguishing modulation lives in that plane, where the maximal magic overhang
is exactly sqrt(2) -- and (1 + b^2/2)/3 then forces Q = 2/3.  The contrast cases
(b=1 stabilizer, b=sqrt(3) full-3D magic) give 1/2 and 5/6, NOT 2/3, so the planar
selection is doing real work.

HONEST STATUS (asserted in the verdict): the number sqrt(2) is forced by the
octahedron's in-plane magic maximum.  What remains ASSUMED (motivated, not
derived) is the identification  sqrt(m)_i = stabilizer-axis baseline + in-plane-
magic modulation  (mass-splitting = magic resource), and that mass enters at the
sqrt level.  This converts "why 2/3?" from a 3-parameter fit into ONE structural
identification with the magic number forced.  It is NOT a complete derivation, and
NOT yet a new-ratio prediction: the phase phi (hence the individual m_e:m_mu:m_tau)
is still free.
"""
import numpy as np


def assert_true(name, value):
    print(f"  {name:<66s} value={value}")
    if not value:
        raise AssertionError(name)


def koide(m):
    m = np.asarray(m, float)
    s = np.sqrt(np.abs(m))
    return float(m.sum() / s.sum() ** 2)


def Q_of_b(b, phi=0.0):
    """Koide ratio for the Z3 spectrum sqrt(m)=1+b cos(phi+2pi i/3) (b<=2 keeps s>=0)."""
    s = np.array([1 + b * np.cos(phi + 2 * np.pi * k / 3) for k in range(3)])
    return koide(s ** 2)


def max_l1_unit(dims, n=2001):
    """Max of |v|_1 over unit vectors in R^dims, by deterministic grid (Koide-blind)."""
    best = 0.0
    if dims == 2:
        for t in np.linspace(0, 2 * np.pi, n):
            best = max(best, abs(np.cos(t)) + abs(np.sin(t)))
    else:
        for th in np.linspace(0, np.pi, n // 4):
            st, ct = np.sin(th), np.cos(th)
            for ph in np.linspace(0, 2 * np.pi, n // 2):
                best = max(best, abs(st * np.cos(ph)) + abs(st * np.sin(ph)) + abs(ct))
    return best


def main():
    print("=" * 80)
    print("PART A  Koide reduces to ONE amplitude:  Q = (1 + b^2/2)/3  =>  2/3 iff b=sqrt2")
    print("=" * 80)
    for b in (1.0, np.sqrt(2), np.sqrt(3)):
        q_closed = (1 + b ** 2 / 2) / 3
        assert_true(f"Q(b={b:.4f}) closed form == numeric spectrum", abs(Q_of_b(b) - q_closed) < 1e-12)
    assert_true("Koide 2/3 <=> b = sqrt(2) exactly", abs(Q_of_b(np.sqrt(2)) - 2 / 3) < 1e-12)

    print("\n" + "=" * 80)
    print("PART B  Octahedron geometry, computed Koide-blind (no Koide input)")
    print("=" * 80)
    axis_stab = 1.0  # stabilizer octahedron vertex |0>=(0,0,1): |v|_1 = 1
    plane_magic = max_l1_unit(2)
    full_magic = max_l1_unit(3)
    print(f"    on-axis stabilizer extent  = {axis_stab:.6f}")
    print(f"    max magic L1 in a 2-plane  = {plane_magic:.6f}   (closed form sqrt2={np.sqrt(2):.6f})")
    print(f"    max magic L1 in full 3D    = {full_magic:.6f}   (closed form sqrt3={np.sqrt(3):.6f})")
    assert_true("in-plane magic maximum == sqrt(2)", abs(plane_magic - np.sqrt(2)) < 2e-3)
    assert_true("full-3D magic maximum  == sqrt(3)", abs(full_magic - np.sqrt(3)) < 2e-3)
    # the equatorial T-state realising the in-plane maximum:
    t_state = np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 0.0])
    assert_true("equatorial T-state (1/sqrt2,1/sqrt2,0) is a unit Bloch vector",
                abs(np.linalg.norm(t_state) - 1) < 1e-12)
    assert_true("...and its L1 norm is sqrt(2) (the in-plane magic overhang)",
                abs(np.abs(t_state).sum() - np.sqrt(2)) < 1e-12)

    print("\n" + "=" * 80)
    print("PART C  Planar Z3 modulation selects b=sqrt(2) -> Koide; the alternatives do not")
    print("=" * 80)
    cases = {
        "stabilizer baseline only (b=1)": (1.0, 1 / 2),
        "IN-PLANE magic (b=sqrt2)  [Z3 is a planar rotation]": (np.sqrt(2), 2 / 3),
        "full-3D magic (b=sqrt3)": (np.sqrt(3), 5 / 6),
    }
    for name, (b, q_exp) in cases.items():
        q = Q_of_b(b)
        print(f"    {name:<52s} Q = {q:.4f}")
        assert_true(f"  -> Q == {q_exp:.4f}", abs(q - q_exp) < 1e-12)
    assert_true("ONLY the in-plane magic amplitude gives Koide 2/3",
                abs(Q_of_b(plane_magic) - 2 / 3) < 2e-3
                and abs(Q_of_b(1.0) - 2 / 3) > 0.1
                and abs(Q_of_b(full_magic) - 2 / 3) > 0.1)

    print(
        """
VERDICT:
  PASS.  The Koide amplitude b=sqrt(2) is FORCED by the record cell's octahedron
  geometry: it is the maximal magic (L1 overhang) of a qubit state confined to the
  Z3 generation-rotation plane -- a Koide-blind extremum -- and (1+b^2/2)/3 then
  pins Q=2/3.  The full-3D magic (sqrt3 -> 5/6) and the bare stabilizer
  (1 -> 1/2) miss, so the planar (Z3-rotation) selection of sqrt(2) is load-bearing.

  What remains ASSUMED (isolated, motivated, NOT derived):
    (i) sqrt(m)_i = stabilizer-axis baseline + in-plane-magic modulation
        (i.e. generation mass-splitting = the magic resource), at the sqrt level.
  This reduces "why 2/3?" from a 3-parameter fit to ONE structural identification
  with the magic number forced.  It is NOT a complete first-principles derivation
  of Koide via THIS geometric route, and this route does not by itself fix the
  phase phi (the octahedron's 4-fold Pauli structure is incommensurate with the
  Z3 3-fold, so no special equatorial angle gives the empirical phase).

  RECONCILIATION with canon (ANCHOR.md sec 5.8, 2026-06-30): the framework already
  derives the FULL charged-lepton spectrum -- R=sqrt(2) from Dirac-operator
  quadrature AND the phase delta=2/9 from defect counting d/N (2-site nu_R defect
  / 9-site plaquette), NOT from an equatorial angle -- reproducing m_e,m_mu to
  0.007%,0.006% and predicting m_tau=1776.969 MeV (Belle II ~0sigma), plus a
  falsifiable NEUTRINO new-ratio prediction (Q_nu=1/2, delta_nu=1/3, Dm2 ratio
  0.030 vs NuFIT 0.0295, Normal Ordering).  So: (a) this magic-octahedron result
  is a SECOND, independent derivation of the same canonical R=sqrt(2); (b) the
  phase and the new-ratio prediction DO exist in the framework, via the defect-
  counting mechanism, not the octahedron-angle one.  Honest open gaps remain
  (canon's own tiering): the neutrino phase OPERATOR is an ansatz (phase-lift
  blocked), the base mass mu_l has a residual /3, the quark sector is partly open.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
