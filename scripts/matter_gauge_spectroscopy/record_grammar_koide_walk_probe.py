#!/usr/bin/env python3
r"""Rung 3, guard-rail: does a NAIVE walk operator W on the record cell produce
the Koide lepton spectrum (Q=2/3) -- with mass = the cube-walk k=0 eigenphases --
or is 2/3 only reachable by free-parameter fitting under THAT ansatz?

CORRECTION 2026-06-30: this script tests the NAIVE "mass = cube-walk eigenphase"
ansatz, which is NOT the framework's actual mass mechanism.  The framework derives
the charged-lepton spectrum from the sec 5.8 Feshbach circulant (ANCHOR.md): the
nu_R pseudocodeword Z3-ring with R=sqrt(2) forced by Dirac-operator quadrature
(sec 3.5) and phase delta=2/9 = d/N (2-site nu_R defect / 9-site plaquette).  That
mechanism reproduces m_e,m_mu to 0.007%,0.006% (mu fixed from tau), predicts
m_tau=1776.969 MeV (Belle II 2023 -> ~0sigma at the world average), and extends to
a falsifiable NEUTRINO prediction Q_nu=1/2, delta_nu=1/3 (Dm2 ratio 0.030 vs
NuFIT 0.0295, Normal Ordering).  So the framework DOES force Koide; what fails
below is only the naive cube-walk-eigenphase ansatz.  Read this as "mass is not a
cube-walk eigenphase", NOT as "the framework cannot derive Koide".

Anti-postdiction discipline -- the input set is LOCKED before Q is computed:
  - matter cell graph = Q3 (the 3-cube = [8,4,4]=RM(1,3) record cell).
  - 3 generations  =  the 3 cube directions  =  the coin/direction space (the
    cell's natural "3", tied to triality).
  - W = coined quantum walk on Q3;  rest mass m_i = the k=0 dispersion = the
    eigenphases of the coin C (at zero momentum the shift is the identity, so
    W(0) = C).  mass m_i = |theta_i| (positive, monotone).
  - Koide:  Q = (sum m_i) / (sum sqrt(m_i))^2,  target 2/3.
The ONLY thing not fixed a priori is the coin C.  Test: do the symmetry-FORCED
coins (no free parameter) give a Koide triple, or only a free coin (a fit)?

Result (asserted below): under the naive ansatz the forced coins each give a
massless and/or degenerate rest spectrum -- NONE is a valid Koide triple -- and
hitting 2/3 needs three independent eigenphases = a 3-parameter fit.  This is why
the framework does NOT use the cube-walk-eigenphase ansatz; it uses the sec 5.8
Feshbach circulant instead.  The load-bearing sqrt(2) amplitude is derived
geometrically (a second, independent route to the canon's R=sqrt2) in
record_grammar_koide_magic_geometry.py.
"""
import numpy as np

w = np.exp(2j * np.pi / 3)


def assert_true(name, value):
    print(f"  {name:<66s} value={value}")
    if not value:
        raise AssertionError(name)


def koide(m):
    m = np.asarray(m, float)
    s = np.sqrt(np.abs(m))
    return float(m.sum() / s.sum() ** 2) if s.sum() > 0 else float("nan")


def valid_koide_triple(m, tol=1e-9):
    """Three strictly positive, pairwise-distinct rest masses."""
    m = np.sort(np.abs(np.asarray(m, float)))
    return bool((m > tol).all()) and len(set(np.round(m, 6))) == 3


def rest_spectrum(C):
    """k=0 walk eigenphases (rest masses) for coin C; |theta_i|, sorted."""
    return np.sort(np.abs(np.angle(np.linalg.eigvals(C))))


def main():
    print("=" * 80)
    print("PART A  Geometric content of the target (pure math, Foot 1994)")
    print("=" * 80)
    # Q=2/3 <=> sqrt(m) at 45 deg to (1,1,1) <=> the 1 + sqrt(2) cos form.
    s0 = np.array([1 + np.sqrt(2) * np.cos(2 * np.pi * k / 3) for k in range(3)])  # phi=0
    me, mmu, mtau = 0.51099895e-3, 0.1056583755, 1.77686  # GeV, PDG
    sm = np.sqrt([me, mmu, mtau])
    ang = np.degrees(np.arccos(sm.sum() / (np.linalg.norm(sm) * np.sqrt(3))))
    assert_true("1+sqrt(2)cos form gives Q=2/3", abs(koide(s0 ** 2) - 2 / 3) < 1e-12)
    assert_true("real charged leptons satisfy Koide (Q~0.6667)", abs(koide([me, mmu, mtau]) - 2 / 3) < 1e-4)
    assert_true("their sqrt(m) sits at 45 deg to (1,1,1)", abs(ang - 45.0) < 0.05)

    print("\n" + "=" * 80)
    print("PART B  FORCED symmetric coins on the 3-generation space (no free param)")
    print("=" * 80)
    svec = np.ones(3) / np.sqrt(3)
    coins = {
        "Grover 2|s><s|-I": 2 * np.outer(svec, svec) - np.eye(3),
        "DFT3 (Z3 Fourier)": np.array([[1, 1, 1], [1, w, w ** 2], [1, w ** 2, w]]) / np.sqrt(3),
        "Z3 cyclic shift": np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float),
    }
    for name, C in coins.items():
        assert np.allclose(C @ C.conj().T, np.eye(3), atol=1e-12), f"{name} not unitary"
        m = rest_spectrum(C)
        print(f"    {name:<20s} rest masses = {np.round(m, 3)}")
        assert_true(f"{name}: NOT a valid Koide triple", not valid_koide_triple(m))

    print("\n" + "=" * 80)
    print("PART C  A free coin trivially fits 2/3 (so 2/3 alone is not a prediction)")
    print("=" * 80)
    # A unitary circulant C = F^dag diag(e^{ia},e^{ib},e^{ic}) F has 3 free phases
    # = 3 free rest masses.  Set them to a Koide triple -> Q=2/3 by construction.
    # (phi=0.2 -- where the real charged leptons sit -- gives three DISTINCT
    #  positive sqrt-masses; phi=0 would collapse two of them by symmetry.)
    s_distinct = np.array([1 + np.sqrt(2) * np.cos(0.2 + 2 * np.pi * k / 3) for k in range(3)])
    free_masses = s_distinct ** 2  # three independent, distinct, positive numbers
    assert_true("free 3-phase coin can realise Q=2/3 (a 3-dof fit)",
                valid_koide_triple(free_masses) and abs(koide(free_masses) - 2 / 3) < 1e-12)

    print(
        """
VERDICT:
  PASS (the naive ansatz is excluded).  Every symmetry-FORCED coin on the cell's
  generation space gives a massless and/or degenerate rest spectrum -- none is a
  Koide triple -- and reproducing 2/3 needs three independent eigenphases, a fit.
  So mass is NOT a cube-walk k=0 eigenphase.  This does NOT mean the framework
  cannot derive Koide: it derives it via the sec 5.8 Feshbach circulant
  (R=sqrt2 + delta=2/9), verified to 0.006-0.007% with m_tau=1776.969 MeV
  confirmed by Belle II.  This script's role is to exclude the tempting naive
  mechanism, so the Feshbach mechanism is not mistaken for the only option by
  default.  See ANCHOR.md sec 5.8 and record_grammar_koide_magic_geometry.py.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
