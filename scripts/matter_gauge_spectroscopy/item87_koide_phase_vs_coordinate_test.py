#!/usr/bin/env python3
r"""Item 87 -- reframing test: is the Koide delta a removable COORDINATE of a real
spectrum, or a forced PHASE?  The proposal (to dissolve the "why d/N rad" residual)
was that the charged-lepton masses are CP-trivial reals, so delta is just the
coordinate you get when you write them in Koide form -- no phase to derive.

This script tests that proposal against the framework's OWN structure: the three
nu_R pseudocodewords form a Z3-CIRCULANT ring on the generation register
(ANCHOR sec 5.8, line 669).  The question is whether a Z3-circulant can give the
real, NON-DEGENERATE lepton hierarchy WITHOUT a complex phase.

Honest outcome (asserted below): it CANNOT.  A real circulant is either
degenerate (symmetric, c1=c2) or has complex/unphysical eigenvalues (c1!=c2).
The non-degenerate real hierarchy forces a HERMITIAN circulant with a complex
off-diagonal phase delta.  So under the framework's Z3 ring, delta is a FORCED,
mass-setting phase -- the coordinate reframing FAILS.  What survives: delta is
CP-TRIVIAL (the circulant eigenvectors are the delta-independent DFT modes, so
delta sets only masses, not mixing), hence it is a mass-SHAPE phase, distinct
from the CP-violating portal phase Phi.  "Why delta=2/9" stands as a real
question; the arc/winding analysis was not mis-posed.
"""
import numpy as np

R2 = np.sqrt(2)
w = np.exp(2j * np.pi / 3)
me, mmu, mtau = 0.51099895, 105.6583755, 1776.86  # MeV, PDG


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def circulant(c0, c1, c2):
    return np.array([[c0, c1, c2], [c2, c0, c1], [c1, c2, c0]])


def main():
    print("ITEM 87 -- is delta a removable coordinate, or a forced phase?")
    print("=" * 78)

    print("\n[1] A real spectrum is reproducible WITHOUT a phase -- if you DROP Z3-circulant")
    M_general_real = np.diag([me, mmu, mtau])  # real, phaseless, but NOT circulant
    ev = np.sort(np.linalg.eigvalsh(M_general_real))
    check(np.allclose(ev, sorted([me, mmu, mtau])), "a general real (non-circulant) matrix gives the 3 masses, phaseless")

    print("\n[2] But a REAL Z3-CIRCULANT cannot give the real non-degenerate hierarchy")
    ev_sym = np.linalg.eigvals(circulant(1.0, 0.3, 0.3))            # real symmetric circulant
    degen = np.min([abs(ev_sym[i] - ev_sym[j]) for i in range(3) for j in range(i + 1, 3)])
    print(f"    real symmetric circulant (c1=c2): eigenvalues {np.round(ev_sym.real,4)} -> degenerate pair")
    check(degen < 1e-9, "real symmetric circulant FORCES a degenerate pair (no 3-way hierarchy)")
    ev_asym = np.linalg.eigvals(circulant(1.0, 0.3, 0.1))           # real non-symmetric circulant
    print(f"    real non-symmetric circulant (c1!=c2): max |Im(eig)| = {np.max(np.abs(ev_asym.imag)):.4f}")
    check(np.max(np.abs(ev_asym.imag)) > 1e-3, "real non-symmetric circulant gives COMPLEX (unphysical) masses")

    print("\n[3] The hierarchy forces a HERMITIAN circulant with a complex phase delta")
    delta = 2.0 / 9.0
    c1 = (R2 / 2.0) * np.exp(1j * delta)                            # |c1|=R2/2 -> 2|c1|=sqrt2
    M_herm = circulant(1.0, c1, np.conj(c1))                        # Hermitian circulant
    check(np.max(np.abs(M_herm - M_herm.conj().T)) < 1e-12, "Hermitian circulant is a valid (real-eigenvalue) operator")
    sqrt_factors = np.sort(np.linalg.eigvalsh(M_herm))             # = 1 + sqrt2 cos(delta+2pi n/3)
    masses = sqrt_factors ** 2
    mass_ratios = masses / masses.min()
    pdg_ratios = np.array(sorted([me, mmu, mtau])) / me
    print(f"    eigenvalues (sqrt m/mu) = {np.round(sqrt_factors,4)}  (real, non-degenerate)")
    print(f"    mass ratios  pred {np.round(mass_ratios,1)}  vs PDG {np.round(pdg_ratios,1)}")
    check(np.max(np.abs(np.linalg.eigvalsh(M_herm).imag if np.iscomplexobj(np.linalg.eigvalsh(M_herm)) else 0)) < 1e-12,
          "eigenvalues are real")
    check(abs(mass_ratios[1] / pdg_ratios[1] - 1) < 0.01 and abs(mass_ratios[2] / pdg_ratios[2] - 1) < 0.01,
          "delta=2/9 Hermitian circulant reproduces the lepton mass ratios (~1%)")

    print("\n[4] ...but that phase is CP-TRIVIAL: it sets masses, not mixing")
    # Circulant eigenvectors are the DFT modes for ANY c0,c1,c2 -> independent of delta.
    dft = np.array([[1, 1, 1], [1, w, w ** 2], [1, w ** 2, w]]) / np.sqrt(3)
    _, vecs = np.linalg.eigh(M_herm)
    # each eigenvector must be (a phase times) a DFT column -> |DFT^dagger vecs| is a permutation of identity
    overlap = np.abs(dft.conj().T @ vecs)
    is_perm = np.allclose(np.sort(overlap, axis=0)[-1], 1.0, atol=1e-9)
    check(is_perm, "eigenvectors are the delta-INDEPENDENT DFT modes -> delta sets only masses, not mixing (CP-trivial)")

    print(
        """
[5] VERDICT (honest -- the reframing FAILS)
    The coordinate reframing does NOT survive the framework's own Z3 structure.
    A real Z3-circulant is either degenerate (c1=c2) or complex/unphysical
    (c1!=c2); the real, non-degenerate lepton hierarchy therefore FORCES a
    Hermitian circulant with a complex off-diagonal phase delta.  So delta is a
    forced, mass-SETTING phase, not a removable coordinate -- "why delta=2/9"
    stands as a genuine question, and the arc/winding analysis was not mis-posed.

    What the test DID clarify (the surviving, useful half):
      * delta is CP-TRIVIAL: the circulant eigenvectors are the delta-independent
        DFT modes, so delta fixes only the real masses, carrying no Jarlskog /
        CP content.  It is a mass-SHAPE phase.
      * It is therefore a DIFFERENT object from the CP-VIOLATING portal phase Phi
        (item87 K3 Majorana portal).  The canon conflates them by testing
        Phi = delta_nu; this test shows they are distinct in kind (CP-trivial
        mass-shape vs CP-odd baryogenesis phase).
    Net: the residual is NOT dissolved.  delta=2/9 is a real, forced, CP-trivial
    mass-shape phase whose value is still fit-selected; the documented Feshbach
    (real, generation-diagonal, degenerate) does not supply it; and the genuine
    CP phase Phi remains the separate open portal.  Item 87 stays open; I was
    wrong that the reframing would dissolve it.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
