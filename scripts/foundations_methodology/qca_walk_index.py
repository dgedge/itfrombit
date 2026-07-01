#!/usr/bin/env python
r"""qca_walk_index.py -- QCA scoping probe for the chiral-gauge frontier.

CONTEXT. The chiral-lattice-gauge frontier (DRIFT chiral cluster) reduces to the field-wide
non-abelian SU(2)_L overlap-measure construction, and the framework proved (the homotopy theorem,
DRIFT "FINAL CLOSURE") that the substrate gives ZERO leverage there. Every route tried realises
chirality via sgn(gamma5 D) -- the overlap/admissibility setup the homotopy theorem covers. The one
frame NOT in the exhausted list: the walk W = S.C is literally a Quantum Cellular Automaton, and the
post-2020 "chiral fermion as the boundary of a (d+1)D anomalous QCA via a NON-ONSITE symmetry" line
evades Nielsen-Ninomiya through a different door (N-N assumes an ONSITE symmetry). This probe asks the
sharpest tractable QCA question for the actual walk:

   Does the framework's walk W = S.C carry a nonzero 1D chiral (GNVW) index?

A free 1D QCA's GNVW index = winding number of det W(k) around the Brillouin zone (= net chiral
information transport per step; = sum of the quasi-energy band windings). Nonzero index = a genuine
anomalous/chiral QCA (the N-N loophole). Zero index = trivial = the doubling is mandatory.

RESULT (honest scoping negative). For the framework's CHIRAL-BALANCED shift S = diag(e^{ik}, e^{-ik})
(one right-mover + one left-mover, the documented walk structure) the index is ZERO for ANY coin C
(det S = 1, so winding(det W) = winding(det C) = 0; the two modes wind +1 and -1, summing to 0) --
i.e. the strict-1D QCA framing reproduces the doubling (consistent with item 97's 4 right-handed
doublers), NOT a loophole. A nonzero index (single chiral mode) requires an UNBALANCED shift
diag(e^{ik}, 1) -- a valid local QCA, but NOT gamma5-symmetric and NOT the framework's walk; it is
the "anomalous QCA" that can only live as the boundary of a higher-D bulk.

CONCLUSION. The QCA escape is real in principle but (i) is genuinely (d+1)D-bulk / non-onsite, not a
property of the framework's 1D walk (index 0 here), and (ii) even granted a single chiral fermion the
BINDING wall is the non-abelian gauge MEASURE, on which the homotopy theorem already proves zero
substrate leverage. So the QCA route is a long-shot that does not obviously crack the actual wall.
Recorded as an open-frontier scoping item, not a closure. See qca_chiral_index_scope.md.
"""
import numpy as np


def winding(vals):
    """Integer winding number of a closed loop of nonzero complex samples vals[0..N-1] (then ->vals[0])."""
    v = vals / np.abs(vals)
    inc = np.angle(v[1:] / v[:-1]).sum() + np.angle(v[0] / v[-1])
    return inc / (2 * np.pi)


def winding_det(Wk, n=2048):
    """The 1D free-QCA GNVW index = winding number of det W(k) around the BZ."""
    ks = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    dets = np.array([np.linalg.det(Wk(k)) for k in ks])
    return winding(dets)


def coin(theta, phi=0.0):
    """A representative 2x2 unitary coin (rotation by theta with overall phase phi). det = e^{2i phi}."""
    c, s = np.cos(theta), np.sin(theta)
    return np.exp(1j * phi) * np.array([[c, -s], [s, c]], dtype=complex)


def chiral_shift(k):
    """Framework-style BALANCED chiral shift: up = right-mover (e^{ik}), down = left-mover (e^{-ik})."""
    return np.diag([np.exp(1j * k), np.exp(-1j * k)])


def walk(k, C):
    return chiral_shift(k) @ C


def anomalous_shift(k):
    """UNBALANCED single-mode shift diag(e^{ik}, 1): a valid local QCA with index 1 -- NOT gamma5-symmetric, NOT the walk."""
    return np.diag([np.exp(1j * k), 1.0 + 0j])


def main():
    print("=== QCA chiral (GNVW) index of the framework's walk W = S.C ===")
    print("    index := winding number of det W(k) around the BZ (= sum of quasi-energy band windings)")
    print()
    print("  framework walk  S=diag(e^{ik},e^{-ik}) . C(theta,phi)  -> expect index 0 (coin-independent):")
    idxs = []
    for th in (0.0, 0.3, 0.7, 1.1, np.pi / 4, 1.5):
        idx = winding_det(lambda k: walk(k, coin(th, phi=0.4)))
        idxs.append(idx)
        print(f"    coin theta={th:.3f}:  index = {idx:+.4f}  -> {int(round(idx))}")

    ks = np.linspace(0.0, 2 * np.pi, 2048, endpoint=False)
    w_up = winding(np.exp(1j * ks))
    w_dn = winding(np.exp(-1j * ks))
    print(f"\n  chiral-shift decoupled modes: right-mover winding {w_up:+.0f}, left-mover winding {w_dn:+.0f}"
          f"  (sum {w_up + w_dn:+.0f} = the doubling)")

    idx_anom = winding_det(anomalous_shift)
    print(f"\n  anomalous single-mode QCA  diag(e^{{ik}},1):  index = {idx_anom:+.4f}  -> {int(round(idx_anom))}"
          f"  (a valid local QCA, NOT gamma5-symmetric, NOT the walk)")

    print("\n[verdict] The framework's 1D walk has chiral (GNVW) index 0 for every coin -> the doubling is")
    print("  mandatory at the strict-1D QCA level (consistent with item 97's 4 right-handed doublers); the")
    print("  QCA framing offers NO 1D loophole. A nonzero index needs the unbalanced single-mode shift,")
    print("  which is a valid local QCA but NOT gamma5-symmetric and NOT the walk -- the genuine QCA escape")
    print("  is therefore the (d+1)D anomalous-BULK / non-onsite-symmetry construction, not a property of")
    print("  the 1D walk. And even granted a single chiral fermion, the BINDING wall is the non-abelian")
    print("  gauge MEASURE, on which the homotopy theorem proves zero substrate leverage.")
    print("  => QCA route = long-shot open frontier; does NOT obviously crack the actual wall.")

    assert all(abs(i - round(i)) < 1e-3 and round(i) == 0 for i in idxs), \
        "framework walk must have index 0 (the doubling), coin-independent"
    assert abs(idx_anom - 1) < 1e-3, "the unbalanced single-mode QCA must have index 1 (the contrast case)"
    assert abs(w_up - 1) < 1e-2 and abs(w_dn + 1) < 1e-2, "decoupled chiral-shift modes must wind +1 / -1"
    print("\nGATES PASSED -- index 0 for the walk (doubling); index 1 for the anomalous single-mode QCA;")
    print("decoupled modes wind +/-1. The 1D QCA framing gives no loophole; the escape is higher-D and")
    print("still faces the gauge-measure wall. Recorded as an open-frontier scoping item.")
    print("exit 0")


if __name__ == "__main__":
    main()
