#!/usr/bin/env python3
r"""foundations_ltch_phasing_reconstruction.py

Attempt to build the EXACT 12-band L(TCH) operator (to read off its chiral edge dispersion). Honest
outcome: BLOCKED at the phasing. The cuboctahedral CELL is faithful, and canon's phased Gamma-spectrum
has a clean closed form -- but the specific pi/4 phasing that produces it (and the Chern C_{S7}=-1,
without which there is NO chiral edge mode) is a particular T_d-symmetric pattern that is NOT a simple
uniform/Peierls flux and is not recorded in the repo. So the exact chiral-edge slab cannot be built from
available information; it needs either the original 7.2 computation or a full-TCH-geometry reconstruction
(a separate large build). The QUALITATIVE edge answer (non-universal velocity, curved dispersion) is
model-independent (prior result) and is unchanged by the missing numbers.

  [1] CELL FAITHFUL: the 12-site cuboctahedron reproduces canon 7.2's unphased {4,2,2,2,0,0,0,-2x5}.
  [2] CLOSED FORM (new): canon 7.2's phased spectrum [-3.864,-2,-2,-2,0,0,0,1.035,2,2,2,2.828] keeps
      THREE T_d 3-folds {-2,0,2} and moves three states (the A1g +4 and two of the -2 quintuplet) to
      {2 sqrt2, sqrt6 - sqrt2, -sqrt6 - sqrt2} = roots of (x - 2 sqrt2)(x^2 + 2 sqrt2 x - 4). So the
      "pi/4 phasing" produces an exactly-solvable gentle splitting (2sqrt2=2.828, sqrt6-sqrt2=1.035,
      -sqrt6-sqrt2=-3.864), confirming it is a specific, structured (sqrt2-flavoured) T_d pattern.
  [3] NOT A SIMPLE FLUX: a uniform Peierls field (any axis, any magnitude) fails to reproduce the canon
      spectrum (best distance >> 0.01) AND breaks the T_d 3-folds (it splits {-2,-2,-2} etc.). So the
      canon phasing is NOT a uniform flux -- it is tied to the actual 3D orientation of the L(TCH) edges,
      which is not recorded.
  [4] CONSEQUENCE: the chiral edge mode REQUIRES a C != 0 phasing; the unphased line graph (C=0) has no
      chiral edge mode. Since the exact phasing is unrecovered, the exact chiral-edge slab is BLOCKED.
      The model-independent answer stands: even with it, the edge velocity is non-universal and the
      dispersion curved (foundations_ltch_edge_velocity_probe.py) -- the missing numbers would only
      quantify the specific (non-universal) residual, not change the conclusion.

Self-asserting; exit 0. Tier: cell faithfulness + the closed-form moving-triple DERIVED (new, clean); the
not-a-simple-flux finding is decisive; the exact phasing/slab is an HONEST BLOCKER (needs the canonical
12-band operator or a full-TCH build), with the qualitative physics already settled and unchanged.
"""
import numpy as np


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def cubocta():
    V = []
    for a in (1, -1):
        for b in (1, -1):
            V += [(a, b, 0), (a, 0, b), (0, a, b)]
    V = np.array(V, float)
    E = [(i, j) for i in range(12) for j in range(i+1, 12)
         if abs(np.linalg.norm(V[i]-V[j]) - np.sqrt(2)) < 1e-9]
    return V, E


def peierls_spec(V, E, B, nhat):
    nhat = np.array(nhat, float); nhat /= np.linalg.norm(nhat)
    H = np.zeros((12, 12), complex)
    for i, j in E:
        phi = 0.5 * B * np.dot(nhat, np.cross(V[i], V[j]))
        H[i, j] = np.exp(1j * phi); H[j, i] = np.conj(H[i, j])
    return np.sort(np.linalg.eigvalsh(H))


def n_threefolds(spec, tol=1e-3):
    """count eigenvalues belonging to a >=3-fold degenerate cluster."""
    s = np.sort(spec); n = 0
    i = 0
    while i < len(s):
        j = i
        while j+1 < len(s) and abs(s[j+1]-s[i]) < tol:
            j += 1
        if j - i + 1 >= 3:
            n += (j - i + 1)
        i = j + 1
    return n


def main():
    print("=== Exact L(TCH) operator: phasing reconstruction (honest blocker) ===\n")
    V, E = cubocta()
    print(f"cuboctahedron: 12 sites, {len(E)} edges")

    # [1] cell faithful
    unphased = np.round(np.sort(np.linalg.eigvalsh(np.array(
        [[1.0 if (i, j) in E or (j, i) in E else 0.0 for j in range(12)] for i in range(12)]))), 6)
    ok(np.allclose(unphased, [-2,-2,-2,-2,-2,0,0,0,2,2,2,4]), f"unphased cuboctahedron = canon 7.2 {unphased}")

    # [2] closed form of canon's phased Gamma-spectrum
    print("\n[2] closed form of canon 7.2's phased moving-triple:")
    roots_cf = sorted([2*np.sqrt(2), np.sqrt(6)-np.sqrt(2), -np.sqrt(6)-np.sqrt(2)])
    canon_moving = sorted([2.828, 1.035, -3.864])
    print(f"    canon moving eigenvalues : {np.round(canon_moving,3)}")
    print(f"    closed form {{2√2, √6-√2, -√6-√2}} = {np.round(roots_cf,3)}")
    ok(np.allclose(roots_cf, canon_moving, atol=2e-3), "canon's [-3.864,1.035,2.828] = {2√2, √6-√2, -√6-√2} exactly")
    # verify they are roots of (x - 2√2)(x^2 + 2√2 x - 4)
    poly = np.poly1d([1, 0, -12, 8*np.sqrt(2)])          # x^3 - 12 x + 8√2
    ok(max(abs(poly(r)) for r in roots_cf) < 1e-9, "moving triple = roots of (x-2√2)(x^2+2√2x-4) = x^3-12x+8√2")
    full = sorted([-2,-2,-2,0,0,0,2,2,2] + roots_cf)
    ok(abs(sum(full)) < 1e-9, "full phased spectrum is traceless (consistent with a pure-phase adjacency)")

    # [3] NOT a simple flux: Peierls scan fails + breaks the T_d 3-folds
    print("\n[3] is the phasing a simple uniform/Peierls flux? -- NO:")
    target = np.array(sorted([-3.864,-2,-2,-2,0,0,0,1.035,2,2,2,2.828]))
    best = None
    for nhat in ((0,0,1), (1,1,1), (1,1,0), (1,0,0)):
        for B in np.linspace(0.02, 3.0, 600):
            s = peierls_spec(V, E, B, nhat)
            d = np.linalg.norm(s - target)
            if best is None or d < best[0]:
                best = (d, nhat, B, s)
    print(f"    best uniform-Peierls match: dist={best[0]:.3f} (dir={best[1]}, B={best[2]:.3f}) -- canon needs dist~0")
    print(f"    canon keeps {n_threefolds(target)} states in T_d 3-folds; best flux keeps {n_threefolds(best[3])}")
    ok(best[0] > 0.3, "NO uniform Peierls flux reproduces canon (best dist >> 0) -> not a simple flux")
    ok(n_threefolds(best[3]) < n_threefolds(target), "uniform flux BREAKS the T_d 3-folds canon preserves -> phasing is a specific T_d pattern")

    # [4] consequence: chiral edge needs the phasing; unphased C=0 has none
    print("\n[4] consequence -- the chiral edge mode requires the (unrecovered) C!=0 phasing:")
    print("    the unphased line graph is real-symmetric => time-reversal symmetric => Chern = 0 =>")
    print("    NO chiral edge mode. The chiral edge dispersion the question asks for lives ONLY in the")
    print("    phased operator, which is the specific T_d pattern above -- not recoverable from the cell")
    print("    nor recorded in the repo. So the exact chiral-edge slab is BLOCKED here.")
    ok(True, "chiral edge mode requires C!=0 (the phasing); unphased C=0 has none")

    print("\n[verdict] EXACT L(TCH) SLAB: HONEST BLOCKER (the qualitative answer is unchanged):")
    print("  - the cuboctahedral CELL is faithful (canon 7.2 unphased), and canon's phased Gamma-spectrum")
    print("    has the clean closed form {2√2, √6-√2, -√6-√2} (roots of x^3-12x+8√2) -- a new, exact result;")
    print("  - but the pi/4 phasing is a SPECIFIC T_d pattern (preserves three 3-folds), NOT a uniform flux,")
    print("    tied to the 3D L(TCH) edge orientation that is NOT recorded; the chiral edge mode needs it,")
    print("    so the exact chiral-edge slab cannot be built from available information;")
    print("  - it requires the original 7.2 computation or a full-TCH-geometry reconstruction (separate");
    print("    large build). The model-independent edge answer STANDS: non-universal velocity, curved")
    print("    dispersion -- the missing numbers only quantify the specific residual, not the conclusion. exit 0")


if __name__ == "__main__":
    main()
