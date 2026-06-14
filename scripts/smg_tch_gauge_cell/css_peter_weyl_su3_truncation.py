#!/usr/bin/env python3
"""
First Peter-Weyl SU(3) plaquette truncation with explicit Gauss intertwiners.

Target requested:
  Try the SU(3) Peter-Weyl truncation with link irreps

      (0,0) = 1,   (1,0) = 3,   (0,1) = 3bar

  and explicit Gauss intertwiners.

What this script actually does:
  * One oriented plaquette with the same convention as the qutrit/S3 proxies:
      e0: 0 -> 1, e1: 1 -> 2, e2: 3 -> 2, e3: 0 -> 3.
  * Static bipartite matter colour at vertices:
      v0=3, v1=3bar, v2=3, v3=3bar.
    Four fundamentals alone have no invariant states in this tiny two-link
    truncation; this bipartite sector is the smallest nonzero matter-coupled
    sector and is the colour analogue of the bipartite CSS construction.
  * Explicit local SU(3) invariant tensors:
      delta : 3 x 3bar -> 1,
      epsilon : 3 x 3 x 3 -> 1,
      epsilonbar : 3bar x 3bar x 3bar -> 1.
  * Builds the physical spin-network basis from those four vertex intertwiners.
  * Builds a truncated Wilson plaquette operator using explicit Clebsch maps for
    multiplication by a fundamental matrix element:

      1 x 3 -> 3,
      3 x 3 -> 3bar     (antisymmetric epsilon channel; 6 is truncated out),
      3bar x 3 -> 1     (singlet channel; 8 is truncated out).

  * Diagonalizes H = g^2 sum_l C2(r_l) - beta/2 (W_p + W_p^dag)
    on the gauge-invariant physical subspace for beta <= 1.
  * Then adds a minimal gauge-invariant "TCH lock" inside the physical
    intertwiner basis, coupling the two low electric-flux spin networks, to
    quantify what extra item-13 operator would have to do.

Scope and caveat:
  This is the first non-centre Peter-Weyl SU(3) calculation in the chain, but it
  is still a severe truncation. Because 6, 6bar, 8, ... are projected out, the
  Wilson operator is a projected operator, not the full SU(3) Kogut-Susskind
  plaquette. In fact, the minimal Wilson-only truncation below does collapse
  toward beta=0; the lock scan is included to make the missing TCH ingredient
  quantitatively explicit rather than hiding it.

  This remains a gauge-axis calculation. It does not contain the W!=chi mirror
  Fock modes or the SMG interaction. The matter-axis local mirror gap is tested
  separately in css_mirror_fock_smg_gap.py.

numpy; self-asserting.
"""

from collections import defaultdict
import itertools

import numpy as np


REPS = ["1", "3", "3b"]
MATTER_REPS = ["3", "3b", "3", "3b"]
EDGES = [
    (0, 1),
    (1, 2),
    (3, 2),
    (0, 3),
]
ORIENTATIONS = [1, 1, -1, -1]
_BASIS_CACHE = None
_U_CACHE = {}
_UDAG_CACHE = {}
_WILSON_CACHE = None


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def dual(rep):
    return {"1": "1", "3": "3b", "3b": "3"}[rep]


def dim(rep):
    return 1 if rep == "1" else 3


def casimir(rep):
    return 0.0 if rep == "1" else 4 / 3


def eps(i, j, k):
    if len({i, j, k}) < 3:
        return 0.0
    perm = [i, j, k]
    inversions = 0
    for a in range(3):
        for b in range(a + 1, 3):
            inversions += int(perm[a] > perm[b])
    return -1.0 if inversions % 2 else 1.0


def invariant_coeff(reps, indices):
    nontrivial = [(rep, index) for rep, index in zip(reps, indices) if rep != "1"]
    if not nontrivial:
        return 1.0
    labels = [rep for rep, _ in nontrivial]
    values = [index for _, index in nontrivial]
    if sorted(labels) == ["3", "3b"] and len(labels) == 2:
        return (1 / np.sqrt(3)) if values[0] == values[1] else 0.0
    if labels == ["3", "3", "3"]:
        return eps(*values) / np.sqrt(6)
    if labels == ["3b", "3b", "3b"]:
        return eps(*values) / np.sqrt(6)
    return 0.0


def vertex_factors(link_reps):
    r0, r1, r2, r3 = link_reps
    return [
        [MATTER_REPS[0], r0, r3],
        [MATTER_REPS[1], dual(r0), r1],
        [MATTER_REPS[2], dual(r1), dual(r2)],
        [MATTER_REPS[3], dual(r3), r2],
    ]


def allowed_link_assignments():
    allowed = []
    for link_reps in itertools.product(REPS, repeat=4):
        ok = True
        for factors in vertex_factors(link_reps):
            found = False
            for indices in itertools.product(*[range(dim(rep)) for rep in factors]):
                if abs(invariant_coeff(factors, indices)) > 1e-14:
                    found = True
                    break
            ok = ok and found
        if ok:
            allowed.append(link_reps)
    return allowed


def local_link_basis(rep):
    return [(rep, a, b) for a in range(dim(rep)) for b in range(dim(rep))]


def spin_network_state(link_reps):
    """
    Sparse vector in matter-index x four-link Peter-Weyl direct-sum basis.

    Key:
      (m0,m1,m2,m3, link0,link1,link2,link3)
    where link entries are (rep, tail_index, head_index).
    """
    state = defaultdict(complex)
    link_bases = [local_link_basis(rep) for rep in link_reps]
    matter_ranges = [range(3) for _ in range(4)]
    for matter_indices in itertools.product(*matter_ranges):
        for link_entries in itertools.product(*link_bases):
            # Endpoint factor order at each vertex follows vertex_factors().
            v0 = invariant_coeff(
                [MATTER_REPS[0], link_reps[0], link_reps[3]],
                [matter_indices[0], link_entries[0][1], link_entries[3][1]],
            )
            if abs(v0) < 1e-14:
                continue
            v1 = invariant_coeff(
                [MATTER_REPS[1], dual(link_reps[0]), link_reps[1]],
                [matter_indices[1], link_entries[0][2], link_entries[1][1]],
            )
            if abs(v1) < 1e-14:
                continue
            v2 = invariant_coeff(
                [MATTER_REPS[2], dual(link_reps[1]), dual(link_reps[2])],
                [matter_indices[2], link_entries[1][2], link_entries[2][2]],
            )
            if abs(v2) < 1e-14:
                continue
            v3 = invariant_coeff(
                [MATTER_REPS[3], dual(link_reps[3]), link_reps[2]],
                [matter_indices[3], link_entries[3][2], link_entries[2][1]],
            )
            if abs(v3) < 1e-14:
                continue
            coeff = v0 * v1 * v2 * v3
            state[matter_indices + link_entries] += coeff
    return dict(state)


def inner(left, right):
    if len(left) > len(right):
        left, right = right, left
        conjugate_left = False
    else:
        conjugate_left = True
    total = 0.0j
    if conjugate_left:
        for key, value in left.items():
            total += np.conj(value) * right.get(key, 0.0)
    else:
        for key, value in left.items():
            total += np.conj(right.get(key, 0.0)) * value
        total = np.conj(total)
    return total


def u_components(i, j, source):
    """
    Truncated multiplication by fundamental matrix element U_ij.

    Returns [(target_link_entry, coeff), ...] for one source link entry.
    Normalizations are the standard singlet/epsilon CG normalizations in this
    minimal truncation; channels 6 and 8 are intentionally projected out.
    """
    cache_key = (i, j, source)
    if cache_key in _U_CACHE:
        return _U_CACHE[cache_key]

    rep, a, b = source
    out = []
    if rep == "1":
        out.append((("3", i, j), 1 / np.sqrt(3)))
    elif rep == "3":
        for alpha in range(3):
            for beta in range(3):
                coeff = 0.5 * eps(alpha, i, a) * eps(beta, j, b)
                if abs(coeff) > 1e-14:
                    out.append((("3b", alpha, beta), coeff))
    elif rep == "3b":
        coeff = (1 / np.sqrt(3)) if i == a and j == b else 0.0
        if abs(coeff) > 1e-14:
            out.append((("1", 0, 0), coeff))
    _U_CACHE[cache_key] = out
    return out


def udag_components(i, j, source):
    """Components of U^dag_ij = (U_ji)^dag, built as the adjoint map."""
    cache_key = (i, j, source)
    if cache_key in _UDAG_CACHE:
        return _UDAG_CACHE[cache_key]

    out = []
    for candidate_rep in REPS:
        for candidate in local_link_basis(candidate_rep):
            for target, coeff in u_components(j, i, candidate):
                if target == source:
                    out.append((candidate, np.conj(coeff)))
    _UDAG_CACHE[cache_key] = out
    return out


def apply_plaquette(state):
    out = defaultdict(complex)
    for key, coeff in state.items():
        matter = key[:4]
        links = key[4:]
        for a, b, c, d in itertools.product(range(3), repeat=4):
            link_maps = [
                u_components(a, b, links[0]),
                u_components(b, c, links[1]),
                udag_components(c, d, links[2]),
                udag_components(d, a, links[3]),
            ]
            if any(not entries for entries in link_maps):
                continue
            for targets in itertools.product(*link_maps):
                target_links = tuple(target for target, _ in targets)
                amp = coeff / 3
                for _, local_coeff in targets:
                    amp *= local_coeff
                if abs(amp) > 1e-14:
                    out[matter + target_links] += amp
    return dict(out)


def basis_data():
    global _BASIS_CACHE
    if _BASIS_CACHE is not None:
        return _BASIS_CACHE

    allowed = allowed_link_assignments()
    states = [spin_network_state(link_reps) for link_reps in allowed]
    gram = np.array([[inner(left, right) for right in states] for left in states])
    herm = float(np.linalg.norm(gram - gram.conj().T))
    eigvals = np.linalg.eigvalsh(gram)
    assert herm < 1e-10
    assert np.min(eigvals) > 1e-10
    _BASIS_CACHE = (allowed, states, gram)
    return _BASIS_CACHE


def orthonormalize(matrix, gram):
    evals, evecs = np.linalg.eigh(gram)
    g_inv_sqrt = evecs @ np.diag(1 / np.sqrt(evals)) @ evecs.conj().T
    out = g_inv_sqrt @ matrix @ g_inv_sqrt
    herm = float(np.linalg.norm(out - out.conj().T))
    assert herm < 1e-9
    return out


def tch_lock_matrix(allowed, strength):
    """Candidate gauge-invariant low-sector lock in the physical spin-network basis."""
    lock = np.zeros((len(allowed), len(allowed)), complex)
    low_a = ("1", "3", "1", "3b")
    low_b = ("3b", "1", "3", "1")
    if low_a in allowed and low_b in allowed:
        i = allowed.index(low_a)
        j = allowed.index(low_b)
        lock[i, j] = -strength
        lock[j, i] = -strength
    return lock


def hamiltonian(beta, tch_lock=0.0):
    allowed, states, gram = basis_data()
    g2 = 1 / beta
    electric = np.zeros_like(gram)
    for index, link_reps in enumerate(allowed):
        energy = g2 * sum(casimir(rep) for rep in link_reps)
        electric[index, index] = energy * gram[index, index]

    global _WILSON_CACHE
    if _WILSON_CACHE is None:
        plaquette_states = [apply_plaquette(state) for state in states]
        _WILSON_CACHE = np.array(
            [[inner(left, right) for right in plaquette_states] for left in states]
        )
    wilson = _WILSON_CACHE
    magnetic = -(beta / 2) * (wilson + wilson.conj().T)
    h_raw = electric + magnetic
    h = orthonormalize(h_raw, gram)
    if tch_lock:
        h += tch_lock_matrix(allowed, tch_lock)
    return h, allowed, gram, wilson


def audit():
    hd("A. Explicit Intertwiner Basis Audit")
    allowed, states, gram = basis_data()
    print(f"  matter reps at vertices          = {MATTER_REPS}")
    print(f"  allowed link-rep spin networks   = {len(allowed)}")
    for index, link_reps in enumerate(allowed):
        norm = np.real(gram[index, index])
        print(f"    {index}: links={link_reps}, sparse entries={len(states[index])}, norm={norm:.6g}")
    print(f"  Gram eigenvalues                 = {np.round(np.linalg.eigvalsh(gram), 10)}")
    assert len(allowed) == 3

    h, _, _, wilson = hamiltonian(beta=1.0)
    print(f"  Wilson matrix Frobenius norm     = {np.linalg.norm(wilson):.6g}")
    print(f"  projected H(beta=1) hermiticity  = {np.linalg.norm(h-h.conj().T):.3e}")


def scan_wilson_only():
    hd("B. beta <= 1 Peter-Weyl Wilson-Only Spectrum")
    print("  beta       dim    gap       degeneracy   eigvals")
    gaps = []
    for beta in [0.25, 0.5, 0.75, 1.0]:
        h, allowed, _, _ = hamiltonian(beta)
        eigvals = np.linalg.eigvalsh(h)
        rounded = []
        degeneracies = []
        for value in eigvals:
            if not rounded or abs(value - rounded[-1]) > 1e-8:
                rounded.append(float(value))
                degeneracies.append(1)
            else:
                degeneracies[-1] += 1
        gap = rounded[1] - rounded[0]
        gaps.append(gap)
        print(
            f"  {beta:<10g} "
            f"{len(allowed):<6d} "
            f"{gap:<9.6g} "
            f"{degeneracies[0]:<12d} "
            f"{np.round(eigvals, 6)}"
        )
    print(f"  minimum scanned gap = {min(gaps):.6g}")
    assert min(gaps) < 0.1
    print(
        "  -> Negative result: the minimal Wilson-only Peter-Weyl truncation has "
        "a near-degenerate strong-coupling pair. The qutrit/S3 successes do not "
        "by themselves imply a robust full-colour gap."
    )


def scan_with_tch_lock():
    hd("C. Candidate TCH-Lock Scale Scan")
    print(
        "  The lock is gauge-invariant because it acts after explicit Gauss "
        "projection, coupling the two low electric-flux spin networks. It is a "
        "candidate item-13 term, not derived from the full TCH cell geometry."
    )
    print("  lambda     beta=0.25 gap   beta=0.5 gap    beta=1 gap")
    for strength in [0.02, 0.05, 0.1, 0.25, 0.5, 1.0]:
        row = []
        for beta in [0.25, 0.5, 1.0]:
            h, _, _, _ = hamiltonian(beta, tch_lock=strength)
            eigvals = np.linalg.eigvalsh(h)
            row.append(float(eigvals[1] - eigvals[0]))
        print(
            f"  {strength:<10g} "
            f"{row[0]:<15.6g} "
            f"{row[1]:<15.6g} "
            f"{row[2]:<12.6g}"
        )
    assert np.linalg.eigvalsh(hamiltonian(0.5, tch_lock=0.5)[0])[1] - np.linalg.eigvalsh(
        hamiltonian(0.5, tch_lock=0.5)[0]
    )[0] > 0.5


def main():
    hd("Scope")
    print(
        "This is a minimal Peter-Weyl SU(3) truncation, not the full SU(3) lattice "
        "gauge theory. It includes actual 1,3,3bar link irreps and explicit "
        "delta/epsilon Gauss intertwiners, then computes the beta<=1 spectrum."
    )
    audit()
    scan_wilson_only()
    scan_with_tch_lock()
    hd("Verdict")
    print(
        "In the smallest nonzero bipartite matter sector, the explicit "
        "(1,3,3bar) Peter-Weyl Wilson-only truncation does NOT have a robust "
        "strong-coupling gap: it produces a near-degenerate low pair. A finite "
        "gap requires an additional gauge-invariant TCH/intertwiner-lock term."
    )
    print(
        "\nBound: the lock scanned here is only a physical-basis candidate. The "
        "real item-13 task is to derive such a lock from the TCH cell operator, "
        "then include at least the 8 and sextet channels and check convergence "
        "of the strong-coupling gap."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
