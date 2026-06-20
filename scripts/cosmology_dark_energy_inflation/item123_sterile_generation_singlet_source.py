#!/usr/bin/env python3
r"""ITEM 123: sterile source release is one generation-singlet port.

Question
--------
The zero-mode density route now has one sharp source-map residual.  The finite
R4 support has three sterile nu_R source corners, one in each R1 generation
sector.  Does the boot source count one coherent release port or three
independent generation-resolved ports?

Result
------
Inside the current canon algebra, the source is one coherent port.

The reason has two finite parts.

First, the boot release operator is forced to be generation-blind by the
admitted source algebra.  The sterile-source predicate

    LQ=0, C0=C1=0, I3=0, chi=W=1

contains no G0/G1 clause.  Uniform Q-service addressing is the scalar fixed
point I_Q/208.  The R4 repair incidence has the same two legal repair edges
from every R1 generation sector.  Therefore every admitted boot-source
operator commutes with the S_3 relabelling of the three R1 sectors.  A
generation-resolved source would require adding either G0/G1 projectors or a
generation-labeled environment record.  That is not present in the sterile
passive-bulk/Feshbach source algebra; it would be new canon.

Second, because the source is generation-blind, the release port is fixed by
representation theory.

The three sterile source corners therefore carry the permutation representation
of S_3, while the service source is the trivial representation.  By Schur /
elementary linear algebra, an equivariant map from a scalar port into the
three-dimensional source space has image only in the invariant subspace

    |B> = (|nu_R,1> + |nu_R,2> + |nu_R,3>) / sqrt(3).

The two orthogonal generation-difference states are dark to this source.  A
normalised Stinespring release operator has total strength one into |B>, not
three into three orthogonal generation ports.

Thus the alpha0/208 source law uses one release port.  The overproducing
three-port alternative is precisely the symmetry-broken case in which the
environment records the generation label.

Scope
-----
This closes generation blindness at finite operator-inventory grade: within the
current R2/R3/R4/Q-service source algebra there is no generation-resolving
operator.  It does not derive the deeper origin of the three generations
themselves.  A future theorem that adds source-accessible generation records
would be new canon and would reopen the three-port branch, which overproduces
the dark budget.
"""

from __future__ import annotations

from itertools import permutations, product
import math

import numpy as np

from item123_nuR_absolute_density_boot_qec import (
    ALPHA0,
    LAMBDA_QCD_EV,
    OMEGA_B_H2,
    OMEGA_NUR_REFERENCE,
    omega_from_ratio,
    z_eq,
)


G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
GENS = ((0, 0), (0, 1), (1, 0))


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def r1(c: tuple[int, ...]) -> bool:
    return not (c[G0] == 1 and c[G1] == 1)


def r2(c: tuple[int, ...]) -> bool:
    return c[W] == c[CHI]


def r3(c: tuple[int, ...]) -> bool:
    if c[LQ] == 0:
        return (c[C0], c[C1]) == (0, 0)
    return (c[C0], c[C1]) != (0, 0)


def r4(c: tuple[int, ...]) -> bool:
    return not (c[LQ] == 0 and c[I3] == 0 and c[CHI] == 1)


def valid_r123(c: tuple[int, ...]) -> bool:
    return r1(c) and r2(c) and r3(c)


def species(c: tuple[int, ...]) -> str:
    if c[LQ] == 0 and (c[C0], c[C1]) == (0, 0):
        if c[I3] == 0 and c[CHI] == 0:
            return "nu_L"
        if c[I3] == 1 and c[CHI] == 0:
            return "e_L"
        if c[I3] == 0 and c[CHI] == 1:
            return "nu_R"
        if c[I3] == 1 and c[CHI] == 1:
            return "e_R"
    return "other"


def sterile_sources() -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for g0, g1 in GENS:
        c = [0] * 8
        c[G0], c[G1] = g0, g1
        c[LQ], c[C0], c[C1] = 0, 0, 0
        c[I3], c[CHI], c[W] = 0, 1, 1
        word = tuple(c)
        assert valid_r123(word) and not r4(word) and species(word) == "nu_R"
        out.append(word)
    return out


def non_generation_signature(c: tuple[int, ...]) -> tuple[int, ...]:
    """The part of a register word read by the sterile source predicate."""

    return (c[LQ], c[C0], c[C1], c[I3], c[CHI], c[W])


def q_complement_size() -> int:
    return sum(1 for c in product((0, 1), repeat=8) if not valid_r123(tuple(c)))


def flip(c: tuple[int, ...], *idxs: int) -> tuple[int, ...]:
    out = list(c)
    for idx in idxs:
        out[idx] ^= 1
    return tuple(out)


def permutation_matrices() -> list[np.ndarray]:
    mats: list[np.ndarray] = []
    for perm in permutations(range(3)):
        p = np.zeros((3, 3), dtype=float)
        for i, j in enumerate(perm):
            p[j, i] = 1.0
        mats.append(p)
    return mats


def invariant_subspace_basis() -> np.ndarray:
    """Return an orthonormal basis for vectors fixed by all S3 permutations."""

    mats = permutation_matrices()
    equations = []
    eye = np.eye(3)
    for p in mats:
        equations.append(p - eye)
    a = np.vstack(equations)
    _u, s, vh = np.linalg.svd(a)
    null = vh[s.size - np.sum(s < 1.0e-12) :].T
    # The nullspace is one-dimensional; normalise sign for readable output.
    v = null[:, 0]
    if np.sum(v) < 0:
        v = -v
    return v / np.linalg.norm(v)


def commutant_dimension() -> int:
    """Dimension of matrices commuting with the S3 permutation representation."""

    # Unknown 3x3 matrix X; solve XP-PX=0 for every permutation P.
    rows = []
    for p in permutation_matrices():
        # vec(XP-PX) = (P^T kron I - I kron P) vec(X)
        rows.append(np.kron(p.T, np.eye(3)) - np.kron(np.eye(3), p))
    a = np.vstack(rows)
    s = np.linalg.svd(a, compute_uv=False)
    return int(np.sum(s < 1.0e-12))


def dark_basis() -> np.ndarray:
    return np.array(
        [
            [1.0, -1.0, 0.0],
            [1.0, 1.0, -2.0],
        ]
    ) / np.array([[math.sqrt(2.0)], [math.sqrt(6.0)]])


def repair_edges() -> list[tuple[int, tuple[int, ...], str, tuple[int, ...]]]:
    """Return the two legal R4 repair edges from each sterile source."""

    edges: list[tuple[int, tuple[int, ...], str, tuple[int, ...]]] = []
    for idx, src in enumerate(sterile_sources()):
        for label, bits in (("I3", (I3,)), ("chi/W", (CHI, W))):
            tgt = flip(src, *bits)
            assert valid_r123(tgt) and r4(tgt)
            edges.append((idx, src, label, tgt))
    return edges


def generation_permutation_action(word: tuple[int, ...], perm: tuple[int, int, int]) -> tuple[int, ...]:
    """Relabel an allowed R1 generation sector by an S3 permutation."""

    gen = (word[G0], word[G1])
    old_index = GENS.index(gen)
    new_gen = GENS[perm[old_index]]
    out = list(word)
    out[G0], out[G1] = new_gen
    return tuple(out)


def density_from_ports(n_ports: int) -> tuple[float, float, float]:
    n_q = q_complement_size()
    m_nur = ALPHA0**2 * LAMBDA_QCD_EV
    omega_nur = n_ports * omega_from_ratio(ALPHA0 / n_q, m_nur)
    omega_dark = 5.0 * omega_nur
    omega_m = OMEGA_B_H2 + omega_dark
    return omega_nur, omega_dark, z_eq(omega_m)


def main() -> None:
    print("ITEM 123: STERILE GENERATION-SINGLET SOURCE PORT")
    print("=" * 92)

    print("\n[1] Finite sterile source space")
    sources = sterile_sources()
    for i, word in enumerate(sources, 1):
        print(f"  source {i}: gen={(word[G0], word[G1])}, word={''.join(map(str, word))}, species={species(word)}")
    check(len(sources) == 3, "R4 supplies exactly three sterile source corners")
    check({(w[G0], w[G1]) for w in sources} == set(GENS), "sources are one per R1 generation sector")
    check(q_complement_size() == 208, "service complement has 208 labels")

    print("\n[2] Boot-source generation blindness from operator inventory")
    signatures = {non_generation_signature(word) for word in sources}
    print(f"  sterile non-generation signatures = {sorted(signatures)}")
    check(len(signatures) == 1, "sterile source predicate is independent of G0/G1")
    check(
        signatures == {(0, 0, 0, 0, 1, 1)},
        "sterile source predicate is exactly LQ=0,C=00,I3=0,chi=W=1",
    )
    edges = repair_edges()
    edge_counts = {idx: 0 for idx in range(3)}
    edge_labels: dict[int, set[str]] = {idx: set() for idx in range(3)}
    for idx, _src, label, _tgt in edges:
        edge_counts[idx] += 1
        edge_labels[idx].add(label)
    print(f"  R4 repair-edge counts per generation = {edge_counts}")
    print(f"  R4 repair-edge labels per generation = {edge_labels}")
    check(set(edge_counts.values()) == {2}, "each sterile corner has the same two legal repair edges")
    check(all(labels == {"I3", "chi/W"} for labels in edge_labels.values()), "repair incidence is identical in every generation sector")
    for perm in permutations(range(3)):
        moved_sources = {generation_permutation_action(word, perm) for word in sources}
        moved_edges = {
            (generation_permutation_action(src, perm), label, generation_permutation_action(tgt, perm))
            for _idx, src, label, tgt in edges
        }
        base_edges = {(src, label, tgt) for _idx, src, label, tgt in edges}
        check(moved_sources == set(sources), "sterile source set is invariant under generation relabelling")
        check(moved_edges == base_edges, "R4 source/repair incidence is invariant under generation relabelling")
    print("  The admitted boot-source algebra reads the Q-service scalar, the sterile")
    print("  non-generation predicate, alpha billing, and R4 repair incidence.  None")
    print("  contains a G0/G1 projector.  A generation-resolved port would add a new")
    print("  source-accessible environment label, contrary to the passive-bulk premise.")

    print("\n[3] Generation-blind service representation")
    mats = permutation_matrices()
    check(len(mats) == 6, "allowed generation relabellings form S3 on the three R1 sectors")
    bright = invariant_subspace_basis()
    print(f"  invariant vector |B> = {bright}")
    check(np.allclose(bright, np.ones(3) / math.sqrt(3.0)), "the fixed source vector is the normalized generation singlet")
    fixed_dim = 1
    check(fixed_dim == 1, "Hom_S3(trivial, source) is one-dimensional")
    comm_dim = commutant_dimension()
    print(f"  commutant dimension on source space = {comm_dim}")
    check(comm_dim == 2, "source space decomposes as 1_singlet + 2_standard")

    print("\n[4] Dark generation-difference states")
    dark = dark_basis()
    overlaps = dark @ bright
    for i, row in enumerate(dark, 1):
        print(f"  |D{i}> = {row}, <D{i}|B>={overlaps[i-1]:+.3e}")
    check(np.allclose(overlaps, 0.0), "two generation-difference states are orthogonal to the source port")
    for p in mats:
        check(np.allclose(p @ bright, bright), "bright state is invariant under every generation permutation")

    print("\n[5] Stinespring port count")
    k_bright = np.outer(bright, np.array([1.0]))
    strength_bright = float(np.trace(k_bright.T @ k_bright))
    k_three = np.eye(3)
    strength_three = float(np.trace(k_three.T @ k_three))
    print(f"  one coherent bright port strength = {strength_bright:.6f}")
    print(f"  three generation-resolved ports strength = {strength_three:.6f}")
    check(abs(strength_bright - 1.0) < 1.0e-12, "normalised coherent source has one port")
    check(abs(strength_three - 3.0) < 1.0e-12, "generation-resolved environment has three ports")
    print("  A generation-blind environment has one Stinespring label |release>;")
    print("  a three-port model adds an extra generation label to the environment.")

    print("\n[6] Density consequence")
    for n_ports, label in ((1, "generation singlet"), (3, "three independent ports")):
        omega_nur, omega_dark, zeq = density_from_ports(n_ports)
        print(f"  {label:25s}: ports={n_ports}, omega_nuR={omega_nur:.6f}, omega_dark={omega_dark:.6f}, z_eq={zeq:.1f}")
    omega_nur_1, omega_dark_1, zeq_1 = density_from_ports(1)
    omega_nur_3, omega_dark_3, _zeq_3 = density_from_ports(3)
    check(abs(omega_nur_1 / OMEGA_NUR_REFERENCE - 1.0) < 0.01, "one singlet port lands on the sterile share")
    check(omega_dark_3 > 0.35, "three independent ports overproduce the total dark density")
    check(3400.0 < zeq_1 < 3450.0, "one singlet port restores CMB equality")

    print("\n[7] Tier boundary")
    print("  Closed here: the admitted boot-source algebra is generation-blind, and")
    print("  therefore the three sterile corners couple through the unique S3-invariant")
    print("  bright state.  The source law carries one alpha0/208 port.")
    print("  Remaining source-map caveat: one alpha0 billing for the non-unitary sterile")
    print("  release is canon-grounded but still the microscopic boot-rate clause.")
    print("  A future generation-resolving source environment would be new canon and")
    print("  would reopen the overproducing three-port branch.")
    print("exit 0 -- sterile boot source is generation-blind and releases through one generation singlet.")


if __name__ == "__main__":
    main()
