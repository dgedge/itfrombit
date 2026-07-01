#!/usr/bin/env python3
r"""ITEM 118: sterile-nu_R active-mixing theorem gate.

Target
------
Promote the item-118 X-ray flux input

    theta = m_s / v_R4

from a dimensional estimate to a finite-algebra + Schur-complement theorem.

The theorem proved here is conditional but sharp:

    If the R4-forbidden neutral sterile-active edge has the service-normalised
    bright Schur spectral moment

        K_B = <q_B|(QHQ)^{-1}|q_B> = 1/v_R4,

    and the two edge insertions are the already-normalised sterile mass
    insertion m_s, then

        m_D = m_s^2 / v_R4,
        tan(2 theta) = 2 m_D / m_s,
        sin^2(2 theta) = 4 (m_s/v_R4)^2 / (1 + 4 (m_s/v_R4)^2).

What is actually proved by the finite register in this file:

  1. the R4 sterile source has exactly two local repair edges per generation;
  2. only one of them is neutral and lands on nu_L: the chi/W edge;
  3. generation-blindness makes the bright sterile port couple to the bright
     active-neutrino port with coefficient exactly 1, not sqrt(3) or 3;
  4. all generation-difference sterile states are dark to this source.

What remains a named premise:

  The finite register and S3 equivariance compress the virtual sector to one
  bright scalar spectral moment.  They do not by themselves derive the value of
  that moment.  The remaining premise is the service normalisation
  K_B=1/v_R4.  This may be realised by a single virtual state of gap v_R4, or
  by unresolved microstates whose weighted spectral sum equals 1/v_R4.  Under
  that premise the flux branch is no longer a fitted mixing angle; without it,
  the script reports the single missing scalar.
"""

from __future__ import annotations

from itertools import permutations
import math

import numpy as np


# Bit convention from item123_nuR_absolute_density_boot_qec.py:
# (G0, G1, LQ, C0, C1, I3, CHI, W)
G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
GENS = ((0, 0), (0, 1), (1, 0))

ALPHA0 = 1.0 / 137.035999
LAMBDA_QCD_EV = 0.332e9
V_R4_EV = 246.0e9


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


def word_for(gen: tuple[int, int], *, i3: int, chi: int) -> tuple[int, ...]:
    c = [0] * 8
    c[G0], c[G1] = gen
    c[LQ], c[C0], c[C1] = 0, 0, 0
    c[I3], c[CHI], c[W] = i3, chi, chi
    return tuple(c)


def sterile_sources() -> list[tuple[int, ...]]:
    out = [word_for(gen, i3=0, chi=1) for gen in GENS]
    for word in out:
        assert valid_r123(word) and not r4(word) and species(word) == "nu_R"
    return out


def active_neutrinos() -> list[tuple[int, ...]]:
    out = [word_for(gen, i3=0, chi=0) for gen in GENS]
    for word in out:
        assert valid_r123(word) and r4(word) and species(word) == "nu_L"
    return out


def flip(c: tuple[int, ...], *idxs: int) -> tuple[int, ...]:
    out = list(c)
    for idx in idxs:
        out[idx] ^= 1
    return tuple(out)


def repair_edges() -> list[tuple[int, str, tuple[int, ...], tuple[int, ...]]]:
    edges: list[tuple[int, str, tuple[int, ...], tuple[int, ...]]] = []
    for idx, src in enumerate(sterile_sources()):
        for label, bits in (("I3", (I3,)), ("chi/W", (CHI, W))):
            tgt = flip(src, *bits)
            assert valid_r123(tgt) and r4(tgt)
            edges.append((idx, label, src, tgt))
    return edges


def permutation_matrices() -> list[np.ndarray]:
    mats: list[np.ndarray] = []
    eye = np.eye(3)
    for perm in permutations(range(3)):
        mats.append(eye[list(perm)])
    return mats


def bright() -> np.ndarray:
    return np.ones(3) / math.sqrt(3.0)


def dark_basis() -> np.ndarray:
    return np.array(
        [
            [1.0, -1.0, 0.0],
            [1.0, 1.0, -2.0],
        ]
    ) / np.array([[math.sqrt(2.0)], [math.sqrt(6.0)]])


def neutral_edge_matrix() -> np.ndarray:
    """Matrix from sterile generation basis to active-neutrino generation basis."""

    sources = sterile_sources()
    targets = active_neutrinos()
    mat = np.zeros((3, 3), dtype=float)
    for src_idx, label, src, tgt in repair_edges():
        if label == "chi/W":
            tgt_idx = targets.index(tgt)
            assert sources[src_idx] == src
            mat[tgt_idx, src_idx] = 1.0
    return mat


def schur_effective_mass(m_s: float, v_r4: float, edge_coeff: float = 1.0) -> float:
    """One virtual R4-enforcement denominator: m_D = c^2 m_s^2 / v_R4."""

    return edge_coeff**2 * m_s * m_s / v_r4


def spectral_moment(gaps: list[float], couplings: list[float]) -> float:
    """Bright Schur moment K_B = sum_a |g_a|^2 / Delta_a."""

    if len(gaps) != len(couplings):
        raise ValueError("gaps and couplings must have the same length")
    return sum((g * g) / gap for gap, g in zip(gaps, couplings, strict=True))


def schur_effective_mass_from_moment(m_s: float, k_b: float) -> float:
    """Feshbach/Schur complement: m_D = m_s^2 K_B."""

    return m_s * m_s * k_b


def exact_sin2_2theta(m_s: float, m_d: float) -> float:
    """For [[0,m_D],[m_D,m_s]], tan(2theta)=2m_D/m_s."""

    return 4.0 * m_d * m_d / (m_s * m_s + 4.0 * m_d * m_d)


def main() -> None:
    print("ITEM 118: sterile-nu_R active-mixing theorem gate")
    print("=" * 92)

    print("\n[1] Local R4 repair incidence")
    edges = repair_edges()
    for idx, label, src, tgt in edges:
        print(
            f"  gen={idx + 1}, edge={label:5s}: "
            f"{species(src):4s} {''.join(map(str, src))} -> "
            f"{species(tgt):4s} {''.join(map(str, tgt))}"
        )
    edge_labels = {i: {label for j, label, _src, _tgt in edges if i == j} for i in range(3)}
    check(all(labels == {"I3", "chi/W"} for labels in edge_labels.values()), "each sterile source has exactly two local R4 repair edges")
    neutral_edges = [(i, src, tgt) for i, label, src, tgt in edges if label == "chi/W"]
    charged_edges = [(i, src, tgt) for i, label, src, tgt in edges if label == "I3"]
    check(all(species(tgt) == "nu_L" for _i, _src, tgt in neutral_edges), "the chi/W repair edge is the unique neutral nu_R -> nu_L edge")
    check(all(species(tgt) == "e_R" for _i, _src, tgt in charged_edges), "the I3 repair edge is charged and cannot supply the X-ray active-neutrino mixing")

    print("\n[2] Generation-singlet normalization")
    M = neutral_edge_matrix()
    b = bright()
    d = dark_basis()
    bright_coeff = float(b @ M @ b)
    dark_to_bright = d @ M @ b
    bright_to_dark = d @ M.T @ b
    print(f"  neutral edge matrix =\n{M}")
    print(f"  <B_nu|M_neutral|B_s> = {bright_coeff:.12f}")
    print(f"  <D_i|M_neutral|B_s>  = {dark_to_bright}")
    print(f"  <D_i|M_neutral^T|B_nu> = {bright_to_dark}")
    check(np.allclose(M, np.eye(3)), "neutral repair is the identity on generation labels")
    check(abs(bright_coeff - 1.0) < 1.0e-12, "bright sterile port couples to bright active neutrino with coefficient exactly 1")
    check(np.allclose(dark_to_bright, 0.0), "generation-difference active states do not couple to the bright sterile source")
    check(np.allclose(bright_to_dark, 0.0), "generation-difference sterile states do not couple to the bright active channel")
    for P in permutation_matrices():
        check(np.allclose(P.T @ M @ P, M), "neutral mixing edge is S3-equivariant")

    print("\n[3] Schur-complement denominator tightening")
    print("  Feshbach form: m_D = m_s^2 K_B, K_B=<q_B|(QHQ)^-1|q_B>")
    print("  The finite/S3 part fixes the bright channel; microstructure can only change kappa=v_R4 K_B.")
    denominator_cases = {
        "single virtual state": ([V_R4_EV], [1.0]),
        "two degenerate, unnormalised": ([V_R4_EV, V_R4_EV], [1.0, 1.0]),
        "two degenerate, service-normalised": ([V_R4_EV, V_R4_EV], [1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0)]),
        "split gaps, service-normalised": ([0.8 * V_R4_EV, 1.2 * V_R4_EV], [1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0)]),
    }
    kappas: dict[str, float] = {}
    for label, (gaps, couplings) in denominator_cases.items():
        k_b = spectral_moment(gaps, couplings)
        kappa = V_R4_EV * k_b
        kappas[label] = kappa
        print(f"  {label:34s}: kappa=v_R4 K_B={kappa:.12f}")
    check(abs(kappas["single virtual state"] - 1.0) < 1.0e-12, "one gap v_R4 gives the one-denominator branch")
    check(abs(kappas["two degenerate, service-normalised"] - 1.0) < 1.0e-12, "unresolved degenerate microstates do not change the branch when service-normalised")
    check(abs(kappas["two degenerate, unnormalised"] - 2.0) < 1.0e-12, "un-normalised hidden multiplicity would double the Schur moment")
    check(abs(kappas["split gaps, service-normalised"] - 1.0) > 1.0e-2, "non-degenerate unresolved gaps are a measurable scalar correction")
    print("  [TIGHTENED] the live premise is K_B=1/v_R4, not literally one microscopic state.")

    print("\n[4] Schur-complement mixing law")
    m_s = ALPHA0**2 * LAMBDA_QCD_EV
    k_b = spectral_moment([V_R4_EV], [bright_coeff])
    m_d = schur_effective_mass_from_moment(m_s, k_b)
    eps = m_s / V_R4_EV
    sin2 = exact_sin2_2theta(m_s, m_d)
    sin2_small = 4.0 * eps * eps
    print(f"  m_s = alpha0^2 Lambda_QCD = {m_s/1e3:.6f} keV")
    print(f"  v_R4 = {V_R4_EV/1e9:.3f} GeV")
    print(f"  epsilon = m_s/v_R4 = {eps:.12e}")
    print(f"  m_D = m_s^2/v_R4 = {m_d*1e3:.6f} meV")
    print(f"  sin^2(2theta) exact two-level = {sin2:.12e}")
    print(f"  small-angle 4 epsilon^2       = {sin2_small:.12e}")
    check(abs(sin2 / sin2_small - 1.0) < 1.0e-12, "small-angle formula is indistinguishable at this scale")
    check(2.0e-14 < sin2 < 2.2e-14, "mixing lands at the item-118 R4-scale value")

    print("\n[5] Alternative denominators / hidden multiplicity controls")
    alternatives = {
        "one bright port (derived)": 1.0,
        "three incoherent ports": math.sqrt(3.0),
        "two repair edges counted": math.sqrt(2.0),
        "charged+neutral and 3 ports": math.sqrt(6.0),
    }
    for label, coeff in alternatives.items():
        md_alt = schur_effective_mass(m_s, V_R4_EV, coeff)
        sin2_alt = exact_sin2_2theta(m_s, md_alt)
        print(f"  {label:28s}: edge coeff={coeff:.6f}, sin^2(2theta)={sin2_alt:.3e}")
    check(exact_sin2_2theta(m_s, schur_effective_mass(m_s, V_R4_EV, math.sqrt(3.0))) > 8.0 * sin2, "three-port hidden multiplicity would be observationally distinct")
    check(exact_sin2_2theta(m_s, schur_effective_mass(m_s, V_R4_EV, math.sqrt(2.0))) > 3.0 * sin2, "counting the charged edge would be observationally distinct")

    print(
        """
[6] VERDICT
  Finite part CLOSED:
    The X-ray-active sterile mixing channel is the unique neutral R4 repair edge
    nu_R -> nu_L (chi/W).  It is the identity on generation labels, so the
    generation-singlet bright port has coefficient exactly 1.  The two
    generation-difference sterile states are dark to this source.  This removes
    the hidden sqrt(2), sqrt(3), and 3-port normalisation freedoms.

  Schur-complement tightening:
    The local generation-blind virtual sector contributes only through the
    bright spectral moment K_B=<q_B|(QHQ)^-1|q_B>.  Hidden microstates do not add
    a matrix or a generation multiplicity; they only change the scalar
    kappa=v_R4 K_B.  The named one-denominator branch is the service-normalised
    condition K_B=1/v_R4.  It can be realised by one virtual state or by a
    normalised unresolved multiplet.

  Remaining named premise:
    The microscopic Hamiltonian must derive kappa=1 for the neutral R4 repair
    moment.  If a future derivation finds kappa != 1, the flux rescales by
    kappa^2.  But there is no longer a free active-sterile mixing matrix,
    generation factor, or hidden denominator count inside the branch.

exit 0 -- mixing theorem closed at conditional Schur-complement grade.
"""
    )


if __name__ == "__main__":
    main()
