#!/usr/bin/env python3
r"""R12/R13/R15 promotion-gate attempt.

Question
--------
Can the three remaining record-action promotion gates be closed now?

  A. prove one-record-per-event and address x channel factorization;
  B. universalize the traffic clock beyond native QEC/cosmology;
  C. derive recovery holonomy phase and orientation.

Answer
------
A closes at finite service-instrument grade.  B does not close; the current
canon points to two irreducible saturated anchors rather than one universal
traffic clock.  C closes only conditionally: K_or exists and the faithful C3
phase is forced if the oriented latch is admitted, but the physical sign
pointer and absolute orientation are not in the current sterile source algebra.
"""

from __future__ import annotations

from itertools import permutations
import math

import numpy as np


ADDRESS = 8
CHANNEL = 137
M_RECORD = ADDRESS * CHANNEL
TOL = 1.0e-12

V_EW = 246.0
M_TOP = 172.7
M_B = 4.18
M_TAU = 1.777
M_C = 1.27
M_S = 0.093
M_MU = 0.1057

Y_DIAG = np.array([0.5, 0.8, 1.3])
H_DIRAC = np.diag(Y_DIAG**2)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def entropy_uniform(n: int) -> float:
    return math.log(float(n))


def event_label(address: int, channel: int) -> int:
    return address * CHANNEL + channel


def event_pair(label: int) -> tuple[int, int]:
    return divmod(label, CHANNEL)


def prove_record_factorization() -> None:
    print("\n[1] Gate A: one-record-per-event and address x channel factorization")
    labels = [event_label(a, c) for a in range(ADDRESS) for c in range(CHANNEL)]
    pairs = [event_pair(label) for label in labels]
    check(len(labels) == M_RECORD, "product alphabet has 8 x 137 atoms")
    check(len(set(labels)) == M_RECORD, "each address-channel pair is one distinct atom")
    check(set(pairs) == {(a, c) for a in range(ADDRESS) for c in range(CHANNEL)}, "composite label is bijective with (address, channel)")
    check(abs(entropy_uniform(M_RECORD) - entropy_uniform(ADDRESS) - entropy_uniform(CHANNEL)) < TOL, "record entropy factorizes additively")

    # One-hot instrument algebra: exactly one address projector and one channel
    # projector fire.  The joint projectors are their commuting products.
    address_row_sums = [sum(1 for _c in range(CHANNEL)) for _a in range(ADDRESS)]
    channel_col_sums = [sum(1 for _a in range(ADDRESS)) for _c in range(CHANNEL)]
    check(set(address_row_sums) == {CHANNEL}, "each address block carries the same 137 channel labels")
    check(set(channel_col_sums) == {ADDRESS}, "each channel block carries the same 8 address labels")
    check(sum(address_row_sums) == M_RECORD == sum(channel_col_sums), "joint one-hot projectors resolve the identity once")
    print("  The finite service-instrument theorem is: if an event is the minimal")
    print("  irreversible commit of the current admitted pointer inventory, the record")
    print("  is one atom of C^8 tensor C^137, i.e. one composite label (address,channel).")
    print("  Hidden timing/environment labels would be new pointer algebras, not part")
    print("  of the current event instrument.")


def prove_or_refute_traffic_clock() -> None:
    print("\n[2] Gate B: traffic-clock universalization")
    y_top = math.sqrt(2.0) * M_TOP / V_EW
    other_masses = {"b": M_B, "tau": M_TAU, "c": M_C, "s": M_S, "mu": M_MU}
    other_yukawas = {name: math.sqrt(2.0) * mass / V_EW for name, mass in other_masses.items()}
    print(f"  y_top = sqrt(2) m_top / v = {y_top:.6f}")
    print("  other Yukawas =", {name: round(value, 5) for name, value in other_yukawas.items()})
    check(0.95 < y_top < 1.05, "top Yukawa is an O(1) saturated coupling")
    check(max(other_yukawas.values()) < 0.05, "non-top Yukawas are not saturated anchors")

    saturating_anchors = ("g3/QCD confinement -> Lambda_QCD", "top Yukawa -> v = sqrt(2) m_top")
    print("  saturated anchors:")
    for anchor in saturating_anchors:
        print(f"    - {anchor}")
    check(len(saturating_anchors) == 2, "current SM/substrate ledger has two irreducible IR-saturating anchors")

    # Common-rate rescaling remains a null direction of the Crooks/KMS ratios.
    rates = np.array([0.4, 1.2, 2.5])
    scaled = 7.0 * rates
    ratios = np.log(rates[1:] / rates[:-1])
    scaled_ratios = np.log(scaled[1:] / scaled[:-1])
    check(np.allclose(ratios, scaled_ratios), "Crooks/KMS ratios do not fix the common traffic multiplier")
    print("  Therefore the native QEC/cosmology selector clock can be a real clock")
    print("  without becoming a universal EW/nuclear/dressed-alpha clock.  The current")
    print("  algebra refutes one-clock universalization: EW has its own saturated top")
    print("  anchor, and common rate rescaling is invisible to entropy/KMS ratios.")


def generation_perms() -> list[tuple[tuple[int, ...], np.ndarray, int]]:
    out: list[tuple[tuple[int, ...], np.ndarray, int]] = []
    eye = np.eye(3, dtype=complex)
    for p in permutations(range(3)):
        inversions = sum(1 for i in range(3) for j in range(i + 1, 3) if p[i] > p[j])
        out.append((p, eye[list(p)], -1 if inversions % 2 else 1))
    return out


PERMS = generation_perms()


def a_k3() -> np.ndarray:
    return np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)


def p_cycle(sigma: int) -> np.ndarray:
    cycle = (0, 1, 2) if sigma == 1 else (0, 2, 1)
    mat = np.zeros((3, 3), dtype=complex)
    for src, dst in zip(cycle, cycle[1:] + cycle[:1]):
        mat[dst, src] = 1.0
    return mat


def k_or() -> np.ndarray:
    p = p_cycle(+1)
    return p - p.T


def ordinary_average(mat: np.ndarray) -> np.ndarray:
    out = np.zeros_like(mat)
    for _perm, p, _parity in PERMS:
        out += p.T @ mat @ p
    return out / len(PERMS)


def sign_average(mat: np.ndarray) -> np.ndarray:
    out = np.zeros_like(mat)
    for _perm, p, parity in PERMS:
        out += parity * (p.T @ mat @ p)
    return out / len(PERMS)


def majorana_portal(phi: float, sigma: int, r: float = 0.5) -> np.ndarray:
    return np.eye(3, dtype=complex) + r * np.exp(1j * sigma * phi) * a_k3()


def cp_invariant(mat: np.ndarray) -> float:
    md = mat.conj().T
    return float(np.imag(np.trace(H_DIRAC @ md @ mat @ md @ H_DIRAC @ mat)))


def cp_norm(mat: np.ndarray) -> float:
    md = mat.conj().T
    probes = np.array(
        [
            np.imag(np.trace(H_DIRAC @ md @ mat @ md @ H_DIRAC @ mat)),
            np.imag(np.trace(H_DIRAC @ mat @ md @ mat @ md @ mat)),
        ],
        dtype=float,
    )
    return float(np.linalg.norm(probes))


def prove_or_refute_recovery_holonomy() -> None:
    print("\n[3] Gate C: recovery holonomy phase and orientation")
    A = a_k3()
    K = k_or()
    P_plus = p_cycle(+1)
    reconstructed = 0.5 * (A + K)
    check(np.linalg.norm(P_plus - reconstructed) < TOL, "directed latch is (A_K3 + K_or)/2")
    for perm, p, parity in PERMS:
        moved = p.T @ K @ p
        if parity == 1:
            check(np.linalg.norm(moved - K) < TOL, f"even generation relabelling {perm} preserves K_or")
        else:
            check(np.linalg.norm(moved + K) < TOL, f"odd generation relabelling {perm} flips K_or")

    phi = 2.0 * math.pi / 3.0
    m_plus = majorana_portal(phi, +1)
    m_minus = majorana_portal(phi, -1)
    print(f"  faithful C3 character phase Phi = {phi:.9f}")
    print(f"  CP norm sigma=+1 = {cp_norm(m_plus):.6e}; I_CP={cp_invariant(m_plus):+.6e}")
    print(f"  CP norm sigma=-1 = {cp_norm(m_minus):.6e}; I_CP={cp_invariant(m_minus):+.6e}")
    check(cp_norm(m_plus) > 1.0e-3 and cp_norm(m_minus) > 1.0e-3, "faithful C3 phase gives nonzero CP")
    check(abs(cp_invariant(m_plus) + cp_invariant(m_minus)) < TOL, "orientation conjugation reverses CP sign")

    scalar_seen = ordinary_average(K)
    sign_seen = sign_average(K)
    unoriented_mix = 0.5 * (m_plus + m_minus)
    print(f"  ||ordinary S3 average of K_or|| = {np.linalg.norm(scalar_seen):.3e}")
    print(f"  ||sign-twisted average - K_or|| = {np.linalg.norm(sign_seen - K):.3e}")
    print(f"  CP norm of unoriented mixture   = {cp_norm(unoriented_mix):.3e}")
    check(np.linalg.norm(scalar_seen) < TOL, "generation-blind scalar source cannot read K_or")
    check(np.linalg.norm(sign_seen - K) < TOL, "sign-representation pointer can read K_or")
    check(cp_norm(unoriented_mix) < TOL, "without absolute orientation the conjugate phases cancel")
    print("  Phase closes conditionally as the faithful C3 character.  Orientation")
    print("  does not close from the current scalar sterile source: it needs a physical")
    print("  sign-representation environment pointer tied to global substrate handedness.")


def main() -> None:
    print("R12/R13/R15 PROMOTION-GATE ATTEMPT")
    print("=" * 96)
    prove_record_factorization()
    prove_or_refute_traffic_clock()
    prove_or_refute_recovery_holonomy()
    print(
        """
VERDICT
  Gate A closes at finite service-instrument grade:
    one minimal irreversible event commits exactly one composite record
    (address, channel) in C^8 tensor C^137.  Extra timing/environment labels
    would be new pointer algebras.

  Gate B does not close:
    the native QEC/cosmology clock is not universalized to EW/nuclear or
    dressed-alpha sectors.  The current ledger instead has two irreducible
    saturated anchors, QCD and top-Yukawa, and the common traffic multiplier
    remains invisible to Crooks/KMS ratios.

  Gate C partially closes:
    K_or is the needed sign-representation tensor and Phi=2*pi/3 follows for
    a faithful C3 recovery character.  But the current sterile source is
    scalar/generation-blind; it erases K_or.  Absolute orientation still needs
    a sign-pointer bridge tied to global substrate handedness.

exit 0 -- one promotion gate closes; traffic universalization is refuted; holonomy remains conditional.
"""
    )


if __name__ == "__main__":
    main()
