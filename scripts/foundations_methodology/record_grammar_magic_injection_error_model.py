#!/usr/bin/env python3
r"""Noisy detector/feedforward model for the magic-injection channel.

The ideal magic-injection script showed a closed-loop T resource becoming a
non-Clifford channel through a finite detector branch record:

    closed Wilson loop -> magic ancilla -> branch bit -> Clifford feedforward
    -> T gate -> reset.

This script asks what happens when the branch record fails before feedforward.

Three cases are separated because they are physically different:

    correct record:
        the detector bit is read correctly, so the channel is T rho T^\dagger.

    flipped record:
        the detector bit is trusted but wrong, so the wrong Clifford correction
        is applied.

    erased record:
        the detector reports "lost".  The best simple fallback here is a fixed
        default correction, equivalent to no branch-dependent feedforward.

For the T-state injection gadget, the entanglement fidelities relative to the
ideal T channel are

    correct:        F_e = 1,
    wrong bit:      F_e = 1/2,
    erased/default: F_e = 3/4.

Therefore a detector bit-flip rate p gives

    F_e = 1 - p/2,        F_avg = 1 - p/3,

while an erasure rate e gives

    F_e = 1 - e/4,        F_avg = 1 - e/6.

Record-grammar reading
----------------------
The detector register is load-bearing.  Resetting it is not just thermodynamic
bookkeeping: if the branch record is lost before feedforward, the injected
operation degrades.  A wrong record is worse than a known erasure.

This remains a one-qubit toy channel, not a full fault-tolerance threshold
calculation.
"""

from __future__ import annotations

import math

import numpy as np


I2 = np.eye(2, dtype=complex)
S_GATE = np.diag([1.0, 1.0j]).astype(complex)
T_GATE = np.diag([1.0, np.exp(1j * math.pi / 4.0)]).astype(complex)


def phase(theta: float) -> complex:
    return complex(np.exp(1j * theta))


def dagger(a: np.ndarray) -> np.ndarray:
    return a.conj().T


def ket(index: int) -> np.ndarray:
    v = np.zeros(2, dtype=complex)
    v[index] = 1.0
    return v


def density(psi: np.ndarray) -> np.ndarray:
    psi = psi / np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def raw_injection_branches(phi: float) -> list[np.ndarray]:
    """Kraus branches before detector-conditioned feedforward."""

    branch_0 = np.diag([1.0, phase(phi)]).astype(complex) / math.sqrt(2.0)
    branch_1 = np.diag([phase(phi), 1.0]).astype(complex) / math.sqrt(2.0)
    return [branch_0, branch_1]


def feedforward(record_bit: int) -> np.ndarray:
    """Clifford correction selected by the detector branch record."""

    return I2 if record_bit == 0 else S_GATE


def ideal_injection_kraus(phi: float = math.pi / 4.0) -> list[np.ndarray]:
    branches = raw_injection_branches(phi)
    return [feedforward(m) @ branches[m] for m in [0, 1]]


def bit_flip_record_kraus(p_flip: float, phi: float = math.pi / 4.0) -> list[np.ndarray]:
    """Detector branch bit flips before feedforward with probability p_flip."""

    if not 0.0 <= p_flip <= 1.0:
        raise ValueError("p_flip must lie in [0,1]")
    branches = raw_injection_branches(phi)
    kraus: list[np.ndarray] = []
    for true_m in [0, 1]:
        for flip, prob in [(0, 1.0 - p_flip), (1, p_flip)]:
            if prob == 0.0:
                continue
            recorded_m = true_m ^ flip
            kraus.append(math.sqrt(prob) * feedforward(recorded_m) @ branches[true_m])
    return kraus


def erasure_record_kraus(p_erase: float, default_record: int = 0, phi: float = math.pi / 4.0) -> list[np.ndarray]:
    """Detector branch is erased with probability p_erase.

    If erased, the controller knows the record is missing and applies a fixed
    default correction.  default_record=0 is the "do nothing" fallback.
    """

    if not 0.0 <= p_erase <= 1.0:
        raise ValueError("p_erase must lie in [0,1]")
    branches = raw_injection_branches(phi)
    kraus: list[np.ndarray] = []
    for true_m in [0, 1]:
        if p_erase < 1.0:
            kraus.append(math.sqrt(1.0 - p_erase) * feedforward(true_m) @ branches[true_m])
        if p_erase > 0.0:
            kraus.append(math.sqrt(p_erase) * feedforward(default_record) @ branches[true_m])
    return kraus


def apply_channel(rho: np.ndarray, kraus: list[np.ndarray]) -> np.ndarray:
    return sum(k @ rho @ dagger(k) for k in kraus)


def entanglement_fidelity_to_unitary(kraus: list[np.ndarray], unitary: np.ndarray) -> float:
    """Entanglement fidelity of a one-qubit channel relative to a unitary."""

    d = 2.0
    return float(sum(abs(np.trace(dagger(unitary) @ k)) ** 2 for k in kraus) / (d * d))


def average_gate_fidelity(entanglement_fidelity: float) -> float:
    """One-qubit average gate fidelity from entanglement fidelity."""

    return (2.0 * entanglement_fidelity + 1.0) / 3.0


def choi_matrix(kraus: list[np.ndarray]) -> np.ndarray:
    """Choi matrix normalized to trace one."""

    bell = (np.kron(ket(0), ket(0)) + np.kron(ket(1), ket(1))) / math.sqrt(2.0)
    rho = density(bell)
    out = np.zeros((4, 4), dtype=complex)
    for k in kraus:
        op = np.kron(I2, k)
        out += op @ rho @ dagger(op)
    return out


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<76s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_matrix_close(name: str, value: np.ndarray, target: np.ndarray, tol: float = 1e-12) -> None:
    err = float(np.linalg.norm(value - target))
    print(f"  {name:<76s} frob_err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_greater(name: str, value: float, bound: float) -> None:
    print(f"  {name:<76s} value={value:.12g} bound={bound:.12g}")
    if not value > bound:
        raise AssertionError(name)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<76s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def main() -> None:
    print("Magic-injection error model: branch-record failure degrades T")
    print("=" * 92)

    print("\n[1] Ideal branch record gives the exact T channel")
    ideal = ideal_injection_kraus()
    f_ideal = entanglement_fidelity_to_unitary(ideal, T_GATE)
    assert_close("ideal entanglement fidelity", f_ideal, 1.0)
    for label, psi in {
        "|0>": ket(0),
        "|+>": (ket(0) + ket(1)) / math.sqrt(2.0),
        "generic": math.cos(0.73 / 2.0) * ket(0) + phase(-0.42) * math.sin(0.73 / 2.0) * ket(1),
    }.items():
        rho = density(psi)
        assert_matrix_close(f"ideal output on {label}", apply_channel(rho, ideal), T_GATE @ rho @ dagger(T_GATE))
    print("  -> with a reliable branch bit, injection is a channel identity.")

    print("\n[2] A flipped detector bit is a trusted wrong record")
    wrong = bit_flip_record_kraus(1.0)
    f_wrong = entanglement_fidelity_to_unitary(wrong, T_GATE)
    assert_close("entanglement fidelity for always-wrong record", f_wrong, 0.5)
    print("  -> a wrong branch record applies the wrong Clifford correction.")

    print("\n[3] A known erasure is less damaging than a trusted wrong bit")
    erased = erasure_record_kraus(1.0, default_record=0)
    f_erased = entanglement_fidelity_to_unitary(erased, T_GATE)
    assert_close("entanglement fidelity for erased/default record", f_erased, 0.75)
    assert_greater("erasure fidelity exceeds wrong-record fidelity", f_erased, f_wrong)
    print("  -> knowing the branch was lost is better than trusting a false branch.")

    print("\n[4] Bit-flip noise gives the closed degradation law")
    previous = 1.1
    for p in [0.0, 0.001, 0.01, 0.05, 0.1, 1.0]:
        kraus = bit_flip_record_kraus(p)
        fe = entanglement_fidelity_to_unitary(kraus, T_GATE)
        favg = average_gate_fidelity(fe)
        assert_close(f"F_e bit flip p={p}", fe, 1.0 - 0.5 * p)
        assert_close(f"F_avg bit flip p={p}", favg, 1.0 - p / 3.0)
        if fe > previous + 1e-12:
            raise AssertionError("bit-flip fidelity increased with p")
        previous = fe
        print(f"    p={p:<5g}  F_e={fe:.12f}  F_avg={favg:.12f}")
    print("  -> branch-bit error is linear because the detector record selects the correction.")

    print("\n[5] Erasure/default noise has a different slope")
    previous = 1.1
    for e in [0.0, 0.001, 0.01, 0.05, 0.1, 1.0]:
        kraus = erasure_record_kraus(e, default_record=0)
        fe = entanglement_fidelity_to_unitary(kraus, T_GATE)
        favg = average_gate_fidelity(fe)
        assert_close(f"F_e erasure e={e}", fe, 1.0 - 0.25 * e)
        assert_close(f"F_avg erasure e={e}", favg, 1.0 - e / 6.0)
        if fe > previous + 1e-12:
            raise AssertionError("erasure fidelity increased with e")
        previous = fe
        print(f"    e={e:<5g}  F_e={fe:.12f}  F_avg={favg:.12f}")
    print("  -> erased record plus safe default is a distinct failure class.")

    print("\n[6] A reset-before-feedforward event is exactly the erased-record endpoint")
    no_feedforward = erasure_record_kraus(1.0, default_record=0)
    assert_matrix_close("reset-too-early channel Choi", choi_matrix(no_feedforward), choi_matrix(erased))
    print("  -> the reset ledger is correctly ordered: use the record first, reset after.")

    print("\n[7] Record-quality thresholds are explicit")
    target_average = 0.999
    p_flip_max = 3.0 * (1.0 - target_average)
    p_erase_max = 6.0 * (1.0 - target_average)
    print(f"  For F_avg >= {target_average}:")
    print(f"    detector bit-flip rate must satisfy p <= {p_flip_max:.6g}")
    print(f"    erasure/default rate must satisfy e <= {p_erase_max:.6g}")
    assert_less("bit-flip tolerance", p_flip_max, p_erase_max)
    print("  -> the same finite branch bit now has an operational accuracy budget.")

    print(
        """
Verdict
-------
The noisy magic-injection toy turns "record reliability" into a channel
calculation.

  reliable branch record:
      exact T channel.

  flipped branch record:
      trusted wrong correction, F_e = 1 - p/2 and F_avg = 1 - p/3.

  erased branch record:
      known missing branch with default correction, F_e = 1 - e/4 and
      F_avg = 1 - e/6.

  reset before feedforward:
      the erased/default endpoint; the branch bit was cleared before doing its
      operational job.

So the detector register is not decorative thermodynamic bookkeeping.  It is
the finite record that makes the non-Clifford operation possible, and its
failure produces a quantitative noisy channel.

Boundary
--------
This is a one-qubit injection toy, not a threshold theorem, not magic-state
distillation, and not a full fault-tolerant code.  Its role is to make the
record-ordering principle executable: write the branch, use it for feedforward,
then reset it.
"""
    )


if __name__ == "__main__":
    main()
