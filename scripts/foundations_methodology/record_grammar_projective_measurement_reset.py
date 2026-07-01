#!/usr/bin/env python3
"""Projective measurement as finite detector record plus reset.

This script adds the explicit measurement/reset example to the record-grammar
ladder.

Setup
-----
There is one system qubit S and one finite detector register R.  The detector
starts in the blank state |0_R>.  A repeatable Z measurement is represented by
the reversible copy/premeasurement unitary

    U_meas |0_S 0_R> = |0_S 0_R>,
    U_meas |1_S 0_R> = |1_S 1_R>,

i.e. a CNOT from S to R.

For an input

    |psi> = sqrt(p)|0> + exp(i phi)sqrt(1-p)|1>,

the copy creates the correlated record

    sqrt(p)|00> + exp(i phi)sqrt(1-p)|11>.

Tracing or reading the detector gives the usual non-selective projective
measurement channel

    rho -> P0 rho P0 + P1 rho P1.

Record-grammar reading
----------------------
The reversible copy is not the thermodynamic arrow.  It is an isometry/unitary
correlation step.  The reusable detector must later be reset to the blank state.
If the detector is reset as a local memory without using side information, the
erased entropy is the Shannon entropy of its outcome distribution:

    H(p) = -p ln p - (1-p) ln(1-p) nats.

If the system-detector correlation is still coherently available, the copy can
be uncomputed reversibly; that is not a committed measurement record.  This is
the precise finite version of "measurement correlation is not the arrow; reset
of a reusable finite record is the arrow."
"""

from __future__ import annotations

import math

import numpy as np


I2 = np.eye(2, dtype=complex)
Z0 = np.array([[1, 0], [0, 0]], dtype=complex)
Z1 = np.array([[0, 0], [0, 1]], dtype=complex)


def ket(index: int, dim: int = 2) -> np.ndarray:
    v = np.zeros(dim, dtype=complex)
    v[index] = 1.0
    return v


def kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.kron(a, b)


def density(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, psi.conj())


def input_state(p: float, phi: float) -> np.ndarray:
    return math.sqrt(p) * ket(0) + np.exp(1j * phi) * math.sqrt(1.0 - p) * ket(1)


def cnot_system_to_detector() -> np.ndarray:
    """CNOT in basis |S R> = |00>, |01>, |10>, |11>."""

    u = np.zeros((4, 4), dtype=complex)
    # S=0: R unchanged.
    u[0, 0] = 1.0  # |00> -> |00>
    u[1, 1] = 1.0  # |01> -> |01>
    # S=1: R flips.
    u[3, 2] = 1.0  # |10> -> |11>
    u[2, 3] = 1.0  # |11> -> |10>
    return u


def partial_trace_two(rho: np.ndarray, keep: str) -> np.ndarray:
    tensor = rho.reshape(2, 2, 2, 2)  # row_S,row_R,col_S,col_R
    if keep == "S":
        return np.einsum("srtr->st", tensor)
    if keep == "R":
        return np.einsum("srsv->rv", tensor)
    raise ValueError("keep must be S or R")


def dephase_z(rho_s: np.ndarray) -> np.ndarray:
    return Z0 @ rho_s @ Z0 + Z1 @ rho_s @ Z1


def von_neumann_entropy_nats(rho: np.ndarray) -> float:
    evals = np.linalg.eigvalsh((rho + rho.conj().T) / 2.0)
    return float(-sum(v * math.log(v) for v in evals if v > 1e-14))


def shannon_nats(probs: list[float]) -> float:
    return -sum(p * math.log(p) for p in probs if p > 0)


def detector_readout_projector(outcome: int) -> np.ndarray:
    p_r = Z0 if outcome == 0 else Z1
    return kron(I2, p_r)


def conditional_system_state(rho_sr: np.ndarray, outcome: int) -> tuple[float, np.ndarray]:
    p_op = detector_readout_projector(outcome)
    prob = float(np.real_if_close(np.trace(p_op @ rho_sr)))
    post = p_op @ rho_sr @ p_op / prob
    return prob, partial_trace_two(post, "S")


def classical_record_state(p: float) -> np.ndarray:
    """p |00><00| + (1-p) |11><11| after record decoherence/commit."""

    rho = np.zeros((4, 4), dtype=complex)
    rho[0, 0] = p
    rho[3, 3] = 1.0 - p
    return rho


def local_reset_detector(rho_sr: np.ndarray) -> np.ndarray:
    """CP reset channel on R: |0><0| and |0><1| Kraus maps."""

    k0 = kron(I2, np.outer(ket(0), ket(0).conj()))
    k1 = kron(I2, np.outer(ket(0), ket(1).conj()))
    return k0 @ rho_sr @ k0.conj().T + k1 @ rho_sr @ k1.conj().T


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<62s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_matrix_close(name: str, value: np.ndarray, target: np.ndarray, tol: float = 1e-12) -> None:
    err = float(np.linalg.norm(value - target))
    print(f"  {name:<62s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def main() -> None:
    print("Projective measurement channel with finite detector register and reset")
    print("=" * 88)

    p = 0.37
    phi = 0.91
    psi_s = input_state(p, phi)
    rho_s = density(psi_s)
    blank_r = ket(0)
    psi_sr0 = np.kron(psi_s, blank_r)
    rho_sr0 = density(psi_sr0)
    u = cnot_system_to_detector()

    print("\n[1] The detector write is reversible premeasurement")
    assert_matrix_close("U_meas^dag U_meas = I", u.conj().T @ u, np.eye(4))
    rho_pre = u @ rho_sr0 @ u.conj().T
    rho_uncomputed = u.conj().T @ rho_pre @ u
    assert_matrix_close("U_meas^dag restores |psi>|blank>", rho_uncomputed, rho_sr0)
    assert_close("S(SR) after coherent copy", von_neumann_entropy_nats(rho_pre), 0.0)
    print("  -> copying the record coherently is reversible; it is not yet the thermodynamic arrow.")

    print("\n[2] Tracing the detector gives the projective measurement channel")
    rho_s_after = partial_trace_two(rho_pre, "S")
    assert_matrix_close("Tr_R[U rho U^dag] = P0 rho P0 + P1 rho P1", rho_s_after, dephase_z(rho_s))
    assert_close("off-diagonal before measurement |rho_01|", abs(rho_s[0, 1]), math.sqrt(p * (1.0 - p)))
    assert_close("off-diagonal after non-selective measurement", abs(rho_s_after[0, 1]), 0.0)

    print("\n[3] Selective detector readout gives the usual conditional records")
    p0, s0 = conditional_system_state(rho_pre, 0)
    p1, s1 = conditional_system_state(rho_pre, 1)
    assert_close("Pr(detector=0)", p0, p)
    assert_close("Pr(detector=1)", p1, 1.0 - p)
    assert_matrix_close("S | detector=0 -> |0><0|", s0, Z0)
    assert_matrix_close("S | detector=1 -> |1><1|", s1, Z1)

    print("\n[4] The committed classical record has Shannon entropy H(p)")
    rho_classical = classical_record_state(p)
    rho_r_classical = partial_trace_two(rho_classical, "R")
    h = shannon_nats([p, 1.0 - p])
    assert_close("H_shannon(p)", h, -p * math.log(p) - (1.0 - p) * math.log(1.0 - p))
    assert_close("S(detector register)", von_neumann_entropy_nats(rho_r_classical), h)
    assert_close("S(classical S,R joint)", von_neumann_entropy_nats(rho_classical), h)
    print("  -> the finite detector register holds exactly the outcome entropy to be cleared.")

    print("\n[5] Explicit reset channel clears the detector and exposes the Landauer ledger")
    rho_reset = local_reset_detector(rho_classical)
    detector_blank = partial_trace_two(rho_reset, "R")
    system_after_reset = partial_trace_two(rho_reset, "S")
    assert_matrix_close("detector after reset = |0><0|", detector_blank, Z0)
    assert_matrix_close("system after local reset remains dephased", system_after_reset, dephase_z(rho_s))
    assert_close("entropy erased from standalone detector register", von_neumann_entropy_nats(rho_r_classical), h)
    print(f"  standalone detector reset exports at least H(p) = {h:.12g} nats")
    print(f"  for p=1/2 this would be ln2 = {math.log(2.0):.12g} nats")

    print("\n[6] Side-information caveat: reset cost belongs to committed reusable memory")
    # If the coherent premeasurement has not been committed/read, U^dag uncomputes it at zero cost.
    assert_matrix_close("coherent uncopy leaves detector blank", partial_trace_two(rho_uncomputed, "R"), Z0)
    # If the classical correlation is still jointly accessible, a controlled uncopy can reset R while
    # preserving total entropy; this is erasure with side information, not standalone detector reset.
    rho_joint_uncopy = u.conj().T @ rho_classical @ u
    assert_matrix_close("classical correlated uncopy blanks detector", partial_trace_two(rho_joint_uncopy, "R"), Z0)
    assert_close("joint entropy preserved by uncopy", von_neumann_entropy_nats(rho_joint_uncopy), h)
    print("  -> a reusable detached detector pays H(p); a coherent or side-informed uncopy is not the same operation.")

    print(
        """
Verdict
-------
The projective-measurement record machine is explicit and finite:

  the detector write is a reversible CNOT/isometry;
  tracing or reading the detector gives rho -> P0 rho P0 + P1 rho P1;
  selective readout gives the usual conditional outcome records;
  the detector register carries H(p) nats of outcome entropy;
  local reset maps the detector back to blank and exports that entropy;
  the arrow is reset of reusable memory, not the mere existence of correlation.

This is the simplest executable version of the measurement service ledger.
"""
    )


if __name__ == "__main__":
    main()
