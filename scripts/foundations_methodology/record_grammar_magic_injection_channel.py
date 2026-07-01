#!/usr/bin/env python3
r"""Magic-injection channel with a finite detector register.

The previous examples showed that a non-Clifford T phase can be carried by a
closed Wilson loop.  This script takes the next step: use that loop phase as a
finite magic resource that acts on a data qubit through a measurement channel.

Protocol
--------
Prepare a route/ancilla qubit from a closed Wilson loop

    |A_phi> = (|0> + exp(i phi)|1>)/sqrt(2).

For phi = pi/4 this is the usual T magic state.  Couple a data qubit to the
ancilla with CNOT(data -> ancilla), measure the ancilla in Z, write the result
to a one-bit detector register, and apply the Clifford correction S on the data
when the detector reads 1.

The two data Kraus branches are

    K_0 = diag(1, exp(i phi)) / sqrt(2),
    K_1 = S diag(exp(i phi), 1) / sqrt(2).

At phi = pi/4 both branches implement the same T gate up to branch weight and
global phase:

    K_0 rho K_0^dagger + K_1 rho K_1^dagger = T rho T^dagger.

Record-grammar reading
----------------------
This is the smallest executable example of a magic resource becoming an
operation:

    closed Wilson loop -> magic ancilla -> detector outcome -> feedforward
    -> reset of a finite detector register.

The detector register is not decorative.  Without the outcome record and the
feedforward, the channel is not the T gate.  A stabilizer replacement for the
magic resource also fails by the familiar 14.6% process-fidelity gap.

This is a toy channel, not a fault-tolerant distillation protocol.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


I2 = np.eye(2, dtype=complex)
S_GATE = np.diag([1.0, 1.0j]).astype(complex)
T_GATE = np.diag([1.0, np.exp(1j * math.pi / 4.0)]).astype(complex)


def ket(index: int) -> np.ndarray:
    v = np.zeros(2, dtype=complex)
    v[index] = 1.0
    return v


def density(psi: np.ndarray) -> np.ndarray:
    psi = psi / np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def phase(theta: float) -> complex:
    return complex(np.exp(1j * theta))


def dagger(a: np.ndarray) -> np.ndarray:
    return a.conj().T


def kron(*ops: np.ndarray) -> np.ndarray:
    out = np.array([[1.0]], dtype=complex)
    for op in ops:
        out = np.kron(out, op)
    return out


def magic_resource(phi: float) -> np.ndarray:
    """Wilson-loop route state (|0> + e^{i phi}|1>)/sqrt(2)."""

    return (ket(0) + phase(phi) * ket(1)) / math.sqrt(2.0)


def arbitrary_state(theta: float, varphi: float) -> np.ndarray:
    """A generic one-qubit state for checking the channel away from axes."""

    return math.cos(theta / 2.0) * ket(0) + phase(varphi) * math.sin(theta / 2.0) * ket(1)


@dataclass(frozen=True)
class DiamondLinks:
    """Four U(1) link phases whose closed loop supplies the resource phase."""

    sa: complex
    ad: complex
    sb: complex
    bd: complex

    @classmethod
    def from_angles(cls, theta_sa: float, theta_ad: float, theta_sb: float, theta_bd: float) -> "DiamondLinks":
        return cls(phase(theta_sa), phase(theta_ad), phase(theta_sb), phase(theta_bd))

    def upper(self) -> complex:
        return self.sa * self.ad

    def lower(self) -> complex:
        return self.sb * self.bd

    def loop(self) -> complex:
        return self.upper() * np.conj(self.lower())

    def flux(self) -> float:
        return float(np.angle(self.loop()))


def transform_link(u_ij: complex, lambda_i: float, lambda_j: float) -> complex:
    return complex(np.exp(1j * (lambda_i - lambda_j)) * u_ij)


def gauge_transform(links: DiamondLinks, lambdas: dict[str, float]) -> DiamondLinks:
    return DiamondLinks(
        sa=transform_link(links.sa, lambdas["S"], lambdas["A"]),
        ad=transform_link(links.ad, lambdas["A"], lambdas["D"]),
        sb=transform_link(links.sb, lambdas["S"], lambdas["B"]),
        bd=transform_link(links.bd, lambdas["B"], lambdas["D"]),
    )


def injection_kraus(phi: float, feedforward: bool = True) -> list[np.ndarray]:
    """Data-only Kraus operators after ancilla Z-measurement.

    If feedforward=True, branch 1 applies the Clifford correction S.  For the
    T resource this makes both branches the same logical T operation up to an
    irrelevant global phase.
    """

    u_phi = np.diag([1.0, phase(phi)]).astype(complex)
    branch_1_raw = np.diag([phase(phi), 1.0]).astype(complex)
    k0 = u_phi / math.sqrt(2.0)
    k1 = branch_1_raw / math.sqrt(2.0)
    if feedforward:
        k1 = S_GATE @ k1
    return [k0, k1]


def apply_channel(rho: np.ndarray, kraus: list[np.ndarray]) -> np.ndarray:
    return sum(k @ rho @ dagger(k) for k in kraus)


def branch_probabilities(rho: np.ndarray, kraus: list[np.ndarray]) -> list[float]:
    return [float(np.real_if_close(np.trace(k @ rho @ dagger(k)))) for k in kraus]


def detector_joint_after_feedforward(rho: np.ndarray, phi: float) -> np.ndarray:
    """Data + one-bit detector register after measurement and feedforward.

    The detector register stores the branch label.  The output is block
    diagonal in that finite record basis.
    """

    p0 = density(ket(0))
    p1 = density(ket(1))
    k0, k1 = injection_kraus(phi, feedforward=True)
    return kron(k0 @ rho @ dagger(k0), p0) + kron(k1 @ rho @ dagger(k1), p1)


def reset_detector_register(joint: np.ndarray) -> np.ndarray:
    """Reset the one-bit detector register to |0>.

    This is a many-to-one operation on the register.  In this protocol the
    branch bit is fair and independent of the corrected data, so the detached
    reset ledger is ln 2 nats.
    """

    r0 = np.outer(ket(0), ket(0).conj())
    r1 = np.outer(ket(0), ket(1).conj())
    k0 = kron(I2, r0)
    k1 = kron(I2, r1)
    return k0 @ joint @ dagger(k0) + k1 @ joint @ dagger(k1)


def shannon_entropy_nats(probs: list[float]) -> float:
    return -sum(p * math.log(p) for p in probs if p > 0.0)


def entanglement_fidelity_to_unitary(kraus: list[np.ndarray], unitary: np.ndarray) -> float:
    """Entanglement fidelity of a channel relative to a one-qubit unitary."""

    d = 2.0
    return float(sum(abs(np.trace(dagger(unitary) @ k)) ** 2 for k in kraus) / (d * d))


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


def main() -> None:
    print("Magic-injection channel: closed holonomy -> T gate -> detector reset")
    print("=" * 92)

    phi_t = math.pi / 4.0
    links = DiamondLinks.from_angles(theta_sa=phi_t, theta_ad=0.0, theta_sb=0.0, theta_bd=0.0)
    redistributed = DiamondLinks.from_angles(phi_t / 4.0, phi_t / 4.0, -phi_t / 4.0, -phi_t / 4.0)
    transformed = gauge_transform(links, {"S": 0.37, "A": -1.11, "B": 2.04, "D": -0.58})

    print("\n[1] The resource is the closed loop, not a particular open-link phase")
    assert_close("arg original W", links.flux(), phi_t)
    assert_close("arg redistributed W", redistributed.flux(), phi_t)
    assert_close("arg gauge-transformed W", transformed.flux(), phi_t)
    print("  -> different link gauges with the same closed holonomy define the same resource.")

    print("\n[2] The T-resource injection channel equals the T unitary")
    test_states = {
        "|0>": ket(0),
        "|1>": ket(1),
        "|+>": (ket(0) + ket(1)) / math.sqrt(2.0),
        "|+i>": (ket(0) + 1.0j * ket(1)) / math.sqrt(2.0),
        "generic": arbitrary_state(theta=0.73, varphi=-0.42),
    }
    kraus_t = injection_kraus(links.flux(), feedforward=True)
    for label, psi in test_states.items():
        rho = density(psi)
        actual = apply_channel(rho, kraus_t)
        target = T_GATE @ rho @ dagger(T_GATE)
        assert_matrix_close(f"T injection on {label}", actual, target)
        probs = branch_probabilities(rho, kraus_t)
        assert_close(f"branch-0 probability for {label}", probs[0], 0.5)
        assert_close(f"branch-1 probability for {label}", probs[1], 0.5)
    print("  -> the detector branch is fair, but after feedforward both branches are the same T gate.")

    print("\n[3] The detector register is explicit and must be reset")
    rho_plus = density(test_states["|+>"])
    joint = detector_joint_after_feedforward(rho_plus, links.flux())
    target_data = T_GATE @ rho_plus @ dagger(T_GATE)
    expected_joint = 0.5 * kron(target_data, density(ket(0))) + 0.5 * kron(target_data, density(ket(1)))
    assert_matrix_close("data + detector branch record", joint, expected_joint)
    entropy = shannon_entropy_nats([0.5, 0.5])
    assert_close("detached branch-record entropy", entropy, math.log(2.0))
    reset = reset_detector_register(joint)
    assert_matrix_close("after reset, detector is blank", reset, kron(target_data, density(ket(0))))
    print("  -> the finite detector bit is the record that makes feedforward possible; clearing it costs ln 2 nats.")

    print("\n[4] Without feedforward the same resource does not implement T")
    no_feed = injection_kraus(phi_t, feedforward=False)
    f_no_feed = entanglement_fidelity_to_unitary(no_feed, T_GATE)
    assert_close("entanglement fidelity without feedforward", f_no_feed, 0.75)
    print("  -> the outcome record is operational, not decorative.")

    print("\n[5] Stabilizer replacements cannot inject the T operation")
    for label, phi in {"|+> resource": 0.0, "|+i> resource": math.pi / 2.0}.items():
        f = entanglement_fidelity_to_unitary(injection_kraus(phi, feedforward=True), T_GATE)
        print(f"  {label:<22s} entanglement fidelity to T = {f:.12g}")
        assert_close(f"{label} fidelity", f, (2.0 + math.sqrt(2.0)) / 4.0)
    f_magic = entanglement_fidelity_to_unitary(kraus_t, T_GATE)
    assert_close("T-resource fidelity", f_magic, 1.0)
    print("  -> replacing the holonomy magic by a stabilizer resource loses 14.6% in process fidelity.")

    print("\n[6] The protocol is a channel identity, not just a state match")
    target_kraus = [T_GATE]
    f_redist = entanglement_fidelity_to_unitary(injection_kraus(redistributed.flux(), feedforward=True), T_GATE)
    f_gauge = entanglement_fidelity_to_unitary(injection_kraus(transformed.flux(), feedforward=True), T_GATE)
    assert_close("redistributed-loop process fidelity", f_redist, 1.0)
    assert_close("gauge-transformed-loop process fidelity", f_gauge, 1.0)
    assert_close("target unitary self-fidelity", entanglement_fidelity_to_unitary(target_kraus, T_GATE), 1.0)
    print("  -> the injected operation depends on the closed loop class, not the open-link placement.")

    print(
        """
Verdict
-------
The magic-injection channel gives the first finite operational use of the
holonomy-carried T resource:

  closed Wilson loop W = exp(i pi/4)
    -> magic ancilla
    -> detector branch record
    -> Clifford feedforward
    -> exact T channel
    -> reset of one finite detector bit.

This is the record-grammar version of "magic is a resource": it is not merely a
phase that can be read; it can be consumed by a finite measurement channel to
produce a non-Clifford operation.  The branch record is necessary, and a
stabilizer replacement leaves the characteristic 14.6% fidelity gap.

Boundary
--------
This is not a fault-tolerant distillation protocol, not a many-qubit magic
resource theory, and not a derivation of gauge dynamics.  It is the smallest
self-checking channel identity linking closed holonomy, detector readout,
feedforward, and reset.
"""
    )


if __name__ == "__main__":
    main()
