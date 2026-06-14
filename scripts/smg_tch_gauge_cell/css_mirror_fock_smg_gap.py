#!/usr/bin/env python3
"""
Matter-axis mirror Fock SMG gap with coarse colour Gauss projection.

Why this exists:
  The gauge-sector ladder (Z3 centre -> S3 -> Peter-Weyl SU(3)) does not by
  itself compute the mirror gap. The mirror gap requires actual charged mirror
  fermions plus the SMG interaction.

This script is the first finite matter-axis check:

  * Use the existing 8-bit cell register as an actual 8-mode fermion Fock
    space. A basis state is an occupation bitstring

        |n_G0 n_G1 n_LQ n_C0 n_C1 n_I3 n_chi n_W>.

  * The mirror sector is W != chi, i.e. the R2-violating modes identified in
    mirror_invalid_subspace.py.

  * Build the CSS SMG interaction as genuine even fermion operators:
        Z stabilizers = products of occupation parities (-1)^n,
        X stabilizers = products of Jordan-Wigner Majoranas gamma_j=c_j+c_j^dag.
    The weight-4 and weight-8 X terms are even fermion interactions, not Pauli
    relabeling of a one-particle colour qutrit.

  * Project to the mirror Fock sector and compute the excitation gap.

  * Add the qutrit colour-orbit compensator/Gauss projection from the previous
    build and compute the gauge-invariant mirror excitation gap.

What this establishes:
  In the single-cell mirror Fock model, the SMG interaction gives a unique
  mirror ground state and a finite excitation gap even after coarse colour-orbit
  Gauss projection.

What this does not establish:
  This is not the full 3+1D chiral lattice gauge theory. It has no spatial
  mirror hopping, no dynamical open-link SU(3) Peter-Weyl fields, no sign/locality
  fermion measure problem, and no thermodynamic-limit test. It is the missing
  local matter-axis check, not the field-wide sufficiency theorem.

numpy; self-asserting.
"""

import itertools

import numpy as np


G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
GENERATORS = [
    [1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 1, 1, 1, 1],
    [0, 0, 1, 1, 0, 0, 1, 1],
    [0, 1, 0, 1, 0, 1, 0, 1],
]
N_MODES = 8
DIM = 2**N_MODES
OMEGA = np.exp(2j * np.pi / 3)
ALL_BITS = [tuple((state >> (N_MODES - 1 - bit)) & 1 for bit in range(N_MODES)) for state in range(DIM)]


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def idx(bits):
    return sum(bits[bit] << (N_MODES - 1 - bit) for bit in range(N_MODES))


def jw_majorana(bit):
    """gamma_bit = c_bit + c_bit^dag in Jordan-Wigner convention."""
    matrix = np.zeros((DIM, DIM), complex)
    for source, bits in enumerate(ALL_BITS):
        target = list(bits)
        target[bit] ^= 1
        sign = (-1) ** sum(bits[:bit])
        matrix[idx(target), source] = sign
    return matrix


def parity(bit):
    return np.diag([1 if bits[bit] == 0 else -1 for bits in ALL_BITS]).astype(complex)


def product(operators, dim=DIM):
    matrix = np.eye(dim, dtype=complex)
    for operator in operators:
        matrix = matrix @ operator
    return matrix


def css_fermion_stabilizers():
    gammas = [jw_majorana(bit) for bit in range(N_MODES)]
    parities = [parity(bit) for bit in range(N_MODES)]
    x_stabs = []
    z_stabs = []
    for row in GENERATORS:
        support = [bit for bit, value in enumerate(row) if value]
        x_stabs.append(product([gammas[bit] for bit in support]))
        z_stabs.append(product([parities[bit] for bit in support]))
    return z_stabs, x_stabs


def mirror_indices():
    return [index for index, bits in enumerate(ALL_BITS) if bits[CHI] != bits[W]]


def physical_indices():
    return [index for index, bits in enumerate(ALL_BITS) if bits[CHI] == bits[W]]


def unique_levels(eigvals, tol=1e-8):
    levels = []
    degeneracies = []
    for value in eigvals:
        if not levels or abs(value - levels[-1]) > tol:
            levels.append(float(value))
            degeneracies.append(1)
        else:
            degeneracies[-1] += 1
    return levels, degeneracies


def color_cycle_unitary():
    cycle = {
        (0, 0): (0, 0),
        (0, 1): (1, 0),
        (1, 0): (1, 1),
        (1, 1): (0, 1),
    }
    matrix = np.zeros((DIM, DIM), complex)
    for source, bits in enumerate(ALL_BITS):
        target = list(bits)
        target[C0], target[C1] = cycle[(bits[C0], bits[C1])]
        matrix[idx(target), source] = 1
    return matrix


def qutrit_ops():
    clock = np.diag([1, OMEGA, OMEGA**2]).astype(complex)
    shift = np.zeros((3, 3), complex)
    for state in range(3):
        shift[(state + 1) % 3, state] = 1
    return clock, shift


def adjoint_components(operator, unitary):
    powers = [np.linalg.matrix_power(unitary, power) for power in range(3)]
    components = []
    for charge in range(3):
        component = np.zeros_like(operator)
        for power in range(3):
            component += (
                OMEGA ** (-charge * power)
                * powers[power]
                @ operator
                @ powers[power].conj().T
            )
        components.append(component / 3)
    return components


def dress_by_color_orbit(operator, unitary):
    clock, _ = qutrit_ops()
    dressed = np.zeros((DIM * 3, DIM * 3), complex)
    for charge, component in enumerate(adjoint_components(operator, unitary)):
        dressed += np.kron(component, np.linalg.matrix_power(clock, charge))
    return dressed


def spectrum_report(label, hamiltonian):
    eigvals = np.linalg.eigvalsh(hamiltonian)
    levels, degeneracies = unique_levels(eigvals)
    gap = levels[1] - levels[0]
    print(
        f"  {label:<36s}: dim={hamiltonian.shape[0]:<4d} "
        f"E0={levels[0]:<8.6g} gap={gap:<8.6g} ground degeneracy={degeneracies[0]}"
    )
    return gap, degeneracies[0], levels


def full_and_mirror_fock_report():
    hd("A. Actual Fermion Occupation CSS/SMG Gap")
    z_stabs, x_stabs = css_fermion_stabilizers()
    all_stabs = z_stabs + x_stabs
    max_comm = max(float(np.linalg.norm(left @ right - right @ left)) for left in all_stabs for right in all_stabs)
    max_herm = max(float(np.linalg.norm(stab - stab.conj().T)) for stab in all_stabs)
    print(f"  max fermionic stabilizer commutator = {max_comm:.3e}")
    print(f"  max fermionic stabilizer Hermiticity = {max_herm:.3e}")
    assert max_comm < 1e-10
    assert max_herm < 1e-10

    h_full = -sum(all_stabs)
    full_gap, full_deg, _ = spectrum_report("full 8-mode Fock CSS", h_full)

    for label, indices in [
        ("mirror sector W!=chi", mirror_indices()),
        ("physical sector W=chi", physical_indices()),
    ]:
        h_sector = h_full[np.ix_(indices, indices)]
        gap, degeneracy, _ = spectrum_report(label, h_sector)
        projected_x_norms = [
            float(np.linalg.norm(x_stab[np.ix_(indices, indices)]))
            for x_stab in x_stabs
        ]
        print(f"    projected X-stabilizer norms: {np.round(projected_x_norms, 6)}")
        assert gap > 1.9
        assert degeneracy == 1

    assert full_gap > 1.9
    assert full_deg == 1
    print(
        "  -> In actual occupation language, the mirror sector is gapped by the "
        "fermionic CSS SMG interaction. The X stabilizer that flips W xor chi "
        "drops out under mirror projection, but the remaining terms still give "
        "a unique mirror ground state with gap 2."
    )
    return z_stabs, x_stabs


def gauge_invariant_mirror_report(z_stabs, x_stabs):
    hd("B. Coarse Colour-Orbit Gauss-Projected Mirror Gap")
    color_cycle = color_cycle_unitary()
    _, shift = qutrit_ops()
    total_cycle = np.kron(color_cycle, shift)

    dressed_stabs = [
        dress_by_color_orbit(stab, color_cycle)
        for stab in z_stabs + x_stabs
    ]
    max_comm = max(
        float(np.linalg.norm(left @ right - right @ left))
        for left in dressed_stabs
        for right in dressed_stabs
    )
    max_gauss_residual = max(
        float(np.linalg.norm(total_cycle @ stab @ total_cycle.conj().T - stab))
        for stab in dressed_stabs
    )
    print(f"  max dressed stabilizer commutator       = {max_comm:.3e}")
    print(f"  max dressed colour-orbit Gauss residual = {max_gauss_residual:.3e}")
    assert max_comm < 1e-10
    assert max_gauss_residual < 1e-10

    h_dressed = -sum(dressed_stabs)
    mirror_base = mirror_indices()
    mirror_qutrit_indices = [
        state * 3 + qutrit
        for state in mirror_base
        for qutrit in range(3)
    ]
    h_mirror = h_dressed[np.ix_(mirror_qutrit_indices, mirror_qutrit_indices)]
    projected_x_norms = [
        float(np.linalg.norm(stab[np.ix_(mirror_qutrit_indices, mirror_qutrit_indices)]))
        for stab in dressed_stabs[4:]
    ]
    print(f"  dressed mirror-sector X norms           = {np.round(projected_x_norms, 6)}")

    identity = np.eye(DIM * 3, dtype=complex)
    gauss_projector = (identity + total_cycle + total_cycle @ total_cycle) / 3
    pg_mirror = gauss_projector[np.ix_(mirror_qutrit_indices, mirror_qutrit_indices)]
    gauss_error = float(np.linalg.norm(pg_mirror @ pg_mirror - pg_mirror))
    comm_error = float(np.linalg.norm(pg_mirror @ h_mirror - h_mirror @ pg_mirror))
    gauss_rank = round(float(np.real(np.trace(pg_mirror))))
    print(f"  mirror colour-Gauss projector rank      = {gauss_rank}")
    print(f"  ||P_G^2-P_G||                           = {gauss_error:.3e}")
    print(f"  ||[P_G,H_mirror]||                      = {comm_error:.3e}")
    assert gauss_error < 1e-10
    assert comm_error < 1e-10

    eigenvalues, eigenvectors = np.linalg.eigh(pg_mirror)
    invariant_basis = eigenvectors[:, eigenvalues > 0.5]
    h_physical = invariant_basis.conj().T @ h_mirror @ invariant_basis
    gap, degeneracy, levels = spectrum_report("Gauss-invariant mirror Fock", h_physical)
    assert invariant_basis.shape[1] == 128
    assert gap > 1.9
    assert degeneracy == 1
    print(
        "  -> The coarse colour-orbit Gauss projection removes the qutrit gauge "
        "redundancy but does not close the local mirror SMG gap."
    )
    return gap, levels


def main():
    hd("Scope")
    print(
        "This is the matter-axis local check: actual 8-mode fermion occupation, "
        "mirror sector W!=chi, fermionic CSS SMG interaction, and coarse "
        "colour-orbit Gauss projection. It is not a full spatial/dynamical "
        "chiral lattice gauge theory."
    )
    z_stabs, x_stabs = full_and_mirror_fock_report()
    gauge_invariant_mirror_report(z_stabs, x_stabs)
    hd("Verdict")
    print(
        "The local mirror Fock gap is computed rather than assumed: the W!=chi "
        "mirror sector has a unique SMG ground state and gap 2, and the same "
        "gap survives the coarse colour-orbit Gauss projection."
    )
    print(
        "\nRemaining wall: put this charged mirror Fock/SMG block on a spatial "
        "lattice with dynamical open-link gauge fields and test the gap in the "
        "thermodynamic/strong-coupling regime. That is the genuine chiral "
        "lattice gauge / SMG problem."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
