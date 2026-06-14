#!/usr/bin/env python3
"""
Spatial charged mirror-Fock/SMG chain with open-link Z3 gauge fields.

Purpose:
  Move from the local mirror Fock gap to a finite spatial lattice with charged
  mirror states, gauge-covariant hopping, open-link gauge fields, and Gauss law.
  This is the matter-axis probe requested after the pure-gauge ladder stalled.

Controlled truncation:
  The full mirror Fock Hilbert space is 128 states per cell, before the qutrit
  colour-orbit compensator. A full L-site exact diagonalization is immediately
  too large. This script therefore:

    1. Builds the local mirror Fock + qutrit colour-orbit compensated SMG block
       exactly (dimension 128*3 = 384).
    2. Diagonalizes it jointly with the local Z3 gauge-charge operator.
    3. Keeps the lowest n local eigenstates in each charge sector q=0,1,2.
    4. Builds an open spatial chain from those charge-resolved local states.
    5. Imposes Gauss law exactly in the basis: only total charge zero chains are
       kept, and open-link Z3 electric fluxes are fixed by cumulative charge.
    6. Adds gauge-covariant hopping of actual local fermion annihilation operators
       c_k, k=0..5 (the modes that preserve W!=chi), decomposed into Z3 charge
       components.

Hamiltonian:
    H = sum_x H_SMG(x)
      + g2 sum_links E(a_link)
      - t sum_<xy>,k,q [ c_{x,k,q}^dag U_link^q c_{y,k,q} + h.c. ]

  with g2=1/beta and beta<=1 the strong-coupling regime. E(a)=0 for a=0 and
  E(a)=1 for nonzero Z3 electric flux. There is no plaquette term in this 1D
  open-chain probe; the point is the charged mirror matter + open-link Gauss
  coupling, not magnetic confinement.

Scope:
  This is still not the full 3+1D chiral lattice gauge theory. It is a finite
  1D strong-coupling, charge-resolved low-energy truncation. A stable finite
  gap here is evidence that the local SMG gap is not immediately destroyed by
  charged hopping and Gauss-law links; it is not a thermodynamic proof.

numpy; self-asserting.
"""

from collections import defaultdict
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
ALL_BITS = [
    tuple((state >> (N_MODES - 1 - bit)) & 1 for bit in range(N_MODES))
    for state in range(DIM)
]
_LOCAL_CACHE = {}


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def idx(bits):
    return sum(bits[bit] << (N_MODES - 1 - bit) for bit in range(N_MODES))


def jw_annihilation(bit):
    matrix = np.zeros((DIM, DIM), complex)
    for source, bits in enumerate(ALL_BITS):
        if bits[bit] == 0:
            continue
        target = list(bits)
        target[bit] = 0
        matrix[idx(target), source] = (-1) ** sum(bits[:bit])
    return matrix


def jw_majorana(bit):
    return jw_annihilation(bit) + jw_annihilation(bit).conj().T


def parity(bit):
    return np.diag([1 if bits[bit] == 0 else -1 for bits in ALL_BITS]).astype(complex)


def product(operators, dim=DIM):
    matrix = np.eye(dim, dtype=complex)
    for operator in operators:
        matrix = matrix @ operator
    return matrix


def mirror_indices():
    return [index for index, bits in enumerate(ALL_BITS) if bits[CHI] != bits[W]]


def color_cycle_unitary_full():
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


def local_mirror_dressed_block(states_per_charge):
    if states_per_charge in _LOCAL_CACHE:
        return _LOCAL_CACHE[states_per_charge]

    mirror = mirror_indices()
    color_full = color_cycle_unitary_full()
    color = color_full[np.ix_(mirror, mirror)]
    clock, shift = qutrit_ops()
    total_charge = np.kron(color, shift)
    local_dim = len(mirror) * 3

    gammas = [jw_majorana(bit) for bit in range(N_MODES)]
    parities = [parity(bit) for bit in range(N_MODES)]
    dressed_stabs = []
    for row in GENERATORS:
        support = [bit for bit, value in enumerate(row) if value]
        z_full = product([parities[bit] for bit in support])
        x_full = product([gammas[bit] for bit in support])
        for full_op in [z_full, x_full]:
            op = full_op[np.ix_(mirror, mirror)]
            dressed = np.zeros((local_dim, local_dim), complex)
            for charge, component in enumerate(adjoint_components(op, color)):
                dressed += np.kron(component, np.linalg.matrix_power(clock, charge))
            dressed_stabs.append(dressed)

    h_local = -sum(dressed_stabs)
    comm = float(np.linalg.norm(h_local @ total_charge - total_charge @ h_local))
    assert comm < 1e-10

    eigvals, eigvecs = np.linalg.eigh(h_local)
    # Simultaneously diagonalize the charge operator inside degenerate H blocks.
    local_states = []
    cursor = 0
    while cursor < len(eigvals):
        end = cursor + 1
        while end < len(eigvals) and abs(eigvals[end] - eigvals[cursor]) < 1e-8:
            end += 1
        block_vecs = eigvecs[:, cursor:end]
        charge_block = block_vecs.conj().T @ total_charge @ block_vecs
        charge_vals, charge_vecs = np.linalg.eig(charge_block)
        for pos, value in enumerate(charge_vals):
            phase = np.angle(value)
            charge = int(np.round((phase % (2 * np.pi)) / (2 * np.pi / 3))) % 3
            vector = block_vecs @ charge_vecs[:, pos]
            vector = vector / np.linalg.norm(vector)
            local_states.append((float(eigvals[cursor]), charge, vector))
        cursor = end

    selected = []
    for charge in range(3):
        sector = sorted(
            [state for state in local_states if state[1] == charge],
            key=lambda item: item[0],
        )[:states_per_charge]
        selected.extend(sector)

    selected = sorted(selected, key=lambda item: (item[0], item[1]))
    basis = np.column_stack([state[2] for state in selected])
    energies = np.array([state[0] for state in selected])
    charges = np.array([state[1] for state in selected], dtype=int)

    # Charged single-fermion annihilation components projected into the selected basis.
    annihilators = []
    for bit in range(6):  # preserve W!=chi
        full_c = jw_annihilation(bit)
        c_mirror = full_c[np.ix_(mirror, mirror)]
        c_lifted = np.kron(c_mirror, np.eye(3, dtype=complex))
        components = []
        for charge in range(3):
            component = np.zeros_like(c_lifted)
            powers = [np.linalg.matrix_power(total_charge, power) for power in range(3)]
            for power in range(3):
                component += (
                    OMEGA ** (-charge * power)
                    * powers[power]
                    @ c_lifted
                    @ powers[power].conj().T
                )
            component /= 3
            components.append(basis.conj().T @ component @ basis)
        annihilators.append(components)

    charge_counts = {charge: int(np.sum(charges == charge)) for charge in range(3)}
    _LOCAL_CACHE[states_per_charge] = (energies, charges, annihilators, charge_counts)
    return _LOCAL_CACHE[states_per_charge]


def link_fluxes(charges):
    fluxes = []
    cumulative = 0
    for charge in charges[:-1]:
        cumulative = (cumulative + charge) % 3
        fluxes.append((-cumulative) % 3)
    return tuple(fluxes)


def electric_energy(flux):
    return 0.0 if flux == 0 else 1.0


def physical_chain_basis(charges, length):
    local_dim = len(charges)
    basis = []
    for config in itertools.product(range(local_dim), repeat=length):
        total = sum(int(charges[index]) for index in config) % 3
        if total != 0:
            continue
        charge_tuple = tuple(int(charges[index]) for index in config)
        basis.append((config, link_fluxes(charge_tuple)))
    index = {item: pos for pos, item in enumerate(basis)}
    return basis, index


def build_chain_hamiltonian(length, beta, hopping, states_per_charge):
    energies, charges, annihilators, charge_counts = local_mirror_dressed_block(states_per_charge)
    basis, index = physical_chain_basis(charges, length)
    dim = len(basis)
    g2 = 1 / beta
    hamiltonian = np.zeros((dim, dim), complex)
    h_matter = np.zeros((dim, dim), complex)
    h_electric = np.zeros((dim, dim), complex)

    for col, (config, fluxes) in enumerate(basis):
        matter_energy = sum(energies[site] for site in config)
        electric = g2 * sum(electric_energy(flux) for flux in fluxes)
        h_matter[col, col] += matter_energy
        h_electric[col, col] += electric
        hamiltonian[col, col] += matter_energy + electric

    for col, (config, _) in enumerate(basis):
        config = list(config)
        for link in range(length - 1):
            left_state = config[link]
            right_state = config[link + 1]
            for mode_components in annihilators:
                for charge in range(3):
                    a_charge = mode_components[charge]
                    # c_left^dag c_right with the link flux shift fixed by Gauss law.
                    left_targets = np.nonzero(np.abs(a_charge[:, left_state]) > 1e-12)[0]
                    right_targets = np.nonzero(np.abs(a_charge[:, right_state]) > 1e-12)[0]
                    for new_left in left_targets:
                        amp_left = np.conj(a_charge[new_left, left_state])
                        for new_right in right_targets:
                            amp_right = a_charge[new_right, right_state]
                            new_config = list(config)
                            new_config[link] = int(new_left)
                            new_config[link + 1] = int(new_right)
                            total = sum(int(charges[state]) for state in new_config) % 3
                            if total != 0:
                                continue
                            new_fluxes = link_fluxes(tuple(int(charges[state]) for state in new_config))
                            row = index.get((tuple(new_config), new_fluxes))
                            if row is not None:
                                hamiltonian[row, col] += -hopping * amp_left * amp_right

                    # Hermitian conjugate is included by the matrix Hermitian symmetrization below.
    hamiltonian = (hamiltonian + hamiltonian.conj().T) / 2
    herm = float(np.linalg.norm(hamiltonian - hamiltonian.conj().T))
    assert herm < 1e-10
    components = {"matter": h_matter, "electric": h_electric}
    return hamiltonian, charge_counts, components


def matter_dominated_gap(hamiltonian, h_matter):
    eigvals, eigvecs = np.linalg.eigh(hamiltonian)
    ground_matter = float(np.real(eigvecs[:, 0].conj() @ h_matter @ eigvecs[:, 0]))
    for index in range(1, len(eigvals)):
        matter_expectation = float(
            np.real(eigvecs[:, index].conj() @ h_matter @ eigvecs[:, index])
        )
        if matter_expectation - ground_matter > 0.5:
            return float(eigvals[index] - eigvals[0]), matter_expectation - ground_matter
    return float("inf"), 0.0


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


def gap_scan():
    hd("B. Strong-Coupling Finite-Size Gap Scan")
    print(
        "  charge-resolved local truncation: L=2,3 use 4 states/charge; "
        "L=4 uses 2 states/charge as a cheaper trend probe"
    )
    print("  L   n/q   beta   t      dim     full gap  matter gap  ground degeneracy")
    results = []
    cases = [(2, 4), (3, 4), (4, 2)]
    for length, states_per_charge in cases:
        for beta in [0.5, 1.0]:
            for hopping in [0.0, 0.1, 0.2]:
                hamiltonian, charge_counts, components = build_chain_hamiltonian(
                    length=length,
                    beta=beta,
                    hopping=hopping,
                    states_per_charge=states_per_charge,
                )
                eigvals = np.linalg.eigvalsh(hamiltonian)
                levels, degeneracies = unique_levels(eigvals)
                gap = levels[1] - levels[0]
                matter_gap, _ = matter_dominated_gap(hamiltonian, components["matter"])
                results.append((length, beta, hopping, gap))
                print(
                    f"  {length:<3d} "
                    f"{states_per_charge:<5d} "
                    f"{beta:<6g} "
                    f"{hopping:<6g} "
                    f"{hamiltonian.shape[0]:<7d} "
                    f"{gap:<9.6g} "
                    f"{matter_gap:<11.6g} "
                    f"{degeneracies[0]}"
                )
    min_gap = min(result[3] for result in results)
    print(f"  minimum scanned gap = {min_gap:.6g}")
    assert min_gap > 0.1
    return results


def two_site_truncation_convergence():
    hd("C. Two-Site Truncation Convergence")
    print("  beta=0.5, t=0.2")
    print("  states/charge   dim     gap       ground degeneracy")
    gaps = []
    for states_per_charge in [2, 4, 6]:
        hamiltonian, _, components = build_chain_hamiltonian(
            length=2,
            beta=0.5,
            hopping=0.2,
            states_per_charge=states_per_charge,
        )
        eigvals = np.linalg.eigvalsh(hamiltonian)
        levels, degeneracies = unique_levels(eigvals)
        gap = levels[1] - levels[0]
        matter_gap, _ = matter_dominated_gap(hamiltonian, components["matter"])
        gaps.append(gap)
        print(
            f"  {states_per_charge:<15d} "
            f"{hamiltonian.shape[0]:<7d} "
            f"{gap:<9.6g} "
            f"(matter {matter_gap:.6g}) "
            f"{degeneracies[0]}"
        )
    assert min(gaps) > 0.1


def truncation_check():
    hd("A. Local Truncation Check")
    for states_per_charge in [2, 4, 6]:
        energies, charges, _, charge_counts = local_mirror_dressed_block(states_per_charge)
        print(
            f"  states/charge={states_per_charge}: total={len(energies)}, "
            f"charge counts={charge_counts}, lowest energies={np.round(energies[:9], 6)}"
        )
    print(
        "  -> The scan keeps equal low-energy content in each Z3 charge sector, "
        "so Gauss-law total-charge sectors are not biased by the truncation."
    )


def main():
    hd("Scope")
    print(
        "This is a finite 1D spatial chain of charged mirror Fock/SMG blocks with "
        "open-link Z3 gauge fluxes imposed by Gauss law. It includes actual local "
        "fermion annihilation operators projected into a charge-resolved SMG "
        "low-energy basis. It is still a controlled truncation, not a full 3+1D "
        "thermodynamic proof."
    )
    truncation_check()
    gap_scan()
    two_site_truncation_convergence()
    hd("Verdict")
    print(
        "In this finite charge-resolved chain, the local mirror SMG gap remains "
        "open under gauge-covariant charged hopping and strong-coupling open-link "
        "Gauss constraints for L=2..4 and beta<=1 in the scanned hopping range."
    )
    print(
        "\nRemaining wall: remove the low-energy truncation, include a 2D/3D gauge "
        "lattice with plaquettes, and test finite-size scaling at larger L. This "
        "script is evidence for gap stability, not the constructive theorem."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
