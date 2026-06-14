#!/usr/bin/env python3
"""
2D plaquette of charged mirror-Fock/SMG cells with Z3 open links and Lanczos.

This is the next step after css_spatial_mirror_fock_gauge_gap.py:

  * four spatial sites on one plaquette,
  * each site is a charged mirror-Fock/SMG block,
  * four open Z3 gauge links with exact Gauss-law basis enumeration,
  * electric energy plus magnetic plaquette shift,
  * gauge-covariant hopping of the projected fermion annihilation operators,
  * matrix-free Lanczos gap computation.

It is still a low-energy local truncation and still Z3 colour-orbit gauge, not
full SU(3). But it is now genuinely 2D and has a plaquette term, so it tests the
next obstruction beyond the 1D open-chain result.

numpy; self-asserting.
"""

import importlib.util
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
LOCAL_BLOCK_PATH = HERE / "css_spatial_mirror_fock_gauge_gap.py"
EDGES = [
    (0, 1),  # bottom
    (1, 2),  # right
    (3, 2),  # top, dagger in plaquette
    (0, 3),  # left, dagger in plaquette
]
PLAQUETTE_SHIFT = (1, 1, -1, -1)


def load_local_block():
    spec = importlib.util.spec_from_file_location("spatial_block", LOCAL_BLOCK_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


block = load_local_block()


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def electric_energy(flux):
    return 0.0 if flux % 3 == 0 else 1.0


def gauss_ok(charges, fluxes):
    for vertex in range(4):
        value = charges[vertex]
        for link, (tail, head) in enumerate(EDGES):
            if tail == vertex:
                value += fluxes[link]
            if head == vertex:
                value -= fluxes[link]
        if value % 3 != 0:
            return False
    return True


def build_basis(charges, n_sites=4):
    local_dim = len(charges)
    basis = []
    for config in np.ndindex(*(local_dim for _ in range(n_sites))):
        charge_tuple = tuple(int(charges[state]) for state in config)
        if sum(charge_tuple) % 3 != 0:
            continue
        for fluxes in np.ndindex(3, 3, 3, 3):
            if gauss_ok(charge_tuple, fluxes):
                basis.append((tuple(config), tuple(int(f) for f in fluxes)))
    index = {item: pos for pos, item in enumerate(basis)}
    return basis, index


def flux_shift(fluxes, link, amount):
    out = list(fluxes)
    out[link] = (out[link] + amount) % 3
    return tuple(out)


def plaquette_shift(fluxes, amount):
    return tuple((value + amount * shift) % 3 for value, shift in zip(fluxes, PLAQUETTE_SHIFT))


def transition_tables(annihilators):
    ann = []
    create = []
    for mode_components in annihilators:
        ann_by_charge = []
        create_by_charge = []
        for charge, matrix in enumerate(mode_components):
            ann_for_source = []
            create_for_source = []
            for source in range(matrix.shape[1]):
                ann_for_source.append(
                    [
                        (int(target), matrix[target, source])
                        for target in np.nonzero(np.abs(matrix[:, source]) > 1e-12)[0]
                    ]
                )
                dagger = matrix.conj().T
                create_for_source.append(
                    [
                        (int(target), dagger[target, source])
                        for target in np.nonzero(np.abs(dagger[:, source]) > 1e-12)[0]
                    ]
                )
            ann_by_charge.append(ann_for_source)
            create_by_charge.append(create_for_source)
        ann.append(ann_by_charge)
        create.append(create_by_charge)
    return ann, create


class PlaquetteHamiltonian:
    def __init__(self, beta, hopping, states_per_charge):
        self.beta = beta
        self.hopping = hopping
        self.g2 = 1 / beta
        self.energies, self.charges, self.annihilators, self.charge_counts = (
            block.local_mirror_dressed_block(states_per_charge)
        )
        self.basis, self.index = build_basis(self.charges)
        self.dim = len(self.basis)
        self.ann, self.create = transition_tables(self.annihilators)
        self.diagonal = np.zeros(self.dim)
        self.matter_diagonal = np.zeros(self.dim)
        for pos, (config, fluxes) in enumerate(self.basis):
            matter = sum(self.energies[state] for state in config)
            electric = self.g2 * sum(electric_energy(flux) for flux in fluxes)
            self.matter_diagonal[pos] = matter
            self.diagonal[pos] = matter + electric
        self.transitions = self.build_transitions()

    def build_transitions(self):
        transitions = []
        for col, (config, fluxes) in enumerate(self.basis):
            col_transitions = []

            # Magnetic plaquette term.
            for amount in [1, -1]:
                target = (config, plaquette_shift(fluxes, amount))
                row = self.index.get(target)
                if row is not None:
                    col_transitions.append((row, -(self.beta / 2)))

            # Gauge-covariant hopping on each oriented link.
            if abs(self.hopping) < 1e-15:
                transitions.append(col_transitions)
                continue
            config_list = list(config)
            for link, (tail, head) in enumerate(EDGES):
                old_tail = config_list[tail]
                old_head = config_list[head]
                for mode in range(len(self.ann)):
                    for charge in range(3):
                        # create at tail, annihilate at head, shift flux +charge.
                        for new_tail, amp_tail in self.create[mode][charge][old_tail]:
                            for new_head, amp_head in self.ann[mode][charge][old_head]:
                                new_config = list(config_list)
                                new_config[tail] = new_tail
                                new_config[head] = new_head
                                new_fluxes = flux_shift(fluxes, link, charge)
                                target = (tuple(new_config), new_fluxes)
                                row = self.index.get(target)
                                if row is not None:
                                    col_transitions.append(
                                        (row, -self.hopping * amp_tail * amp_head)
                                    )

                        # annihilate at tail, create at head, shift flux -charge.
                        for new_tail, amp_tail in self.ann[mode][charge][old_tail]:
                            for new_head, amp_head in self.create[mode][charge][old_head]:
                                new_config = list(config_list)
                                new_config[tail] = new_tail
                                new_config[head] = new_head
                                new_fluxes = flux_shift(fluxes, link, -charge)
                                target = (tuple(new_config), new_fluxes)
                                row = self.index.get(target)
                                if row is not None:
                                    col_transitions.append(
                                        (row, -self.hopping * amp_tail * amp_head)
                                    )
            transitions.append(col_transitions)
        return transitions

    def matvec(self, vector):
        out = self.diagonal * vector
        for col, col_transitions in enumerate(self.transitions):
            amp = vector[col]
            if abs(amp) < 1e-15:
                continue
            for row, coeff in col_transitions:
                out[row] += coeff * amp

        return out

    def matter_expectation(self, vector):
        return float(np.real(np.vdot(vector, self.matter_diagonal * vector)))


def lanczos_lowest(model, n_iter=80, n_eigs=8, seed=143):
    rng = np.random.default_rng(seed)
    q = rng.normal(size=model.dim) + 1j * rng.normal(size=model.dim)
    q = q / np.linalg.norm(q)
    q_prev = np.zeros_like(q)
    beta_prev = 0.0
    q_vectors = []
    alphas = []
    betas = []

    for _ in range(min(n_iter, model.dim)):
        q_vectors.append(q)
        v = model.matvec(q)
        alpha = float(np.real(np.vdot(q, v)))
        v = v - alpha * q - beta_prev * q_prev
        # Full reorthogonalization; dim is moderate and this stabilizes the low gap.
        for basis_vector in q_vectors:
            v -= np.vdot(basis_vector, v) * basis_vector
        beta = float(np.linalg.norm(v))
        alphas.append(alpha)
        if beta < 1e-12:
            break
        betas.append(beta)
        q_prev = q
        q = v / beta
        beta_prev = beta

    k = len(alphas)
    tridiag = np.diag(alphas)
    for index in range(k - 1):
        tridiag[index, index + 1] = betas[index]
        tridiag[index + 1, index] = betas[index]
    evals, evecs = np.linalg.eigh(tridiag)
    q_matrix = np.column_stack(q_vectors)
    ritz_vectors = q_matrix @ evecs[:, :n_eigs]
    return evals[:n_eigs], ritz_vectors


def exact_check_small(model):
    if model.dim > 1600:
        return None
    dense = np.column_stack([model.matvec(np.eye(model.dim, dtype=complex)[:, i]) for i in range(model.dim)])
    herm = float(np.linalg.norm(dense - dense.conj().T))
    eigvals = np.linalg.eigvalsh(dense)
    assert herm < 1e-8
    return eigvals[:6]


def format_gap(value):
    if np.isfinite(value):
        return f"{value:<11.6g}"
    return "not in win "


def scan():
    hd("A. 2D Plaquette Strong-Coupling Lanczos Scan")
    print("  n/q beta   t      dim    full gap  matter gap  exact check")
    results = []
    for states_per_charge in [2, 3]:
        for beta in [0.5, 1.0]:
            for hopping in [0.0, 0.1, 0.2]:
                model = PlaquetteHamiltonian(
                    beta=beta,
                    hopping=hopping,
                    states_per_charge=states_per_charge,
                )
                evals, vectors = lanczos_lowest(model, n_iter=90, n_eigs=10)
                full_gap = float(evals[1] - evals[0])
                ground_matter = model.matter_expectation(vectors[:, 0])
                matter_gap = float("inf")
                for index in range(1, vectors.shape[1]):
                    delta_matter = model.matter_expectation(vectors[:, index]) - ground_matter
                    if delta_matter > 0.5:
                        matter_gap = float(evals[index] - evals[0])
                        break
                exact = exact_check_small(model)
                exact_text = "no"
                if exact is not None:
                    exact_text = f"{exact[1]-exact[0]:.6g}"
                    full_gap = float(exact[1] - exact[0])
                results.append((states_per_charge, beta, hopping, full_gap, matter_gap))
                print(
                    f"  {states_per_charge:<3d} "
                    f"{beta:<6g} "
                    f"{hopping:<6g} "
                    f"{model.dim:<6d} "
                    f"{full_gap:<9.6g} "
                    f"{format_gap(matter_gap)} "
                    f"{exact_text}"
                )
    min_full = min(row[3] for row in results)
    min_matter = min(row[4] for row in results if np.isfinite(row[4]))
    print(f"  minimum full gap   = {min_full:.6g}")
    print(f"  minimum matter gap = {min_matter:.6g}")
    assert min_full > 0.1
    assert min_matter > 1.5


def larger_truncation_probe():
    hd("B. Larger Local-Truncation Probe")
    print("  one plaquette, states/charge=4, t=0.2")
    print("  beta   dim     full gap  matter gap")
    results = []
    for beta in [0.5, 1.0]:
        model = PlaquetteHamiltonian(beta=beta, hopping=0.2, states_per_charge=4)
        evals, vectors = lanczos_lowest(model, n_iter=80, n_eigs=12)
        full_gap = float(evals[1] - evals[0])
        ground_matter = model.matter_expectation(vectors[:, 0])
        matter_gap = float("inf")
        for index in range(1, vectors.shape[1]):
            delta_matter = model.matter_expectation(vectors[:, index]) - ground_matter
            if delta_matter > 0.5:
                matter_gap = float(evals[index] - evals[0])
                break
        results.append((full_gap, matter_gap))
        print(
            f"  {beta:<6g} "
            f"{model.dim:<7d} "
            f"{full_gap:<9.6g} "
            f"{format_gap(matter_gap)}"
        )
    finite_matter = [row[1] for row in results if np.isfinite(row[1])]
    assert min(row[0] for row in results) > 0.1
    assert finite_matter and min(finite_matter) > 1.5


def volume_scaling_note():
    hd("C. Finite-Volume Scaling Target")
    print(
        "  next 2D volume: a 2-plaquette strip has 6 sites and 7 links. "
        "With the first matter-resolving local truncation, states/charge=2, "
        "the exact Gauss-law basis has about 1.4e5 states."
    )
    print(
        "  next 3D volume: a one-cube Z3 proxy has 8 sites, 12 links, and 6 "
        "plaquettes. Even states/charge=1 gives about 5.3e5 physical states, "
        "but that truncation has no local matter excitation. A real matter-gap "
        "cube therefore needs a more optimized Krylov/matvec implementation."
    )


def main():
    hd("Scope")
    print(
        "This is a one-plaquette 2D Z3 gauge calculation with four charged "
        "mirror-Fock/SMG sites, exact Gauss-law basis enumeration, magnetic "
        "plaquette shifts, and matrix-free Lanczos. It is still a local-truncated "
        "proxy, not the full 3D thermodynamic problem."
    )
    scan()
    larger_truncation_probe()
    volume_scaling_note()
    hd("Verdict")
    print(
        "The 2D plaquette scan keeps the low-cutoff gap open for beta<=1 in the "
        "tested one-plaquette range. After the hopping-truncation audit, this "
        "should be read only as a finite truncated-model result: a convincing "
        "mirror-gap test must first show that the projected hopping visibly "
        "moves the gap as t becomes O(1)."
    )
    print(
        "\nNext escalation: implement a streamed transition table for the "
        "2-plaquette states/charge=2 strip, then reuse the same Krylov code for "
        "a minimal cube. At that point the bottleneck is sparse/matrix-free "
        "engineering, not representation bookkeeping."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
