#!/usr/bin/env python3
"""
Beyond-centre open-link colour proxy with vertex Gauss projection.

Why this exists:
  css_openlink_qutrit_su3_tch.py tests the Z3 centre/qutrit sector. That is
  necessary but not sufficient for full SU(3): off-centre colour dynamics could
  still reintroduce a small strong-coupling gap. This script adds the next
  finite, computable layer: the colour-permutation/Weyl subgroup S3.

Scope:
  This is NOT full SU(3). It is a finite non-abelian open-link gauge-Higgs
  proxy with explicit vertex Gauss projection. It includes off-diagonal colour
  permutations beyond the Z3 centre, and it computes the physical gap from the
  Hamiltonian spectrum rather than assigning the CSS value by hand. Passing this
  proxy is a useful sanity check, not a substitute for Peter-Weyl SU(3).

Gauge fixing:
  Four vertex colour qutrits are fixed to a reference colour by local S3 gauge
  transformations (unitary gauge). The residual gauge group at each vertex is
  H=S2, the stabilizer of that reference colour. Physical states are residual
  H^4 orbits of four open S3 link variables.

Hamiltonian on residual-Gauss orbit basis:
  H = H_electric + H_higgs + H_plaquette

  H_electric: class-sum hopping on each link over non-identity S3 elements.
  H_higgs:    -kappa if a link preserves the reference colour.
  H_plaquette:-beta Re Tr(U0 U1 U2^-1 U3^-1)/3.

The beta<=1 scan is the same confining-regime diagnostic used in the centre
proxy, now with non-centre colour permutations present.

numpy; self-asserting.
"""

import itertools

import numpy as np


EDGES = [
    (0, 1),
    (1, 2),
    (3, 2),
    (0, 3),
]


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


PERMS = list(itertools.permutations(range(3)))
P_INDEX = {perm: index for index, perm in enumerate(PERMS)}
IDENTITY = P_INDEX[(0, 1, 2)]
RESIDUAL_SWAP = P_INDEX[(0, 2, 1)]
RESIDUAL = [IDENTITY, RESIDUAL_SWAP]
NONIDENTITY = [index for index in range(len(PERMS)) if index != IDENTITY]


def compose(left, right):
    """Permutation composition left after right."""
    p = PERMS[left]
    q = PERMS[right]
    return P_INDEX[tuple(p[q[i]] for i in range(3))]


def inverse(perm_index):
    perm = PERMS[perm_index]
    inv = [0, 0, 0]
    for source, target in enumerate(perm):
        inv[target] = source
    return P_INDEX[tuple(inv)]


INV = [inverse(index) for index in range(len(PERMS))]


def transform_link_config(config, vertex_residuals):
    out = []
    for link, (tail, head) in enumerate(EDGES):
        transformed = compose(
            vertex_residuals[tail],
            compose(config[link], INV[vertex_residuals[head]]),
        )
        out.append(transformed)
    return tuple(out)


def orbit(seed):
    return {
        transform_link_config(seed, residuals)
        for residuals in itertools.product(RESIDUAL, repeat=4)
    }


def build_orbits():
    remaining = set(itertools.product(range(len(PERMS)), repeat=4))
    orbits = []
    config_to_orbit = {}
    while remaining:
        seed = next(iter(remaining))
        current = orbit(seed)
        orbit_index = len(orbits)
        sorted_current = sorted(current)
        orbits.append(sorted_current)
        for config in current:
            config_to_orbit[config] = orbit_index
        remaining -= current
    return orbits, config_to_orbit


def trace_norm(perm_index):
    perm = PERMS[perm_index]
    return sum(1 for i, value in enumerate(perm) if value == i) / 3


def preserves_reference_colour(perm_index):
    return 1.0 if PERMS[perm_index][0] == 0 else 0.0


def plaquette_perm(config):
    return compose(
        config[0],
        compose(config[1], compose(INV[config[2]], INV[config[3]])),
    )


def diagonal_energy(config, beta, kappa):
    higgs = -kappa * sum(preserves_reference_colour(link) for link in config)
    plaquette = -beta * trace_norm(plaquette_perm(config))
    return higgs + plaquette


def shifted_config(config, link, group_element):
    out = list(config)
    out[link] = compose(group_element, out[link])
    return tuple(out)


def hamiltonian(beta, kappa=1.0):
    if beta <= 0:
        raise ValueError("beta must be positive")
    g2 = 1 / beta
    orbits, config_to_orbit = build_orbits()
    dim = len(orbits)
    hamiltonian = np.zeros((dim, dim), complex)

    for src_index, configs in enumerate(orbits):
        orbit_size = len(configs)
        diag_values = {round(diagonal_energy(config, beta, kappa), 12) for config in configs}
        assert len(diag_values) == 1
        hamiltonian[src_index, src_index] += diagonal_energy(configs[0], beta, kappa)

        transition_counts = {}
        for config in configs:
            for link in range(4):
                for group_element in NONIDENTITY:
                    target = shifted_config(config, link, group_element)
                    target_index = config_to_orbit[target]
                    transition_counts[target_index] = transition_counts.get(target_index, 0) + 1

        for target_index, count in transition_counts.items():
            target_size = len(orbits[target_index])
            # Class-sum hopping over non-identity group elements on each link.
            coeff = -g2 / len(NONIDENTITY)
            hamiltonian[target_index, src_index] += (
                coeff * count / np.sqrt(orbit_size * target_size)
            )

    hermiticity_error = float(np.linalg.norm(hamiltonian - hamiltonian.conj().T))
    assert hermiticity_error < 1e-10
    return hamiltonian, orbits


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


def orbit_audit():
    hd("A. Residual Vertex-Gauss Orbit Audit")
    orbits, _ = build_orbits()
    sizes = sorted(len(current) for current in orbits)
    print(f"  raw unitary-gauge link configurations = {len(PERMS) ** 4}")
    print(f"  residual H^4 physical orbits          = {len(orbits)}")
    print(f"  orbit size distribution               = {dict((size, sizes.count(size)) for size in sorted(set(sizes)))}")
    assert sum(sizes) == len(PERMS) ** 4
    assert len(orbits) < len(PERMS) ** 4


def scan():
    hd("B. Beyond-Centre Strong-Coupling Spectrum")
    print("  beta       dim    gap       degeneracy   E0")
    gaps = []
    for beta in [0.25, 0.5, 0.75, 1.0]:
        h, orbits = hamiltonian(beta=beta, kappa=1.0)
        eigvals = np.linalg.eigvalsh(h)
        levels, degeneracies = unique_levels(eigvals)
        gap = levels[1] - levels[0]
        gaps.append(gap)
        print(
            f"  {beta:<10g} "
            f"{len(orbits):<6d} "
            f"{gap:<9.6g} "
            f"{degeneracies[0]:<12d} "
            f"{levels[0]:.6g}"
        )

    print(f"  minimum scanned gap = {min(gaps):.6g}")
    assert min(gaps) > 0.1
    return gaps


def main():
    hd("Scope")
    print(
        "This is a finite S3/Weyl open-link colour proxy, not full SU(3). "
        "It is meant to test whether adding non-centre colour permutations "
        "immediately reintroduces the beta<=1 gap collapse."
    )
    orbit_audit()
    scan()
    hd("Verdict")
    print(
        "The non-centre S3 colour-permutation proxy keeps a finite physical gap "
        "through beta<=1 after residual vertex-Gauss projection. This is a real "
        "computed spectrum, not a hardcoded CSS value."
    )
    print(
        "\nBound: S3 is only the Weyl/permutation subgroup. The decisive remaining "
        "test is still the Peter-Weyl SU(3) truncation with irreps such as "
        "(0,0), (1,0), and (0,1), including explicit intertwiners."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
