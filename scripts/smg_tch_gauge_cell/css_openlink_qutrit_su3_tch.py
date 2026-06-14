#!/usr/bin/env python3
"""
Open-link qutrit-colour Wilson/TCH Hamiltonian with vertex Gauss projection.

This is the next executable probe after css_qutrit_color_build.py.

Important scope:
  A literal Kogut-Susskind SU(3) link has the infinite Hilbert space L^2(SU(3)).
  This script does not claim to construct that full object. It constructs the
  first open-link truncation that directly targets the obstruction we found:
  the SU(3) Z3 centre / qutrit-colour sector, with explicit vertex Gauss
  projection and a Wilson-like plaquette Hamiltonian.

Lattice:
  One oriented plaquette with vertices 0,1,2,3 and links

      e0: 0 -> 1
      e1: 1 -> 2
      e2: 3 -> 2        (plaquette uses e2^dag)
      e3: 0 -> 3        (plaquette uses e3^dag)

Degrees of freedom:
  * Matter-colour qutrit m_v in Z3 at each vertex. This is the qutrit colour
    orbit layer from css_qutrit_color_build.py, not a full quark field.
  * Open link qutrit a_l in Z3 on each link.

Vertex Gauss action:
      m_v -> m_v + h_v,
      a_l -> a_l + h_tail(l) - h_head(l).

Gauge-invariant open-link variables:
      b_l = m_head(l) - m_tail(l) + a_l  mod 3.

Hamiltonian in the Gauss-projected basis:
      H = -g2/2 sum_l (T_l + T_l^dag)
          -kappa sum_l cos(2pi b_l/3)
          -beta cos(2pi (b0+b1-b2-b3)/3)
          -eta/2 sum_v (M_v + M_v^dag)

  T_l shifts b_l. M_v is the gauge-invariant remnant of matter-colour orbit
  kinetic motion; it shifts the covariant b variables adjacent to v. eta=0 is
  the most conservative CSS-like choice, while eta>0 tests whether colour-orbit
  mobility destabilizes the strong-coupling gap.

What it tests:
  The old Pauli SU(3) shadow had a no-flip matter gap 0.52 at beta=0.5 and
  1.25 at beta=1. This script asks whether the qutrit-colour open-link,
  Gauss-projected Hamiltonian has a finite physical gap in that same confining
  beta<=1 regime. That physical gap is the qutrit colour-centre/gauge-Higgs
  gap, not yet an independent full mirror-fermion CSS gap.

What it still does not test:
  Full SU(3) non-centre colour generators, Peter-Weyl irrep mixing, fermionic
  mirror hopping, or locality of a non-abelian fermion measure. Passing this
  script is a necessary centre-sector sanity check, not the full chiral gauge
  construction.

numpy; self-asserting.
"""

import itertools

import numpy as np


OMEGA = np.exp(2j * np.pi / 3)
EDGES = [
    (0, 1),  # e0
    (1, 2),  # e1
    (3, 2),  # e2, dagger in plaquette
    (0, 3),  # e3, dagger in plaquette
]
ORIENTATION = [1, 1, -1, -1]


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def cos_z3(value):
    return float(np.cos(2 * np.pi * (value % 3) / 3))


def all_configs():
    return list(itertools.product(range(3), repeat=8))


def split_config(config):
    return tuple(config[:4]), tuple(config[4:])


def covariant_links(matter, links):
    out = []
    for link, (tail, head) in enumerate(EDGES):
        out.append((matter[head] - matter[tail] + links[link]) % 3)
    return tuple(out)


def plaquette_flux(b_links):
    return sum(sign * value for sign, value in zip(ORIENTATION, b_links)) % 3


def gauss_transform(config, shifts):
    matter, links = split_config(config)
    new_matter = tuple((matter[v] + shifts[v]) % 3 for v in range(4))
    new_links = []
    for value, (tail, head) in zip(links, EDGES):
        new_links.append((value + shifts[tail] - shifts[head]) % 3)
    return new_matter + tuple(new_links)


def canonical_orbit_label(config):
    orbit = [
        gauss_transform(config, shifts)
        for shifts in itertools.product(range(3), repeat=4)
    ]
    return min(orbit)


def orbit_audit():
    hd("A. Explicit Vertex Gauss-Orbit Audit")
    configs = all_configs()
    labels = {canonical_orbit_label(config) for config in configs}
    b_labels = {covariant_links(*split_config(config)) for config in configs}
    print(f"  raw matter+link configurations = {len(configs)}")
    print(f"  Gauss orbits                  = {len(labels)}")
    print(f"  covariant b_l labels          = {len(b_labels)}")
    print(f"  expected qutrit-link basis    = {3**4}")

    label_to_b = {}
    for config in configs:
        label = canonical_orbit_label(config)
        b_links = covariant_links(*split_config(config))
        if label in label_to_b:
            assert label_to_b[label] == b_links
        else:
            label_to_b[label] = b_links

    assert len(labels) == 3**4
    assert len(b_labels) == 3**4
    print(
        "  -> The physical Gauss-projected basis is exactly the four open-link "
        "covariant qutrit variables b_l=m_head-m_tail+a_l."
    )


def b_basis():
    basis = list(itertools.product(range(3), repeat=4))
    index = {state: i for i, state in enumerate(basis)}
    return basis, index


def shift_b(state, link, amount):
    out = list(state)
    out[link] = (out[link] + amount) % 3
    return tuple(out)


def matter_shift_b(state, vertex, amount):
    out = list(state)
    for link, (tail, head) in enumerate(EDGES):
        if head == vertex:
            out[link] = (out[link] + amount) % 3
        if tail == vertex:
            out[link] = (out[link] - amount) % 3
    return tuple(out)


def hamiltonian(beta, kappa=1.0, eta=0.0):
    if beta <= 0:
        raise ValueError("Use beta > 0 so g2=1/beta is finite.")
    g2 = 1 / beta
    basis, index = b_basis()
    dim = len(basis)
    h = np.zeros((dim, dim), complex)

    for col, state in enumerate(basis):
        higgs_energy = -kappa * sum(cos_z3(value) for value in state)
        magnetic_energy = -beta * cos_z3(plaquette_flux(state))
        h[col, col] += higgs_energy + magnetic_energy

        for link in range(4):
            for amount in [1, -1]:
                row = index[shift_b(state, link, amount)]
                h[row, col] += -g2 / 2

        for vertex in range(4):
            for amount in [1, -1]:
                row = index[matter_shift_b(state, vertex, amount)]
                h[row, col] += -eta / 2

    hermiticity_error = float(np.linalg.norm(h - h.conj().T))
    assert hermiticity_error < 1e-10
    return h


def spectrum(beta, kappa=1.0, eta=0.0):
    eigvals = np.linalg.eigvalsh(hamiltonian(beta, kappa=kappa, eta=eta))
    return eigvals


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


def old_pauli_shadow_gap(beta):
    # Values from css_matter_coupled_gauss_proxy.py, rounded enough for comparison.
    table = {
        0.5: 0.520538,
        1.0: 1.25499,
    }
    return table.get(float(beta), None)


def scan_confining_regime():
    hd("B. Confining-Regime Open-Link Qutrit Scan")
    print(
        "  beta       eta   g2=1/beta   gap      degeneracy  "
        "old Pauli no-flip"
    )
    worst_by_eta = {}
    for eta in [0.0, 0.25, 0.5, 1.0]:
        gaps = []
        for beta in [0.25, 0.5, 0.75, 1.0]:
            eigvals = spectrum(beta=beta, kappa=1.0, eta=eta)
            levels, degeneracies = unique_levels(eigvals)
            gap = levels[1] - levels[0]
            gaps.append(gap)
            old_gap = old_pauli_shadow_gap(beta)
            old_text = "n/a" if old_gap is None else f"{old_gap:.6g}"
            print(
                f"  {beta:<10g} "
                f"{eta:<5g} "
                f"{1 / beta:<12.6g} "
                f"{gap:<8.6g} "
                f"{degeneracies[0]:<11d} "
                f"{old_text}"
            )
        worst_by_eta[eta] = min(gaps)
        print(f"  eta={eta:<4g} worst scanned physical gap = {worst_by_eta[eta]:.6g}\n")

    assert worst_by_eta[0.0] > 0.1
    assert worst_by_eta[0.5] > 0.1
    return worst_by_eta


def beta_grid_stability():
    hd("C. Strong-Coupling Grid Stability")
    print("  eta=0.5 with beta grid 0.2..1.0")
    gaps = []
    for beta in np.linspace(0.2, 1.0, 9):
        eigvals = spectrum(beta=float(beta), kappa=1.0, eta=0.5)
        levels, degeneracies = unique_levels(eigvals)
        gap = levels[1] - levels[0]
        gaps.append(gap)
        print(
            f"  beta={beta:.2f}: gap={gap:.6g}, "
            f"ground degeneracy={degeneracies[0]}, E0={levels[0]:.6g}"
        )
    print(f"  minimum grid gap = {min(gaps):.6g}")
    assert min(gaps) > 0.1


def main():
    hd("Scope")
    print(
        "This is an open-link Z3/qutrit-colour gauge-Higgs truncation of the "
        "requested SU(3) Wilson/TCH problem. It has explicit vertex Gauss "
        "projection and scans the confining beta<=1 regime, but it does not "
        "yet include non-centre SU(3) Peter-Weyl irreps."
    )
    orbit_audit()
    scan_confining_regime()
    beta_grid_stability()
    hd("Verdict")
    print(
        "The open-link qutrit-colour centre/TCH truncation passes the first "
        "strong-coupling test: after explicit vertex Gauss projection, the "
        "qutrit colour-centre physical gap stays finite for beta<=1 across the "
        "scanned eta values. This removes the old Pauli Z2 strong-coupling "
        "collapse in the sector that diagnostic actually tested: the "
        "colour-centre/qutrit sector."
    )
    print(
        "\nThis is still one tier below the full problem. The next escalation is a "
        "Peter-Weyl SU(3) open-link truncation, e.g. link irreps "
        "(0,0), (1,0), (0,1), with explicit Clebsch-Gordan/Gauss intertwiners, "
        "then the same beta<=1 matter-gap scan."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
