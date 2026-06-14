#!/usr/bin/env python3
"""
Background-twist gap diagnostic for the cell-level CSS SMG Hamiltonian.

Question:
  Does the [[8,0,4]] CSS symmetric gap survive simple background gauge/Wilson
  twists?

Model used here:
  This is a finite-cell "shadow" of a Wilson/holonomy test. For a diagonal gauge
  charge q(c) on codeword c and CSS generator g, replace the bare flip X_g by

      X_g(theta)|c> =
          exp(i theta [q(c xor g) - q(c)] / 2) |c xor g>.

  This is the usual covariant-shift phase rule for a background diagonal gauge
  connection. The operator is Hermitian/unitary because the phase reverses under
  c <-> c xor g. If all four X stabilizers use the same theta, the Hamiltonian
  is only a unitary conjugate of the zero-background CSS Hamiltonian, so the
  gap must remain exactly 2. Nonuniform theta_j are the finite-cell analogue of
  threading relative holonomy/flux through the stabilizer terms.

What this can and cannot establish:
  It can falsify naive robustness of the undynamical cell CSS gap under simple
  background twists. It cannot prove the full TCH dynamically gauged construction,
  because the actual gauge-dressed Wilson/plaquette degrees of freedom are not
  yet identified.

numpy; self-asserting.
"""

import itertools
import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1, -1]).astype(complex)

G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
ALL_WORDS = [tuple(bits) for bits in itertools.product([0, 1], repeat=8)]

GENERATORS = [
    [1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 1, 1, 1, 1],
    [0, 0, 1, 1, 0, 0, 1, 1],
    [0, 1, 0, 1, 0, 1, 0, 1],
]


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def idx(bits):
    return int("".join(map(str, bits)), 2)


def op(singles):
    out = np.array([[1]], complex)
    for bit in range(8):
        out = np.kron(out, singles.get(bit, I2))
    return out


def z_codeword(codeword):
    return op({bit: Z for bit, value in enumerate(codeword) if value})


def charge_zf(bits):
    return 1 if bits[I3] == 0 else -1


def sum_colour_projectors(bits):
    return -3 if (bits[C0], bits[C1]) == (0, 0) else -1


def electric_charge(bits):
    return 0.5 * charge_zf(bits) + (1 / 3) * sum_colour_projectors(bits) + 0.5


def weak_t3(bits):
    if bits[CHI] != 0:
        return 0.0
    return 0.5 if bits[I3] == 0 else -0.5


def hypercharge(bits):
    return electric_charge(bits) - weak_t3(bits)


def colour_cartan_shadow(bits):
    # Minimal diagonal colour shadow for the three nonzero colour vertices.
    # Leptons are colour-neutral. This is not a full SU(3) treatment.
    colour = (bits[C0], bits[C1])
    if colour == (0, 0):
        return 0.0
    if colour == (0, 1):
        return 1.0
    if colour == (1, 0):
        return -1.0
    return 0.0


def twisted_x(codeword, theta, charge_fn):
    matrix = np.zeros((256, 256), complex)
    row = tuple(codeword)
    for bits in ALL_WORDS:
        target = tuple(bit ^ flip for bit, flip in zip(bits, row))
        phase = np.exp(0.5j * theta * (charge_fn(target) - charge_fn(bits)))
        matrix[idx(target), idx(bits)] = phase
    return matrix


def css_hamiltonian(thetas, charge_fn):
    hamiltonian = -sum(z_codeword(row) for row in GENERATORS)
    for row, theta in zip(GENERATORS, thetas):
        hamiltonian -= twisted_x(row, theta, charge_fn)
    hermiticity_error = float(np.linalg.norm(hamiltonian - hamiltonian.conj().T))
    return hamiltonian, hermiticity_error


def gap(thetas, charge_fn):
    hamiltonian, hermiticity_error = css_hamiltonian(thetas, charge_fn)
    eigvals = np.linalg.eigvalsh(hamiltonian)
    return float(eigvals[1] - eigvals[0]), eigvals[:8], hermiticity_error


def scan_one_stabilizer(charge_fn, n_grid=65):
    best = (float("inf"), None, None)
    for theta in np.linspace(0, 2 * np.pi, n_grid):
        value, eigvals, hermiticity_error = gap([theta, 0, 0, 0], charge_fn)
        if value < best[0]:
            best = (value, theta, eigvals)
        assert hermiticity_error < 1e-10
    return best


def scan_uniform(charge_fn, n_grid=65):
    best = (float("inf"), None, None)
    for theta in np.linspace(0, 2 * np.pi, n_grid):
        value, eigvals, hermiticity_error = gap([theta, theta, theta, theta], charge_fn)
        if value < best[0]:
            best = (value, theta, eigvals)
        assert hermiticity_error < 1e-10
    return best


def scan_random(charge_fn, n_samples=500, seed=143):
    rng = np.random.default_rng(seed)
    best = (float("inf"), None, None)
    for _ in range(n_samples):
        thetas = rng.uniform(0, 2 * np.pi, 4)
        value, eigvals, hermiticity_error = gap(thetas, charge_fn)
        if value < best[0]:
            best = (value, thetas, eigvals)
        assert hermiticity_error < 1e-10
    return best


def report_charge_sector(name, charge_fn):
    hd(name)
    zero_gap, zero_eigvals, zero_herm = gap([0, 0, 0, 0], charge_fn)
    uniform_gap, uniform_theta, _ = scan_uniform(charge_fn)
    one_gap, one_theta, one_eigvals = scan_one_stabilizer(charge_fn)
    random_gap, random_thetas, random_eigvals = scan_random(charge_fn)

    print(f"  zero-background gap                 : {zero_gap:.12g}")
    print(f"  min uniform-twist gap               : {uniform_gap:.12g} at theta={uniform_theta:.6f}")
    print(f"  min one-stabilizer twist gap        : {one_gap:.12g} at theta={one_theta:.6f}")
    print(f"  one-stabilizer lowest eigvals       : {np.round(one_eigvals[:6], 8)}")
    print(f"  min random nonuniform gap           : {random_gap:.12g}")
    print(f"  random argmin thetas                : {np.round(random_thetas, 6)}")
    print(f"  random lowest eigvals               : {np.round(random_eigvals[:6], 8)}")

    assert zero_herm < 1e-10
    assert abs(zero_gap - 2.0) < 1e-10
    assert uniform_gap > 1.99
    return one_gap, random_gap


def main():
    hd("Scope")
    print(
        "This probes the finite-cell CSS gap under diagonal background holonomies. "
        "Uniform twists are pure conjugations and should preserve the gap. "
        "Nonuniform twists are the first discrete shadow of relative Wilson flux."
    )

    results = {
        "U(1)_Y hypercharge twist": report_charge_sector("U(1)_Y hypercharge twist", hypercharge),
        "SU(2)_L T3 Cartan twist": report_charge_sector("SU(2)_L T3 Cartan twist", weak_t3),
        "electric Q twist": report_charge_sector("electric Q twist", electric_charge),
        "SU(3)_c diagonal shadow twist": report_charge_sector(
            "SU(3)_c diagonal shadow twist", colour_cartan_shadow
        ),
    }

    hd("Verdict")
    for name, (one_gap, random_gap) in results.items():
        print(f"  {name:32s}: one-stab min={one_gap:.3e}, random min={random_gap:.3e}")

    min_one = min(value[0] for value in results.values())
    min_random = min(value[1] for value in results.values())
    assert min_one < 1e-8
    assert min_random < 0.1

    print(
        "\nThe undynamical cell CSS gap is not robust against simple nonuniform "
        "background twists: one-stabilizer or random relative holonomies can drive "
        "the first gap to zero or near-zero, while uniform pure-gauge twists leave "
        "the gap at 2. This is a negative/diagnostic result, not a no-go theorem."
    )
    print(
        "\nNext required step: repeat the same spectral scan after replacing the "
        "ad hoc diagonal twist by the actual TCH gauge-dressed Wilson/plaquette "
        "operators. A surviving finite gap there would be real evidence for "
        "locality-in-gauge-fields; this bare/twisted-cell result does not provide it."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
