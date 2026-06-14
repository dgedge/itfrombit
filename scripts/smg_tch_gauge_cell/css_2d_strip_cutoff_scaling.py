#!/usr/bin/env python3
"""
2D strip cutoff scaling after the hopping-truncation audit.

The n/q=3 two-plaquette strip showed a visible t-response, but the 1D audit
shows that this cutoff is still too small to establish gap stability. This
script pushes the same 2D charge-block/Krylov representation one cutoff step
harder and makes the computational wall explicit.

What is computed:

  * n/q=3 reference using the checked charge-block strip audit.
  * n/q=4 streaming Lanczos estimates without storing all Krylov vectors.

What is not honestly computable in this representation:

  * n/q=6 and n/q=8 full-vector Krylov scans. Their dimensions are about
    1.02e8 and 5.73e8 respectively; the matvec cost, not the algebra, becomes
    the blocker. Those need a more compressed/tensor-network matvec.

Acceptance logic:

  "stays gapped" is not accepted from one low cutoff. The meaningful diagnostic
  here is the cutoff trend of the t=4 gap and the t-response.
"""

import importlib.util
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
AUDIT_PATH = HERE / "css_2d_strip_chargeblock_krylov_audit.py"


def load_audit():
    spec = importlib.util.spec_from_file_location("strip_audit", AUDIT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


audit = load_audit()


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def estimate_dimension(states_per_charge):
    return 2187 * states_per_charge**6


def memory_gib(elements, bytes_per_element):
    return elements * bytes_per_element / 1024**3


def streaming_lanczos_evals(model, n_iter=32, seed=20260602):
    rng = np.random.default_rng(seed)
    q = rng.normal(size=model.dim) + 1j * rng.normal(size=model.dim)
    q = q / np.linalg.norm(q)
    q_prev = np.zeros_like(q)
    beta_prev = 0.0
    alphas = []
    betas = []

    for _ in range(min(n_iter, model.dim)):
        v = model.matvec(q)
        alpha = float(np.real(np.vdot(q, v)))
        v = v - alpha * q - beta_prev * q_prev
        beta = float(np.linalg.norm(v))
        alphas.append(alpha)
        if beta < 1e-12:
            break
        betas.append(beta)
        q_prev = q
        q = v / beta
        beta_prev = beta

    tridiag = np.diag(alphas)
    for index in range(len(alphas) - 1):
        tridiag[index, index + 1] = betas[index]
        tridiag[index + 1, index] = betas[index]
    evals, evecs = np.linalg.eigh(tridiag)
    residuals = np.zeros_like(evals)
    if betas:
        last_beta = betas[-1]
        residuals = np.abs(last_beta * evecs[-1, :])
    return evals, residuals


def gap_estimate(states_per_charge, beta, hopping, n_iter):
    model = audit.ChargeBlockStripHamiltonian(
        states_per_charge=states_per_charge,
        beta=beta,
        hopping=hopping,
    )
    herm = audit.hermiticity_probe(model)
    assert herm < 1e-10
    evals, residuals = streaming_lanczos_evals(model, n_iter=n_iter)
    return {
        "dim": model.dim,
        "gap": float(evals[1] - evals[0]),
        "e0_resid": float(residuals[0]),
        "e1_resid": float(residuals[1]),
        "herm": herm,
    }


def resource_table():
    hd("A. Full-Vector Resource Scale")
    print("  n/q   dim          one complex vector   4 work vectors   diagonal pair")
    for states_per_charge in [3, 4, 6, 8]:
        dim = estimate_dimension(states_per_charge)
        one_vector = memory_gib(dim, 16)
        work = 4 * one_vector
        diagonals = 2 * memory_gib(dim, 8)
        print(
            f"  {states_per_charge:<5d} "
            f"{dim:<12d} "
            f"{one_vector:<20.3f} "
            f"{work:<16.3f} "
            f"{diagonals:.3f} GiB"
        )


def cutoff_scan():
    hd("B. 2D Strip Cutoff Scan")
    print("  beta=0.5, two plaquettes, streaming Lanczos")
    print("  n/q   t      iter   dim       full gap   e0 resid   e1 resid")
    rows = []
    cases = [
        (3, 24, [0.2, 4.0]),
        (4, 32, [0.2, 4.0]),
    ]
    for states_per_charge, n_iter, hoppings in cases:
        for hopping in hoppings:
            result = gap_estimate(states_per_charge, beta=0.5, hopping=hopping, n_iter=n_iter)
            rows.append((states_per_charge, hopping, result))
            print(
                f"  {states_per_charge:<5d} "
                f"{hopping:<6g} "
                f"{n_iter:<6d} "
                f"{result['dim']:<9d} "
                f"{result['gap']:<10.6g} "
                f"{result['e0_resid']:<10.3g} "
                f"{result['e1_resid']:.3g}"
            )

    by_cutoff = {}
    for states_per_charge, hopping, result in rows:
        by_cutoff.setdefault(states_per_charge, {})[hopping] = result["gap"]
    for states_per_charge, gaps in by_cutoff.items():
        response = abs(gaps[0.2] - gaps[4.0])
        print(
            f"  response n/q={states_per_charge}: "
            f"|gap(0.2)-gap(4)|={response:.6g}, "
            f"relative={response / max(abs(gaps[0.2]), 1e-12):.3f}"
        )
    assert by_cutoff[4][4.0] < by_cutoff[3][4.0]
    return rows


def main():
    hd("Scope")
    print(
        "This is a cutoff-scaling audit, not a stability claim. The question is "
        "whether the 2D t=4 gap keeps falling as the local mirror-Fock cutoff is "
        "increased beyond n/q=3."
    )
    resource_table()
    rows = cutoff_scan()
    hd("Verdict")
    print(
        "The requested n/q=3->4 direction is now directly tested. n/q=6 and "
        "n/q=8 are not practical with this full-vector Krylov representation; "
        "they require a compressed/tensor-network matvec or another reduction."
    )
    print(f"\nminimum scanned full gap: {min(row[2]['gap'] for row in rows):.6g}")
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
