#!/usr/bin/env python3
"""
HPC runner for the two-plaquette charge-block mirror-Fock/SMG strip.

This is the "bring bigger hardware" extension of css_2d_strip_cutoff_scaling.py.
It targets the cutoff ladder that was too expensive for the current session:

    n/q = 6, 8

The calculation is still finite-volume and still the Z3 colour-orbit proxy. It
does not prove the constructive chiral-gauge theorem. It tests the concrete
numerical question that survived the audits:

    Does the two-plaquette t=4 gap keep falling as the local mirror-Fock cutoff
    is increased, or does it saturate?

Recommended first runs:

    python css_2d_strip_hpc_lanczos.py \
      --states-per-charge 6 \
      --hoppings 0.2 1.0 2.0 4.0 \
      --iterations 96 \
      --checkpoint-dir checkpoints_n6 \
      --out strip_n6.csv

    python css_2d_strip_hpc_lanczos.py \
      --states-per-charge 8 \
      --hoppings 0.2 4.0 \
      --iterations 128 \
      --checkpoint-dir checkpoints_n8 \
      --out strip_n8.csv

If SciPy is available and the machine has enough RAM, ARPACK can give a more
orthogonal low-spectrum check:

    python css_2d_strip_hpc_lanczos.py \
      --solver eigsh --states-per-charge 6 --hoppings 0.2 4.0 \
      --k 4 --ncv 24 --tol 1e-8 --out strip_n6_eigsh.csv

Memory warning:
  Full-vector Krylov scales as dim = 2187*(n/q)^6. Streaming Lanczos stores only
  a few full complex vectors plus model diagonals. eigsh stores many more basis
  vectors internally.
"""

import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import time

import numpy as np


HERE = Path(__file__).resolve().parent
AUDIT_PATH = HERE / "css_2d_strip_chargeblock_krylov_audit.py"


def load_audit():
    spec = importlib.util.spec_from_file_location("strip_audit", AUDIT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


audit = load_audit()


def estimate_dimension(states_per_charge):
    return 2187 * states_per_charge**6


def gib(bytes_count):
    return bytes_count / 1024**3


def memory_estimate(states_per_charge, solver, ncv):
    dim = estimate_dimension(states_per_charge)
    complex_vector = dim * np.dtype(np.complex128).itemsize
    real_vector = dim * np.dtype(np.float64).itemsize
    if solver == "streaming":
        # q, q_prev, v/out plus diagonal. Temporaries can double this briefly.
        steady = 3 * complex_vector + real_vector
        conservative = 5 * complex_vector + 2 * real_vector
    else:
        # ARPACK/Scipy stores roughly ncv basis vectors, plus work vectors/model.
        steady = (ncv + 3) * complex_vector + real_vector
        conservative = (ncv + 6) * complex_vector + 2 * real_vector
    return {
        "dim": dim,
        "one_complex_vector_gib": gib(complex_vector),
        "steady_gib": gib(steady),
        "conservative_gib": gib(conservative),
    }


def checkpoint_path(checkpoint_dir, states_per_charge, beta, hopping, seed):
    if checkpoint_dir is None:
        return None
    safe_beta = str(beta).replace(".", "p")
    safe_t = str(hopping).replace(".", "p")
    name = f"strip_spc{states_per_charge}_beta{safe_beta}_t{safe_t}_seed{seed}.npz"
    return Path(checkpoint_dir) / name


def save_checkpoint(path, iteration, q, q_prev, beta_prev, alphas, betas):
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as handle:
        np.savez(
            handle,
            iteration=np.array([iteration], dtype=np.int64),
            beta_prev=np.array([beta_prev], dtype=np.float64),
            q=q,
            q_prev=q_prev,
            alphas=np.array(alphas, dtype=np.float64),
            betas=np.array(betas, dtype=np.float64),
        )
    os.replace(tmp, path)


def load_checkpoint(path):
    if path is None or not path.exists():
        return None
    data = np.load(path)
    return {
        "iteration": int(data["iteration"][0]),
        "beta_prev": float(data["beta_prev"][0]),
        "q": data["q"],
        "q_prev": data["q_prev"],
        "alphas": list(map(float, data["alphas"])),
        "betas": list(map(float, data["betas"])),
    }


def tridiagonal_result(alphas, betas):
    tridiag = np.diag(alphas)
    for index in range(len(alphas) - 1):
        tridiag[index, index + 1] = betas[index]
        tridiag[index + 1, index] = betas[index]
    evals, evecs = np.linalg.eigh(tridiag)
    residuals = np.zeros_like(evals)
    if betas:
        residuals = np.abs(betas[-1] * evecs[-1, :])
    return evals, residuals


def streaming_lanczos(model, iterations, seed, checkpoint=None, checkpoint_every=0, resume=False):
    state = load_checkpoint(checkpoint) if resume else None
    if state is None:
        rng = np.random.default_rng(seed)
        q = rng.normal(size=model.dim) + 1j * rng.normal(size=model.dim)
        q = q / np.linalg.norm(q)
        q_prev = np.zeros_like(q)
        beta_prev = 0.0
        alphas = []
        betas = []
        start = 0
    else:
        q = state["q"]
        q_prev = state["q_prev"]
        beta_prev = state["beta_prev"]
        alphas = state["alphas"]
        betas = state["betas"]
        start = state["iteration"]

    last_save = time.time()
    for iteration in range(start, iterations):
        step_start = time.time()
        v = model.matvec(q)
        alpha = float(np.real(np.vdot(q, v)))
        v = v - alpha * q - beta_prev * q_prev
        beta = float(np.linalg.norm(v))
        alphas.append(alpha)
        if beta < 1e-12:
            save_checkpoint(checkpoint, iteration + 1, q, q_prev, beta_prev, alphas, betas)
            break
        betas.append(beta)
        q_prev = q
        q = v / beta
        beta_prev = beta

        if checkpoint_every and (iteration + 1) % checkpoint_every == 0:
            save_checkpoint(checkpoint, iteration + 1, q, q_prev, beta_prev, alphas, betas)
            last_save = time.time()

        evals, residuals = tridiagonal_result(alphas, betas)
        gap = float(evals[1] - evals[0]) if len(evals) > 1 else float("nan")
        print(
            json.dumps(
                {
                    "event": "iteration",
                    "iteration": iteration + 1,
                    "seconds": round(time.time() - step_start, 3),
                    "seconds_since_save": round(time.time() - last_save, 3),
                    "ritz_e0": float(evals[0]),
                    "ritz_gap": gap,
                    "resid_e0": float(residuals[0]),
                    "resid_e1": float(residuals[1]) if len(residuals) > 1 else None,
                }
            ),
            flush=True,
        )

    evals, residuals = tridiagonal_result(alphas, betas)
    return evals, residuals, len(alphas)


def eigsh_solver(model, k, ncv, tol, maxiter, return_vectors):
    try:
        from scipy.sparse.linalg import LinearOperator, eigsh
    except Exception as exc:
        raise RuntimeError("SciPy is required for --solver eigsh") from exc

    operator = LinearOperator(
        shape=(model.dim, model.dim),
        dtype=np.complex128,
        matvec=model.matvec,
    )
    evals, evecs = eigsh(
        operator,
        k=k,
        which="SA",
        ncv=ncv,
        tol=tol,
        maxiter=maxiter,
        return_eigenvectors=return_vectors,
    )
    if return_vectors:
        order = np.argsort(evals)
        return evals[order], evecs[:, order]
    return np.sort(evals), None


def append_csv(path, row):
    if path is None:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    with open(path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def run_case(args, states_per_charge, hopping):
    mem = memory_estimate(states_per_charge, args.solver, args.ncv)
    print(json.dumps({"event": "case_start", "spc": states_per_charge, "t": hopping, **mem}), flush=True)
    start = time.time()
    model = audit.ChargeBlockStripHamiltonian(
        states_per_charge=states_per_charge,
        beta=args.beta,
        hopping=hopping,
        build_matter=args.compute_matter,
    )
    build_seconds = time.time() - start

    herm = None
    if args.hermiticity_probe:
        herm = audit.hermiticity_probe(model)
        if herm > args.hermiticity_tol:
            raise RuntimeError(f"Hermiticity probe failed: {herm}")

    solve_start = time.time()
    if args.solver == "streaming":
        ckpt = checkpoint_path(args.checkpoint_dir, states_per_charge, args.beta, hopping, args.seed)
        evals, residuals, used_iterations = streaming_lanczos(
            model,
            iterations=args.iterations,
            seed=args.seed,
            checkpoint=ckpt,
            checkpoint_every=args.checkpoint_every,
            resume=args.resume,
        )
        vectors = None
    else:
        evals, vectors = eigsh_solver(
            model,
            k=args.k,
            ncv=args.ncv,
            tol=args.tol,
            maxiter=args.maxiter,
            return_vectors=args.compute_matter,
        )
        residuals = np.full_like(evals, np.nan, dtype=float)
        used_iterations = None

    solve_seconds = time.time() - solve_start
    gap = float(evals[1] - evals[0])
    matter_delta_1 = None
    if args.compute_matter and vectors is not None:
        matter = [model.matter_expectation(vectors[:, index]) for index in range(vectors.shape[1])]
        matter_delta_1 = float(matter[1] - matter[0])

    row = {
        "solver": args.solver,
        "states_per_charge": states_per_charge,
        "beta": args.beta,
        "hopping": hopping,
        "dim": model.dim,
        "gap": gap,
        "e0": float(evals[0]),
        "e1": float(evals[1]),
        "resid_e0": float(residuals[0]) if len(residuals) else None,
        "resid_e1": float(residuals[1]) if len(residuals) > 1 else None,
        "iterations": used_iterations,
        "build_seconds": round(build_seconds, 3),
        "solve_seconds": round(solve_seconds, 3),
        "hermiticity": herm,
        "matter_delta_1": matter_delta_1,
    }
    print(json.dumps({"event": "case_done", **row}), flush=True)
    append_csv(args.out, row)
    return row


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--states-per-charge", type=int, nargs="+", default=[6])
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--hoppings", type=float, nargs="+", default=[0.2, 1.0, 2.0, 4.0])
    parser.add_argument("--solver", choices=["streaming", "eigsh"], default="streaming")
    parser.add_argument("--iterations", type=int, default=96)
    parser.add_argument("--seed", type=int, default=20260602)
    parser.add_argument("--checkpoint-dir", type=str, default=None)
    parser.add_argument("--checkpoint-every", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--out", type=str, default="strip_hpc_results.csv")
    parser.add_argument("--hermiticity-probe", action="store_true")
    parser.add_argument("--hermiticity-tol", type=float, default=1e-10)
    parser.add_argument("--compute-matter", action="store_true")
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--ncv", type=int, default=24)
    parser.add_argument("--tol", type=float, default=1e-8)
    parser.add_argument("--maxiter", type=int, default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    all_rows = []
    for states_per_charge in args.states_per_charge:
        for hopping in args.hoppings:
            all_rows.append(run_case(args, states_per_charge, hopping))

    by_cutoff = {}
    for row in all_rows:
        by_cutoff.setdefault(row["states_per_charge"], {})[row["hopping"]] = row["gap"]
    for states_per_charge, gaps in by_cutoff.items():
        if 0.2 in gaps and 4.0 in gaps:
            response = abs(gaps[0.2] - gaps[4.0])
            print(
                json.dumps(
                    {
                        "event": "response",
                        "states_per_charge": states_per_charge,
                        "gap_t0p2": gaps[0.2],
                        "gap_t4": gaps[4.0],
                        "absolute_response": response,
                        "relative_response": response / max(abs(gaps[0.2]), 1e-12),
                    }
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
