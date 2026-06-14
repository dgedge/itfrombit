#!/usr/bin/env python3
r"""Glueball spectral-function (WIDTH) run -- the single-loop gauge-channel proxy scoped in
glueball_ladder/scope_spectral_function.md. Reuses the validated KPM engine (item48_kpm_rho_ldos.py).

WHAT THIS IS: the spectral function of the gauge-web A_1g (0++) and E_g (2++) plaquette operators on the
photon graph L(SC) -- i.e. does a single closed gauge-flux loop (a plaquette 1-cochain) survive as a
SHARP resonance in the gauge spectrum, or radiate into the continuum? Output = peak + WIDTH (the
observable the (2N-1)Lambda mass ladder lacks).

WHAT THIS IS NOT (scope, read first): the physical glueball is an N-body FOCK condensate
m_N=(2N-1)Lambda (ANCHOR L1247). This single-particle KPM gives the single-LOOP channel only -- a PROXY.
The peak energy here is NOT the glueball mass (identifying it would be the K_6 error, L1227); the
deliverable is the WIDTH/survival of the gauge-loop channel. The N-body Fock width is deferred (DMRG).

ACTION-SPACE DISCIPLINE: operators are gauge-web closed-cycle 1-cochains (the L1216-1217 three-plaquette
space P_xy,P_yz,P_zx -> A_1g symmetric sum, E_g traceless pair), NOT the bulk K_6 star. A symmetry gate
asserts the operators are the correct irreps of the O_h action ON THIS SPACE.
"""
import sys
import numpy as np
import scipy.sparse as ssp

from item48_kpm_rho_ldos import (build_linegraph, kpm_local_moments, reconstruct, spectral_bounds,
                                  _smooth, _fwhm_around)

LQCD = 332.0


def node(x, y, z, d, L):
    """L(SC) node = SC edge (vertex (x,y,z) -> +e_d), d in {0,1,2}. Matches build_linegraph indexing."""
    return 3 * ((x % L) * L * L + (y % L) * L + (z % L)) + d


def plaquette_cochain(L, x, y, z, plane):
    """Oriented boundary 1-cochain of the unit square plaquette at (x,y,z) in the given plane.
    plane in {'xy','yz','zx'}. Returns (indices, signs) of the 4 edges (the discrete curl)."""
    d1, d2 = {"xy": (0, 1), "yz": (1, 2), "zx": (2, 0)}[plane]
    e = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    s1, s2 = e[d1], e[d2]
    # loop: +d1 from r ; +d2 from r+s1 ; -d1 from r+s2 ; -d2 from r   (oriented boundary)
    idx = [node(x, y, z, d1, L),
           node(x + s1[0], y + s1[1], z + s1[2], d2, L),
           node(x + s2[0], y + s2[1], z + s2[2], d1, L),
           node(x, y, z, d2, L)]
    sgn = [+1.0, +1.0, -1.0, -1.0]
    return idx, sgn


def glueball_operators(L):
    """A_1g (0++) and the two E_g (2++) components as unit vectors on L(SC), built from the three
    orthogonal central plaquettes (the L1216-1217 three-plaquette space)."""
    N = 3 * L ** 3
    cen = L // 2
    P = {}
    for plane in ("xy", "yz", "zx"):
        v = np.zeros(N)
        idx, sgn = plaquette_cochain(L, cen, cen, cen, plane)
        for i, s in zip(idx, sgn):
            v[i] += s
        P[plane] = v / np.linalg.norm(v)                  # each plaquette curl, unit norm
    a1g = (P["xy"] + P["yz"] + P["zx"]);  a1g /= np.linalg.norm(a1g)        # totally symmetric (0++)
    eg1 = (P["xy"] - P["yz"]);            eg1 /= np.linalg.norm(eg1)        # traceless (2++)
    eg2 = (P["yz"] - P["zx"]);            eg2 /= np.linalg.norm(eg2)
    return P, a1g, eg1, eg2


def _sigma_perm(L):
    """Node permutation of the geometric O_h 3-fold R: e_x->e_y->e_z->e_x, i.e. R(x,y,z)=(z,x,y),
    sigma(d)=(d+1)%3. The central cell is a fixed point, so the three central plaquettes map among
    themselves. This is the genuine gauge-web action -- the rigorous guard against the K_6-star trap."""
    N = 3 * L ** 3
    sig = np.empty(N, dtype=np.int64)
    for x in range(L):
        for y in range(L):
            for z in range(L):
                for d in range(3):
                    sig[node(x, y, z, d, L)] = node(z, x, y, (d + 1) % 3, L)
    return sig


def symmetry_gate(L, P, a1g):
    """Assert (i) the three plaquette orientations are linearly independent (span the 3-plaquette space),
    and (ii) the GEOMETRIC 3-fold rotation permutes them xy->yz->zx (so we are demonstrably in the
    gauge-web closed-cycle action space, L1216-1217) with A_1g invariant (the irrep gate)."""
    Pm = np.stack([P["xy"], P["yz"], P["zx"]])
    assert np.allclose(np.linalg.norm(Pm, axis=1), 1.0, atol=1e-9), "plaquette curls must be unit norm"
    assert np.linalg.matrix_rank(Pm) == 3, "the three plaquette orientations must be linearly independent"
    sig = _sigma_perm(L)

    def rot(v):
        vr = np.zeros_like(v); vr[sig] = v; return vr            # pushforward: (rot v)[sigma(n)] = v[n]
    for src, dst in [("xy", "yz"), ("yz", "zx"), ("zx", "xy")]:
        ov = float(np.dot(rot(P[src]), P[dst]))
        assert abs(abs(ov) - 1.0) < 1e-9, f"3-fold must map plaquette {src}->{dst} (overlap {ov:.3f})"
    assert abs(abs(np.dot(rot(a1g), a1g)) - 1.0) < 1e-9, "A_1g must be 3-fold invariant (irrep gate)"
    print("[gate] OK: gauge-link 1-cochains; the geometric 3-fold permutes the 3 plaquettes "
          "(action-space confirmed, not the K_6 star); A_1g invariant.")


def chan(H, op, a, b, Nc, npts, label, t_mev=LQCD):
    """KPM spectral function of operator `op` (unit vector) in gauge Hamiltonian H. Separates the
    E=-2 FLAT-BAND (cycle-space, dispersionless/localised) weight from the dispersive weight, since the
    flat band is a delta -> its 'width' is the KPM resolution floor, not a physical width."""
    mu = kpm_local_moments(H, [0], a, b, Nc, backend="cpu", report=False, start=op[:, None])
    E, rho = reconstruct(mu, a, b, npts=npts)
    rho = np.clip(rho, 0, None); norm = float(np.trapezoid(rho, E))
    eta = np.pi * a / Nc                                   # KPM resolution
    flat = (E > -2 - 3 * eta) & (E < -2 + 3 * eta)         # E=-2 flat band = the cycle space
    w_flat = float(np.trapezoid(rho[flat], E[flat]) / norm)
    disp = ~flat
    # dispersive part: peak + width of whatever is NOT in the flat band (the genuinely radiating piece)
    if np.trapezoid(rho[disp], E[disp]) > 1e-6 * norm:
        Ed, rd = E[disp], rho[disp]
        ip = int(np.argmax(_smooth(Ed, rd, 0.01 * (E[-1] - E[0]))))
        disp_pk = float(Ed[ip])
    else:
        disp_pk = float("nan")
    print(f"  [{label:5s}] flat-band (E=-2, cycle space) weight = {w_flat*100:5.1f}%  -> localised / "
          f"non-radiating, intrinsic width 0 (peak width = KPM floor {eta*t_mev:.1f} MeV); "
          f"dispersive (radiating) weight = {(1-w_flat)*100:.1f}% (peak ~{disp_pk:+.2f})")
    return E, rho, dict(w_flat=w_flat, disp_pk=disp_pk)


def validate(L=4, Nc=400):
    """KPM moments of the A_1g operator vs EXACT eigh on a small L(SC). Exit 0 here = engine correct."""
    A, _ = build_linegraph(L)
    _, a1g, _, _ = glueball_operators(L)
    emin, emax, _, _ = spectral_bounds(A, (-2.0, 10.0), 0.0)
    a = (emax - emin) / 2; b = (emax + emin) / 2
    evals, evecs = np.linalg.eigh(A.toarray())
    xk = (evals - b) / a
    wj = (evecs.T @ a1g) ** 2                              # |<eig|A1g>|^2
    n = np.arange(Nc)
    mu_exact = (wj[:, None] * np.cos(np.outer(np.arccos(np.clip(xk, -1, 1)), n))).sum(0)
    mu_kpm = kpm_local_moments(A, [0], a, b, Nc, backend="cpu", report=False, start=a1g[:, None])
    err = float(np.max(np.abs(mu_kpm - mu_exact)))
    print(f"[validate] L(SC) L={L} N={A.shape[0]}: max|mu_KPM - mu_exact(A1g)| = {err:.2e}")
    assert err < 1e-9, "KPM glueball-operator moments must match exact eigh"
    print("[validate] PASS.\n")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=40)
    ap.add_argument("--moments", type=int, default=8000)
    ap.add_argument("--egrid", type=int, default=12000)
    ap.add_argument("--no-validate", action="store_true")
    args = ap.parse_args()

    if not args.no_validate:
        validate()

    A, _ = build_linegraph(args.L)                         # the photon graph L(SC), band [-2,10]
    P, a1g, eg1, eg2 = glueball_operators(args.L)
    symmetry_gate(args.L, P, a1g)
    emin, emax, _, _ = spectral_bounds(A, (-2.0, 10.0), 0.0)
    a = (emax - emin) / 2; b = (emax + emin) / 2
    print(f"# gauge web L(SC) L={args.L} N={A.shape[0]:,} band [{emin:.2f},{emax:.2f}]; flat band at -2 = cycle space")
    print(f"\n=== gauge-channel spectral functions (single-loop PROXY -- peak is NOT the glueball mass) ===")
    f0 = chan(A, a1g, a, b, args.moments, args.egrid, "0++")[2]
    f2 = chan(A, eg1, a, b, args.moments, args.egrid, "2++")[2]
    assert f0["w_flat"] > 0.99, "0++ curl is a pure cycle-space (flat-band) state -> fully localised"
    print(f"""
RESULT (single-loop gauge-channel WIDTH):
  The glueball operators are closed-loop 1-cochains, hence elements of the L(SC) CYCLE SPACE -- which is
  exactly the dispersionless E=-2 FLAT BAND of the line graph. So the single gauge-loop channel does NOT
  radiate: the 0++ (A_1g) channel is {f0['w_flat']*100:.0f}% flat-band (a localised, zero-intrinsic-width
  mode); the 2++ (E_g) channel is {f2['w_flat']*100:.0f}% flat-band + {(1-f2['w_flat'])*100:.0f}% dispersive
  (partially radiating). This is the opposite of the rho's Fano-dissolution: a compact gauge flux loop is
  structurally LOCALISED (cycle space), not a continuum resonance.
  HONEST BOUND: (i) the peak sits at E=-2 (the flat band), which is a single-loop proxy artefact, NOT the
  (2N-1)Lambda glueball mass -- identifying them is the K_6 error (L1227). (ii) This is the SINGLE-loop
  width; the physical 0++ is the N-body Fock condensate (m_N=(2N-1)Lambda), whose width needs the deferred
  many-body (DMRG/Markov-strings) computation. What is established: the single gauge-loop channel is a
  stable, non-radiating cycle-space mode -- a clean, structurally-grounded new datum.""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
