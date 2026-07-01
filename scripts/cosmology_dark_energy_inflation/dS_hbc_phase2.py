#!/usr/bin/env python
r"""dS_hbc_phase2.py -- HBC dynamics spec, Phase 2: de Sitter thermality (Gibbons-Hawking T ~ 1/R_H).

Phase 1 established the horizon AREA law and that adiabatic printing sustains it. Phase 2 tests the
THERMAL content: a de Sitter horizon radiates at T_GH = H/2pi = 1/(2 pi R_H), i.e. the horizon
temperature scales as 1/R_H. We test this via the ENTANGLEMENT FIRST LAW: for a causal-patch ball
of radius R, a small uniform excitation changes the entanglement by dS and the ball energy by dE,
and the entanglement temperature T_ent = dE/dS. For a CFT, T_ent ~ 1/R would be the Gibbons-Hawking
scaling. RESULT (honest negative): at lattice radii the naive first-law T_ent is UV-area-law-
dominated and RISES as ~R, NOT 1/R -- the Gibbons-Hawking temperature is not cleanly extractable
this way; it needs UV-subtraction (quartic fit) or the Bisognano-Wichmann modular-weight profile.
The horizon THERMALITY itself is already established in Phase 1 (dS_horizon_modular: modular gap
closing); this script documents that the clean TEMPERATURE SCALING needs the proper UV-safe method.

Perturbation = a small global thermal population (only the IR modes excited); dX,dP carry the
Bose factor n_k. Free Gaussian; bosonic region entropy from symplectic eigenvalues.
"""
import numpy as np

def correlators(L, m2, tau):
    k = 2*np.pi*np.arange(L)/L
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    Kk = m2 + 6 - 2*(np.cos(KX)+np.cos(KY)+np.cos(KZ))
    om = np.sqrt(Kk)
    if tau is None:                       # vacuum
        coth = np.ones_like(om)
    else:
        coth = 1.0/np.tanh(om/(2*tau))
    Xk = 0.5*coth/om
    Pk = 0.5*coth*om
    Xr = np.fft.ifftn(Xk).real
    Pr = np.fft.ifftn(Pk).real
    return Xr, Pr, Kk

def Kreal(L, m2):
    k = 2*np.pi*np.arange(L)/L
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    Kk = m2 + 6 - 2*(np.cos(KX)+np.cos(KY)+np.cos(KZ))
    return np.fft.ifftn(Kk).real          # real-space kernel row K(r)

def sub(M, sites, L):
    n = len(sites); A = np.empty((n, n))
    for a, sa in enumerate(sites):
        for b, sb in enumerate(sites):
            A[a, b] = M[(sa[0]-sb[0]) % L, (sa[1]-sb[1]) % L, (sa[2]-sb[2]) % L]
    return A

def entropy(XA, PA):
    nu = np.sqrt(np.clip(np.linalg.eigvals(XA@PA).real, 0, None))
    nu = np.clip(nu, 0.5+1e-12, None)
    return float(np.sum((nu+0.5)*np.log(nu+0.5) - (nu-0.5)*np.log(nu-0.5)))

def ball(R, L):
    c = L//2
    return [(c+x, c+y, c+z) for x in range(-R, R+1) for y in range(-R, R+1) for z in range(-R, R+1)
            if x*x+y*y+z*z <= R*R]

def ball_energy(Xr, Pr, Kr, sites, L):
    # E = sum_{i in B} [ 1/2 P_ii + 1/2 (K X)_ii ]; local since translation-invariant
    P0 = Pr[0, 0, 0]
    # (K X)_ii = sum_r K(r) X(-r) = sum_r K(r) X(r) (both even) -> a scalar per site
    KX_ii = float(np.sum(Kr * Xr))        # = (K X)_00, same for every site by translation inv.
    return len(sites) * 0.5 * (P0 + KX_ii)

def linfit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None); yh = A@c
    ss, st = np.sum((y-yh)**2), np.sum((y-y.mean())**2)
    return float(c[0]), float(1-ss/st if st > 0 else 0)

def main():
    L, m2, tau = 24, 2e-3, 0.15
    Xv, Pv, _ = correlators(L, m2, None)      # vacuum
    Xt, Pt, _ = correlators(L, m2, tau)       # small thermal excitation
    Kr = Kreal(L, m2)
    Rs = [3, 4, 5, 6, 7, 8]
    print(f"=== Phase 2: de Sitter thermality, entanglement temperature T_ent=dE/dS (L={L}, tau={tau}) ===")
    rows = []
    for R in Rs:
        b = ball(R, L)
        Sv = entropy(sub(Xv, b, L), sub(Pv, b, L))
        St = entropy(sub(Xt, b, L), sub(Pt, b, L))
        dS = St - Sv
        dE = ball_energy(Xt, Pt, Kr, b, L) - ball_energy(Xv, Pv, Kr, b, L)
        Tent = dE/dS if dS > 1e-12 else float("nan")
        rows.append((R, dS, dE, Tent))
        print(f"  R={R}: dS={dS:.4f}  dE={dE:.4f}  T_ent=dE/dS={Tent:.4f}   (1/R={1/R:.4f})")
    Rv = np.array([r for r, *_ in rows], float)
    Tent = np.array([t for *_, t in rows])
    slope_ll, r2_ll = linfit(np.log(Rv), np.log(Tent))
    dEv = np.array([e for _, _, e, _ in rows]); dSv = np.array([s for _, s, _, _ in rows])
    slope_E, _ = linfit(np.log(Rv), np.log(dEv))
    slope_S, _ = linfit(np.log(Rv), np.log(dSv))
    print(f"\n  scaling: dS ~ R^{slope_S:.2f},  dE ~ R^{slope_E:.2f}  ->  T_ent = dE/dS ~ R^{slope_ll:+.2f}  (R2={r2_ll:.3f})")
    print("\n[verdict] HONEST NEGATIVE: the naive first-law route is UV-DOMINATED, NOT yet Gibbons-Hawking.")
    print("  dE ~ R^3 (volume thermal energy) but dS ~ R^2 (the lattice UV AREA-law term, not the CHM")
    print("  first-law R^4 piece) -> T_ent RISES as ~R: the thermodynamic E/S of the populated subsystem,")
    print("  not the entanglement temperature. At these radii the UV area-law swamps the CHM 1/R signal, so")
    print("  T_GH ~ 1/R_H is NOT cleanly extractable this way. (The horizon THERMALITY itself is not in")
    print("  doubt -- Phase 1 / dS_horizon_modular shows the modular gap closing.) Proper UV-safe routes:")
    print("  (i) UV-subtracted fit dS = a R^2 + b R^4, T_ent from the R^4 piece (fragile at R<=8, needs")
    print("  larger lattices); (ii) the Bisognano-Wichmann modular-weight profile beta(r) ~ (R^2-r^2)/R")
    print("  -- the rigorous UV-safe route. This slice needs the proper method or larger L.")
    assert np.all(dSv > 0) and np.all(Tent > 0), "first law sane (dS, T_ent > 0)"
    assert slope_E > 2.5, "dE ~ volume (R^3), as expected for a thermal energy"
    assert slope_ll > 0.5, "documents the UV-domination: T_ent RISES with R (the naive route does NOT give 1/R)"
    print("\nASSERTIONS PASSED -- documents the UV-dominated NEGATIVE; Gibbons-Hawking 1/R needs the proper route.")
    print("exit 0")

if __name__ == "__main__":
    main()
