#!/usr/bin/env python3
r"""foundations_frame_coarse_graining_pilot.py

PHASE 1 of FRAME_COARSE_GRAINING_PLAN.md (the micro->macro frame coarse-graining, dynamical-frame residual (i)).
Establishes that the MACRO strain field exists and is recovered by coarse-graining the GAUGE-INVARIANT metric
g_{munu}=e^a_mu e^b_nu eta_ab -- NOT the frame e itself (which is corrupted by the local-Lorentz gauge).

Setup: a substrate patch (L^3 cells). Each cell carries a 4D vierbein e = Lambda_i (I + S(x_i) + n_i):
  - S(x) = a smooth macro conformal strain (the metric perturbation we want to recover);
  - n_i  = random anisotropic micro fluctuation (scale sigma);
  - Lambda_i = a random LOCAL-LORENTZ gauge rotation (each cell's frame independently oriented).

Results (the Phase 1 gate):
  (1) GAUGE INVARIANCE: g = e^T eta e is identical whether computed from the gauged or ungauged frame
      (Lambda^T eta Lambda = eta) -- the local-Lorentz noise DROPS. So g, not e, is the coarse-graining object.
  (2) MACRO FIELD EMERGES: block-averaging g recovers the smooth macro strain -- the residual (coarse g minus
      the noise-free macro) falls as 1/sqrt(N_block), and the correlation of the coarse conformal factor with
      the macro -> 1. The random micro fluctuations average out; the macro field survives.
  (3) RANDOM ANISOTROPY WASHES OUT: the anisotropic (spatial-traceless) part of the coarse g falls as
      1/sqrt(N_block) -> emergent isotropy from coarse-graining (for the RANDOM part).
  (4) CONTROL: averaging the FRAME e and then forming g=<e>^T eta <e> is corrupted by the gauge noise (error
      stays O(1), ratio to the g-route grows with block size) -- confirming you MUST coarse-grain g.

HONEST SCOPE: this is the KINEMATIC coarse-graining -- it establishes the macro strain field is well-defined
and gauge-invariant, and that RANDOM micro anisotropy + gauge noise wash out as 1/sqrt(N). It does NOT
establish emergent LI for SYSTEMATIC lattice anisotropy (that is the separate Route C dispersion result,
anisotropy ~ (ka)^2; block-averaging alone does not kill a systematic anisotropy), NOR the DYNAMICS (Phase 2:
is the frame kinetic term induced-dominated? -- the crux, the gravitational sister of the failed T-R2 photon).
Phase 1 supplies the macro field that Phase 2 then tests.

Self-asserting; exit 0.
"""
import numpy as np
from scipy.linalg import expm


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    rng = np.random.default_rng(0); eta = np.diag([1, -1, -1, -1.])
    gens = []
    for i, j in [(1, 2), (1, 3), (2, 3)]:
        M = np.zeros((4, 4)); M[i, j] = -1; M[j, i] = 1; gens.append(M)
    for i in (1, 2, 3):
        M = np.zeros((4, 4)); M[0, i] = 1; M[i, 0] = 1; gens.append(M)
    ok(all(np.allclose(M.T @ eta + eta @ M, 0) for M in gens), "Lorentz generators satisfy M^T eta + eta M = 0")
    pool = np.array([expm(sum(t * M for t, M in zip(rng.normal(size=6) * 0.5, gens))) for _ in range(2000)])
    ok(max(np.linalg.norm(Lam.T @ eta @ Lam - eta) for Lam in pool[:50]) < 1e-10, "gauge pool elements are Lorentz (Lambda^T eta Lambda = eta)")

    L = 24
    phi = lambda x, y, z: 0.04 * np.exp(-(((x - 12) ** 2 + (y - 12) ** 2 + (z - 12) ** 2) / (2 * 6.0 ** 2)))
    sig = 0.05
    e = np.zeros((L, L, L, 4, 4)); g = np.zeros((L, L, L, 4, 4)); g0 = np.zeros((L, L, L, 4, 4)); gM = np.zeros((L, L, L, 4, 4))
    for x in range(L):
        for y in range(L):
            for z in range(L):
                p = phi(x, y, z); S = np.diag([0, p, p, p])
                n = rng.normal(size=(4, 4)) * sig; n = (n + n.T) / 2
                e0 = np.eye(4) + S + n; Lam = pool[rng.integers(2000)]
                e[x, y, z] = Lam @ e0
                g[x, y, z] = e[x, y, z].T @ eta @ e[x, y, z]
                g0[x, y, z] = e0.T @ eta @ e0
                gM[x, y, z] = (np.eye(4) + S).T @ eta @ (np.eye(4) + S)

    ok(np.max(np.abs(g - g0)) < 1e-10, f"[1] g gauge-invariant: g(gauged frame) == g(ungauged), max diff {np.max(np.abs(g-g0)):.1e} -> local-Lorentz noise drops")

    conf = lambda M: -(M[..., 1, 1] + M[..., 2, 2] + M[..., 3, 3]) / 3

    def aniso(M):
        sp = M[..., 1:, 1:]; tr = (sp[..., 0, 0] + sp[..., 1, 1] + sp[..., 2, 2]) / 3
        a = sp.copy()
        for k in range(3):
            a[..., k, k] -= tr
        return np.sqrt((a ** 2).sum(axis=(-1, -2)))

    def blocks(field, b):
        nb = L // b; out = np.zeros((nb, nb, nb) + field.shape[3:])
        for i in range(nb):
            for j in range(nb):
                for k in range(nb):
                    out[i, j, k] = field[i * b:(i + 1) * b, j * b:(j + 1) * b, k * b:(k + 1) * b].mean(axis=(0, 1, 2))
        return out

    print("\n[2] coarse-grain g: macro strain emerges; residual + random anisotropy fall as 1/sqrt(N_block)")
    print("    b   N=b^3   resid-std   corr(conformal)   <aniso>")
    Ns, resids, anisos, corrs = [], [], [], []
    for b in (2, 3, 4, 6, 8):
        cg = blocks(g, b); cM = blocks(gM, b)
        resid = np.sqrt(((cg - cM) ** 2).sum(axis=(-1, -2))).std()
        corr = np.corrcoef(conf(cg).ravel(), conf(cM).ravel())[0, 1]; an = aniso(cg).mean()
        print(f"    {b}   {b**3:4d}    {resid:.4f}      {corr:.3f}             {an:.4f}")
        Ns.append(b ** 3); resids.append(resid); anisos.append(an); corrs.append(corr)
    sresid = np.polyfit(np.log(Ns), np.log(resids), 1)[0]; saniso = np.polyfit(np.log(Ns), np.log(anisos), 1)[0]
    ok(-0.6 < sresid < -0.4, f"[2] residual ~ N^{sresid:.2f} ~ 1/sqrt(N): the macro field emerges as fluctuations average out")
    ok(corrs[-1] > 0.95, f"[2] coarse conformal factor recovers the macro strain (corr={corrs[-1]:.3f} at b=8)")
    ok(-0.6 < saniso < -0.4, f"[3] random anisotropy ~ N^{saniso:.2f} ~ 1/sqrt(N): washes out under coarse-graining")

    print("\n[4] control: coarse-graining the FRAME e then forming g=<e>^T eta <e> is corrupted by gauge noise")
    ratios = []
    for b in (2, 8):
        ce = blocks(e, b); ge = np.einsum('...ai,ab,...bj->...ij', ce, eta, ce); cM = blocks(gM, b)
        err_e = np.sqrt(((ge - cM) ** 2).sum(axis=(-1, -2))).mean()
        err_g = np.sqrt(((blocks(g, b) - cM) ** 2).sum(axis=(-1, -2))).mean()
        print(f"    b={b}: <e>-route error={err_e:.3f}  g-route error={err_g:.3f}  ratio {err_e/err_g:.0f}x")
        ratios.append(err_e / err_g)
    ok(ratios[1] > ratios[0] > 5, "[4] <e>-route stays corrupted (ratio grows with block size) -> MUST coarse-grain g, not e")

    print("\n[verdict] Phase 1 -- the macro strain field exists and is recovered by gauge-invariant coarse-graining:")
    print("  - g = e^T eta e drops the local-Lorentz gauge noise (the key Phase 1 insight; <e> fails);")
    print("  - block-averaging g recovers the smooth macro strain (corr -> 0.98) with noise + random anisotropy")
    print("    falling as 1/sqrt(N_block) -> the macro field is well-defined in the continuum limit.")
    print("  HONEST SCOPE: kinematic only. Systematic lattice anisotropy washout = the separate Route C dispersion")
    print("  result (~(ka)^2); the DYNAMICS (induced-dominance, the crux) is Phase 2 -- the graviton <TT> self-")
    print("  energy, the gravitational sister of the failed T-R2 photon. Phase 1 supplies the field Phase 2 tests. exit 0")


if __name__ == "__main__":
    main()
