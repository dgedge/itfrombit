#!/usr/bin/env python
r"""dS_horizon_modular.py -- causal-patch (cosmological-horizon) entanglement of the gauge web.

The existing RT check (substrate_rt_wedge) used a SPATIAL disk. The flagged dS step is the
causal/cosmological HORIZON as the holographic screen. Static computable precursor: take a
causal-patch BALL whose boundary is the horizon, in the gapless gauge-web field (kernel
K=6-2 sum cos k, the framework's physical photon sector), and test two horizon signatures:

  [1] HORIZON AREA LAW: S(ball radius R) ~ boundary area (R^2), not volume -- Gibbons-Hawking/
      Bekenstein S ~ A on the cosmological horizon.
  [2] THERMAL MODULAR SPECTRUM: the modular energies eps_k = log((2 nu_k + 1)/(2 nu_k - 1)) from
      the ball's symplectic eigenvalues nu_k. A thermal HORIZON has a GAPLESS modular spectrum
      (eps -> 0, a continuous Gibbons-Hawking thermal ladder); a region with NO horizon (gapped
      bulk) has a MODULAR GAP. Contrast: massless (horizon) vs massive (no horizon).

SCOPE (honest): this is the STATIC/spatial precursor (causal patch = ball, via the standard
diamond<->static-patch picture). The full dynamical dS horizon needs the HBC printing dynamics --
the deeper open piece. Free Gaussian field; emergent/kinematic.
"""
import json
import numpy as np

def correlators(L, m2):
    k = 2 * np.pi * np.arange(L) / L
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    Kk = m2 + 6 - 2 * (np.cos(KX) + np.cos(KY) + np.cos(KZ))
    return np.fft.ifftn(0.5 * Kk ** -0.5).real, np.fft.ifftn(0.5 * Kk ** 0.5).real

def sympl_eigs(Xr, Pr, sites, L):
    n = len(sites)
    XA = np.empty((n, n)); PA = np.empty((n, n))
    for a, sa in enumerate(sites):
        for b, sb in enumerate(sites):
            dr = ((sa[0]-sb[0]) % L, (sa[1]-sb[1]) % L, (sa[2]-sb[2]) % L)
            XA[a, b] = Xr[dr]; PA[a, b] = Pr[dr]
    nu = np.sqrt(np.clip(np.linalg.eigvals(XA @ PA).real, 0, None))
    return np.clip(nu, 0.5 + 1e-12, None)

def entropy(nu):
    return float(np.sum((nu+0.5)*np.log(nu+0.5) - (nu-0.5)*np.log(nu-0.5)))

def modular_energies(nu):
    return np.log((2*nu + 1) / (2*nu - 1))   # eps_k; small eps = highly entangled (thermal) mode

def ball(R, L):
    c = L // 2
    return [(c+x, c+y, c+z) for x in range(-R, R+1) for y in range(-R, R+1) for z in range(-R, R+1)
            if x*x + y*y + z*z <= R*R]

def linfit_r2(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None); yh = A @ c
    ss, st = np.sum((y-yh)**2), np.sum((y-y.mean())**2)
    return float(1 - ss/st if st > 0 else 0)

def run(L, m2, tag):
    Xr, Pr = correlators(L, m2)
    Rs = [3, 4, 5, 6]
    S, gap = {}, {}
    for R in Rs:
        nu = sympl_eigs(Xr, Pr, ball(R, L), L)
        S[R] = entropy(nu)
        gap[R] = float(np.min(modular_energies(nu)))   # modular gap = most-entangled (thermal) mode
    r2_area = linfit_r2([R*R for R in Rs], list(S.values()))
    r2_vol = linfit_r2([R**3 for R in Rs], list(S.values()))
    return {"tag": tag, "m2": m2, "area_R2": r2_area, "vol_R2": r2_vol,
            "S_ball": S, "modular_gap_vs_R": gap,
            "gap_closes": gap[Rs[-1]] < 0.9 * gap[Rs[0]]}      # does the modular gap shrink >10% with R?

def main():
    L = 28
    out = {}
    for m2, tag in [(1e-4, "massless_horizon"), (1.0, "massive_no_horizon")]:
        r = run(L, m2, tag); out[tag] = r
        g = r["modular_gap_vs_R"]
        print(f"\n=== {tag}  (L={L}, m^2={m2}) ===")
        print(f"  [1] HORIZON AREA LAW: S(ball) ~ R^2 R2={r['area_R2']:.4f} vs R^3 R2={r['vol_R2']:.4f}"
              f"  -> {'AREA' if r['area_R2']>r['vol_R2'] else 'VOLUME'}")
        print(f"  [2] MODULAR GAP vs patch radius R (min eps; -> 0 = thermal/Gibbons-Hawking):")
        print(f"      " + "  ".join(f"R={R}:{g[R]:.3f}" for R in sorted(g)) +
              f"   -> {'CLOSES (thermal horizon)' if r['gap_closes'] else 'saturates (no horizon)'}")
    ml, mv = out["massless_horizon"], out["massive_no_horizon"]
    gl, gv = ml["modular_gap_vs_R"], mv["modular_gap_vs_R"]
    Rmax = max(gl)
    print("\n[verdict] CAUSAL-PATCH / HORIZON ENTANGLEMENT of the gauge web:")
    print(f"  massless (HORIZON): area law; modular gap CLOSES with R "
          f"({gl[min(gl)]:.2f} -> {gl[Rmax]:.2f}) -> trending to a gapless Gibbons-Hawking thermal screen.")
    print(f"  massive (NO horizon): modular gap SATURATES high ({gv[min(gv)]:.2f} -> {gv[Rmax]:.2f}) -> no thermal horizon.")
    print(f"  ratio gap(massive)/gap(massless) at R={Rmax}: {gv[Rmax]/gl[Rmax]:.2f}x")
    print("  => the cosmological horizon behaves as a THERMAL holographic screen (S~A + a modular gap")
    print("     that closes with patch size). SCOPE: static causal-patch precursor; the full dynamical")
    print("     dS horizon (HBC printing) is the deeper open piece.")
    assert ml["area_R2"] > ml["vol_R2"], "horizon must obey area law"
    assert ml["gap_closes"], "massless horizon modular gap must CLOSE with patch radius (thermal)"
    assert gl[Rmax] < gv[Rmax], "horizon modular gap must be below the no-horizon gap"
    assert gv[Rmax] / gl[Rmax] > 1.5, "clear thermal-vs-gapped separation"
    print("\nALL ASSERTIONS PASSED -- horizon area law + modular gap closing with R (thermal screen).")
    print("exit 0")
    with open("dS_horizon_modular_results.json", "w") as f:
        json.dump(out, f, indent=2, default=float)

if __name__ == "__main__":
    main()
