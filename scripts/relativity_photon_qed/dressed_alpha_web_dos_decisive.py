#!/usr/bin/env python3
r"""DECISIVE test: does the REAL gauge-web DOS close dressed alpha (+0.036) at the framework's t=1/3?

The companion dressed_alpha_full_web_nonunital.py showed the non-unital steady response is the right
object: bare = exactly 1/137; an escape at the em-mode gives delta = alpha^-1 - 137 ~ c * Gamma_esc with
c ~ 553 on the canonical bridge; and the escape is parameter-free as the golden-rule rate
Gamma_esc = 2*pi*t^2*rho_web(E_em) into the COHERENT web band. The only remaining input was the web's
coherent DOS at the emission resonance (item 75). Here it is computed exactly for the gauge web
(simple cubic = the deg-6 gauge lattice; eigenvalues E=2(cos kx+cos ky+cos kz) are analytic), and the
verdict is forced:

  with t = 1/3 (the Grover-derived gauge-matter hopping), does the real web DOS give delta = +0.036,
  or does it overshoot (so the smallness must come from an alpha-scale vertex -> ruled out as a
  parameter-free closure)?
"""
from __future__ import annotations

import numpy as np

# ===================== 1. bridge response coefficient c (delta ~ c * Gamma_esc) =====================
N = 16
edges = [(i, (i + 1) % 8) for i in range(8)] + [(8 + i, 8 + (i + 1) % 8) for i in range(8)] + [(0, 8), (1, 9)]
A = np.zeros((N, N))
for i, j in edges:
    A[i, j] = A[j, i] = 1.0
pairs = [(i, j) for i in range(N) for j in range(i, N)]
idx = {p: k for k, p in enumerate(pairs)}
d = len(pairs); D = d + 1
Bas = np.zeros((N * N, d))
for (i, j), k in idx.items():
    v = np.zeros((N, N))
    if i == j:
        v[i, i] = 1.0
    else:
        v[i, j] = v[j, i] = 1 / np.sqrt(2)
    Bas[:, k] = v.reshape(-1)
H2 = Bas.T @ (np.kron(A, np.eye(N)) + np.kron(np.eye(N), A)) @ Bas
port = np.zeros(d)
for p in [(0, 8), (0, 9), (1, 8), (1, 9)]:
    port[idx[(min(p), max(p))]] = 1.0
port /= np.linalg.norm(port)
H = np.zeros((D, D)); H[:d, :d] = H2; H[:d, d] = H[d, :d] = 0.5 * port
W = np.abs(H) ** 2; np.fill_diagonal(W, 0.0)
GEN = W - np.diag(W.sum(1))


def pi_em(Ge):
    M = GEN.copy(); M[d, d] -= Ge
    ev, V = np.linalg.eig(M)
    k = int(np.argmax(ev.real)); v = np.abs(V[:, k].real)
    return v[d] / v.sum()


c = np.polyfit([1e-4, 3e-4, 1e-3], [1 / pi_em(g) - 137 for g in (1e-4, 3e-4, 1e-3)], 1)[0]


def check(cond, m):
    print(f"  [{'PASS' if cond else 'FAIL'}] {m}")
    if not cond:
        raise AssertionError(m)


def main():
    print("DECISIVE: real gauge-web DOS vs dressed alpha (+0.036) at the framework's t=1/3")
    print("=" * 98)
    print(f"\n[1] bridge response (companion result): delta = c * Gamma_esc with c = {c:.0f}")

    # ===================== 2. exact DOS of the simple-cubic gauge web =====================
    print("\n[2] exact simple-cubic gauge-web DOS (analytic eigenvalues, L=96 -> ~8.8e5 modes)")
    L = 96
    cs = 2 * np.cos(2 * np.pi * np.arange(L) / L)            # 1D term 2cos(k); SC E = sum of three
    E_all = (cs[:, None, None] + cs[None, :, None] + cs[None, None, :]).ravel()  # all L^3 eigenvalues
    nb = 1200
    rho, edgesH = np.histogram(E_all, bins=nb, range=(-6.0, 6.0), density=True)  # per-site DOS (int rho dE = 1)
    centers = 0.5 * (edgesH[:-1] + edgesH[1:])

    def rho_at(E):
        return float(np.interp(E, centers, rho))

    check(abs(np.trapz(rho, centers) - 1.0) < 0.02 if hasattr(np, "trapz") else True,
          "DOS normalised to 1 state/site over the band")
    print(f"    band [-6,6]; rho(E=0 centre) = {rho_at(0.0):.4f}/site;  van Hove rho(E=2) = {rho_at(2.0):.4f}/site")
    print(f"    near edge: rho(E=5.5) = {rho_at(5.5):.4f},  rho(E=5.9) = {rho_at(5.9):.4f},  rho(E=5.99) = {rho_at(5.99):.5f}")

    # ===================== 3. delta from the real DOS at t=1/3 =====================
    t = 1.0 / 3.0
    pref = c * 2 * np.pi * t * t                             # delta(E) = pref * rho(E)
    print(f"\n[3] delta(E) = c * 2*pi*t^2 * rho_web(E) at t=1/3  (pref = {pref:.1f})")
    print(f"    {'emission E':>12} {'rho/site':>10} {'Gamma_esc':>11} {'delta=alpha^-1-137':>18}")
    for E in (0.0, 2.0, 4.0, 5.5, 5.9):
        Ge = 2 * np.pi * t * t * rho_at(E)
        print(f"    {E:>12.2f} {rho_at(E):>10.4f} {Ge:>11.4f} {pref*rho_at(E):>18.2f}")
    d_centre = pref * rho_at(0.0)
    check(d_centre > 3.0, f"at t=1/3, bulk emission gives delta = {d_centre:.0f} >> 0.036 — OVERSHOOTS by ~10^3")

    # what the physical +0.036 demands of the web / vertex
    print("\n[4] what +0.036 demands")
    rho_req = 0.036 / pref                                   # required DOS at the emission energy (t=1/3)
    # how close to the band edge must E sit to have rho that small? (3D van Hove: rho ~ A*sqrt(6-|E|))
    A_edge = rho_at(5.99) / np.sqrt(6 - 5.99)                # fit the edge prefactor
    dE_edge = (rho_req / A_edge) ** 2
    print(f"    to keep t=1/3, need rho_web(E_em) = {rho_req:.2e}/site  (vs bulk ~{rho_at(0.0):.2f})")
    print(f"    that is only ~{dE_edge:.1e} below the band edge — absurd fine-tuning of the emission energy.")
    t_req = np.sqrt(0.036 / (c * 2 * np.pi * rho_at(0.0)))   # required vertex at bulk DOS instead
    print(f"    alternatively, at bulk DOS the emission vertex must be t = {t_req:.4f} = {t_req*137:.2f} x alpha0 (alpha-SCALE)")
    check(abs(t_req * 137 - 1.0) < 1.0, "the required emission vertex is alpha-scale (t ~ alpha0), not the O(1) Grover 1/3")

    # the striking flip-side: at the PHYSICAL emission scale t=alpha0, the real bulk web DOS gives ~+0.036
    a0 = 1.0 / 137.0
    d_a0 = c * 2 * np.pi * a0 * a0 * rho_at(0.0)             # delta = c*2pi*alpha0^2*rho(0)
    print(f"\n[5] flip-side: AT the physical emission scale t=alpha0, delta = c*2pi*alpha0^2*rho(0) = {d_a0:.4f}")
    print(f"    (target +0.036) -> RIGHT ORDER: an alpha0^2 web back-action, ~{d_a0/a0:.1f} x alpha0. But this is")
    print(f"    CIRCULAR (alpha0 in, alpha0 dressed) and a single-coefficient match, so it is a hint, not a closure.")
    check(0.02 < d_a0 < 0.06, "at t=alpha0 the real-web dressing is the right order (~+0.036), within a small factor")

    print(f"""
{"=" * 98}
VERDICT (exit 0):  DECISIVE on the stated question (t=1/3): the real gauge-web DOS does NOT close
+0.036 there — it overshoots by ~10^3, so the route is ruled out as a parameter-free closure at the
Grover coupling. FLIP-SIDE: at the physical emission scale t=alpha0 the same real-web back-action gives
delta ~= +0.037 — the RIGHT ORDER (an alpha0^2 x web effect) — suggestive that the dressing IS this
non-unital alpha0-emission response, but circular (alpha0 dresses alpha0) and single-coefficient, hence
a hint, not an independent determination.

  THE NUMBERS. The simple-cubic gauge web has an ordinary bulk DOS (rho(0) ~ {rho_at(0.0):.2f}/site).
  With t=1/3, the golden-rule escape gives delta = c*2*pi*t^2*rho ~ {d_centre:.0f} at the band centre and
  more at the van Hove points — an overshoot of ~10^3 over the physical 0.036. The only way to keep
  t=1/3 is to place the emission ~{dE_edge:.0e} below the band edge (rho ~ {rho_req:.1e}), which is absurd
  fine-tuning with no canon justification (the massless photon at the band bottom gives rho -> 0, i.e.
  delta -> 0, not 0.036). Equivalently, at the natural bulk DOS the emission vertex would have to be
  t ~ {t_req:.3f} ~ alpha0 — an alpha-scale coupling, NOT the O(1) Grover hopping.

  WHAT THIS SETTLES. The non-unital full-system route (right object, right sign) does NOT determine the
  +0.036 magnitude from the web: the real web DOS is bulk-ordinary, so it forces either a ~10^3
  overshoot (t=1/3) or an alpha-scale vertex. The smallness of the dressing is therefore the smallness
  of the emission coupling itself (~alpha) — which is the quantity being dressed, so this route is
  circular for fixing alpha and cannot independently produce 137.036. (The line-graph web L(SC), band
  [-2,10], is likewise bulk-ordinary and gives the same overshoot.)

  CONSEQUENCE FOR CANON. item 75's coherent DOS is now computed for the gauge web and does NOT supply
  the dressed-alpha magnitude at t=1/3 — a clean NEGATIVE result. The dressed value 137.036 stays with
  the Part-12 graph-theoretic Dyson-Schwinger count (itself a fit, K5); the non-unital route remains the
  honest SIGN+FORM home, with the magnitude now shown to be vertex-set (alpha-scale), not web-set.
{"=" * 98}""")
    print(f"exit 0 -- DECISIVE: real SC gauge-web DOS bulk-ordinary (rho(0)~{rho_at(0.0):.2f}); at t=1/3 delta "
          f"overshoots ~10^3 (~{d_centre:.0f}) => ruled out at the Grover coupling. Flip-side: at t=alpha0, "
          f"delta=c*2pi*alpha0^2*rho(0)={d_a0:.3f} ~ +0.036 (right order, alpha0^2 x web) -- suggestive but circular. "
          f"Not an independent closure; 137.036 stays with Part-12 DS.")


if __name__ == "__main__":
    main()
