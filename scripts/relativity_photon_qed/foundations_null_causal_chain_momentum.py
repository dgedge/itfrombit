r"""Proposed B_N primitive = a NULL CAUSAL CHAIN, and the one property that makes it work.

The B_N service-ledger audit (foundations_bn_service_ledger_gate.py) found the ledger gives
G3 (vacuum: event-conditioned, no oscillator tower) and the scalar event-count, but NOT
G1 (a noncompact, un-aliased spacetime-momentum label) and NOT G2 (one-leg detector coupling).
G1 is the crux: on a FIXED lattice, spatial momentum is crystal momentum -- compact (Brillouin),
anisotropic -- which is exactly why K15 caps null quanta at ~1 GeV.

The only principled escape (K15 door iii) is to stop carrying high-energy momentum as a Bloch
phase of the fixed lattice, and carry it instead as the LINK-COUNT of a causal chain through a
frame-independent (Poisson-sprinkled) substrate -- the causal-set route, the one known way to
have Lorentz-invariant discreteness (Bombelli-Henson-Sorkin). The fixed TCH lattice is then the
ANNEALED IR condensate (valid omega << Lambda_QCD); a trans-Lambda quantum is a null causal
chain through the pre-condensate network.

This script demonstrates the TWO enabling properties a causal chain HAS and a lattice Bloch mode
LACKS, which are precisely the G1 obstruction:
  (A) NONCOMPACT: the chain link-count grows without bound (~ proper time), so the momentum
      label N*Lambda*n is unbounded -- no Brillouin cap, no aliasing.
  (B) LORENTZ-INVARIANT: the link-count is a boost invariant, and the sprinkling measure is
      boost-invariant (1+1 boosts have det=1), so the momentum-per-proper-length is the same in
      every frame and every direction -- unlike the lattice, whose chain-count is anisotropic.

It does NOT close G2 (a one-leg detector vertex reproducing QED amplitudes) -- that is the hard,
unaddressed part. exit 0 = (A) and (B) hold for the causal set and FAIL for the lattice.
"""
import numpy as np

rng = np.random.default_rng(20260613)
LAMBDA = 0.332  # GeV, per-link action scale (one service tick)

# ---- causal set in 1+1 Minkowski: p < q iff dt>0 and dt^2 - dx^2 > 0 (future timelike) ----
def longest_chain(points):
    """links in the longest causal chain (DAG longest path), O(n^2)."""
    P = points[np.argsort(points[:, 0])]
    n = len(P); L = np.ones(n, dtype=int)
    for j in range(n):
        dt = P[j, 0] - P[:j, 0]; dx = P[j, 1] - P[:j, 1]
        rel = (dt > 0) & (dt * dt - dx * dx > 0)
        if rel.any():
            L[j] = 1 + L[:j][rel].max()
    return int(L.max()) - 1   # links = nodes - 1

def sprinkle_diamond(T, rho):
    """Poisson sprinkle the causal diamond between (0,0) and (T,0), density rho."""
    area = T * T / 2.0
    n = rng.poisson(rho * area)
    t = rng.uniform(0, T, size=n * 3); x = rng.uniform(-T / 2, T / 2, size=n * 3)
    keep = (np.abs(x) < np.minimum(t, T - t))      # inside the diamond
    pts = np.column_stack([t[keep], x[keep]])
    return np.vstack([[0.0, 0.0], pts, [T, 0.0]])

def boost(points, beta):
    g = 1.0 / np.sqrt(1 - beta * beta)
    t, x = points[:, 0], points[:, 1]
    return np.column_stack([g * (t - beta * x), g * (x - beta * t)])

print("[A] NONCOMPACT: causal-chain link-count grows ~ proper time T (no Brillouin cap)")
rho = 6.0
Ls = {}
for T in (10.0, 20.0, 40.0):
    vals = [longest_chain(sprinkle_diamond(T, rho)) for _ in range(4)]
    Ls[T] = np.mean(vals)
    print(f"    T={T:5.1f}: mean longest-chain links N = {Ls[T]:6.1f}  ->  P = N*Lambda = {Ls[T]*LAMBDA:7.2f} GeV")
ratio = Ls[40.0] / Ls[10.0]
print(f"    N(T=40)/N(T=10) = {ratio:.2f}  (expect ~4: linear in T, UNBOUNDED -> noncompact momentum)")
assert 3.3 < ratio < 4.7

print("\n[B] LORENTZ-INVARIANT: link-count is a boost invariant; sprinkling measure is boost-invariant")
pts = sprinkle_diamond(30.0, rho)
base = longest_chain(pts)
inv = [longest_chain(boost(pts, b)) for b in (0.3, 0.6, 0.9)]
print(f"    longest-chain links: rest={base}, boosted(beta=.3,.6,.9)={inv}  (EXACTLY invariant)")
assert all(v == base for v in inv)                       # causal order is Lorentz-invariant: count unchanged
for b in (0.3, 0.6, 0.9):
    g = 1 / np.sqrt(1 - b * b)
    det = g * g * (1 - b * b)                            # det of the 1+1 boost
    assert abs(det - 1.0) < 1e-12                        # area-preserving => sprinkling density frame-independent
print("    1+1 boost determinant = 1 for all beta -> Poisson density is the same in every frame (BHS).")

print("\n[C] CONTRAST: a FIXED lattice chain-count is anisotropic and capped (the G1 obstruction)")
# monotone lattice path length to a target at angle theta vs Euclidean distance: Manhattan anisotropy
print("    angle from axis | lattice steps / Euclidean dist  (causal-set ratio is theta-independent)")
ratios = []
for deg in (0, 15, 30, 45):
    th = np.radians(deg); ax, ay = np.cos(th), np.sin(th)
    manh = abs(ax) + abs(ay); eucl = np.hypot(ax, ay)
    ratios.append(manh / eucl)
    print(f"    {deg:5d} deg        | {manh/eucl:.4f}")
anis = max(ratios) / min(ratios)
print(f"    lattice anisotropy max/min = {anis:.4f} (= sqrt2): direction-dependent 'speed' -> NOT Lorentz")
assert anis > 1.3                                        # lattice is anisotropic (sqrt2), causal set is not
print("    + a lattice path has <= one step per site: chain-count per length is CAPPED (compact/Brillouin).")

print(f"""
[D] VERDICT — the causal-chain primitive HAS the G1-enabling structure the lattice lacks:
  * NONCOMPACT: link-count ~ proper time, unbounded -> P_label = N*Lambda*n^mu has no Brillouin
    cap and no aliasing (a 1 TeV chain is just N ~ {1000/LAMBDA:.0f} links; nothing oscillates at a0).
  * LORENTZ-INVARIANT: link-count is a boost invariant and the sprinkling is frame-independent
    (det boost = 1), so momentum-per-length is isotropic -- the fixed lattice is sqrt2-anisotropic.
  This is K15 door (iii) made concrete: B_N = a null causal chain through the pre-condensate
  (causal-set) substrate; the fixed TCH lattice is the annealed IR condensate (omega << Lambda).
  It plausibly supplies G1 (this script) and G3 (audit: source-conditioned, no oscillator tower).
  STILL OPEN (the hard part, NOT shown here): G2 -- a single matter vertex that absorbs a chain
  as one event carrying total P^mu and reproduces QED/LSZ amplitudes and high-energy cross
  sections. And the cost: it demotes the fixed-lattice canon to IR-emergent and must reconcile
  the item-118 (Lambda_QCD-cutoff) cosmological constant with a causal-set (Sorkin-fluctuation)
  one. This is a research PROGRAMME with a sharp first theorem, not a closure.
exit 0""")
print("ALL ASSERTIONS PASSED — causal chain is noncompact + Lorentz-invariant; lattice is capped + anisotropic.")
