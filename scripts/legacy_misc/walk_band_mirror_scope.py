#!/usr/bin/env python3
"""
DRIFT K4 route (b): the walk Bloch band structure -> mirror kinetic scale.

Question (the "t vs t_c" decider for SMG mirror-gap survival under dynamical gauge):
  Does the TCH walk W = S.C supply a physical mirror hopping t below the SMG gap t_c~2
  (gap survives), or comparable/above (gap closes)?

Method: build the walk's effective Bloch Hamiltonian on the chi(x)I3 Dirac spinor
(ANCHOR §3.5; Euclidean gamma5 = sigma_y^chi (x) I, K4), with the documented body-diagonal
kinetic operator (gw_nogo: V(k) = 8(sin kx cos ky cos kz, +cyc)) and the walk-sourced Wilson
scalar (walk_kernel_overlap: W_body = 8(1 - cos kx cos ky cos kz)). Guardrail against the
documented results, then read off the kinetic scale and what it does/does not decide.

Result (honest, negative-leaning):
  - GUARDRAILS reproduce canon: {gamma5, H_kin} = 0 exactly; the walk's own kernel gives
    4 degenerate light species at the EVEN corners (the documented walk-kernel result).
  - The body-diagonal kinetic scale is 8x ENHANCED vs a face-direction lattice; the mirror
    modes carry the SAME (unsuppressed) bandwidth -> NO t << t_c hierarchy (negative signal).
  - CRUCIALLY: the free Bloch band structure CANNOT decide gap survival, which is the
    non-perturbative PMS-vs-symmetric question. So route (b) is not the decider: t is
    computable (enhanced, no suppression); t_c (gap survival) is non-perturbative.
  => continuum SMG closure stays negative-leaning; the decider is the TN/MC computation.

Needs only numpy. Self-asserting on the guardrails.
"""
import numpy as np, itertools
np.set_printoptions(suppress=True)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)
kron = np.kron

# documented Clifford set on chi (x) I3 (ANCHOR §3.5 + K4 Euclidean gamma5)
g5 = kron(sy, I2)
beta = kron(sz, I2)
alpha = [kron(sx, sx), kron(sx, sy), kron(sx, sz)]


def V(k):  # body-diagonal kinetic (gw_nogo documented closed form)
    kx, ky, kz = k
    return 8 * np.array([np.sin(kx) * np.cos(ky) * np.cos(kz),
                         np.cos(kx) * np.sin(ky) * np.cos(kz),
                         np.cos(kx) * np.cos(ky) * np.sin(kz)])


def Wbody(k):  # walk-sourced finite-range Wilson scalar (walk_kernel_overlap)
    kx, ky, kz = k
    return 8 * (1 - np.cos(kx) * np.cos(ky) * np.cos(kz))


def Hkin(k):
    return sum(V(k)[i] * alpha[i] for i in range(3))


def Hfull(k, m, r):
    return Hkin(k) + (m + r * Wbody(k)) * beta


def band(k, m, r):
    return np.sort(np.linalg.eigvalsh(Hfull(k, m, r)))


def hd(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


hd("GUARDRAILS (reproduce documented results)")
rng = np.random.default_rng(0)
acom = max(np.linalg.norm(g5 @ Hkin(p) + Hkin(p) @ g5) for p in rng.uniform(-np.pi, np.pi, (100, 3)))
print(f"  {{gamma5, H_kin}} = {acom:.1e}  (massless chiral kinetic, expect ~0)")
assert acom < 1e-12
corners = list(itertools.product([0.0, np.pi], repeat=3))
m, r = 0.3, 1.0
print(f"  corner masses (min|E|) at m={m}, r={r}:")
even, odd = [], []
for c in corners:
    g = min(abs(band(c, m, r)))
    par = int(round(sum(np.array(c) / np.pi))) % 2
    (even if par == 0 else odd).append(g)
    print(f"    k/pi={tuple(int(round(x/np.pi)) for x in c)} parity={par}: min|E|={g:.3f}")
print(f"  -> EVEN corners (4): min|E|={np.mean(even):.3f}=m (light, degenerate species);"
      f" ODD corners (4): min|E|={np.mean(odd):.3f}=m+16r (lifted)")
assert abs(np.mean(even) - m) < 1e-9 and len(even) == 4 and min(odd) > 10
print("  CONFIRMED: walk's own kernel -> 4 degenerate light species (the documented result).")

hd("(A) the walk's KINETIC scale vs a standard face-direction lattice")


def Vface(k):
    return np.array([np.sin(k[0]), np.sin(k[1]), np.sin(k[2])])


def gradnorm(Vfun, k0=(0.0, 0.0, 0.0)):
    h = 1e-5
    kp = np.array(k0, float); kp[0] += h
    km = np.array(k0, float); km[0] -= h
    return np.linalg.norm((Vfun(kp) - Vfun(km)) / (2 * h))


gb, gf = gradnorm(V), gradnorm(Vface)
print(f"  body-diagonal kinetic coefficient |dV/dk|_0   = {gb:.3f}")
print(f"  face-direction kinetic coefficient |dV/dk|_0   = {gf:.3f}")
print(f"  ratio (body-diagonal / face) = {gb/gf:.2f}  -> body-diagonal kinetics ENHANCED")
ks = rng.uniform(-np.pi, np.pi, (4000, 3))
bw = max(max(abs(band(k, 0.0, 0.0))) for k in ks)
print(f"  massless kinetic bandwidth max|E| over BZ = {bw:.2f}   (vs SMG local gap = 2)")
assert abs(gb - 8.0) < 1e-6 and abs(gf - 1.0) < 1e-6

hd("(B) is there a NATURAL suppression of the mirror's kinetics?  (the t << t_c hope)")
print("  The mirror modes (W!=chi partners) ride the SAME dispersion as the physical mode:")
print("  the 4 even-corner species are DEGENERATE (all mass m), none kinetically suppressed.")
print(f"  -> no k-dependent suppression at the mirror; the walk supplies NO t << gap hierarchy.")
print(f"  -> if anything the body-diagonal geometry ENHANCES the kinetic scale (x{gb/gf:.0f},")
print(f"     plus nodal LINES not points, gw_nogo) -- a NEGATIVE-direction signal for gapping.")

hd("(C) what the band structure CANNOT decide (the honest limit of route (b))")
print("  Gap survival under the SMG interaction is whether a symmetric 4-fermion term gaps the")
print("  mirror WITHOUT condensing (symmetric phase) vs flowing to PMS/symmetry-breaking at the")
print("  physical coupling. That is a NON-PERTURBATIVE many-body question. H_full(k) here is a")
print("  FREE (quadratic) Bloch Hamiltonian -- it fixes the kinetic scale (t) but cannot evaluate")
print("  the interacting gap. So 't vs t_c' conflated a free-field scale (computable: enhanced,")
print("  no mirror suppression) with a non-perturbative gap-survival question (NOT computable here).")

hd("VERDICT - route (b)")
print("""  Route (b) is COMPUTED and GUARDRAILED, and it does NOT close to positive:
   - the walk band structure is reproduced (4 even-corner degenerate species, body-diagonal
     kinetics, {gamma5,H_kin}=0);
   - the walk supplies a relativistic kinetic term with NO suppression at the mirror modes and
     an ENHANCED (x8, nodal-line) body-diagonal scale -> a negative-direction signal;
   - crucially, the band structure (free-field) CANNOT decide gap survival, which is the
     non-perturbative PMS-vs-symmetric question. So the 't vs t_c' framing is refined: t (the
     kinetic scale) is computable and shows no t<<t_c hierarchy; t_c (gap survival) is
     non-perturbative and beyond band structure.
  NET: (b) strengthens the negative lean (enhanced kinetics, no mirror suppression) and proves
  the genuine decider is the NON-PERTURBATIVE computation (route 2: tensor-network / Monte Carlo),
  which the cutoff-scaling already found falling. Continuum SMG closure stays NEGATIVE-LEANING;
  band structure is not the decider.""")
print("\nALL ASSERTS PASSED.")
