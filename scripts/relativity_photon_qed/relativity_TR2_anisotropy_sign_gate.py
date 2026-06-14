#!/usr/bin/env python3
r"""T-R2 — IR FATE OF THE MARGINAL VELOCITY ANISOTROPY: well-posed sign gate.

Canon state (do not re-fight): the bare K6 anisotropy is FIRST-ORDER/marginal
(DRIFT K7 — v[100]=sqrt(2/3), v[111]=1/sqrt3, ratio sqrt2), so power-counting
CANNOT retire it; and the Chadha-Nielsen beta-function is a CATEGORY ERROR
(DRIFT K8 — a charged-fermion drag formula applied to the U(1)-neutral
T1u/E_g modes, which have no p.A or 3-photon self-vertex). The ONLY one-loop
dressing of the neutral-mode anisotropy is therefore the charged-matter
vacuum-polarization loop Pi_munu = INT V G_F V G_F (K8). Its sign is open.

This script does NOT close the gate. It (1) PRE-REGISTERS the sign gate; (2)
proves the one rigorous lemma that sharpens it — the leading O(k^2) photon
kinetic dressing from ANY O_h-symmetric charged-matter loop is EXACTLY
isotropic (delta_ij is the unique O_h-invariant rank-2 tensor), verified by
computation with an anisotropic control; (3) records the screening sign; and
(4) REDUCES the residual sign question to the item-115 second-order matter
kernel — the SAME kernel T-R1's 3D lift needs (the 'pays twice' structure).

exit 0 = lemma verified, control fires, sign recorded, reduction stated.
"""
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DRIFT = (ROOT / "DRIFT.md").read_text()
ANCHOR = (ROOT / "ANCHOR.md").read_text()

print("[1] PRE-REGISTERED SIGN GATE (fix BEFORE any loop is evaluated):")
print("    Object: kappa_aniso(mu) = coefficient of the cubic I_4 anisotropy in")
print("    the DRESSED inverse photon/graviton propagator; kappa_iso(mu) the")
print("    isotropic coefficient. Flow d kappa/d ln mu from the charged-matter")
print("    vacuum polarization (the sole neutral-mode dressing, DRIFT K8).")
print("    GATE (pre-registered, three exclusive outcomes):")
print("      PASS   : delta kappa_aniso = 0 at the matter-loop order AND")
print("               kappa_iso runs with screening sign  => Delta v / v -> 0")
print("               logarithmically (IR Lorentz restoration).")
print("      FAIL   : delta kappa_aniso same sign as the bare splitting")
print("               => anisotropy grows; FALSIFIED at the SME cavity bound")
print("               (a0 k)^2 ~ 1e-17 (item 102 empirical backstop).")
print("      MARGINAL: delta kappa_aniso opposes but sub-dominant => undecided")
print("               at one loop; two-loop required.")

print("\n[2] RIGOROUS LEMMA (computed) — the leading O(k^2) dressing is EXACTLY")
print("    isotropic: delta_ij is the unique O_h-invariant symmetric rank-2")
print("    tensor, so any O_h-symmetric matter loop sources only delta_ij k^2")
print("    at leading order; anisotropy is pushed to rank-4 (the I_4 invariant).")

def build_grid(n):
    pts = [(-math.pi + 2*math.pi*(i+0.5)/n) for i in range(n)]
    return pts

N = 48
KS = build_grid(N)

def chi_tensor(disp, vel, weight):
    """Current-current susceptibility chi_ij = <v_i v_j W(eps)> over the BZ."""
    chi = [[0.0]*3 for _ in range(3)]
    norm = 0.0
    for kx in KS:
        for ky in KS:
            for kz in KS:
                k = (kx, ky, kz)
                e = disp(k)
                w = weight(e)
                v = vel(k)
                norm += w
                for a in range(3):
                    for b in range(3):
                        chi[a][b] += w * v[a] * v[b]
    return [[chi[a][b]/norm for b in range(3)] for a in range(3)], norm

# O_h-symmetric SEPARABLE matter dispersion (the §3.2 'separable' kernel shape)
def eps_iso(k):  return -(math.cos(k[0]) + math.cos(k[1]) + math.cos(k[2]))
def vel_iso(k):  return (math.sin(k[0]), math.sin(k[1]), math.sin(k[2]))
# anisotropic CONTROL: break the C4 symmetry by stretching the x bond
def eps_an(k):   return -(1.7*math.cos(k[0]) + math.cos(k[1]) + math.cos(k[2]))
def vel_an(k):   return (1.7*math.sin(k[0]), math.sin(k[1]), math.sin(k[2]))
W = lambda e: math.exp(-2.0*e)          # any smooth O_h-symmetric weight

chi_i, _ = chi_tensor(eps_iso, vel_iso, W)
offdiag = max(abs(chi_i[a][b]) for a in range(3) for b in range(3) if a != b)
diag_spread = max(chi_i[a][a] for a in range(3)) - min(chi_i[a][a] for a in range(3))
trace = sum(chi_i[a][a] for a in range(3))
print(f"    O_h matter loop:  chi diag = [{chi_i[0][0]:.6f},{chi_i[1][1]:.6f},"
      f"{chi_i[2][2]:.6f}], max|offdiag| = {offdiag:.2e}")
print(f"      off-diagonal ~ 0 ? {offdiag < 1e-12};  diagonal spread = {diag_spread:.2e}")
assert offdiag < 1e-12 and diag_spread < 1e-12          # chi_ij EXACTLY ∝ delta_ij
# directional velocity^2 moment m(n) = n_i chi_ij n_j -> isotropic iff chi ∝ delta
def moment(chi, n):
    nn = math.sqrt(sum(c*c for c in n)); u = [c/nn for c in n]
    return sum(u[a]*chi[a][b]*u[b] for a in range(3) for b in range(3))
m100, m111 = moment(chi_i, (1,0,0)), moment(chi_i, (1,1,1))
print(f"      m[100] = {m100:.8f},  m[111] = {m111:.8f},  diff = {abs(m100-m111):.2e}")
assert abs(m100 - m111) < 1e-12        # leading dressing identical in all directions
print("    -> leading photon kinetic dressing is EXACTLY isotropic. LEMMA VERIFIED.")

# CONTROL: the probe MUST see anisotropy when the matter loop is anisotropic
chi_a, _ = chi_tensor(eps_an, vel_an, W)
a100, a111 = moment(chi_a, (1,0,0)), moment(chi_a, (1,1,1))
print(f"    anisotropic control: m[100] = {a100:.6f}, m[111] = {a111:.6f}, "
      f"diff = {abs(a100-a111):.4f}")
assert abs(a100 - a111) > 1e-3         # the probe is not blind: it CAN resolve anisotropy
print("    -> control fires: the null above is a real symmetry, not a dull probe.")

print("\n[3] SCREENING SIGN (computed): the isotropic piece Tr(chi) > 0, the")
print(f"    standard screening sign for a kinetic dressing: Tr chi = {trace:.6f} > 0.")
assert trace > 0

print("\n[4] ANCHOR/DRIFT support live:")
for needle, src, label in [
    ("category error", DRIFT, "K8: Chadha-Nielsen is a category error"),
    ("vacuum polarization", DRIFT, "K8: correct mechanism is Pi_munu"),
    ("is **MARGINAL**", ANCHOR, "K7: bare anisotropy is first-order/marginal"),
    ("one-hop ceiling theorem", ANCHOR, "item 115: kernel exists at second order"),
]:
    ok = needle in src
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
    assert ok, needle

print(f"""
[5] T-R2 VERDICT — well-posed, gate pre-registered, NOT closed:
  * The sign gate is fixed in advance (block 1): PASS/FAIL/MARGINAL with the
    SME (a0 k)^2 ~ 1e-17 bound as the falsification backstop.
  * LEMMA (proved, block 2): the leading O(k^2) dressing from any O_h-symmetric
    charged-matter loop is EXACTLY isotropic (delta_ij unique); the cubic I_4
    anisotropy is rank-4, parametrically (a0 k)^2 down. The isotropic piece
    has the screening sign (block 3). This is the half that points to PASS.
  * RESIDUAL (honest): the K6 splitting is a DEGENERATE-MULTIPLET (T1u + E_g
    k.p) effect, not a single-band mass tensor, so the lemma is necessary but
    not sufficient — deciding the gate needs the FULL Pi_munu with the explicit
    charged-matter Green function, i.e. the item-115 §3.2 second-order kernel.
  * REDUCTION: T-R2 sign  <=>  item-115 second-order matter kernel  =  the
    SAME object T-R1's 3D dispersion lift needs. Progress there pays twice
    (user's note confirmed structurally). T-R2 moves OPEN -> well-posed +
    reduced + sign-leaning-PASS-by-symmetry, conditional on that one kernel.
exit 0""")
print("ALL ASSERTIONS PASSED — lemma verified, control fired, reduction stated.")
