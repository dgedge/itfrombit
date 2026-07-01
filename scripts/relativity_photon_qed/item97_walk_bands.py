#!/usr/bin/env python3
"""
DRIFT K4 / item-97 dig: the real O_h Bloch bands of the walk W = S.C.

Goal: test whether the FULL unitary walk (not route-(b)'s simplified body-diagonal kinetic
kernel) resolves the doubling / gaps the mirror -> a positive for continuum SMG closure.

Finding (honest, negative-reinforcing):
  1. The LITERAL body-diagonal W=S.C is not cleanly computable: the shift's 8 body-diagonal
     directions (3 bridge bits) vs the chi(x)I3 Dirac spinor (4-dim) map is the OPEN item-97
     question itself -- not a shortcut around it.
  2. Two grounded RECONSTRUCTIONS exist; both lean negative:
       (b) body-diagonal kinetic sum (route b)         -> 4 massive even-corner species;
       (this) 3D axis-product unitary walk, real coin  -> 4 even-corner GAPLESS RIGHT-handed
              doublers (the zero-controlled-CNOT coin masses chi=0 only).
  So the walk-kernel "4 species" are pinned as right-handed doublers, ungapped by the chiral
  coin -- strictly worse than route (b)'s massive mirror. Continuum SMG closure stays negative.

Construction (documented): ANCHOR §3.5 alpha_i = sx^chi (x) sigma_i; coin = zero-controlled CNOT
(fire U(m) on I3 iff chi=0); the 3D walk is the axis-product extension of gw_walk_doubling's
2D W2 = Sy Sx C. Needs only numpy. Self-asserting on the guardrails + key findings.
"""
import numpy as np, itertools, statistics
np.set_printoptions(suppress=True)
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)
kron = np.kron

a1, a2, a3 = kron(sx, sx), kron(sx, sy), kron(sx, sz)   # alpha_i = sx^chi (x) sigma_i
chi_op = kron(sz, I2)                                    # chi = sz^chi  (+1 left, -1 right)
P0 = kron((I2 + sz) / 2, I2)
P1 = kron((I2 - sz) / 2, I2)


def C(m):  # zero-controlled CNOT coin: fire U(m) on I3 iff chi=0
    return P0 @ (np.cos(m) * np.eye(4) - 1j * np.sin(m) * kron(I2, sx)) + P1


def Si(ki, a):  # = exp(-i ki a), a^2 = I
    return np.cos(ki) * np.eye(4) - 1j * np.sin(ki) * a


def W3(k, m):   # 3D axis-product walk Sz Sy Sx C
    return Si(k[2], a3) @ Si(k[1], a2) @ Si(k[0], a1) @ C(m)


def qe(W):      # |omega| in [0, pi]
    return np.sort(np.abs(np.angle(np.linalg.eigvals(W))))


def hd(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


m = 0.30
THRESH = 0.15

hd("GUARDRAILS vs gw_walk_doubling (documented)")
def W1(k, mm):
    return Si(k, a1) @ C(mm)
w0 = qe(W1(0.0, m))
print(f"  1D W=Sx C, k=0: |omega| = {np.array2string(w0, precision=3)}")
print("  -> chi=0 (left) massed at omega~m; chi=1 (right) UNMASSED by the zero-ctrl coin (omega~0).")
assert sum(x < THRESH for x in w0) == 2 and sum(abs(x - m) < 1e-6 for x in w0) == 2

hd("(A) 3D walk W = Sz Sy Sx C : omega~0 Dirac points over the 8 BZ corners")
corners = list(itertools.product([0.0, np.pi], repeat=3))
total = 0
even_counts, odd_counts = [], []
for c in corners:
    w = qe(W3(c, m))
    nz = int(np.sum(w < THRESH))
    total += nz
    parity = int(round(sum(np.array(c) / np.pi))) % 2
    (even_counts if parity == 0 else odd_counts).append(nz)
    print(f"   k/pi={tuple(int(round(x/np.pi)) for x in c)} parity={parity}: "
          f"|omega|={np.array2string(w, precision=3)}  near-0 modes={nz}")
print(f"  total omega~0 modes = {total}  (even corners: {even_counts}, odd corners: {odd_counts})")
assert total == 8 and even_counts == [2, 2, 2, 2] and odd_counts == [0, 0, 0, 0]

hd("(B) chirality of the gapless modes")
gapless_chi = []
for c in corners:
    val, vec = np.linalg.eig(W3(c, m))
    om = np.abs(np.angle(val))
    for j in range(4):
        if om[j] < THRESH:
            v = vec[:, j]
            gapless_chi.append(np.real(v.conj() @ chi_op @ v / (v.conj() @ v)))
print(f"  <chi> of the {len(gapless_chi)} gapless modes: mean = {statistics.mean(gapless_chi):+.3f}")
print(f"  -> all gapless modes are RIGHT-handed (chi=-1); the coin masses LEFT (chi=+1) only.")
assert len(gapless_chi) == 8 and all(x < -0.99 for x in gapless_chi)

hd("(C) real asymmetric coin vs route-(b) symmetric mass")
print("  route (b) (H = D_kin + m*beta, SYMMETRIC): 4 even-corner species, each massive (gap 2m).")
print("  real walk (W=S.C, zero-ctrl coin, ASYMMETRIC): 4 even-corner RIGHT-handed doublers, GAPLESS.")
print("  => the real chiral coin is strictly WORSE for mirror-gapping: it masses one chirality only,")
print("     leaving the mirror chirality fully gapless. The SMG would have to gap a gapless sector.")

hd("VERDICT - item 97 dig")
print("""  item 97 is the OPEN STRUCTURAL QUESTION (the 8-body-diagonal-direction <-> chi(x)I3 spinor
  map), NOT a tractable shortcut to a positive. The most faithful full-unitary reconstruction
  (3D axis-product walk, exact zero-controlled-CNOT coin) gives 4 even-corner GAPLESS RIGHT-handed
  doublers -- the walk-kernel '4 species' pinned as right-handed, ungapped by the chiral coin,
  strictly worse than route (b)'s massive mirror. Both grounded reconstructions reinforce the
  negative. Continuum SMG closure stays NEGATIVE-LEANING; bare theta_UV=0 + local SMG gap stand;
  the genuine decider remains the non-perturbative (TN/MC) computation.""")
print("\nALL ASSERTS PASSED.")
