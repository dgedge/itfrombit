#!/usr/bin/env python3
"""
SMG operator check: is the proposed 10-channel quartic H = U sum_a |psi^T C Gamma_a psi|^2
a viable gapping operator for a SINGLE SO(10) 16, or does it vanish?

Rep theory: 16 (x) 16 = 10_s (+) 120_a (+) 126_s  (no singlet -> the user's correct point:
no 2-body symmetric mass exists). But the 10 and 126 are in the SYMMETRIC part and 120 in
the ANTISYMMETRIC part. A fermion bilinear psi_i psi_j is ANTISYMMETRIC (anticommuting).
So:
  * if the 16x16 intertwiner C*Gamma_a (the 10) is SYMMETRIC  -> psi^T C Gamma_a psi == 0
    -> the proposed H_SMG is the ZERO operator (cannot gap anything).
  * the nonzero bilinear-squared singlet channel is the ANTISYMMETRIC 120 (C*Gamma_abc).

This script builds the SO(10) Clifford algebra on 5 flavour qubits (Dirac spinor 2^5=32),
restricts to the even-chirality 16 (= the P3 even-minus weight basis), and checks the
symmetry of the 10-channel (CGamma_a) and 120-channel (CGamma_abc) intertwiners directly.
Self-validating: Clifford algebra, C^T=-C, and C Gamma_a C^{-1} = +-Gamma_a^T are asserted.

numpy only.
"""
import numpy as np
import itertools

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
N = 5  # flavour qubits -> SO(10)


def op(local):  # local: dict qubit->2x2, tensor over N qubits (qubit 0 = leftmost)
    M = np.array([[1]], complex)
    for k in range(N):
        M = np.kron(M, local.get(k, I2))
    return M


# SO(10) gamma matrices via Jordan-Wigner: Gamma_{2k-1}=Z..Z X_k, Gamma_{2k}=Z..Z Y_k
gammas = []
for k in range(N):
    pre = {j: Z for j in range(k)}
    gammas.append(op({**pre, k: X}))
    gammas.append(op({**pre, k: Y}))
G = gammas  # G[0..9]

# (1) Clifford algebra {Gamma_a,Gamma_b}=2 delta_ab
for a in range(10):
    for b in range(10):
        anti = G[a] @ G[b] + G[b] @ G[a]
        assert np.allclose(anti, 2 * (a == b) * np.eye(32)), f"Clifford fails at {a},{b}"
print("(1) Clifford algebra {Gamma_a,Gamma_b}=2 delta_ab  ........ OK (10 gammas on 2^5=32)")

# chirality = product of all gammas (up to phase) == Z^{(x)5}; even-chirality = +1 eigenspace
chir = op({k: Z for k in range(N)})
for a in range(10):
    assert np.allclose(G[a] @ chir + chir @ G[a], 0), "gamma should flip chirality"
even_idx = [s for s in range(32) if bin(s).count("1") % 2 == 0]  # even # of |1> = even chirality
assert len(even_idx) == 16
Vp = np.zeros((32, 16), complex)
for col, s in enumerate(even_idx):
    Vp[s, col] = 1.0
print("(2) chirality Omega=Z^5; Gamma_a flip it; even-chirality subspace = 16  .. OK")

# charge conjugation C = product of the 5 'Y' gammas (Gamma_2,4,6,8,10)
C = np.eye(32, dtype=complex)
for k in range(N):
    C = C @ G[2 * k + 1]
assert np.allclose(C.T, -C), "expected C^T = -C"
eta = None
for a in range(10):
    if np.allclose(C @ G[a] @ np.linalg.inv(C), G[a].T):
        e = +1
    elif np.allclose(C @ G[a] @ np.linalg.inv(C), -G[a].T):
        e = -1
    else:
        raise AssertionError("C Gamma C^-1 != +- Gamma^T")
    eta = e if eta is None else eta
    assert e == eta, "inconsistent eta"
print(f"(3) C = prod(Y-gammas): C^T=-C, and C Gamma_a C^-1 = {eta:+d} Gamma_a^T  .... OK")


def block16(M):  # even-chirality 16x16 block
    return Vp.conj().T @ M @ Vp


# (4) 10-channel intertwiner M^a = (C Gamma_a)|_16 : symmetric?  -> bilinear vanishes
sym10 = []
for a in range(10):
    Ma = block16(C @ G[a])
    s = np.allclose(Ma, Ma.T)
    asym = np.allclose(Ma, -Ma.T)
    sym10.append((np.linalg.norm(Ma), s, asym))
norms10 = [n for n, _, _ in sym10]
allsym10 = all(s for _, s, _ in sym10)
print("\n(4) 10-channel  C*Gamma_a  restricted to the 16:")
print(f"      nonzero blocks: {sum(n>1e-9 for n in norms10)}/10;  all SYMMETRIC: {allsym10}")
print("      => psi^T (C Gamma_a) psi  with antisymmetric psi_i psi_j  ==> IDENTICALLY ZERO.")
assert allsym10 and min(norms10) > 1e-9

# (5) 120-channel C Gamma_a Gamma_b Gamma_c (a<b<c): antisymmetric (nonzero bilinear)
asym120 = []
for a, b, c in itertools.combinations(range(10), 3):
    Mabc = block16(C @ G[a] @ G[b] @ G[c])
    if np.linalg.norm(Mabc) > 1e-9:
        asym120.append(np.allclose(Mabc, -Mabc.T))
print("\n(5) 120-channel  C*Gamma_a Gamma_b Gamma_c  restricted to the 16:")
print(f"      nonzero blocks: {len(asym120)};  all ANTISYMMETRIC: {all(asym120)}")
assert len(asym120) > 0 and all(asym120)

print("""
CONCLUSION
==========
 * Confirmed (the user's deep point): 16 (x) 16 has NO singlet -> no 2-body operator can
   gap the SO(10) mirror symmetrically. The earlier CSS X-stabiliser test used a 2-body,
   gauge-non-covariant operator and therefore did NOT validly refute SMG. The dynamical
   SMG question, with the CORRECT quartic, is OPEN (not negative).
 * BUT the SPECIFIC operator proposed -- the 10-channel  U sum_a |psi^T C Gamma_a psi|^2 --
   is IDENTICALLY ZERO for a single 16: the 10 sits in the SYMMETRIC part of 16(x)16
   (C*Gamma_a is symmetric, verified), while the fermion bilinear psi_i psi_j is
   antisymmetric. So it cannot gap anything.
 * The nonzero bilinear-squared SO(10)-singlet channel is the ANTISYMMETRIC 120
   (C*Gamma_abc, verified antisymmetric). The correct 'build and compute' must use the
   120-channel quartic  U sum_{120} |psi^T C Gamma_abc psi|^2  (a real-rep sum -> singlet),
   or a genuine 4-fermion / mirror-coupled operator -- NOT the 10-channel.
ALL ASSERTS PASSED.
""")
