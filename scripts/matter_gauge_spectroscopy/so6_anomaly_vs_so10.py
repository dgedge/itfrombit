#!/usr/bin/env python3
"""
Why the SO(6) "SMG" test was necessarily SSB, and why SO(6) is NOT a valid SMG analogue of
the framework's SO(10) 16.

SMG (a symmetric gapped vacuum) is possible ONLY for anomaly-free content. Two checks on the
SO(6) ~= SU(4) chiral spinor (the 4 of SU(4)), then the contrast with SO(10)/SM:

  (1) CUBIC ANOMALY of the 4 of SU(4): d^{abc} = Tr(T^a {T^b,T^c}). For SU(N>=3) the
      fundamental is COMPLEX and anomalous (d != 0). The 4+4bar (vector-like) cancels.
  (2) SU(4) SINGLETS in Lambda^N(4): a symmetric gapped vacuum needs an SU(4)-singlet at a
      NONTRIVIAL filling. Counted via the SU(4) Casimir H_C = sum_a G_a^2.

Result (printed): the 4 is anomalous AND has SU(4)-singlets only at N=0,4 (trivial) => NO
symmetric gapped vacuum => any gap is SSB (matches so6_smg_ssb_check.py). So SO(6) cannot show
SMG. CONTRAST: the framework gauges the SM (a SUBGROUP of SO(10)); the 16 is SM-anomaly-free
and HAS SM-singlets at nontrivial filling (1 at N=1, 8 at N=8; sm_singlet_structure.py) -> SMG
not obstructed. The SO(6) analogy is doubly wrong: full (anomalous) group vs the SM subgroup.

numpy only.
"""
import numpy as np
import itertools

# ---- SU(4) fundamental generators (15 Hermitian traceless 4x4) ----
def su4():
    n = 4; T = []
    for i in range(n):
        for j in range(i + 1, n):
            S = np.zeros((n, n), complex); S[i, j] = S[j, i] = 1; T.append(S / 2)
            A = np.zeros((n, n), complex); A[i, j] = -1j; A[j, i] = 1j; T.append(A / 2)
    for k in range(1, n):
        d = np.zeros(n); d[:k] = 1; d[k] = -k
        T.append((np.diag(d) / np.sqrt(2 * k * (k + 1))).astype(complex))
    return T

T = su4(); assert len(T) == 15
# normalisation: Tr(T^a T^b) = 1/2 delta_ab
G = np.array([[np.trace(T[a] @ T[b]).real for b in range(15)] for a in range(15)])
assert np.allclose(G, 0.5 * np.eye(15)), "generator normalisation off"
print("SU(4): 15 generators, Tr(T^a T^b)=1/2 delta_ab  (OK)")

# ---- (1) cubic anomaly ----
def danom(gens):
    return np.array([[[np.trace(gens[a] @ (gens[b] @ gens[c] + gens[c] @ gens[b])).real
                       for c in range(15)] for b in range(15)] for a in range(15)])
d4 = danom(T)
d4bar = danom([-t.conj() for t in T])
print(f"\n(1) cubic anomaly d^abc = Tr(T^a{{T^b,T^c}}):")
print(f"    4 of SU(4):     max|d| = {np.abs(d4).max():.4f}   -> {'ANOMALOUS' if np.abs(d4).max()>1e-9 else 'anomaly-free'}")
print(f"    4 + 4bar:       max|d| = {np.abs(d4 + d4bar).max():.4e}   -> vector-like, cancels")
assert np.abs(d4).max() > 1e-9 and np.abs(d4 + d4bar).max() < 1e-9

# ---- (2) SU(4) singlets in Lambda^N(4) via the Casimir ----
M = 4
NZ = [[(i, j, T[a][i, j]) for i in range(M) for j in range(M) if abs(T[a][i, j]) > 1e-12]
      for a in range(15)]
def pcb(m, k): return bin(m & ((1 << k) - 1)).count("1")
def singlets(NF):
    basis = [sum(1 << x for x in c) for c in itertools.combinations(range(M), NF)]
    idx = {mk: i for i, mk in enumerate(basis)}; D = len(basis)
    HC = np.zeros((D, D), complex)
    for nz in NZ:
        Ga = np.zeros((D, D), complex)
        for col, m in enumerate(basis):
            for i, j, v in nz:
                if not (m >> j) & 1: continue
                s = (-1) ** pcb(m, j); m2 = m ^ (1 << j)
                if i != j and (m2 >> i) & 1: continue
                s *= (-1) ** pcb(m2, i); Ga[idx[m2 | (1 << i)], col] += s * v
        HC += Ga.conj().T @ Ga
    ev = np.linalg.eigvalsh(HC).real
    return D, int(np.sum(ev < 1e-9))

print(f"\n(2) SU(4) singlets in Lambda^N(4)  (gauge group of SO(6) is the FULL SU(4)):")
print(f"    {'N':>3} {'dim':>4} {'#SU(4) singlets':>16}")
for NF in range(0, 5):
    D, ns = singlets(NF)
    print(f"    {NF:>3} {D:>4} {ns:>16}")

print("""
CONCLUSION
==========
 * The SO(6)~=SU(4) chiral spinor (the 4) is ANOMALOUS (d^abc != 0) and has SU(4)-singlets
   ONLY at N=0,4 (trivial). So it has NO symmetric gapped vacuum at nontrivial filling ->
   any gap MUST break the symmetry (SSB). This is the deeper reason so6_smg_ssb_check.py
   found a condensate gap: SMG is *impossible* for SO(6), so SSB is the only option.
 * Therefore SO(6) is NOT a valid SMG analogue of the framework. The framework gauges the
   SM (a SUBGROUP of SO(10)), not SO(10)/SU(4): the 16 is SM-anomaly-free (anomaly_cobordism_
   class.py) and HAS SM-singlets at nontrivial filling (1 at N=1 = nu^c, 8 at N=8;
   sm_singlet_structure.py) -> a symmetric gapped vacuum EXISTS and SMG is not obstructed.
 * The SO(6) test was doubly inappropriate: it gauged the FULL (anomalous) group with no
   nontrivial singlet, whereas the SMG-relevant case is an anomaly-free SUBGROUP with singlets.
   A genuine SMG mechanism test needs ANOMALY-FREE chiral content (the SO(10) 16 with SM
   gauge group: the 2-cell 2^32 frontier; or a tractable anomaly-free chiral toy).
ALL ASSERTS PASSED.
""")
