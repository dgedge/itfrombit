#!/usr/bin/env python3
"""
Item 93/143 — the SMG (symmetric-mass-generation) route to continuum closure.

walk_kernel_overlap.py closed the EXTERNAL route with an honest negative: the body-
diagonal walk's own Wilson term gives 4 species, and a single-species term is barred
in finite-range form by Poincare-Hopf chi(T^3)=0. That route asked the lattice-QCD
question "lift 15 doublers to 1 chiral fermion ultralocally" -- impossible.

This script tests a DIFFERENT escape (the user's "internal walk within the 8-face
octahedron"): symmetric mass generation. A chiral generation is obtained not by lifting
doublers with a Wilson term, but by GAPPING THE MIRROR sector with an interaction that
preserves the chiral gauge symmetry -- legitimate iff the matter content is anomaly-free
INCLUDING the global Z_16 (Wang-Wen 2018; Garcia-Etxebarria-Montero 2018). The 3+1D
magic number is 16 Weyl/gen INCLUDING nu_R (= the SO(10) spinor). N-N is evaded by
interactions, not by the overlap, so the ULTRALOCAL ideal can in principle be restored.

What is RIGOROUSLY decidable here (and tested):
  (A) the framework's [8,4,4]/Q3 code (R1-R4, §2.1-2.2) produces exactly 16 Weyl/gen
      incl nu_R, with the four canonical charges (§2.8).
  (B) the two §2.11 sum rules  sum Q = 0, sum Q^2 = 16  (grounds the code is the real one).
  (C) SMG PRECONDITION: that 16 = the SO(10) spinor; its perturbative anomalies cancel;
      and the 16-vs-15 distinction is the GLOBAL Z_16 -- which the framework's OWN
      §2.11 result  sum(B-L) = +1/gen ("nu_R structural necessity")  is, restated.
  (D) the internal Q3 walk (the 8-face octahedron) is a FINITE space: its walk operator
      is GAPPED (no zero modes), bipartite-chiral ({g5,A}=0). N-N/Poincare-Hopf is a
      theorem about the CONTINUOUS external BZ torus; it does NOT bind a finite internal
      space. So the specific obstruction that capped the external route is ABSENT here.

What is NOT decidable here (flagged, NOT faked):
  (E) whether the code projection P (256->48) actually gaps the lattice MIRROR doublers
      while sparing the physical generation and preserving the gauge symmetry. That needs
      the external-BZ-corner -> codeword map (items 77/97/102, OPEN; NB item 98 is the
      INTERNAL C_8/Q_3 dual, a different question) -- see mirror_invalid_subspace.py. This script
      establishes the PRECONDITION + the absence of the internal no-go, not the dynamical
      closure. It does not, and does not claim to, move item 93/143 to continuum-Locked.

numpy + exact rational (fractions); self-asserting.
"""
import numpy as np, itertools
from fractions import Fraction as Q

def Hd(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

# ----------------------------------------------------------------------------------
Hd("(A) Build the framework's actual [8,4,4]/Q3 code (ANCHOR §2.1-2.2, §2.8)")
# 8-bit register, octant-indexed:  c = (G0,G1,LQ,C0,C1,I3,chi,W)
G0, G1, LQ, C0, C1, I3, CHI, W = range(8)

def R1(c): return not (c[G0] == 1 and c[G1] == 1)          # 3 generations (no 4th)
def R2(c): return c[W] == c[CHI]                           # chirality locks to weak
def R3(c):                                                 # colour separates q/lepton
    if c[LQ] == 0: return (c[C0], c[C1]) == (0, 0)
    return (c[C0], c[C1]) != (0, 0)
def R4(c): return not (c[LQ] == 0 and c[I3] == 0 and c[CHI] == 1)   # nu_R exclusion

allwords = list(itertools.product([0, 1], repeat=8))
valid48  = [c for c in allwords if R1(c) and R2(c) and R3(c)]            # 48
active45 = [c for c in valid48 if R4(c)]                                 # 45
print(f"  full internal Hilbert space  : 2^8 = {len(allwords)}  (Q3 = (C^2)^8)")
print(f"  R1&R2&R3 valid codewords     : {len(valid48)}  (= 3 generations x 16)")
print(f"  + R4 (nu_R excluded) active   : {len(active45)}  (45 SM-fermion codewords)")
assert len(allwords) == 256 and len(valid48) == 48 and len(active45) == 45

# charge via §2.8 :  Q = 1/2 Z_f + 1/3 sum_i Z_ci + 1/2 Z_p   (matter: Z_p=+1)
def Zf(c):  return 1 if c[I3] == 0 else -1                 # up-type(nu,u)=+1, down(e,d)=-1
def sumZc(c):
    if (c[C0], c[C1]) == (0, 0): return -3                 # lepton singlet
    return -1                                              # quark: one +1, two -1
def charge(c): return Q(1, 2) * Zf(c) + Q(1, 3) * sumZc(c) + Q(1, 2) * 1

def species(c):
    lep = (c[LQ] == 0)
    up  = (c[I3] == 0)
    return ("nu" if up else "e") if lep else ("u" if up else "d")

# one generation = fix (G0,G1)=(0,0)
gen0 = [c for c in valid48 if (c[G0], c[G1]) == (0, 0)]
print(f"\n  one generation (G0G1=00): {len(gen0)} Weyl states = "
      f"4 colour-LQ x 2 I3 x 2 chi")
canon = {}
for c in gen0:
    canon.setdefault(species(c), charge(c))
print("  canonical charges from §2.8:")
for s in ("nu", "e", "u", "d"):
    print(f"     {s:>3}: Q = {str(canon[s]):>5}   (SM {'0' if s=='nu' else {'e':'-1','u':'+2/3','d':'-1/3'}[s]})")
assert canon == {"nu": Q(0), "e": Q(-1), "u": Q(2, 3), "d": Q(-1, 3)}
assert len(gen0) == 16
# the single excluded state per generation is exactly nu_R (I3=0,chi=1,lepton)
nuR = [c for c in gen0 if not R4(c)]
assert len(nuR) == 1 and species(nuR[0]) == "nu" and nuR[0][CHI] == 1
print(f"  the 1 R4-excluded state/gen   : nu_R  {''.join(map(str,nuR[0]))}  (I3=0, chi=1) ✓")

# ----------------------------------------------------------------------------------
Hd("(B) §2.11 two-identity sum rules over the 45 active states (grounds the code)")
sumQ  = sum(charge(c) for c in active45)
sumQ2 = sum(charge(c) ** 2 for c in active45)
print(f"  sum_45 Q   = {sumQ}    (U(1)_Y gauge + mixed grav-U(1)_Y anomaly)  [expect 0]")
print(f"  sum_45 Q^2 = {sumQ2}   (QED beta_1 coefficient)                    [expect 16]")
assert sumQ == Q(0) and sumQ2 == Q(16)
print("  -> reproduces ANCHOR §2.11 exactly: the constructed code IS the framework's.")

# ----------------------------------------------------------------------------------
Hd("(C) SMG PRECONDITION: the 16 is the SO(10) spinor; anomaly-free incl global Z_16")
# the 16 Weyl content per generation as LEFT-handed Weyl (Y_eff = -Y for RH conjugates):
#   field   mult  Y          (Y = Q - T3 ; standard, = what the code yields per §2.11)
so10 = [("Q",  6, Q(1, 6)), ("u^c", 3, Q(-2, 3)), ("d^c", 3, Q(1, 3)),
        ("L",  2, Q(-1, 2)), ("e^c", 1, Q(1)),     ("nu^c", 1, Q(0))]
nW = sum(m for _, m, _ in so10)
sY  = sum(m * y    for _, m, y in so10)
sY3 = sum(m * y**3 for _, m, y in so10)
print(f"  16-Weyl content = SO(10) spinor: {nW} Weyl  (= 4 lepton + 12 quark, matches A)")
print(f"  perturbative:  sum Y   = {sY}   [grav x U(1)_Y]        [expect 0]")
print(f"                 sum Y^3 = {sY3}   [U(1)_Y^3]             [expect 0]")
assert nW == 16 and sY == Q(0) and sY3 == Q(0)
# the SU(2)^2 U(1) and SU(3) parts (vector-like / even doublets) -- cited from §2.11 Part 7
print("  (A1 [SU3]^2 U1, A2 [SU2]^2 U1, A5 [SU3]^3, A6 Witten: all 0 -- §2.11 Part 7, exact)")

print("\n  THE 16-vs-15 DISTINCTION IS NOT PERTURBATIVE:")
# nu^c is a total singlet (Y=0): dropping it leaves sum Y, sum Y^3 unchanged.
so10_no_nuR = [t for t in so10 if t[0] != "nu^c"]
sY_15  = sum(m * y    for _, m, y in so10_no_nuR)
sY3_15 = sum(m * y**3 for _, m, y in so10_no_nuR)
print(f"     drop nu_R (Y=0 singlet):  sum Y = {sY_15}, sum Y^3 = {sY3_15}  -- STILL 0.")
print("     so perturbative anomalies cancel with OR without nu_R (15 or 16).")
assert sY_15 == Q(0) and sY3_15 == Q(0)
print("     The 16th state is forced by the GLOBAL Z_16 anomaly (Wang-Wen 2018):")
print("       - 16 Weyl incl nu_R : Z_16 trivial  -> mirror is SYMMETRICALLY GAPPABLE (SMG-able)")
print("       - 15 Weyl, no nu_R  : Z_16 = -1/gen -> NOT SMG-able on its own.")
print("     The framework's §2.11  sum(B-L) = +1/gen  ('nu_R required for gauge closure')")
print("     is the (B-L)-anomaly face of the SAME requirement. The Z_16 is the discrete")
print("     spin-bordism remnant of the (B-L)-related symmetry -- deeply linked, not a")
print("     literal identity -- but both single out exactly the 16th state (nu_R).")
print("     => the framework's matter content is FORCED (by the code) to be the SMG-able 16.")
print("        Precondition for symmetric mirror-gapping: MET, and met uniquely/structurally.")

# ----------------------------------------------------------------------------------
Hd("(D) The internal Q3 walk (8-face octahedron) is FINITE: N-N does not bind it")
# vertices = 3-bit octant addresses; hypercube adjacency (Hamming distance 1)
verts = list(itertools.product([0, 1], repeat=3))
def ham(a, b): return sum(x != y for x, y in zip(a, b))
A = np.array([[1.0 if ham(a, b) == 1 else 0.0 for b in verts] for a in verts])
g5 = np.diag([1.0 if (sum(v) % 2 == 0) else -1.0 for v in verts])   # item 116: g5=diag(+even,-odd)
ev = np.sort(np.linalg.eigvalsh(A))
gap = float(np.min(np.abs(ev)))
anti = float(np.linalg.norm(g5 @ A + A @ g5))
even = [v for v in verts if sum(v) % 2 == 0]; odd = [v for v in verts if sum(v) % 2 == 1]
print(f"  Q3 = 3-cube hypercube graph: 8 vertices, bipartite {len(even)} even + {len(odd)} odd")
print(f"  walk (adjacency) spectrum   : {[int(round(x)) for x in ev]}  (3,1,-1,-3 ; mult 1,3,3,1)")
print(f"  internal chirality g5=diag(+even,-odd): {{g5, A}} = {anti:.1e}  (exactly chiral)")
print(f"  |smallest eigenvalue|        : {gap:.3f}   -> GAPPED, NO zero modes")
assert anti < 1e-12 and gap > 0.5

# contrast: the EXTERNAL body-diagonal walk has zeros on the continuous BZ (walk_kernel_overlap)
def Wbody(k): return 8 * (1 - np.cos(k[0]) * np.cos(k[1]) * np.cos(k[2]))
ext_zeros = [tuple(int(round(x / np.pi)) for x in p)
             for p in itertools.product([0.0, np.pi], repeat=3) if abs(Wbody(np.array(p))) < 1e-9]
print(f"\n  EXTERNAL body-diagonal walk W_body on continuous BZ T^3: zeros at {len(ext_zeros)} corners")
print(f"     {ext_zeros}   -> the 4 surviving species (walk_kernel_overlap.py negative)")
print("  INTERNAL Q3 walk: a FINITE 8-site space, no continuous BZ, no vector field,")
print("     so Poincare-Hopf chi(T^3)=0 has NO jurisdiction -- and the operator is GAPPED.")
print("  => the SPECIFIC no-go that capped the external route does NOT extend to the")
print("     internal walk. A mirror-gapping term sourced INTERNALLY is not forbidden.")

# ----------------------------------------------------------------------------------
Hd("(E) VERDICT — what this does and does not establish")
print("""ESTABLISHED (rigorous, grounded in §2.1-2.2/2.8/2.11):
  * The [8,4,4]/Q3 code produces exactly 16 Weyl/gen incl nu_R = the SO(10) spinor,
    with the four canonical charges and the §2.11 identities (sum Q=0, sum Q^2=16).  (A,B)
  * That 16 is the 3+1D SMG-able / Z_16-trivial set; the 16-vs-15 distinction is the
    GLOBAL anomaly -- the (B-L)-anomaly face of which is the framework's own §2.11
    sum(B-L)=+1/gen nu_R-necessity (deeply linked invariants, same required 16th state).
    So the SMG PRECONDITION (anomaly-free incl global => mirror symmetrically gappable)
    is MET, and met structurally -- the code FORCES the SMG-able content.            (C)
  * The internal Q3 walk is a finite, gapped, bipartite-chiral operator; N-N/Poincare-
    Hopf (a theorem about the continuous external BZ) does not bind it. The specific
    obstruction that gave walk_kernel_overlap.py its 4 species is ABSENT internally.  (D)

NOT ESTABLISHED (open; NOT faked here):
  * That the code projection P (256->48) actually gaps the lattice MIRROR doublers while
    sparing the physical generation and preserving the chiral gauge symmetry. This needs
    the external-BZ-corner -> codeword map (items 77/97/102, OPEN; item 98 is the INTERNAL
    C_8/Q_3 dual, a different question). FOLLOW-UP: mirror_invalid_subspace.py runs this
    under the minimal standard candidate map and finds the mirrors ARE the W!=chi states,
    all landing in the 208-dim invalid subspace, targeted by R2 ultralocally -- BUT a hard
    projection is N-N-forbidden, so legitimacy needs R2 as a symmetric INTERACTION (the
    dynamical SMG step), still unconstructed. So: mechanism identified, dynamics open.

STATUS MOVEMENT (honest):
  Before: the EXTERNAL route is a closed negative -- continuum closure caps at the
          'overlap (exp-local), external-Wilson' tier (standard lattice-QCD, not unique).
  After : a SECOND, uniquely-TCH route exists -- symmetric mass generation via the
          internal Q3 code -- whose PRECONDITION the framework satisfies structurally and
          whose N-N obstruction is ABSENT. Its dynamical closure is OPEN, but it
          is no longer barred by a no-go. The prize (ultralocal continuum chirality) is
          back on the table as a well-posed computation, not a proven result.
  This does NOT move item 93/143 to continuum-Locked. It reopens the path that
  walk_kernel_overlap.py appeared to close, and pins the remaining input: the dynamical
  SMG realisation of R2 (mirror_invalid_subspace.py identifies R2 as the operator).
""")
print("ALL ASSERTS PASSED (computed claims verified; verdict states the open piece honestly).")
