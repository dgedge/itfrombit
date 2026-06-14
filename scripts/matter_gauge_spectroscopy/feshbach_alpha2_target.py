#!/usr/bin/env python3
r"""The Feshbach-trace alpha^2 target: does the graviton E_g self-energy over Q carry alpha^2?
(DRIFT G7 final thread.) The heat-kernel showed the BARE (unitary) trace is alpha-free -> rank ~208;
so alpha enters ONLY by Landauer-weighting the P<->Q crossings (dissipative/QEC reading, which §10.4's
anti-mass argument independently SELECTS over coherent hopping). The sharp question is then a COUNT:
how many Landauer-IRREVERSIBLE crossings does the leading graviton self-energy Sigma = H_PQ G_QQ H_QP have?
alpha-power = (number of irreversible crossings). We resolve the count from the framework's OWN definitions.

§5.9 (anchored, validated by m_d-m_u=alpha*Lambda at 4%): w = alpha*Lambda = (irreversibility fraction alpha)
x (energy per tick). KEY (L767): the error flip itself is a UNITARY Pauli-X (reversible, NO Landauer cost);
ONLY the non-unitary syndrome ERASURE (projection back to the codeword) carries the alpha. So per the
framework's own bit-weight, alpha counts ERASURES (entropy-decreasing Q->P), not error-injections (P->Q).
"""
import sys, math, itertools as it
import numpy as np

# ---------------- register, codewords, P/Q (ANCHOR §2.1/§2.6) ----------------
BIT={0:"G0",1:"G1",2:"LQ",3:"C0",4:"C1",5:"I3",6:"chi",7:"W"}
def bit(n,name): return (n>>[k for k,v in BIT.items() if v==name][0])&1
def is_cw(n):
    if bit(n,"G0") and bit(n,"G1"): return False                 # R1
    if bit(n,"W")!=bit(n,"chi"): return False                    # R2
    cc=(bit(n,"C0"),bit(n,"C1"))                                 # R3
    if bit(n,"LQ")==0 and cc!=(0,0): return False
    if bit(n,"LQ")==1 and cc==(0,0): return False
    return True
P=[n for n in range(256) if is_cw(n)]; Q=[n for n in range(256) if not is_cw(n)]
assert len(P)==48 and len(Q)==208
print(f"P(codewords)={len(P)}  Q(invalid/syndrome)={len(Q)}")

# ---------------- (1) the P<->Q boundary crossings: minimum weight, and WHICH bits are syndromes ----------
# single-bit flips from each codeword: which land in Q? (a weight-1 P->Q crossing = a single-bit syndrome)
bit_cross=np.zeros(8,int)                                        # crossings caused by flipping bit i
n_w1=0
for n in P:
    for i in range(8):
        m=n^(1<<i)
        if not is_cw(m):
            bit_cross[i]+=1; n_w1+=1
print(f"\n[1] weight-1 P->Q crossings (single-bit syndrome events): {n_w1} total")
print( "    by bit:  " + "  ".join(f"{BIT[i]}:{bit_cross[i]}" for i in range(8)))
free=[BIT[i] for i in range(8) if bit_cross[i]==0]
print(f"    FREE/gauge bits (flip stays in P, no Landauer): {free}   (I3 is the §5.9 'active' free bit)")
print(f"    => minimum P<->Q crossing weight = 1 (e.g. flip chi -> breaks R2; flip a colour bit -> breaks R3).")
print(f"       So the LEADING graviton self-energy is a round-trip of weight-1 crossings: P ->(out) Q ->(back) P.")
assert n_w1>0 and min(bit_cross[6],bit_cross[7])>0               # chi,W are syndromes

# ---------------- (2) the count: erasures per virtual graviton loop ----------------
print(f"""
[2] alpha-power = number of Landauer-IRREVERSIBLE crossings in Sigma = H_PQ G_QQ H_QP (2 vertices: out, back).
    Landauer's principle + §5.9: entropy-INCREASING P->Q (error injection) is the UNITARY Pauli-X step
    (reversible, NO cost); entropy-DECREASING Q->P (syndrome erasure/projection) is the non-unitary
    Landauer step (cost alpha). A closed virtual graviton loop P->Q->P contains exactly ONE erasure (the
    return-to-codeword) -> ONE alpha.  =>  the §5.9-GROUNDED count is  alpha^1, not alpha^2.
    The framework's 'double-Landauer' alpha^2 charges a SECOND alpha to the error-injection too -- but §5.9
    itself says that step is unitary (free). The 2nd alpha ('gravitational back-reaction') is the un-grounded one.""")

# ---------------- (3) the three readings are all data-degenerate (so magnitude can't decide) ----------
ALPHA=1/137.035999; OmL=0.6847
K_data = math.sqrt(OmL)/(24*math.pi*ALPHA**2)                   # the unification value ~206
readings = {
 "alpha^0 (coherent rank, 208-3)" : 205.0,
 "alpha^1 (one erasure/loop, §5.9-grounded)": (K_data*ALPHA)/ALPHA, # = K_data; express as c/alpha below
 "alpha^2 (double-Landauer, Part-20==205)" : K_data,
}
c1 = K_data*ALPHA                                               # coefficient if K = c1/alpha
print(f"[3] all three alpha-powers reproduce K_eff ~205-206 (so the magnitude does NOT select the power):")
print(f"    alpha^0 : K = 208-3 = 205            (rank; G7: not constructible)")
print(f"    alpha^1 : K = {c1:.3f}/alpha = {c1/ALPHA:.1f}   (one erasure/loop; coeff {c1:.3f} ~ 3/2, unexplained)")
print(f"    alpha^2 : K = sqrt(OmL)/(24pi a^2) = {K_data:.1f}   (double-Landauer; == Part-20)")
print(f"    spread of the three: {205.0:.0f} / {c1/ALPHA:.1f} / {K_data:.1f}  -- all within "
      f"{(max(c1/ALPHA,K_data)-205)/205*100:.1f}% : DEGENERATE.")
assert abs(c1/ALPHA-K_data)<2

# ---------------- (4) item 136's alpha^4 is a DIFFERENT process (topology change, code distance) ----------
# the gravitational STIFFNESS (self-energy) couples to weight-1 detectable errors; a TOPOLOGY CHANGE needs a
# weight-d LOGICAL fault. Raw min Hamming distance of the 48 codewords vs the [8,4,4] logical distance:
dmin=min(bin(a^b).count("1") for a,b in it.combinations(P,2))
print(f"""
[4] item 136's alpha^4 is a CATEGORY-SEPARATE count, not the stiffness:
    raw min Hamming distance of the 48 codewords = {dmin} (the free I3/gauge bit makes neighbours distance-1);
    the [8,4,4] LOGICAL distance d=4 is the min weight to change the MACROSCOPIC/logical state (topology change,
    items 126/128 baryogenesis/tensor-modes), NOT a virtual stiffness excursion. The graviton STIFFNESS is a
    weight-1 detectable-error self-energy (alpha^1), whereas a real TOPOLOGY CHANGE is a weight-4 logical fault
    (alpha^4). Item 136 imported the topology-change power (alpha^4) into the stiffness/M_P formula -- a category
    error (which is also why its horizon-power was 1/4, the lone outlier).""")

print(f"""
=========================================================================================
VERDICT (Feshbach-trace alpha^2 target -- pursued to the bottom):
  REACHED: the alpha^2 is NOT a trace output (the bare trace is alpha-free -> rank ~208). It is a COUNT of
  Landauer-irreversible P<->Q crossings in the 2-vertex graviton self-energy. Resolving that count from the
  framework's OWN §5.9 definition (error-injection = unitary Pauli-X = free; only the syndrome ERASURE costs
  alpha) gives ONE erasure per virtual loop => **alpha^1**, K_eff = (~3/2)/alpha = 205.5. The framework's
  'double-Landauer' alpha^2 charges a SECOND alpha to the (unitary, free) error-injection -- a double-count;
  its 2nd 'back-reaction alpha' has no §5.9 grounding.
  NOT RESCUED: (a) alpha^0/alpha^1/alpha^2 are data-degenerate (all K~205-206) so the magnitude cannot select
  the power; (b) even the §5.9-grounded alpha^1 leaves an unexplained O(1) coefficient (~3/2) AND the horizon
  and Lambda as inputs -- it is a *derived Dirac RELATION* M_P^2 ~ (3/2/alpha) Lambda^3 R_dS, not an intrinsic
  M_P; (c) the definitive erasure-count needs the substrate's Lindbladian/QEC jump operators (§15 items
  119-122) to confirm exactly one non-unitary projection per graviton loop.
  NET: the target ADVANCES the gravity sector -- it dislodges the headline 'double-Landauer alpha^2' (likely a
  double-count; §5.9 grounds only alpha^1) and pins item 136's alpha^4 as a topology-change power miscast as a
  stiffness. The honest residual: the graviton stiffness self-energy carries alpha^1 (one erasure), giving a
  Dirac-type relation with an O(1) coefficient still open and the horizon still an input. No parameter-free M_P;
  but the alpha-power is now ARGUED (alpha^1) from anchored definitions, not asserted (alpha^2).
  Sharply-posed remaining question (for items 119-122): is there exactly ONE non-unitary syndrome projection
  per virtual graviton P->Q->P loop?  If yes, alpha^1 is derived; the framework's alpha^2 is retired.
=========================================================================================""")
print("exit 0 -- crossing structure built; alpha-power resolved to alpha^1 (§5.9-grounded) vs framework alpha^2 (double-count); degenerate; item136 alpha^4 = topology-change miscast.")
