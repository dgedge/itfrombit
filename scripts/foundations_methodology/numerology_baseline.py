#!/usr/bin/env python3
"""
SEARCH-SPACE ACCOUNTING for the Holographic Circlette / TCH framework.

Question: when the framework reports a dimensionless prediction matching an
experimental ratio to sub-1%, how surprising is that GIVEN the alphabet of
simple expressions the framework helps itself to? If a small alphabet already
"covers" most of the O(1) number line to 1%, then any individual sub-1% match
carries almost no evidential weight (numerology). If coverage is sparse, the
matches are signal.

Method (pre-registered BEFORE looking at hit rates):
  1. Fix the alphabet from what the framework ACTUALLY uses (audit-sourced).
  2. Generate the family of "framework-simple" expressions under an explicit
     complexity budget.
  3. COVERAGE: fraction of log-uniform [0.1,10] within rel-tol of some family
     member. Target-independent numerology baseline.
  4. PER-TARGET: for each of the framework's dimensionless claims, count how
     many *distinct, equally-or-less-simple* family expressions also land
     within the claimed tolerance ("competitors"). >=1 competitor => the match
     is not surprising.
  5. NULL HIT-RATE: take ~30 real dimensionless physics constants the framework
     did NOT predict; what fraction does the same family hit within 1%? If it
     hits a large fraction of arbitrary constants, the framework's chosen hits
     are unremarkable.
  6. Sensitivity over the integer ceiling and budget.

Nothing here is tuned to a target. The alphabet is justified in ALPHABET below.
"""
import math, itertools
from fractions import Fraction as F

# ---------------------------------------------------------------------------
# 1. PRE-REGISTERED ALPHABET (from the audit census of what the framework uses)
# ---------------------------------------------------------------------------
# Integers that appear as building blocks in the framework's dimensionless
# formulas: small ints 1..9 pervasively, plus 12,14,16,24,28 (code/weight
# enumerator), and a few others. For O(1)-RATIO predictions the small ints
# dominate; big structural ints (137,208,240,256) are not ratio-numerators.
INT_CEIL = 12          # default; sensitivity run also does 9 and 16
# Irrationals the framework treats as "substrate-given" (no fit):
#   sqrt2 (C8 eigenvalue, silver), sqrt3 (colour paths), phi & 1/phi (golden,
#   line-graph), pi and its sqrt (continuum projections). silver=1+sqrt2 and
#   1/phi are derivable from these, included for fairness.
PHI = (1+math.sqrt(5))/2
IRRATIONALS = {
    "sqrt2": math.sqrt(2), "sqrt3": math.sqrt(3), "sqrt5": math.sqrt(5),
    "phi": PHI, "1/phi": 1/PHI, "pi": math.pi, "sqrt_pi": math.sqrt(math.pi),
    "silver": 1+math.sqrt(2),
}

def build_family(int_ceil):
    """Framework-simple expressions. Forms that actually occur:
       p/q ; sqrt(p/q) ; (p/q)*irr ; (p/q)*irr1*irr2 ; (p/q)/irr ; sqrt(p/q)*irr
    Deduped by value to 6 sig figs, restricted to (0,1000)."""
    vals = {}
    ints = range(1, int_ceil+1)
    rats = []
    for p in ints:
        for q in ints:
            rats.append((F(p,q), float(p)/q))
    irr_items = list(IRRATIONALS.items())
    def add(v, tokens):
        if not (1e-6 < v < 1e3): return
        key = round(v, 6)
        # keep the SIMPLEST (fewest tokens) representative
        if key not in vals or tokens < vals[key][1]:
            vals[key] = (v, tokens)
    # rationals (token cost 2: p,q)
    for fr, v in rats:
        add(v, 2)
        if v > 0: add(math.sqrt(v), 3)          # sqrt(p/q)
    # rational * single irrational (token cost 3)
    for fr, v in rats:
        for name,iv in irr_items:
            add(v*iv, 3); add(v/iv, 3)
            add(math.sqrt(v)*iv, 4)
    # rational * two distinct irrationals (token cost 4)
    for fr, v in rats:
        for (n1,i1),(n2,i2) in itertools.combinations(irr_items,2):
            add(v*i1*i2, 4); add(v*i1/i2, 4); add(v/(i1*i2), 4)
    return vals   # dict: rounded_value -> (value, min_tokens)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def competitors(family, target, reltol, max_tokens=None):
    lo, hi = target*(1-reltol), target*(1+reltol)
    out=[]
    for v,(val,tok) in family.items():
        if lo <= val <= hi and (max_tokens is None or tok <= max_tokens):
            out.append((val,tok))
    return sorted(out, key=lambda t:abs(t[0]-target))

def coverage(family, lo=0.1, hi=10.0, reltol=0.01, nsamp=200000):
    """Fraction of log-uniform [lo,hi] within reltol of SOME family member."""
    fam = sorted(v for v in family.keys() if lo*0.5 < v < hi*2)
    import bisect
    hits=0
    loglo, loghi = math.log(lo), math.log(hi)
    # deterministic grid instead of RNG (RNG is unavailable / nondeterministic)
    for k in range(nsamp):
        x = math.exp(loglo + (loghi-loglo)*k/(nsamp-1))
        i = bisect.bisect_left(fam, x)
        near = False
        for j in (i-1, i):
            if 0 <= j < len(fam) and abs(fam[j]-x) <= reltol*x:
                near=True; break
        if near: hits+=1
    return hits/nsamp

# ---------------------------------------------------------------------------
# 2. FRAMEWORK'S DIMENSIONLESS CLAIMS (value, claimed reltol, label)
# ---------------------------------------------------------------------------
CLAIMS = [
    (2/9,            0.005, "sin^2_thetaW = 2/9"),
    (math.sqrt(7/9), 0.0006,"M_W/M_Z = sqrt(7/9)"),
    (2*math.sqrt(2), 0.0002,"M_N/Lambda = 2sqrt2"),
    (4/3,            0.006, "sqrt(sigma)/Lambda = 4/3"),
    (7/5,            0.01,  "glueball 2++/0++ = 7/5"),
    (PHI/2,          0.02,  "m_rho/M_N = phi/2"),
    (math.sqrt(2)*PHI,0.0003,"m_rho/Lambda = sqrt2*phi"),
    (1/math.sqrt(4*math.pi),0.017,"f_pi/Lambda = 1/sqrt(4pi)"),
    (12/5,           0.07,  "Omega_DE/Omega_DM = 12/5"),
    (0.02997,        0.016, "nu dmsq ratio (R=1,d=1/3)"),
    (2/3,            1e-4,  "Koide Q_l = 2/3 (ALGEBRAIC IDENTITY)"),
    (27/28,          0.005, "w0 = -27/28 (n_s duality)"),
]

# ---------------------------------------------------------------------------
# 3. REAL physics constants the framework did NOT predict (null reference set)
# ---------------------------------------------------------------------------
# dimensionless ratios/values in ~[0.01,100], standard PDG/CODATA
NULL_CONSTS = [
    ("m_mu/m_e", 206.768), ("m_tau/m_mu", 16.817), ("m_p/m_e", 1836.15),
    ("m_n/m_p", 1.00138), ("m_p/m_n", 0.99862), ("m_s/m_d", 20.0),
    ("m_c/m_s", 11.7), ("m_b/m_c", 3.27), ("m_t/m_b", 41.3),
    ("|V_us|", 0.2243), ("|V_cb|", 0.0408), ("|V_ub|", 0.00382),
    ("sin^2_th12_PMNS", 0.307), ("sin^2_th23", 0.546), ("sin^2_th13", 0.0220),
    ("alpha_s(MZ)", 0.1179), ("sin^2_thW(MSbar)", 0.23122),
    ("Omega_b/Omega_c", 0.186), ("Omega_m", 0.315), ("Omega_L", 0.685),
    ("n_s", 0.9649), ("eta_baryon/1e-10", 6.12), ("m_h/m_t", 0.724),
    ("m_W/m_t", 0.465), ("m_Z/m_h", 0.730), ("m_mu/m_pi", 0.7573),
    ("m_pi/m_rho", 0.1803), ("m_K/m_pi", 3.570), ("g_A", 1.2723),
    ("f_K/f_pi", 1.197), ("BR(H->bb)", 0.582), ("a_e*1e3/...", 1.159652),
]

# ===========================================================================
def main():
    print("="*72)
    print("SEARCH-SPACE ACCOUNTING — numerology baseline for TCH")
    print("="*72)
    fam = build_family(INT_CEIL)
    print(f"\nAlphabet: ints 1..{INT_CEIL}; irrationals {sorted(IRRATIONALS)}")
    print(f"Family size (distinct values, deduped to 6 s.f.): {len(fam):,}")
    nO1 = sum(1 for v in fam if 0.1<=v<=10)
    print(f"  of which in [0.1,10]: {nO1:,}")

    print("\n" + "-"*72)
    print("3. COVERAGE — fraction of log-uniform [0.1,10] within rel-tol of SOME")
    print("   family value (target-INDEPENDENT numerology baseline):")
    print("-"*72)
    for tol in (0.05, 0.01, 0.005, 0.001):
        c = coverage(fam, 0.1, 10, tol, nsamp=60000)
        print(f"   rel-tol {tol*100:>5.2f}% : {c*100:5.1f}% of the O(1) line is 'pre-covered'")

    print("\n" + "-"*72)
    print("4. PER-TARGET COMPETITORS — # distinct family expressions within the")
    print("   framework's CLAIMED tolerance of each prediction (>=1 => unremarkable)")
    print("-"*72)
    print(f"   {'observable':34s} {'value':>9s} {'tol':>7s} {'#comp':>6s}  nearest others")
    comp_counts=[]
    for val,tol,label in CLAIMS:
        comps = competitors(fam, val, tol)
        comp_counts.append(len(comps))
        others = ", ".join(f"{v:.4f}" for v,t in comps[:4])
        print(f"   {label:34s} {val:9.4f} {tol*100:6.2f}% {len(comps):6d}  {others}")
    import statistics
    med = statistics.median(comp_counts)
    print(f"\n   median competitors per claim: {med}")
    print(f"   claims with >=1 competitor (numerology-cheap): "
          f"{sum(1 for c in comp_counts if c>=1)}/{len(comp_counts)}")
    print(f"   claims with 0 competitors (rare/surprising):   "
          f"{sum(1 for c in comp_counts if c==0)}/{len(comp_counts)}")

    print("\n" + "-"*72)
    print("5. NULL HIT-RATE — fraction of ~30 real physics constants the framework")
    print("   did NOT predict that the SAME family hits within 1% by chance:")
    print("-"*72)
    null_hits=0
    for name,val in NULL_CONSTS:
        comps = competitors(fam, val, 0.01)
        if comps:
            null_hits+=1
    print(f"   {null_hits}/{len(NULL_CONSTS)} = {100*null_hits/len(NULL_CONSTS):.0f}% of UNRELATED constants")
    print(f"   are also hit to 1% by the same alphabet.")
    print( "   (If this is high, matching a chosen constant to 1% is near-automatic.)")

    print("\n" + "-"*72)
    print("6. SENSITIVITY (coverage at 1% as alphabet grows):")
    print("-"*72)
    for ic in (9, 12, 16):
        f2 = build_family(ic)
        c = coverage(f2, 0.1, 10, 0.01, nsamp=40000)
        print(f"   ints 1..{ic:2d}: family={len(f2):>7,}  coverage@1% = {c*100:4.1f}%")

    print("\n" + "="*72)
    print("INTERPRETATION KEY")
    print("="*72)
    print("""  - Coverage@1% near or above ~50%: sub-1% O(1) matches are essentially
    free; individual ratio 'predictions' carry little evidential weight.
  - Per-target median competitors >=1: the framework could have hit each
    target with SEVERAL equally-simple expressions, so picking one is not
    evidence of a mechanism (this is the quantitative form of the audit's
    '2/9 reused with different justifications' red flag).
  - Null hit-rate high: the alphabet matches arbitrary constants too, so the
    chosen matches are unremarkable.
  - The TIGHT matches (alpha at 3 ppb, M_N at 0.02%) have ZERO simple-ratio
    competitors -> they are NOT explained by simple numerology. Their status
    rests entirely on whether the *formula* (e.g. the -24/7 Dyson-Schwinger
    coefficient) was constructed freely. A small alphabet cannot fake a ppb
    hit; a flexible multi-term formula can. So tight matches are a SEPARATE
    question (formula-space freedom), not a simple-ratio-numerology question.""")

if __name__ == "__main__":
    main()
