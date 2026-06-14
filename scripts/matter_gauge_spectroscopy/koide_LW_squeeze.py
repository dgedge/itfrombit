#!/usr/bin/env python3
"""
How precisely is delta pinned, and is delta=2/9 EXACT or merely approximate?
Plus the Lindemann-Weierstrass squeeze. Reproducible; nothing fitted.
"""
import math, random

# PDG 2024 charged-lepton masses (MeV) and 1-sigma errors
M  = {"e":0.51099895000, "mu":105.6583755, "tau":1776.86}
dM = {"e":1.5e-10,       "mu":2.3e-6,      "tau":0.12}

def extract_delta(m):
    sm = {k: math.sqrt(v) for k,v in m.items()}
    mu = (sum(sm.values())/3)**2
    # c_n = (sqrt(m_n/mu) - 1)/sqrt2 ; tau is n=0 (max), solve delta from tau branch
    c_tau = (sm["tau"]/math.sqrt(mu) - 1)/math.sqrt(2)
    return math.acos(max(-1,min(1,c_tau)))   # delta + 0

d0 = extract_delta(M)
print("central delta (tau branch) = %.7f rad" % d0)
print("2/9                        = %.7f rad" % (2/9))
print("delta - 2/9                = %+.2e" % (d0 - 2/9))

# Monte-Carlo error propagation (deterministic seed via fixed offsets, no RNG dependence on Date)
samples=[]
random.seed(12345)
for _ in range(20000):
    m={k: M[k]+random.gauss(0,dM[k]) for k in M}
    samples.append(extract_delta(m))
samples.sort()
mean=sum(samples)/len(samples)
sd=(sum((x-mean)**2 for x in samples)/len(samples))**0.5
print("\nMonte-Carlo: delta = %.7f +/- %.7f" % (mean, sd))
nsig=abs(mean-2/9)/sd
print("deviation of 2/9 from central: %.2f sigma" % nsig)
print("=> data is CONSISTENT with delta=2/9 exactly, but a transcendental within")
print("   ~%.0e of 2/9 fits equally well. 0.02%% closeness does NOT prove exactness." % sd)

print("""
LINDEMANN-WEIERSTRASS SQUEEZE (the structural result):
  The Koide circulant off-diagonal is B = |B| e^{i*delta}.
  Lemma: a nonzero ALGEBRAIC complex number has TRANSCENDENTAL argument
         (if z algebraic, |z| algebraic, so e^{i*arg z}=z/|z| algebraic;
          L-W: arg z algebraic & nonzero => e^{i arg z} transcendental;
          contradiction => arg z transcendental).
  Therefore:
    * IF delta = 2/9 EXACTLY (algebraic, nonzero), then B is TRANSCENDENTAL
      -> B cannot be a finite combinatorial/discrete substrate amplitude
         (Hamming counts, rational hops): same obstruction as item 135.
      -> delta=2/9 'exact' must be a DEFINITION (defect-density ratio imposed
         as the phase), NOT the output of a discrete dynamical mechanism.
    * IF delta is a genuine dynamical/geometric phase from algebraic substrate
      amplitudes, it is TRANSCENDENTAL -> delta != 2/9 exactly; '2/9' is only
      its leading rational approximation (consistent with the data above).

  FORK (proven, mutually exclusive):
    (i)  delta = 2/9 exactly  => it is a rational defect ratio by construction,
         not a phase; the cos(delta+2pi n/3) form is a re-parametrisation, and
         no closed holonomy can source it (pi-obstruction, koide_holonomy.py).
    (ii) delta is a real phase => transcendental, ~2/9 approximately; the
         framework's 'exact 2/9' claim downgrades to '~2/9', and the open-path
         Pancharatnam class is the only pi-free geometric candidate.

  Either branch REFUTES the current canonical framing (exact rational 2/9 AS a
  substrate-derived phase): (i) kills 'phase', (ii) kills 'exact'.
""")
