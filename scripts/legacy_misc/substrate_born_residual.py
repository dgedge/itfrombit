#!/usr/bin/env python3
r"""THE BORN RESIDUAL, CANDIDLY — does the microstate-count route close the measure, or relocate it?
And: does the substrate force COMPLEX Hilbert space? (Rung iii + the why-quantum question.)

Rungs i/ii derived the Born STRUCTURE: pointer basis = syndrome basis, non-contextuality from the
self-dual [8,4,4] -> Gleason -> P=|<v|psi>|^2. The residual is the MEASURE step: why does an
equal-amplitude branch carry unit weight? The proposal: ground it in an equal-a-priori-probability
count over substrate (syndrome) microstates. This script tests, honestly, whether that closes the
gap or merely moves it -- and is explicit that it does the latter.

exit 0 = the demonstrations run and the candid verdicts hold:
  [1] equal a-priori counting over microstates gives 1/2, NOT |c|^2 -- the naive branch/microstate
      measure is the WRONG answer;
  [2] recovering Born requires weighting microstates by |c|^2 itself -- Born-in, Born-out (circular);
  [3] the Born-typical set has vanishing fraction under equal counting (the circularity, quantified);
  [4] complex Hilbert space is the locally-tomographic one (real is not), so 'why complex' relocates
      to 'why local tomography', a substrate property -- not an independent derivation.
"""
import math

p = 0.8                 # = |c_0|^2 for a test qubit; deliberately != 1/2
N = 100                 # copies / records

# ================= [1] equal a-priori count over microstates -> 1/2, not |c|^2 =================
print("[1] EQUAL-A-PRIORI COUNT over outcome microstates (the naive 'each branch equally likely'):")
# each of the 2^N outcome sequences weighted equally; mean frequency of outcome 0:
mean_equal = sum(math.comb(N, k) * (k / N) for k in range(N + 1)) / 2 ** N
print(f"    mean frequency of outcome 0 (equal weight over 2^N sequences) = {mean_equal:.4f}")
print(f"    but |c_0|^2 = p = {p}.  Equal counting gives 1/2 regardless of the amplitude.")
assert abs(mean_equal - 0.5) < 1e-9              # equal microstate count -> 1/2, NOT Born
print("    -> the equal-a-priori measure over branches/records is the WRONG measure: it ignores |c|^2.")

# ================= [2] recovering Born needs Born-weighting -> circular =================
print("\n[2] RECOVERING |c|^2 REQUIRES WEIGHTING MICROSTATES BY |c|^2 (Born-in, Born-out):")
# weight each sequence by Prod_i |c_{s_i}|^2 = p^{#0}(1-p)^{#1}; mean frequency of 0:
mean_born = sum(math.comb(N, k) * p ** k * (1 - p) ** (N - k) * (k / N) for k in range(N + 1))
print(f"    mean frequency of outcome 0 (Born weight p^#0 (1-p)^#1) = {mean_born:.4f} = p")
assert abs(mean_born - p) < 1e-9
print("    -> Born is recovered ONLY by inserting the |c|^2 weight on microstates: that IS the Born")
print("       measure. So 'count microstates' closes nothing; it relocates the assumption to the")
print("       a-priori MEASURE on the microstates -- which must already be the Born one.")

# ================= [3] the circularity, quantified: the Born-typical set is negligible under equal count =
print("\n[3] THE CIRCULARITY QUANTIFIED -- equal counting puts ~0 weight on the Born-typical set:")
H = -p * math.log2(p) - (1 - p) * math.log2(1 - p)         # binary entropy
# fraction of sequences with frequency within 0.02 of p, under EQUAL counting:
lo, hi = int((p - 0.02) * N), int((p + 0.02) * N)
frac_equal = sum(math.comb(N, k) for k in range(lo, hi + 1)) / 2 ** N
print(f"    H(p) = {H:.3f} bits; typical-set size ~ 2^(N H) = 2^{N*H:.0f}, but as a FRACTION of all")
print(f"    2^N sequences it is ~ 2^(N(H-1)) -> 0.  Measured fraction (freq within 0.02 of p) = {frac_equal:.2e}")
assert frac_equal < 1e-6                          # the Born-typical frequencies are ~never seen under equal count
print("    -> 'typicality' only concentrates on freq=p in the Born-WEIGHTED measure; under the")
print("       equal-a-priori measure it has vanishing weight. The measure assumption is irreducible.")

# ================= [4] WHY COMPLEX HILBERT SPACE? -> relocates to local tomography =================
print("\n[4] WHY COMPLEX HILBERT SPACE -- it is the LOCALLY-TOMOGRAPHIC one (real is not):")
# state-space dimension K (real parameters for an unnormalised state): complex d^2 ; real d(d+1)/2.
def K_complex(d): return d * d
def K_real(d):    return d * (d + 1) // 2
dA = dB = 2
for name, K in (("complex", K_complex), ("real", K_real)):
    local, glob = K(dA) * K(dB), K(dA * dB)
    lt = (local == glob)
    print(f"    {name:>8s} d=2:  K_A*K_B = {local:>3d}  vs  K_AB = {glob:>3d}  -> locally tomographic: {lt}")
    if name == "complex":
        assert local == glob                      # complex: K_AB = K_A K_B  (locally tomographic)
    else:
        assert local != glob                      # real: fails local tomography
print("    -> GIVEN that the substrate is locally tomographic (records are locally describable),")
print("       complex Hilbert space is singled out over real/quaternionic. This RELOCATES 'why")
print("       complex' to 'why local tomography' -- a substrate premise, not an independent derivation.")

print(f"""
[verdict] CANDID -- rung (iii) RELOCATES the measure; it does not close it; and 'why complex' is
presupposed.
  * The microstate-count route does NOT derive Born. Equal a-priori counting over branches/records
    gives 1/2 [1]; Born is recovered only by weighting microstates by |c|^2 [2], which is the Born
    measure itself; under equal counting the Born-typical frequencies have vanishing weight [3]. So
    the route relocates the irreducible assumption from 'Born probabilities' to 'the a-priori
    measure over syndrome microstates is the Born one' -- a MORE contested premise than the one the
    Gleason route already used.
  * Therefore the framework's TIGHTEST Born argument remains rungs i+ii: non-contextuality (DERIVED
    from self-duality) + Gleason, whose only residual is the near-definitional frame premise that a
    complete measurement yields outcomes with non-contextual probabilities summing to 1. Rung (iii)
    is a step sideways (into the measure debate), not forward.
  * 'Why complex Hilbert space' is not derived: the substrate is a qubit (complex) QEC system by
    construction, so the formalism is presupposed. Complex is singled out over real by LOCAL
    TOMOGRAPHY [4], which relocates the question to a substrate property; deriving the quantum
    formalism from a genuinely pre-quantum (classical/combinatorial) substrate is the open
    why-quantum problem, outside the current framework.
  NET: the Born STRUCTURE is derived (rungs i/ii); the Born MEASURE and the complex formalism are
  the honest, irreducible residuals -- the same ones every reconstruction faces. Better to state
  this plainly than to dress relocation up as closure.
exit 0""")
print("ALL ASSERTIONS PASSED -- equal count=1/2; Born needs Born-weight (circular); typical set negligible; complex<-local tomography.")
