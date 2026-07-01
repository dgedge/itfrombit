#!/usr/bin/env python3
r"""K04 SHADOW-DEPTH NORMALIZATION: honest status of the surviving gravitating fraction.

The K04 frozen wall network overcloses (bare, w=-2/3, a^-1) by ~1e31. Only the
wall-shadow depth-truncation can rescue it: a wall interior has no registers, so
only the recorded boundary strain in the surrounding registered crystal
gravitates. The question: what is the SURVIVING gravitating fraction?

Verified from the island-floor surface (k04_island_floor_surface_results.jsonl,
L=18,20, 16 reps each), three policies (rho_DM/rho_B; observed = 5.36):
  - strict same-level pairing : ~2700  -> overcloses ~500x  -> EXCLUDED
  - boundary-local rescue     : O(few) -> right ballpark, but NOISY/non-converging
                                (median 1.5 @L18, 9.9 @L20)
  - island floor (policy-indep): ~0.002, ZERO-INFLATED (30-50% of runs = 0 exactly)
The boundary-local surface is dominated by RARE depth-one orphan islands
(heavy-tailed); deeper/fully-adaptive readout is excluded by locality (it would
read through the register-free wall interior). Self-asserting on the readings.
"""
def ok(c,m): print(("  PASS " if c else "  FAIL ")+m); assert c,m
OBS = 0.1200/0.02237   # rho_DM/rho_B observed ~ 5.36
# measured medians (k04_island_floor_surface_results.jsonl)
strict = {18:2700.91, 20:2649.76}
rescue = {18:1.525,  20:9.88}
floor  = {18:0.001,  20:0.002}     # island_floor_wall_frac median
zfrac  = {18:0.50,   20:0.31}      # fraction of runs with floor == 0

print("="*72); print("K04 SHADOW-DEPTH STATUS"); print("="*72)
print(f"  observed rho_DM/rho_B = {OBS:.2f}")
print(f"  strict pairing  : {strict}  -> /obs = {strict[18]/OBS:.0f}x , {strict[20]/OBS:.0f}x")
print(f"  boundary-local  : {rescue}  -> /obs = {rescue[18]/OBS:.2f} , {rescue[20]/OBS:.2f}")
print(f"  island floor    : {floor}  (zero-frac {zfrac})")

ok(min(strict.values())/OBS > 100, "strict pairing overcloses by >100x -> excluded (shadow MUST act)")
ok(0.1 < rescue[18]/OBS < 3 and 0.1 < rescue[20]/OBS < 3, "boundary-local rescue lands within O(1) of observed -- shadow rescues in DIRECTION + ROUGH MAGNITUDE")
ok(abs(rescue[20]-rescue[18]) > 0.5*min(rescue.values()), "but rescue is NOT converging across L (1.5 vs 9.9) -- noisy, not a stable number")
ok(max(floor.values()) < 0.01 and min(zfrac.values()) >= 0.3, "island floor is tiny + zero-inflated (30-50% exactly zero)")
ok(not (floor[20] < floor[18]), "floor does NOT decrease over L=18->20 -> cosmological (L->inf) limit NOT established by current data")

print("\n"+"="*72); print("VERDICT")
print("  QUALITATIVE: the wall-shadow truncation DOES rescue K04 from the bare")
print("  ~1e31 catastrophe -- strict 2700x collapses to O(few) under the only")
print("  locality-allowed (boundary-local) readout, the right ballpark for DM.")
print("  QUANTITATIVE: the surviving fraction is NOT a clean derived number. It is")
print("  heavy-tailed (rare depth-one orphan islands), zero-inflated, and does not")
print("  converge across L; the policy-independent island floor is tiny but its")
print("  L->inf (percolated) limit is not established. Deeper readout is excluded")
print("  by locality (no registers inside the wall).")
print("  CONCLUSION: shadow-depth normalization remains OPEN as a precise Omega_K04.")
print("  K04's gravitating fraction is BOUNDED (O(few)x baryons at most, often ~0)")
print("  but not predictable -> K04 is correctly NOT the clean dark-matter component;")
print("  the DM budget stays on the nu_R+R4 zero-mode. The bare 1e31 'crisis' is an")
print("  artifact of ignoring the shadow; the shadowed value is small-but-noisy. exit 0")
