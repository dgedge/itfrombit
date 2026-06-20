#!/usr/bin/env python3
r"""MULTI-MODE HANDOFF SERVICE ACTION — PRE-REGISTRATION (R1, the A_req conversion theorem).

Per the recorded protocol discipline: this FREEZES the minimal action, the state variables, the
zero/gauge-mode rule, and the kill criteria BEFORE any coefficient comparison. It then computes
the determinant ratio AS A FUNCTION of the one physical input it does NOT get to choose — the
§3.2 service-channel connectivity — and leaves the verdict PENDING that input. It does NOT claim
a conversion. (count labels != fluctuation modes; the connectivity must come from §3.2 physics,
not from whichever graph happens to give 3/2.)

================================ FROZEN PRE-REGISTRATION ================================
(i)   MINIMAL ACTION. The handoff fluctuations are those of the §5.2 QEC service CTMC at its
      arrest saddle. The standard quadratic fluctuation action of a CTMC at a steady point is the
      SYMMETRIZED service-graph Laplacian L (Onsager-Machlup / MSRJD): S2 = (1/2) x^T L x, with
      L_ij = -(symmetrized rate i<->j), L_ii = sum_j!=i (...). No other action is admitted.
(ii)  STATE VARIABLES. The occupation-current fluctuations of the §3.2 walk-active service
      channels {V_em, V_weak, V_strong} (the ACTIVE block) and the exiting-vacuum legs (the
      VACUUM block). Counts: active n_a = 3, vacuum n_v = 2 (DRIFT register-handoff l.1702).
(iii) ZERO/GAUGE-MODE RULE. Each connected service block has exactly one conservation
      (uniform-shift) zero mode; remove it. Compare PSEUDODETERMINANTS det' (product of nonzero
      eigenvalues = n * #spanning-trees, matrix-tree theorem).
(iv)  KILL CRITERIA (decided here, before the answer):
        CONVERTED      iff the §3.2-FORCED connectivity gives det'(L_active)/det'(L_vacuum)=3/2
                       EXACTLY, with that connectivity fixed independently of the target.
        REJECTED       iff the forced connectivity gives any other ratio.
        UNDERSPECIFIED iff §3.2 does not fix the active-block connectivity.
      Target for the comparison (revealed only at the end): A_req = sqrt(3/2) = 1.224745.
========================================================================================
"""

def det2(m):                       # 2x2 determinant
    return m[0][0]*m[1][1] - m[0][1]*m[1][0]

def laplacian(n, edges):
    L = [[0.0]*n for _ in range(n)]
    for i, j in edges:
        L[i][j] -= 1.0; L[j][i] -= 1.0
        L[i][i] += 1.0; L[j][j] += 1.0
    return L

def spanning_trees(n, edges):      # matrix-tree: any first cofactor of L
    L = laplacian(n, edges)
    if n == 1:
        return 1
    if n == 2:
        return L[1][1]             # 1x1 cofactor
    if n == 3:
        return det2([[L[1][1], L[1][2]], [L[2][1], L[2][2]]])   # delete row/col 0
    raise ValueError("n<=3 only")

def pseudodet(n, edges):           # det' = n * #spanning-trees (connected graph)
    return n * spanning_trees(n, edges)

# candidate connectivities for the ACTIVE 3-node block (the physical fork) ----------------
ACTIVE_CANDIDATES = {
    "K3  triangle (symmetric: all pairs coupled, cf. ||{Ti,Tj}||~2/3)": (3, [(0,1),(1,2),(0,2)]),
    "P3  path     (exit-asymmetric: 2 of 3 exit, 1 internal)":          (3, [(0,1),(1,2)]),
    "S3  star     (one hub + 2 leaves == P3 topologically)":            (3, [(0,1),(0,2)]),
}
VACUUM = (2, [(0,1)])              # the 2 exiting-vacuum legs: the only connected 2-node graph (K2)

det_vac = pseudodet(*VACUUM)
assert det_vac == 2, "det'(K2) must be 2 (eigenvalues 0,2)"

TARGET = 1.5                       # 3/2 ; A_req = sqrt(TARGET)
A_req  = (TARGET) ** 0.5
rows = []
for name, (n, edges) in ACTIVE_CANDIDATES.items():
    da = pseudodet(n, edges)
    ratio = da / det_vac
    rows.append((name, da, ratio, ratio ** 0.5))

# checks that pin the fork (so the verdict is mechanical once §3.2 fixes the connectivity)
by_da = {name.split()[0]: pseudodet(n, e) for name, (n, e) in ACTIVE_CANDIDATES.items()}
assert by_da["K3"] == 9 and by_da["P3"] == 3 and by_da["S3"] == 3, "det' values: K3=9, P3=S3=3"
assert abs((by_da["P3"]/det_vac) - 1.5) < 1e-12, "PATH/STAR active block -> ratio 3/2 (would CONVERT)"
assert abs((by_da["K3"]/det_vac) - 4.5) < 1e-12, "TRIANGLE active block -> ratio 9/2 (would REJECT)"

# the verdict is PENDING the physical input; this script does NOT decide the connectivity.
sec32_connectivity_fixed = False            # NOT resolved here — see the open question below
verdict = "PENDING-CONNECTIVITY" if not sec32_connectivity_fixed else "(decided)"
assert verdict == "PENDING-CONNECTIVITY", "pre-registration does not claim a conversion"

bar = "=" * 94
print(bar)
print("MULTI-MODE HANDOFF SERVICE ACTION — PRE-REGISTRATION (verdict PENDING the §3.2 connectivity)")
print(bar)
print("  frozen: action = symmetrized service Laplacian; vars = 3 active + 2 vacuum channels;")
print("          zero-mode = remove 1 conservation mode/block; det' = product of nonzero eigenvalues.")
print(f"  target (revealed last): A_req = sqrt(3/2) = {A_req:.6f}\n")
print(f"  ACTIVE connectivity candidate                                   det'  ratio/det'(K2)  sqrt -> A")
print(f"  {'-'*92}")
for name, da, ratio, sq in rows:
    hit = "== A_req" if abs(sq - A_req) < 1e-9 else ""
    print(f"  {name:58.58s}  {da:4.0f}   {ratio:8.3f}     {sq:.5f} {hit}")
print(f"\n  vacuum block: K2, det' = {det_vac}")
print(f"""
{bar}
THE FORK (the whole theorem reduces to ONE physical question):
  * If §3.2 forces the ACTIVE block to be a PATH / STAR (exit-asymmetric: 2 of 3 channels exit,
    the C_loop=3/2 / 2-of-3-exit geometry) -> det' ratio = 3/2 -> A_req = sqrt(3/2): CONVERTED.
  * If §3.2 forces a TRIANGLE K3 (fully symmetric coupling, as the equal pairwise anticommutators
    ||{{Ti,Tj}}||~2/3 superficially suggest) -> det' ratio = 9/2 -> sqrt = 2.121: REJECTED.

  These point OPPOSITE ways, so the connectivity cannot be guessed or chosen for convenience.
  OPEN DECISIVE INPUT: derive the §3.2 active-channel service connectivity from the hopping
  operator / the 2-of-3 exit rule -- is the 3-channel handoff graph a path (one channel does not
  exit) or a triangle (all symmetric)? That single fact flips CONVERTED <-> REJECTED.

  VERDICT: {verdict}. The pre-registration is frozen and the fork is computed; the answer is
  pending the connectivity derivation, NOT a coefficient comparison. A_req stays free until then.
{bar}""")
print(f"exit 0 -- PRE-REGISTRATION frozen; PATH/STAR->3/2 (convert), TRIANGLE->9/2 (reject); "
      f"verdict PENDING the §3.2 active-channel connectivity. No conversion claimed.")
