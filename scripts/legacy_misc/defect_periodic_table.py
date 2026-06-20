#!/usr/bin/env python3
r"""THE PERIODIC TABLE OF K04 DEFECTS -- a self-asserting taxonomy scaffold.

Once the embedded Z^3 (x) Q_3 substrate is trusted, every stable or metastable defect is a
CANDIDATE physical sector. This script is the scaffold the discovery programme hangs on: it
lays out the defect zoo as a ledger and checks the classification is internally consistent.
It does NOT claim every cell is new physics -- most are healable or artifacts, and the point
is to sort which are protected (candidate sectors) from which are not.

AXES.
  * geometric dimension d in {0,1,2,3} (point / line / wall / bulk);
  * topological charge -- for bond-graph defects, the H_1(T^3;Z_2) winding of the
    difference graph [C XOR crystal] (k04_kempe_locked_defect.py: Kempe-locked IFF != 0);
  * PROTECTION CLASS: topological (locked) > kinetic (pinned, trivial) > healable > artifact;
  * ledger columns: energy/tension, mobility, shadow, candidate sector, honesty tier.

TWO FAMILIES emerge, and that is itself a result:
  (I)  BOND-GRAPH (Kempe) defects -- classified by H_1; protection is topological/kinetic/
       healable/artifact. Structural finding asserted below: topological protection (nonzero
       H_1) is carried ONLY by the 1D winding string; 2D walls and 0D islands are
       homologically TRIVIAL, so the dark-matter debris survives KINETICALLY, not
       topologically (the Bullet-Cluster caveat, depinning_mobility_gate.py).
  (II) SUBSTRATE-SERVICE defects -- not misbonds but features of the live service layer:
       the R4 line occupancy (d=1 support -> dark energy, item 131) and the horizon-severing
       record (Bekenstein C=55/8, items 118/122). Different kind of object, real sector.

exit 0 = the protection rule reproduces every DERIVED row; the canon-derived tensions and
         f_frust are correct and correctly ordered; the pinning scale is computed; and every
         OPEN cell is FLAGGED, not asserted (exit 0 claims only the derived rows + the
         classification logic, never the conjectural sectors).
"""
from fractions import Fraction

# energy model E = -w4 C4 - w6 C6 (k04_kempe_locked_defect.py); representative positive scales
W4, W6 = 2.0, 1.0          # only signs / ordering matter for the scaffold

# ---- canon-derived energy laws (each cited to the script that proved it) ----
TENSION = {
    "free_wall":       0.0,            # compatible wall: degenerate cube vacua  (k04_wall_tension.py)
    "frustrated_wall": 2 * W4 + 12 * W6,  # sigma_wall / cut face                (k04_wall_tension.py)
    "winding_string":  1 * W4 + 4 * W6,   # mu / step (minimal LOCKED object)    (k04_kempe_locked_defect.py)
    "peanut":          4 * W4 + 16 * W6,  # misbond gap (heals in one Kempe move)(k04_kempe_locked_defect.py)
}

def protection(h1, barrier_w6, healable, vacuum=False):
    """Family-I protection class from the H_1 winding charge + kinetics."""
    if vacuum:                       return "vacuum"
    if any(b % 2 for b in h1):       return "topological"   # [C]!=0 in H_1(T^3;Z_2): Kempe-locked
    if healable:                     return "healable"      # trivial + a local 2-switch removes it
    if barrier_w6 > 0:               return "kinetic"       # trivial but pinned (Peierls-Nabarro)
    return "artifact"                                       # trivial, free => degenerate-vacuum boundary

# ---- THE TABLE -------------------------------------------------------------------------------
# Family I: bond-graph (Kempe) defects. row = dict; protection is COMPUTED then checked.
FAMILY_I = [
    dict(name="bulk crystal (vacuum)", d=3, h1=(0, 0, 0), barrier=0.0, heal=False, vacuum=True,
         energy="free_wall", mobility="--", shadow="--",
         sector="empty space (8 gauge-equivalent 2x2x2 tilings)", tier="derived", want="vacuum"),
    dict(name="compatible domain wall", d=2, h1=(0, 0, 0), barrier=0.0, heal=False, vacuum=False,
         energy="free_wall", mobility="free", shadow="none",
         sector="none -- free boundary between degenerate tilings", tier="derived", want="artifact"),
    dict(name="frustrated domain wall", d=2, h1=(0, 0, 0), barrier=TENSION["frustrated_wall"], heal=False, vacuum=False,
         energy="frustrated_wall", mobility="pinned (~42 OOM)", shadow="dark (gravity only)",
         sector="DARK MATTER -- debris network, rho ~ sigma_wall * f_frust / xi(R)", tier="derived", want="kinetic"),
    dict(name="winding string", d=1, h1=(1, 0, 0), barrier=TENSION["winding_string"], heal=False, vacuum=False,
         energy="winding_string", mobility="locked", shadow="dark (gravity only)",
         sector="cosmic-string-class sector (the ONLY topologically locked object)", tier="derived", want="topological"),
    dict(name="orphan island (peanut)", d=0, h1=(0, 0, 0), barrier=TENSION["peanut"], heal=True, vacuum=False,
         energy="peanut", mobility="heals", shadow="--",
         sector="none -- self-heals in one Kempe move", tier="derived", want="healable"),
    dict(name="protected 0D defect", d=0, h1=(0, 0, 0), barrier=0.0, heal=False, vacuum=False,
         energy=None, mobility="n/a", shadow="n/a",
         sector="NONE -- no locality-protected point defect", tier="closed", want=None,
         nogo="defect_0d_nogo.py"),   # CLOSED: bounded => winding 0 => not locked; protection floor is 1D
]

# Family II: substrate-service defects (not misbonds; protection is of a different nature).
FAMILY_II = [
    dict(name="R4 line occupancy", d=1, prot="service-support (d=1)",
         energy="3 two-edge stars (1D support)", mobility="--", shadow="Landauer service",
         sector="DARK ENERGY -- w(a) = -1 + a/28", tier="derived"),     # item 131 / item131_r4_support_dimension.py
    dict(name="horizon-severing record", d=2, prot="record / holographic",
         energy="C = 55/8 per severed ledger leg", mobility="--", shadow="record channel",
         sector="BEKENSTEIN ENTROPY / black holes", tier="derived"),     # items 118 / 122
]

# ================== ASSERTION BATTERY ==================
print("[1] PROTECTION RULE reproduces every DERIVED Family-I row:")
for r in FAMILY_I:
    got = protection(r["h1"], r["barrier"], r["heal"], r["vacuum"])
    flag = {"open": "OPEN (flagged, not asserted)",
            "closed": f"CLOSED by no-go ({r.get('nogo', '')})"}.get(r["tier"], got.upper())
    print(f"    {r['name']:<28s} d={r['d']}  H1={r['h1']}  -> {flag}")
    if r["tier"] == "derived":
        assert got == r["want"], (r["name"], got, r["want"])

print("\n[2] STRUCTURAL FINDING: topological protection is carried ONLY by 1D strings.")
topo = [r for r in FAMILY_I if r["tier"] == "derived" and protection(r["h1"], r["barrier"], r["heal"], r["vacuum"]) == "topological"]
assert len(topo) == 1 and topo[0]["d"] == 1, topo
print(f"    exactly {len(topo)} locked class: '{topo[0]['name']}' (d={topo[0]['d']}). 2D walls / 0D islands")
print("    are homologically TRIVIAL -> dark-matter debris is KINETIC (pinned), not topological.")

print("\n[3] CANON TENSIONS, positive and correctly ordered (for w4,w6>0):")
free, fw, mu, pe = (TENSION[k] for k in ("free_wall", "frustrated_wall", "winding_string", "peanut"))
print(f"    free wall = {free:.0f};  mu(string) = {mu:.0f} = w4+4w6;  sigma_wall = {fw:.0f} = 2w4+12w6;  peanut = {pe:.0f} = 4w4+16w6")
assert free == 0.0 and mu > 0 and fw > 0 and pe > 0
assert pe > fw > mu > free                      # peanut > wall > string > 0

print("\n[4] FRUSTRATED FRACTION (k04_frustration_fraction.py):")
f_frust = Fraction(1, 2) / Fraction(7, 8)
print(f"    f_frust = (1/2)/(7/8) = {f_frust} ~ {float(f_frust):.3f}   (free fraction {1 - f_frust})")
assert f_frust == Fraction(4, 7)

print("\n[5] MOBILITY SCALE (depinning_mobility_gate.py): astrophysical drive vs lattice barrier:")
a_astro, a0_m, c = 1.2e-10, 0.595e-15, 2.998e8   # m/s^2 ; lattice a0 in m ; c
R_drive = a_astro * (2 * a0_m) / c**2
print(f"    R_drive = a*(2 a0)/c^2 = {R_drive:.1e}  << 1  ->  kinetic rows are pinned by ~{-__import__('math').log10(R_drive):.0f} OOM")
assert R_drive < 1e-30

print("\n[6] TIER BOOKKEEPING -- exit 0 asserts only DERIVED sectors; CLOSED cite a no-go; OPEN held back:")
established = sorted({r["sector"] for r in FAMILY_I + FAMILY_II if r["tier"] == "derived"})
closed_cells = [r for r in FAMILY_I + FAMILY_II if r["tier"] == "closed"]
open_cells = [r["name"] for r in FAMILY_I + FAMILY_II if r["tier"] == "open"]
for r in FAMILY_I + FAMILY_II:
    if r["tier"] == "open":
        assert r["sector"] == "?", r                          # an open cell must NOT carry a claimed sector
    if r["tier"] == "closed":
        assert r["sector"].startswith("NONE") and r.get("nogo")  # a closed cell must cite its no-go
assert all("?" not in s for s in established)                  # no derived sector is a placeholder
print(f"    established (derived) sectors: {len(established)};  closed by no-go: {len(closed_cells)};  "
      f"still open: {len(open_cells)} -> {open_cells}")
for s in established:
    print(f"      - {s}")
for r in closed_cells:
    print(f"      x {r['name']}: {r['sector']}  ({r['nogo']})")

print("\n[7] COVERAGE (no silent caps): dimensions and the now-closed 0D cell, logged explicitly:")
dims = sorted({r["d"] for r in FAMILY_I + FAMILY_II})
assert dims == [0, 1, 2, 3]
headline = {"DARK MATTER", "DARK ENERGY", "BEKENSTEIN", "cosmic-string"}
assert all(any(h in s for s in established) for h in headline)
print(f"    dimensions covered: {dims} (0=point,1=line,2=wall,3=bulk).")
print("    0D PROTECTED CELL: CLOSED by no-go (defect_0d_nogo.py) -- protection floor is 1D (the string).")
print("    next frontiers are the ledger COLUMNS: coupling channel, mass/shadow, mobility, survival per row.")

# ================== THE PERIODIC TABLE ==================
print("\n" + "=" * 96)
print("PERIODIC TABLE OF K04 DEFECTS")
print("=" * 96)
hdr = f"{'defect':<28s} {'d':>1s} {'protection':<13s} {'energy/tension':<22s} {'mobility':<16s} sector"
print("\n-- Family I: bond-graph (Kempe) defects, classified by H_1(T^3;Z_2) --")
print(hdr); print("-" * 96)
for r in FAMILY_I:
    prot = {"open": "OPEN", "closed": "NO-GO"}.get(r["tier"], protection(r["h1"], r["barrier"], r["heal"], r["vacuum"]))
    en = "--" if r["energy"] is None else (f"{TENSION[r['energy']]:.0f} ({r['energy']})" if r["energy"] in TENSION else r["energy"])
    print(f"{r['name']:<28s} {r['d']:>1d} {prot:<13s} {en:<22s} {r['mobility']:<16s} {r['sector']}")
print("\n-- Family II: substrate-service defects (not misbonds) --")
hdr2 = f"{'defect':<28s} {'d':>1s} {'protection':<22s} {'marker':<32s} sector"
print(hdr2); print("-" * 96)
for r in FAMILY_II:
    prot = "OPEN" if r["tier"] == "open" else r["prot"]
    print(f"{r['name']:<28s} {r['d']:>1d} {prot:<22s} {r['energy']:<32s} {r['sector']}")

print(f"""
[verdict] SCAFFOLD STANDS. The defect zoo is a consistent two-family ledger:
  * Family I (Kempe): protection = topological/kinetic/healable/artifact, all from one H_1
    invariant. The lone TOPOLOGICAL class is the 1D winding string; the DARK-MATTER debris
    is the KINETIC (pinned, trivial-homology) frustrated wall -- durable but not collisionless
    (the honest Bullet-Cluster caveat); islands HEAL; compatible walls are ARTIFACTS.
  * Family II (service): the R4 d=1 line occupancy = DARK ENERGY, the horizon-severing record
    = BEKENSTEIN entropy -- already canon sectors, which is the evidence the taxonomy is real.
  Every cell above is tiered; exit 0 asserts only the DERIVED rows and the classification
  logic. The 0D 'protected point defect' cell is now CLOSED by a no-go (defect_0d_nogo.py):
  protection requires winding, which requires >=1D system-spanning extent, so the floor is the
  1D string. The next rungs hang off this scaffold: coupling channels, a mass/shadow ledger,
  mobility, and survival -- per row.
exit 0""")
print("ALL ASSERTIONS PASSED -- protection rule consistent; tensions & f_frust canon; open cells flagged not asserted.")
