"""
circlette_sim.py  –  Tiny Universe Simulation
Holographic Circlette framework (D.G. Elliman, Part I, 2026)

Bit layout per Table 1 (index 0 = MSB, ring order):
  0: G0   Generation bit 0
  1: G1   Generation bit 1
  2: LQ   Lepton(0) / Quark(1)  ← CNOT CONTROL
  3: C0   Colour bit 0
  4: C1   Colour bit 1
  5: I3   Isospin up(0)/down(1) ← CNOT TARGET
  6: χ    Chirality left(0)/right(1)
  7: W    Weak: doublet(0)/singlet(1)

Ring order: G0 – G1 – LQ – C0 – C1 – I3 – χ – W – (back to G0)

Four constraints (R1–R4) select exactly 45 valid states from 256.
These are LOGICAL constraints, NOT XOR parity sums — that is why the
first version of this code gave 32 instead of 45.

The 9-qubit plaquette = 8 boundary bits + 1 bulk syndrome bit (δ = 2/9).
Vacuum = 0b00000000
CNOT: I3(t+1) = I3(t) XOR LQ(t)
"""

import random
from dataclasses import dataclass, field
from typing import List, Optional, Dict
from collections import Counter

# ---------------------------------------------------------------------------
# Bit indices  (MSB = position 0)
# ---------------------------------------------------------------------------
G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
N_BITS = 8
VACUUM: int = 0

# ---------------------------------------------------------------------------
# Bit helpers
# ---------------------------------------------------------------------------
def get_bit(state: int, idx: int) -> int:
    return (state >> (7 - idx)) & 1

def flip_bit(state: int, idx: int) -> int:
    return state ^ (1 << (7 - idx))

def bits_str(state: int) -> str:
    return format(state, '08b')

# ---------------------------------------------------------------------------
# The four constraints  R1–R4  (Section 2.2)
# ---------------------------------------------------------------------------
def check_R1(state: int) -> bool:
    """Generation bound: (G0,G1) != (1,1) — only 3 generations"""
    return not (get_bit(state, G0) == 1 and get_bit(state, G1) == 1)

def check_R2(state: int) -> bool:
    """Chirality-Weak coupling: χ == W"""
    return get_bit(state, CHI) == get_bit(state, W)

def check_R3(state: int) -> bool:
    """Colour-Lepton exclusion: leptons colourless, quarks coloured"""
    lq = get_bit(state, LQ)
    c0, c1 = get_bit(state, C0), get_bit(state, C1)
    if lq == 0:
        return (c0 == 0 and c1 == 0)
    else:
        return not (c0 == 0 and c1 == 0)

def check_R4(state: int) -> bool:
    """No right-handed neutrino: forbid LQ=0, I3=0, χ=1"""
    return not (get_bit(state, LQ) == 0 and
                get_bit(state, I3) == 0 and
                get_bit(state, CHI) == 1)

RULES      = [check_R1, check_R2, check_R3, check_R4]
RULE_NAMES = ['R1(gen≠11)', 'R2(χ=W)', 'R3(colour-lepton)', 'R4(no νR)']

def is_valid(state: int) -> bool:
    return all(r(state) for r in RULES)

def n_violations(state: int) -> int:
    return sum(1 for r in RULES if not r(state))

def violated_rules(state: int) -> List[str]:
    return [name for r, name in zip(RULES, RULE_NAMES) if not r(state)]

def is_pseudocodeword(state: int) -> bool:
    """Satisfies R1-R3 but violates only R4 → sterile neutrino νR candidate"""
    return (check_R1(state) and check_R2(state) and
            check_R3(state) and not check_R4(state))

# ---------------------------------------------------------------------------
# Valid state sets
# ---------------------------------------------------------------------------
VALID_STATES:  List[int] = [s for s in range(256) if is_valid(s)]
VALID_SET:     set        = set(VALID_STATES)
PSEUDO_STATES: List[int] = [s for s in range(256) if is_pseudocodeword(s)]

# ---------------------------------------------------------------------------
# CNOT gate  (Eq. 2)
# ---------------------------------------------------------------------------
def cnot(state: int) -> int:
    if get_bit(state, LQ) == 1:
        return flip_bit(state, I3)
    return state

# ---------------------------------------------------------------------------
# Particle label
# ---------------------------------------------------------------------------
COLOUR_MAP = {(0,0):'W', (0,1):'R', (1,0):'G', (1,1):'B'}
GEN_MAP    = {(0,0):'1', (0,1):'2', (1,0):'3'}

def particle_label(state: int) -> str:
    if state == VACUUM:
        return 'vacuum'
    lq  = get_bit(state, LQ)
    gen = GEN_MAP.get((get_bit(state,G0), get_bit(state,G1)), '?')
    col = COLOUR_MAP.get((get_bit(state,C0), get_bit(state,C1)), '?')
    i3  = get_bit(state, I3)
    chi = get_bit(state, CHI)
    hand = 'R' if chi == 1 else 'L'
    if lq == 0:
        kind = 'ν' if i3 == 0 else 'e'
        return f"g{gen}-{hand}{kind}"
    else:
        kind = 'u' if i3 == 0 else 'd'
        return f"g{gen}-{hand}{kind}({col})"

# ---------------------------------------------------------------------------
# 4.8.8 lattice builder
# ---------------------------------------------------------------------------
def build_488_lattice(rows: int = 4, cols: int = 4):
    node_at: Dict[tuple, int] = {}
    nid = 0
    for r in range(2 * rows):
        for c in range(2 * cols):
            if (r + c) % 2 == 0:
                node_at[(r, c)] = nid
                nid += 1
    N = nid
    adj: Dict[int, List[int]] = {i: [] for i in range(N)}
    for (r, c), u in node_at.items():
        for dr, dc in [(-1,0),(1,0),(0,-1),(0,1)]:
            nr, nc = (r+dr) % (2*rows), (c+dc) % (2*cols)
            if (nr, nc) in node_at:
                v = node_at[(nr, nc)]
                if v not in adj[u]:
                    adj[u].append(v)
    return list(range(N)), adj, node_at

# ---------------------------------------------------------------------------
# Node
# ---------------------------------------------------------------------------
@dataclass
class Node:
    nid: int
    state: int = VACUUM
    neighbours: List[int] = field(default_factory=list)

    @property
    def is_vacuum(self):   return self.state == VACUUM
    @property
    def is_valid(self):    return self.state in VALID_SET
    @property
    def is_pseudo(self):   return is_pseudocodeword(self.state)
    @property
    def violations(self):  return n_violations(self.state)

# ---------------------------------------------------------------------------
# Tiny Universe
# ---------------------------------------------------------------------------
class TinyUniverse:
    def __init__(self, rows=4, cols=4, seed=None):
        if seed is not None:
            random.seed(seed)
        node_ids, adj, positions = build_488_lattice(rows, cols)
        self.N = len(node_ids)
        self.nodes = [Node(nid=i, state=VACUUM, neighbours=adj[i])
                      for i in node_ids]
        self.step_count = 0
        print(f"TinyUniverse: {self.N} nodes on {rows}×{cols} 4.8.8 torus")
        print(f"Valid matter states: {len(VALID_STATES)}  "
              f"(expected 45, {len(VALID_STATES)-1} non-vacuum)")
        print(f"Sterile ν (νR pseudo): {len(PSEUDO_STATES)}  (expected 3)")
        assert len(VALID_STATES) == 45, \
            f"Constraint error: expected 45, got {len(VALID_STATES)}"

    def inject(self, nid, state):
        self.nodes[nid].state = state
        v = violated_rules(state)
        s = 'valid' if not v else f'INVALID {v}'
        print(f"  Node {nid:3d} ← {bits_str(state)}  "
              f"{particle_label(state):18s}  {s}")

    def inject_random_particles(self, count):
        pool = [s for s in VALID_STATES if s != VACUUM]
        nids = random.sample(range(self.N), min(count, self.N))
        for nid in nids:
            self.nodes[nid].state = random.choice(pool)
        print(f"Injected {len(nids)} random valid particles")

    def inject_anomaly(self, nid, state=0b11100100):
        v = violated_rules(state)
        print(f"  Anomaly node {nid}: {bits_str(state)}  violated={v}")
        self.nodes[nid].state = state

    def _relax(self, node):
        """Greedy single-bit flip to reduce constraint violations."""
        if node.violations == 0:
            return
        idxs = list(range(N_BITS))
        random.shuffle(idxs)
        for bit in idxs:
            candidate = flip_bit(node.state, bit)
            if n_violations(candidate) < node.violations:
                node.state = candidate
                return

    def _propagate(self, node):
        if node.is_vacuum or not node.neighbours:
            return
        vac_nbrs = [self.nodes[nid] for nid in node.neighbours
                    if self.nodes[nid].is_vacuum]
        if vac_nbrs:
            tgt = random.choice(vac_nbrs)
            tgt.state = cnot(node.state)
            node.state = VACUUM

    def step(self, mode="cnot_relax"):
        snapshot = [n.state for n in self.nodes]
        for i, node in enumerate(self.nodes):
            node.state = snapshot[i]
            if mode in ("cnot_only", "cnot_relax"):
                node.state = cnot(node.state)
            if mode == "cnot_relax":
                self._relax(node)
            elif mode == "propagate":
                self._propagate(node)
                self._relax(node)
        self.step_count += 1

    def census(self):
        return dict(
            step      = self.step_count,
            vacuum    = sum(1 for n in self.nodes if n.is_vacuum),
            valid     = sum(1 for n in self.nodes if not n.is_vacuum and n.is_valid),
            pseudo    = sum(1 for n in self.nodes if n.is_pseudo),
            invalid   = sum(1 for n in self.nodes if not n.is_valid and not n.is_pseudo),
            violations= sum(n.violations for n in self.nodes),
        )

    def print_census(self):
        c = self.census()
        print(f"Step {c['step']:4d} | vac={c['vacuum']:3d}  "
              f"valid={c['valid']:3d}  νR={c['pseudo']}  "
              f"invalid={c['invalid']:3d}  viol={c['violations']}")

    def print_particles(self):
        print(f"\n--- Non-vacuum at step {self.step_count} ---")
        for n in self.nodes:
            if not n.is_vacuum:
                tag = 'OK' if n.is_valid else ('νR' if n.is_pseudo else 'INVALID')
                print(f"  Node {n.nid:3d}: {bits_str(n.state)}  "
                      f"{particle_label(n.state):20s}  {tag}")

    def histogram(self):
        return Counter(n.state for n in self.nodes if not n.is_vacuum)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 60)
    print("  Holographic Circlette — Tiny Universe")
    print("  D.G. Elliman, Part I, 2026")
    print("=" * 60)
    print()

    # ---- Constraint validation ----
    print("--- Constraint check ---")
    print(f"  Valid states: {len(VALID_STATES)} (expected 45)")
    print(f"  νR pseudo:    {len(PSEUDO_STATES)} (expected 3)")
    print()
    print("  νR pseudocodewords:")
    for s in PSEUDO_STATES:
        print(f"    {bits_str(s)}  {particle_label(s)}")
    print()
    print("  Sample of valid matter states (first 10):")
    for s in [x for x in VALID_STATES if x != VACUUM][:10]:
        print(f"    {bits_str(s)}  {particle_label(s)}")
    print()

    # ---- Simulation ----
    universe = TinyUniverse(rows=4, cols=4, seed=42)
    print()
    universe.inject_random_particles(count=6)
    print()
    universe.inject_anomaly(nid=0, state=0b11100100)  # G0=G1=1 → R1 violated
    print()

    print("Running 20 steps, mode='cnot_relax'")
    print()
    universe.print_census()
    for _ in range(20):
        universe.step(mode="cnot_relax")
        universe.print_census()

    universe.print_particles()

    print("\nHistogram:")
    for state, cnt in universe.histogram().most_common():
        tag = 'valid' if state in VALID_SET else ('νR' if is_pseudocodeword(state) else 'INVALID')
        print(f"  {bits_str(state)}  {particle_label(state):22s}  ×{cnt}  {tag}")
