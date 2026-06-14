#!/usr/bin/env python3
r"""Measure the boundary-local-rescue island-floor surface rho_D/rho_B(d,L).

Registered target:
    L = 12..16, R = 496..1600, gamma-driver density included.

This worker intentionally separates three readings:

* strict growth pairing:
    adjacent same-level blocks merge in nucleation order; orphaned late wall
    cells remain shallow. This is the diagnostic counterfactual already shown to
    overshoot.

* boundary-local rescue:
    a stranded block can be padded into an ADJACENT registered higher-depth host.
    It cannot read through missing/no-register wall interiors, and disconnected
    crystal islands cannot be rescued.

* island floor:
    the residual shallow wall cells after boundary-local rescue. These are the
    policy-independent floor that can only be reduced by larger-volume/percolated
    crystallisation, not by changing the rescue rule.

Usage:
    # one job, one JSON line
    python k04_island_floor_surface.py single --L 12 --R 496 --rep 1

    # serial grid, append JSONL
    python k04_island_floor_surface.py grid --Ls 12,14,16 --Rs 496,800,1600 \
        --reps 1:4 --out k04_island_floor_surface_results.jsonl

    # write commands for deep/xargs
    python k04_island_floor_surface.py jobs --preset crossing --python /tenpy-env/bin/python \
        --out k04_island_floor_jobs.txt --results k04_island_floor_surface_results.jsonl

    # summarize an existing JSONL surface
    python k04_island_floor_surface.py analyze --input k04_island_floor_surface_results.jsonl \
        --summary-json k04_island_floor_surface_summary.json

Exit 0 means the local invariants checked. The Monte Carlo is intentionally
expensive at L=16; use jobs mode on deep for the full registered surface.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path
from statistics import fmean, median, pstdev

import numpy as np


AX = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
CUBE_SPEC = sorted([3, 1, 1, 1, -1, -1, -1, -3])

ALPHA0 = 1 / 137.0
ETA = 6.1e-10
NGAM = 2 * 1.2020569 / math.pi**2
RHO_B = ETA * NGAM * 2 * math.sqrt(2)
RDB_OBS = 0.1200 / 0.02237

DEFAULT_LS = (12, 14, 16)
DEFAULT_RS = (496, 800, 1200, 1600)
PRESETS = {
    "registered": {
        "Ls": DEFAULT_LS,
        "Rs": DEFAULT_RS,
        "reps": tuple(range(1, 5)),
        "description": "original registered surface",
    },
    "crossing": {
        "Ls": (18, 20, 22),
        "Rs": (1800, 2200, 2400, 2600, 2800, 3200),
        "reps": tuple(range(1, 17)),
        "description": "zero-inflated island-floor crossing band",
    },
}


def q1_of(p: float) -> float:
    return sum(
        math.comb(8, k)
        * p**k
        * (1 - p) ** (8 - k)
        * (0.0 if k <= 3 else (0.5 if k == 4 else 1.0))
        for k in range(9)
    )


Q1 = q1_of(0.0972)


def resid(level: int) -> float:
    return (21 * Q1) ** (2 ** (level - 1)) / 21


def parse_csv_ints(text: str) -> list[int]:
    return [int(part) for part in text.split(",") if part.strip()]


def parse_reps(text: str) -> list[int]:
    if ":" in text:
        a, b = [int(part) for part in text.split(":", 1)]
        return list(range(a, b + 1))
    return parse_csv_ints(text)


def resolve_sweep(args: argparse.Namespace) -> tuple[list[int], list[int], list[int]]:
    if args.preset:
        preset = PRESETS[args.preset]
        ls = list(preset["Ls"]) if args.Ls is None else parse_csv_ints(args.Ls)
        rs = list(preset["Rs"]) if args.Rs is None else parse_csv_ints(args.Rs)
        reps = list(preset["reps"]) if args.reps is None else parse_reps(args.reps)
        return ls, rs, reps
    ls = parse_csv_ints(args.Ls or ",".join(map(str, DEFAULT_LS)))
    rs = parse_csv_ints(args.Rs or ",".join(map(str, DEFAULT_RS)))
    reps = parse_reps(args.reps or "1:4")
    return ls, rs, reps


def stable_seed(L: int, R: int, rep: int) -> int:
    # Keep the earlier pairing_from_dynamics convention so old/new rows compare.
    return 1000 * R + rep + 1000003 * L


class EmbeddedK04:
    def __init__(self, L: int, seed: int):
        assert L % 2 == 0 and L >= 4
        self.L = L
        self.N = L**3
        self.rng = random.Random(seed)
        self.sites = [
            (x, y, z) for x in range(L) for y in range(L) for z in range(L)
        ]
        self.squares = [
            (p, a, b)
            for p in self.sites
            for a in range(3)
            for b in range(3)
            if a < b
        ]

    def wrap(self, p: tuple[int, int, int]) -> tuple[int, int, int]:
        L = self.L
        return (p[0] % L, p[1] % L, p[2] % L)

    @staticmethod
    def add(
        p: tuple[int, int, int], s: tuple[int, int, int]
    ) -> tuple[int, int, int]:
        return (p[0] + s[0], p[1] + s[1], p[2] + s[2])

    def tiling(self) -> set[tuple[tuple[int, int, int], int]]:
        return {(p, a) for p in self.sites for a in range(3) if p[a] % 2 == 0}

    def square_bonds(
        self, p: tuple[int, int, int], a: int, b: int
    ) -> tuple[
        tuple[tuple[tuple[int, int, int], int], tuple[tuple[int, int, int], int]],
        tuple[tuple[tuple[int, int, int], int], tuple[tuple[int, int, int], int]],
    ]:
        return (
            ((p, a), (self.wrap(self.add(p, AX[b])), a)),
            ((p, b), (self.wrap(self.add(p, AX[a])), b)),
        )

    def propose(
        self, bonds: set[tuple[tuple[int, int, int], int]]
    ) -> tuple[
        tuple[tuple[tuple[int, int, int], int], ...],
        tuple[tuple[tuple[int, int, int], int], ...],
    ] | None:
        p, a, b = self.rng.choice(self.squares)
        pair_a, pair_b = self.square_bonds(p, a, b)
        in_a = sum(1 for edge in pair_a if edge in bonds)
        in_b = sum(1 for edge in pair_b if edge in bonds)
        if in_a == 2 and in_b == 0:
            return pair_a, pair_b
        if in_b == 2 and in_a == 0:
            return pair_b, pair_a
        return None

    def nbrs(
        self, bonds: set[tuple[tuple[int, int, int], int]], p: tuple[int, int, int]
    ) -> list[tuple[tuple[int, int, int], tuple[int, int, int]]]:
        out = []
        for a in range(3):
            if (p, a) in bonds:
                out.append((self.wrap(self.add(p, AX[a])), AX[a]))
            neg = (-AX[a][0], -AX[a][1], -AX[a][2])
            m = self.wrap(self.add(p, neg))
            if (m, a) in bonds:
                out.append((m, neg))
        return out

    def cycles_through(
        self,
        bonds: set[tuple[tuple[int, int, int], int]],
        changed: tuple[tuple[tuple[int, int, int], int], ...],
        k: int,
    ) -> set[tuple[tuple[int, int, int], ...]]:
        cycles = set()
        for p, a in changed:
            if (p, a) not in bonds:
                continue
            q = self.wrap(self.add(p, AX[a]))
            stack = [(q, AX[a], (p, q))]
            while stack:
                v, disp, path = stack.pop()
                for w, step in self.nbrs(bonds, v):
                    nd = self.add(disp, step)
                    if len(path) == k:
                        if w == p and nd == (0, 0, 0):
                            best = None
                            for i in range(len(path)):
                                rot = path[i:] + path[:i]
                                for cand in (rot, tuple(reversed(rot))):
                                    if best is None or cand < best:
                                        best = cand
                            cycles.add(best)
                    elif (
                        w not in path
                        and abs(nd[0]) + abs(nd[1]) + abs(nd[2])
                        <= k - len(path)
                    ):
                        stack.append((w, nd, path + (w,)))
        return cycles

    def cells_of(
        self, bonds: set[tuple[tuple[int, int, int], int]]
    ) -> list[frozenset[tuple[int, int, int]]]:
        adj: dict[tuple[int, int, int], list[tuple[int, int, int]]] = {}
        for p, a in bonds:
            q = self.wrap(self.add(p, AX[a]))
            adj.setdefault(p, []).append(q)
            adj.setdefault(q, []).append(p)

        seen: set[tuple[int, int, int]] = set()
        cells = []
        for s0 in adj:
            if s0 in seen:
                continue
            comp = {s0}
            stack = [s0]
            while stack:
                v = stack.pop()
                for w in adj[v]:
                    if w not in comp:
                        comp.add(w)
                        stack.append(w)
            seen |= comp
            if len(comp) != 8:
                continue
            idx = {v: i for i, v in enumerate(sorted(comp))}
            matrix = np.zeros((8, 8))
            for v in comp:
                for w in adj[v]:
                    matrix[idx[v], idx[w]] = 1.0
            spec = sorted(np.round(np.linalg.eigvalsh(matrix)).astype(int).tolist())
            if spec == CUBE_SPEC:
                cells.append(frozenset(comp))
        return cells

    def randomize_hot(
        self, bonds: set[tuple[tuple[int, int, int], int]], sweeps: int = 60
    ) -> set[tuple[tuple[int, int, int], int]]:
        for _ in range(sweeps):
            for _ in range(self.N):
                move = self.propose(bonds)
                if move:
                    remove, add = move
                    bonds = (bonds - set(remove)) | set(add)
        return bonds

    def run_with_history(
        self, R: int, hold: int | None = None, sample_every: int = 10
    ) -> tuple[
        list[frozenset[tuple[int, int, int]]],
        dict[frozenset[tuple[int, int, int]], int],
        set[tuple[int, int, int]],
        float,
    ]:
        bonds = self.randomize_hot(self.tiling())
        w4, w6 = 2.0, 1.0
        hold_sweeps = max(20, R // 10) if hold is None else hold
        total = R + hold_sweeps
        snapshots: list[set[frozenset[tuple[int, int, int]]]] = []

        for sweep in range(total):
            ramp_step = min(sweep, R)
            T = 6.0 * (0.5 / 6.0) ** (ramp_step / max(R - 1, 1))
            for _ in range(self.N):
                move = self.propose(bonds)
                if move is None:
                    continue
                remove, add = move
                trial = (bonds - set(remove)) | set(add)
                d4 = len(self.cycles_through(trial, add, 4)) - len(
                    self.cycles_through(bonds, remove, 4)
                )
                d6 = len(self.cycles_through(trial, add, 6)) - len(
                    self.cycles_through(bonds, remove, 6)
                )
                dE = -w4 * d4 - w6 * d6
                if dE <= 0 or self.rng.random() < math.exp(-dE / T):
                    bonds = trial
            if sweep % sample_every == 0:
                snapshots.append(set(self.cells_of(bonds)))

        final = self.cells_of(bonds)
        snapshots.append(set(final))
        t_nuc = {}
        for cell in final:
            first = len(snapshots) - 1
            for idx in range(len(snapshots) - 1, -1, -1):
                if cell in snapshots[idx]:
                    first = idx
                else:
                    break
            t_nuc[cell] = first

        good = set().union(*final) if final else set()
        d = 1 - len(good) / self.N
        return final, t_nuc, good, d


def cells_adjacent(
    cell_a: frozenset[tuple[int, int, int]],
    cell_b: frozenset[tuple[int, int, int]],
    model: EmbeddedK04,
) -> bool:
    return any(
        model.wrap(model.add(v, step)) in cell_b
        for v in cell_a
        for axis in range(3)
        for step in (AX[axis], (-AX[axis][0], -AX[axis][1], -AX[axis][2]))
    )


def verts_adjacent(
    verts_a: set[tuple[int, int, int]],
    verts_b: set[tuple[int, int, int]],
    model: EmbeddedK04,
) -> bool:
    return any(
        model.wrap(model.add(v, step)) in verts_b
        for v in verts_a
        for axis in range(3)
        for step in (AX[axis], (-AX[axis][0], -AX[axis][1], -AX[axis][2]))
    )


def agglomerate(
    cells: list[frozenset[tuple[int, int, int]]],
    t_nuc: dict[frozenset[tuple[int, int, int]], int],
    model: EmbeddedK04,
) -> tuple[dict[frozenset[tuple[int, int, int]], int], list[dict]]:
    blocks = [
        {"level": 1, "cells": {cell}, "t": t_nuc[cell], "verts": set(cell)}
        for cell in cells
    ]
    changed = True
    while changed:
        changed = False
        blocks.sort(key=lambda b: (b["t"], min(sorted(b["verts"]))))
        for i, bi in enumerate(blocks):
            for j in range(i + 1, len(blocks)):
                bj = blocks[j]
                if bi["level"] != bj["level"] or bi["level"] >= 6:
                    continue
                if not verts_adjacent(bi["verts"], bj["verts"], model):
                    continue
                merged = {
                    "level": bi["level"] + 1,
                    "cells": bi["cells"] | bj["cells"],
                    "t": max(bi["t"], bj["t"]),
                    "verts": bi["verts"] | bj["verts"],
                }
                blocks = [b for k, b in enumerate(blocks) if k not in (i, j)]
                blocks.append(merged)
                changed = True
                break
            if changed:
                break

    depth = {}
    for block in blocks:
        for cell in block["cells"]:
            depth[cell] = block["level"]
    return depth, blocks


def boundary_local_rescue(
    blocks: list[dict],
    depth: dict[frozenset[tuple[int, int, int]], int],
    model: EmbeddedK04,
) -> dict[frozenset[tuple[int, int, int]], int]:
    """Absorb each stranded block into the deepest adjacent registered host.

    One pass over the original block adjacency. This deliberately does not let a
    newly padded orphan become a bridge to rescue another orphan; that would read
    through unregistered padding rather than an actually formed higher block.
    """
    depth2 = dict(depth)
    for block in sorted(blocks, key=lambda b: b["level"]):
        hosts = [
            host
            for host in blocks
            if host is not block
            and host["level"] > block["level"]
            and verts_adjacent(block["verts"], host["verts"], model)
        ]
        if not hosts:
            continue
        host = max(hosts, key=lambda h: (h["level"], len(h["cells"])))
        for cell in block["cells"]:
            depth2[cell] = host["level"]
    return depth2


def wall_cells(
    cells: list[frozenset[tuple[int, int, int]]],
    good: set[tuple[int, int, int]],
    model: EmbeddedK04,
) -> list[frozenset[tuple[int, int, int]]]:
    out = []
    for cell in cells:
        touches_bad = any(
            model.wrap(model.add(v, step)) not in good
            for v in cell
            for axis in range(3)
            for step in (AX[axis], (-AX[axis][0], -AX[axis][1], -AX[axis][2]))
        )
        if touches_bad:
            out.append(cell)
    return out


def cell_components(
    cells: list[frozenset[tuple[int, int, int]]], model: EmbeddedK04
) -> list[list[frozenset[tuple[int, int, int]]]]:
    remaining = set(cells)
    components = []
    while remaining:
        start = remaining.pop()
        comp = [start]
        queue = deque([start])
        while queue:
            cell = queue.popleft()
            neighbours = [other for other in list(remaining) if cells_adjacent(cell, other, model)]
            for other in neighbours:
                remaining.remove(other)
                comp.append(other)
                queue.append(other)
        components.append(comp)
    return components


def depth_hist(
    cells: list[frozenset[tuple[int, int, int]]],
    depth: dict[frozenset[tuple[int, int, int]], int],
) -> dict[str, int]:
    hist = Counter(depth[cell] for cell in cells)
    return {str(k): hist[k] for k in sorted(hist)}


def shadow_readout(
    wall: list[frozenset[tuple[int, int, int]]],
    depth: dict[frozenset[tuple[int, int, int]], int],
    L: int,
) -> tuple[float, float, float]:
    if not wall:
        return 0.0, 0.0, 0.0
    r_eff = fmean(resid(depth[cell]) for cell in wall)
    s1 = 8 * len(wall) / L**3
    rdb = s1 * r_eff * ALPHA0 / RHO_B
    stall_frac = sum(1 for cell in wall if depth[cell] <= 2) / len(wall)
    return r_eff, rdb, stall_frac


def measure_single(L: int, R: int, rep: int, hold: int | None = None) -> dict:
    model = EmbeddedK04(L, stable_seed(L, R, rep))
    cells, t_nuc, good, d = model.run_with_history(R, hold=hold)
    wall = wall_cells(cells, good, model)
    depth_strict, blocks = agglomerate(cells, t_nuc, model)
    depth_rescue = boundary_local_rescue(blocks, depth_strict, model)

    r_eff, rdb, stall_frac = shadow_readout(wall, depth_strict, L)
    r_eff_rescue, rdb_rescue, stall_frac_rescue = shadow_readout(
        wall, depth_rescue, L
    )

    components = cell_components(cells, model)
    component_records = []
    island_floor_wall = 0
    island_floor_cells = 0
    for comp in components:
        comp_wall = [cell for cell in comp if cell in wall]
        comp_depths = [depth_rescue[cell] for cell in comp]
        comp_shallow_wall = [cell for cell in comp_wall if depth_rescue[cell] <= 2]
        if comp_shallow_wall:
            island_floor_wall += len(comp_shallow_wall)
        if max(comp_depths, default=0) <= 2:
            island_floor_cells += len(comp)
        component_records.append(
            {
                "size": len(comp),
                "wall": len(comp_wall),
                "max_depth_rescue": max(comp_depths, default=0),
                "shallow_wall": len(comp_shallow_wall),
            }
        )

    component_size_hist = Counter(record["size"] for record in component_records)
    island_component_records = [
        record for record in component_records if record["shallow_wall"] > 0
    ]
    shallow_component_size_hist = Counter(record["size"] for record in island_component_records)
    island_floor_components = len(island_component_records)
    island_floor_component_cells = sum(record["size"] for record in island_component_records)
    island_floor_component_max_size = max(
        (record["size"] for record in island_component_records), default=0
    )
    island_floor_shallow_wall_max = max(
        (record["shallow_wall"] for record in island_component_records), default=0
    )
    t_stall = [t_nuc[cell] for cell in wall if depth_strict[cell] <= 2]
    t_deep = [t_nuc[cell] for cell in wall if depth_strict[cell] >= 3]
    t_stall_rescue = [t_nuc[cell] for cell in wall if depth_rescue[cell] <= 2]
    t_deep_rescue = [t_nuc[cell] for cell in wall if depth_rescue[cell] >= 3]

    return {
        "L": L,
        "R": R,
        "rep": rep,
        "seed": stable_seed(L, R, rep),
        "hold": max(20, R // 10) if hold is None else hold,
        "d": d,
        "ncells": len(cells),
        "nwall": len(wall),
        "s_wall": 8 * len(wall) / L**3,
        "n_components": len(components),
        "component_size_hist": {str(k): v for k, v in sorted(component_size_hist.items())},
        "island_component_size_hist": {
            str(k): v for k, v in sorted(shallow_component_size_hist.items())
        },
        "hist_strict": depth_hist(wall, depth_strict),
        "r_eff_strict": r_eff,
        "rdb_strict": rdb,
        "overshoot_strict": rdb / RDB_OBS if RDB_OBS else float("inf"),
        "stall_frac_strict": stall_frac,
        "hist_rescue": depth_hist(wall, depth_rescue),
        "r_eff_rescue": r_eff_rescue,
        "rdb_rescue": rdb_rescue,
        "overshoot_rescue": rdb_rescue / RDB_OBS if RDB_OBS else float("inf"),
        "stall_frac_rescue": stall_frac_rescue,
        "island_floor_wall": island_floor_wall,
        "island_floor_cells": island_floor_cells,
        "island_floor_components": island_floor_components,
        "island_floor_component_cells": island_floor_component_cells,
        "island_floor_component_max_size": island_floor_component_max_size,
        "island_floor_shallow_wall_max": island_floor_shallow_wall_max,
        "island_floor_wall_frac": island_floor_wall / max(len(wall), 1),
        "island_floor_cell_frac": island_floor_cells / max(len(cells), 1),
        "t_nuc_stall_strict": fmean(t_stall) if t_stall else None,
        "t_nuc_deep_strict": fmean(t_deep) if t_deep else None,
        "t_nuc_stall_rescue": fmean(t_stall_rescue) if t_stall_rescue else None,
        "t_nuc_deep_rescue": fmean(t_deep_rescue) if t_deep_rescue else None,
    }


def run_grid(args: argparse.Namespace) -> int:
    ls, rs, reps = resolve_sweep(args)
    out_path = Path(args.out)
    if args.preset:
        print(
            f"preset={args.preset} ({PRESETS[args.preset]['description']}): "
            f"Ls={ls}, Rs={rs}, reps={reps}",
            file=sys.stderr,
        )
    with out_path.open("a", encoding="utf-8") as fh:
        for L in ls:
            for R in rs:
                for rep in reps:
                    row = measure_single(L, R, rep, hold=args.hold)
                    fh.write(json.dumps(row, sort_keys=True) + "\n")
                    fh.flush()
                    print(
                        f"L={L} R={R} rep={rep}: d={row['d']:.3f}, "
                        f"rdb_rescue={row['rdb_rescue']:.3g}, "
                        f"island_floor={row['island_floor_wall_frac']:.3f}",
                        file=sys.stderr,
                    )
    return 0


def write_jobs(args: argparse.Namespace) -> int:
    ls, rs, reps = resolve_sweep(args)
    script = Path(__file__).name
    out_path = Path(args.out)
    result_path = args.results
    commands = []
    for L in ls:
        for R in rs:
            for rep in reps:
                hold_part = f" --hold {args.hold}" if args.hold is not None else ""
                cmd = (
                    f"{args.python} {script} single --L {L} --R {R} --rep {rep}"
                    f"{hold_part} >> {result_path}"
                )
                commands.append(cmd)
    out_path.write_text("\n".join(commands) + "\n", encoding="utf-8")
    print(f"wrote {len(commands)} jobs -> {out_path}")
    if args.preset:
        print(f"preset={args.preset}: {PRESETS[args.preset]['description']}")
    print(f"Ls={ls}")
    print(f"Rs={rs}")
    print(f"reps={reps}")
    print(f"run on deep from python_code/: xargs -a {out_path.name} -P 22 -I CMD bash -c 'CMD'")
    return 0


def load_jsonl(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def summarize(rows: list[dict]) -> list[dict]:
    grouped: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(int(row["L"]), int(row["R"]))].append(row)
    out = []
    for (L, R), rs in sorted(grouped.items()):
        def vals(key: str, default: float | None = None) -> list[float]:
            out_vals = []
            for row in rs:
                value = row.get(key, default)
                if value is not None:
                    out_vals.append(float(value))
            return out_vals

        def component_hist(row: dict) -> dict[int, int]:
            return {
                int(size): int(count)
                for size, count in row.get("island_component_size_hist", {}).items()
            }

        def island_component_count(row: dict) -> int:
            if "island_floor_components" in row:
                return int(row["island_floor_components"])
            return sum(component_hist(row).values())

        def island_component_cells(row: dict) -> int:
            if "island_floor_component_cells" in row:
                return int(row["island_floor_component_cells"])
            return sum(size * count for size, count in component_hist(row).items())

        def island_component_max_size(row: dict) -> int:
            if "island_floor_component_max_size" in row:
                return int(row["island_floor_component_max_size"])
            hist = component_hist(row)
            return max(hist, default=0)

        def mean_or_zero(numbers: list[float]) -> float:
            return fmean(numbers) if numbers else 0.0

        def sd_or_zero(numbers: list[float]) -> float:
            return pstdev(numbers) if len(numbers) > 1 else 0.0

        def se_or_zero(numbers: list[float]) -> float:
            return sd_or_zero(numbers) / math.sqrt(len(numbers)) if numbers else 0.0

        def med_or_zero(numbers: list[float]) -> float:
            return median(numbers) if numbers else 0.0

        def nonzero(numbers: list[float]) -> list[float]:
            return [x for x in numbers if x > 0.0]

        d_values = vals("d")
        rdb_values = vals("rdb_rescue")
        overshoot_values = vals("overshoot_rescue")
        island_wall_frac_values = vals("island_floor_wall_frac", 0.0)
        stall_values = vals("stall_frac_rescue", 0.0)
        nwall_values = vals("nwall", 0.0)
        island_wall_counts = vals("island_floor_wall", 0.0)
        island_cell_counts = vals("island_floor_cells", 0.0)
        island_component_counts = [float(island_component_count(row)) for row in rs]
        island_component_cell_counts = [float(island_component_cells(row)) for row in rs]
        island_component_max_sizes = [float(island_component_max_size(row)) for row in rs]
        rdb_nonzero = nonzero(rdb_values)
        island_wall_nonzero = nonzero(island_wall_counts)
        island_component_nonzero = nonzero(island_component_counts)

        record = {
            "L": L,
            "R": R,
            "n": len(rs),
            "d_mean": mean_or_zero(d_values),
            "d_sd": sd_or_zero(d_values),
            "rdb_rescue_mean": mean_or_zero(rdb_values),
            "rdb_rescue_sd": sd_or_zero(rdb_values),
            "rdb_rescue_se": se_or_zero(rdb_values),
            "rdb_rescue_zero_frac": 1.0 - len(rdb_nonzero) / max(len(rdb_values), 1),
            "rdb_rescue_nonzero_n": len(rdb_nonzero),
            "rdb_rescue_nonzero_mean": mean_or_zero(rdb_nonzero),
            "rdb_rescue_nonzero_median": med_or_zero(rdb_nonzero),
            "rdb_rescue_nonzero_max": max(rdb_nonzero, default=0.0),
            "overshoot_rescue_mean": mean_or_zero(overshoot_values),
            "island_floor_wall_frac_mean": mean_or_zero(island_wall_frac_values),
            "stall_frac_rescue_mean": mean_or_zero(stall_values),
            "nwall_mean": mean_or_zero(nwall_values),
            "island_floor_hit_frac": len(island_wall_nonzero)
            / max(len(island_wall_counts), 1),
            "island_floor_zero_frac": 1.0
            - len(island_wall_nonzero) / max(len(island_wall_counts), 1),
            "island_floor_wall_mean": mean_or_zero(island_wall_counts),
            "island_floor_wall_sd": sd_or_zero(island_wall_counts),
            "island_floor_wall_se": se_or_zero(island_wall_counts),
            "island_floor_wall_nonzero_mean": mean_or_zero(island_wall_nonzero),
            "island_floor_wall_nonzero_median": med_or_zero(island_wall_nonzero),
            "island_floor_wall_max": max(island_wall_counts, default=0.0),
            "island_floor_cell_mean": mean_or_zero(island_cell_counts),
            "island_floor_cell_max": max(island_cell_counts, default=0.0),
            "island_floor_components_mean": mean_or_zero(island_component_counts),
            "island_floor_components_nonzero_mean": mean_or_zero(
                island_component_nonzero
            ),
            "island_floor_components_max": max(island_component_counts, default=0.0),
            "island_floor_component_cells_mean": mean_or_zero(
                island_component_cell_counts
            ),
            "island_floor_component_max_size_mean": mean_or_zero(
                island_component_max_sizes
            ),
            "island_floor_component_max_size_max": max(
                island_component_max_sizes, default=0.0
            ),
        }
        out.append(record)
    return out


def analyze(args: argparse.Namespace) -> int:
    rows = load_jsonl(Path(args.input))
    summary = summarize(rows)
    if args.summary_json:
        summary_path = Path(args.summary_json)
        summary_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print("[0] island-floor surface summary")
    print(f"    rows: {len(rows)} from {args.input}")
    print(f"    target rho_D/rho_B: {RDB_OBS:.3f}")
    if args.summary_json:
        print(f"    wrote machine-readable summary: {args.summary_json}")
    print()
    print("L   R     n   d_mean  rdb_rescue  x_obs  island_floor  stall_rescue")
    for rec in summary:
        print(
            f"{rec['L']:<3d} {rec['R']:<5d} {rec['n']:<3d} "
            f"{rec['d_mean']:.3f}   {rec['rdb_rescue_mean']:.3g}      "
            f"{rec['overshoot_rescue_mean']:.2f}   "
            f"{rec['island_floor_wall_frac_mean']:.3f}         "
            f"{rec['stall_frac_rescue_mean']:.3f}"
        )

    print()
    print("[1] zero-inflated rescue-shadow statistics")
    print("    rdb_mean is unconditional; rdb_nz_mean/median condition on nonzero rows.")
    print("L   R     n   p0_rdb  nz  rdb_mean±se       rdb_nz_mean  rdb_nz_med  rdb_max")
    for rec in summary:
        print(
            f"{rec['L']:<3d} {rec['R']:<5d} {rec['n']:<3d} "
            f"{rec['rdb_rescue_zero_frac']:.2f}    "
            f"{rec['rdb_rescue_nonzero_n']:<3d} "
            f"{rec['rdb_rescue_mean']:.3g}±{rec['rdb_rescue_se']:.2g}   "
            f"{rec['rdb_rescue_nonzero_mean']:.3g}        "
            f"{rec['rdb_rescue_nonzero_median']:.3g}       "
            f"{rec['rdb_rescue_nonzero_max']:.3g}"
        )

    print()
    print("[2] island-count statistics")
    print("    hit = fraction of rows with island_floor_wall > 0.")
    print("L   R     n   hit  wall_mean±se  wall_nz_mean  wall_nz_med  wall_max  comp_mean  comp_max")
    for rec in summary:
        print(
            f"{rec['L']:<3d} {rec['R']:<5d} {rec['n']:<3d} "
            f"{rec['island_floor_hit_frac']:.2f} "
            f"{rec['island_floor_wall_mean']:.2f}±{rec['island_floor_wall_se']:.2f}      "
            f"{rec['island_floor_wall_nonzero_mean']:.2f}          "
            f"{rec['island_floor_wall_nonzero_median']:.2f}         "
            f"{rec['island_floor_wall_max']:.0f}        "
            f"{rec['island_floor_components_mean']:.2f}       "
            f"{rec['island_floor_components_max']:.0f}"
        )

    gamma_rows = [rec for rec in summary if rec["R"] == 496]
    if gamma_rows:
        print()
        print("[3] gamma-driver slice R=496")
        for rec in gamma_rows:
            print(
                f"    L={rec['L']}: d={rec['d_mean']:.3f}, "
                f"rho_D/rho_B={rec['rdb_rescue_mean']:.3g} "
                f"({rec['overshoot_rescue_mean']:.2f}x obs), "
                f"island floor={rec['island_floor_wall_frac_mean']:.3f}"
            )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    single = sub.add_parser("single", help="run one L/R/rep measurement")
    single.add_argument("--L", type=int, required=True)
    single.add_argument("--R", type=int, required=True)
    single.add_argument("--rep", type=int, required=True)
    single.add_argument("--hold", type=int, default=None)

    grid = sub.add_parser("grid", help="run a serial grid and append JSONL")
    grid.add_argument("--preset", choices=sorted(PRESETS), default=None)
    grid.add_argument("--Ls", default=None, help="comma list; overrides preset/default")
    grid.add_argument("--Rs", default=None, help="comma list; overrides preset/default")
    grid.add_argument("--reps", default=None, help="N:M or comma list; overrides preset/default")
    grid.add_argument("--hold", type=int, default=None)
    grid.add_argument("--out", default="k04_island_floor_surface_results.jsonl")

    jobs = sub.add_parser("jobs", help="write xargs-friendly deep job commands")
    jobs.add_argument("--preset", choices=sorted(PRESETS), default=None)
    jobs.add_argument("--Ls", default=None, help="comma list; overrides preset/default")
    jobs.add_argument("--Rs", default=None, help="comma list; overrides preset/default")
    jobs.add_argument("--reps", default=None, help="N:M or comma list; overrides preset/default")
    jobs.add_argument("--hold", type=int, default=None)
    jobs.add_argument("--python", default="/tenpy-env/bin/python")
    jobs.add_argument("--out", default="k04_island_floor_jobs.txt")
    jobs.add_argument("--results", default="k04_island_floor_surface_results.jsonl")

    analyse = sub.add_parser("analyze", aliases=["analyse"], help="summarize JSONL")
    analyse.add_argument("--input", default="k04_island_floor_surface_results.jsonl")
    analyse.add_argument("--summary-json", default=None)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd == "single":
        print(json.dumps(measure_single(args.L, args.R, args.rep, args.hold), sort_keys=True))
        return 0
    if args.cmd == "grid":
        return run_grid(args)
    if args.cmd == "jobs":
        return write_jobs(args)
    if args.cmd in {"analyze", "analyse"}:
        return analyze(args)
    raise AssertionError(args.cmd)


if __name__ == "__main__":
    raise SystemExit(main())
