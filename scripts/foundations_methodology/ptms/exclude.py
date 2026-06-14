"""exclude.py — the single source of truth for scripts PTMS must NEVER execute.

Two layers (defence-in-depth): a directory-prefix guard (matches at any depth) and an exact
top-level set. `is_excluded` also rejects any path containing `..` (path-escape guard). Both the
manifest builder (must never *spawn* an excluded script) and the evidence-regression check (must
*explain* why a cited-but-excluded script wasn't run) import this — one list, no drift.
"""
from __future__ import annotations

# repo-relative path prefixes that are LIVE long-running computations — never run, at any depth
DIR_PREFIXES = (
    "python_code/smg_dmrg/",   # live DMRG package (its own cli.py)
    "strong_cp/",              # live; foreign git checkout outside python_code
)

# exact repo-relative scripts excluded by name, with the reason
EXACT = {
    "python_code/css_2d_strip_hpc_lanczos.py": "needs argparse args + writes checkpoints + runs hours",
    "python_code/tenpy_fermion_smoke.py":      "needs the optional TeNPy package + a DMRG ground-state search",
}


def is_excluded(rel_path: str):
    """Return (excluded: bool, reason: str) for a repo-relative POSIX path."""
    p = (rel_path or "").replace("\\", "/")
    if ".." in p.split("/"):
        return True, "path-escape ('..') rejected"
    for pre in DIR_PREFIXES:
        if p.startswith(pre):
            return True, f"under live computation dir {pre}"
    if p in EXACT:
        return True, EXACT[p]
    return False, ""
