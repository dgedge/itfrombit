"""manifest.py — measure + classify the self-asserting scripts, for the evidence runner.

Each runnable script is probed ONCE via an isolated subprocess (NOT runpy — a hang or a stray
sys.exit must not affect PTMS itself; runpy stays only in numeric_tieouts for its one trusted,
fast, __main__-guarded script). The result lands in `evidence_manifest.json`, which doubles as
the incremental cache (skip scripts whose mtime is unchanged) and the known-good evidence
baseline. Excluded scripts are classified statically but NEVER executed. Pure stdlib.
"""
from __future__ import annotations

import platform
import re
import subprocess
import sys
import time
from pathlib import Path

from .exclude import is_excluded

SCHEMA_VERSION = "0.2"
FAST_THRESHOLD_S = 10.0
DEFAULT_TIMEOUT_S = 30

_DEF_CHECK = re.compile(r"(?m)^\s*def\s+check\s*\(")
_ASSERT = re.compile(r"(?m)^\s*assert\s")
_RAISE = re.compile(r"(?m)^\s*raise\s+[A-Za-z_][\w.]*Error")


def classify_assertions(src: str) -> str:
    """check | assert | none — how (and whether) the script gates its exit code.
    'none' = pure-print: exit 0 is meaningless as evidence. raise-on-failure and
    np.testing.assert* fold into 'assert' (exit-affecting)."""
    has_check = bool(_DEF_CHECK.search(src))
    has_sysexit = "sys.exit(" in src
    has_assert = bool(_ASSERT.search(src))
    has_raise = bool(_RAISE.search(src))
    has_nptest = "np.testing.assert" in src or ".assert_allclose" in src
    if has_check and has_sysexit:
        return "check"
    if has_sysexit or has_assert or has_raise or has_nptest:
        return "assert"
    return "none"


def probe_one(abs_path: Path, repo_root: Path, timeout_s=DEFAULT_TIMEOUT_S,
              fast_threshold=FAST_THRESHOLD_S) -> dict:
    rel = abs_path.relative_to(repo_root).as_posix()
    excl, reason = is_excluded(rel)
    assert not excl, f"refusing to execute excluded script {rel}: {reason}"   # call-site guard
    src = abs_path.read_text("utf-8", "replace")
    t0 = time.monotonic()
    timed_out, exit_code, stderr_tail = False, None, ""
    try:
        proc = subprocess.run([sys.executable, str(abs_path)], cwd=str(repo_root),
                              stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                              timeout=timeout_s, check=False)
        exit_code = proc.returncode
        stderr_tail = proc.stderr.decode("utf-8", "replace")[-400:]
    except subprocess.TimeoutExpired:
        timed_out = True
    dur = round(time.monotonic() - t0, 3)
    return {"path": rel, "mtime": round(abs_path.stat().st_mtime, 3), "duration_s": dur,
            "exit_code": exit_code, "timed_out": timed_out,
            "assertion_class": classify_assertions(src),
            "tier": "fast" if (not timed_out and dur < fast_threshold) else "slow",
            "stderr_tail": stderr_tail}


def _stub_excluded(abs_path: Path, repo_root: Path, reason: str) -> dict:
    rel = abs_path.relative_to(repo_root).as_posix()
    src = abs_path.read_text("utf-8", "replace") if abs_path.exists() else ""
    return {"path": rel, "mtime": round(abs_path.stat().st_mtime, 3) if abs_path.exists() else None,
            "excluded": True, "reason": reason, "exit_code": None, "timed_out": False,
            "assertion_class": classify_assertions(src), "tier": "excluded", "stderr_tail": ""}


def cited_scripts(registry) -> list:
    return sorted({sc for c in registry for sc in c.cited_scripts})


def all_scripts(corpus) -> list:
    d = corpus.root / "python_code"
    return sorted(f"python_code/{p.name}" for p in d.glob("*.py")) if d.exists() else []


def build(corpus, registry, scope="cited", force=False, timeout=DEFAULT_TIMEOUT_S,
          fast_threshold=FAST_THRESHOLD_S, prior=None) -> dict:
    prior_scripts = (prior or {}).get("scripts", {})
    targets = cited_scripts(registry) if scope == "cited" else all_scripts(corpus)
    out = {}
    for rel in targets:
        abs_path = corpus.script_path(rel)
        if not abs_path.exists():
            continue                                  # missing files are citation_integrity's job
        excl, reason = is_excluded(rel)
        if excl:
            out[rel] = _stub_excluded(abs_path, corpus.root, reason)
            continue
        cur_mtime = round(abs_path.stat().st_mtime, 3)
        cached = prior_scripts.get(rel)
        if (not force) and cached and not cached.get("excluded") \
                and cached.get("mtime") == cur_mtime and "exit_code" in cached:
            rec = dict(cached); rec["cached"] = True            # unchanged -> reuse, skip run
        else:
            rec = probe_one(abs_path, corpus.root, timeout, fast_threshold); rec["cached"] = False
        out[rel] = rec
    return {"schema_version": SCHEMA_VERSION, "generated_by": "ptms manifest",
            "python": platform.python_version(), "fast_threshold_s": fast_threshold,
            "build_timeout_s": timeout, "scope": scope, "scripts": out}


def load_manifest(path) -> dict:
    p = Path(path)
    if not p.exists():
        return {}
    import json
    return json.loads(p.read_text(encoding="utf-8"))


def save_manifest(path, data) -> None:
    import json
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def lookup(manifest, rel_path):
    return manifest.get("scripts", {}).get(rel_path)


def summary(manifest) -> str:
    sc = manifest.get("scripts", {})
    tiers, classes = {}, {}
    nonzero, timedout, excluded, cached = [], [], [], 0
    for rel, r in sc.items():
        tiers[r.get("tier")] = tiers.get(r.get("tier"), 0) + 1
        classes[r.get("assertion_class")] = classes.get(r.get("assertion_class"), 0) + 1
        if r.get("excluded"):
            excluded.append(rel)
        elif r.get("timed_out"):
            timedout.append(rel)
        elif r.get("exit_code") not in (0, None):
            nonzero.append(rel)
        if r.get("cached"):
            cached += 1
    L = ["=" * 70, " PTMS evidence manifest", "=" * 70,
         f"  scope={manifest.get('scope')}  scripts={len(sc)}  cached(skip-run)={cached}",
         f"  tiers: {tiers}", f"  assertion_class: {classes}",
         f"  nonzero-exit: {len(nonzero)} {nonzero or ''}",
         f"  timed-out: {len(timedout)} {timedout or ''}",
         f"  excluded (never run): {len(excluded)} {excluded or ''}"]
    return "\n".join(L)
