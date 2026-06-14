"""checks/numeric_tieouts.py — absorb methodology_metrics.py (unmodified) via runpy.

HARD if the absorbed script's internal asserts break (a source moved / counts no longer tie
out → exit != 0). The paper-prose discrepancies it detects but only *warns* about (the §6
75-vs-78, the 29-vs-40) are surfaced as SOFT — matching the script's own severity (it does not
fail on them; they are about the methodology PAPER, not the authoritative sources).
"""
from __future__ import annotations

import contextlib
import io
import runpy
from pathlib import Path

from ..findings import Finding, HARD, SOFT

CHECK_ID = "numeric_tieouts"
SCRIPT = "ai_methodology/methodology_metrics.py"


def run(corpus, registry) -> list:
    path = Path(corpus.root) / SCRIPT
    if not path.exists():
        return [Finding(CHECK_ID, HARD, SCRIPT, 0, "absorbed tie-out script missing")]
    buf = io.StringIO()
    code = 0
    try:
        with contextlib.redirect_stdout(buf):
            runpy.run_path(str(path), run_name="__main__")
    except SystemExit as e:
        code = int(e.code or 0)
    except Exception as ex:  # the script raised — itself a failure of the evidence binding
        return [Finding(CHECK_ID, HARD, SCRIPT, 0,
                        f"tie-out script raised {type(ex).__name__}: {ex}")]
    out = buf.getvalue()
    findings = []
    if code != 0:
        fails = [ln.strip() for ln in out.splitlines() if "[FAIL]" in ln]
        findings.append(Finding(CHECK_ID, HARD, SCRIPT, 0,
            f"numeric tie-outs FAILED (exit {code}) — a source moved or counts broke",
            "; ".join(fails[:4]) or f"exit {code}"))
    for ln in out.splitlines():
        s = ln.strip()
        if s.startswith(">> ERROR") or s.startswith(">> INCONSISTENCY"):
            findings.append(Finding(CHECK_ID, SOFT, SCRIPT, 0, s.lstrip(">").strip()))
    return findings
