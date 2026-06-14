"""report.py — itemised markdown dossier. NO aggregate score (by design)."""
from __future__ import annotations

from .findings import HARD, SOFT, order

COLLAPSE_SOFT_OVER = 8          # soft checks with more findings than this are collapsed


def _group(findings: list):
    g: dict = {}
    for f in findings:
        g.setdefault((f.severity, f.check_id), []).append(f)
    return g


def render(findings: list, meta: dict) -> str:
    findings = order(findings)
    g = _group(findings)
    hard_n = sum(1 for f in findings if f.severity == HARD)
    soft_n = sum(1 for f in findings if f.severity == SOFT)

    L = []
    L.append(f"# PTMS canon-linter report — {meta.get('date', 'on-demand')}")
    L.append("")
    L.append(f"sources: ANCHOR.md ({meta.get('anchor_lines', '?')}L), "
             f"DRIFT.md ({meta.get('drift_lines', '?')}L); registry {meta.get('registry_n', 0)} records")
    if meta.get("checks_run"):
        L.append(f"checks run: {', '.join(meta['checks_run'])}")
    L.append("")

    def emit(severity, tag):
        keys = sorted(k for k in g if k[0] == severity)
        if not keys:
            return
        for _, check_id in keys:
            fs = g[(severity, check_id)]
            head = f"## [{tag}] {check_id} ({len(fs)})"
            if severity == SOFT and len(fs) > COLLAPSE_SOFT_OVER:
                L.append(head + "  ▸ collapsed")
                for f in fs[:3]:
                    L.append(f"- {f.file}:{f.line} — {f.message}")
                L.append(f"- … and {len(fs) - 3} more (run `--only {check_id}` to expand)")
            else:
                L.append(head)
                for f in fs:
                    loc = f"{f.file}:{f.line}" if f.line else (f.file or "—")
                    ev = f"  · {f.evidence}" if f.evidence else ""
                    L.append(f"- {loc} — {f.message}{ev}")
            L.append("")

    emit(HARD, "HARD")
    emit(SOFT, "SOFT")
    if not findings:
        L.append("✓ no findings.")
        L.append("")
    L.append("---")
    L.append(f"HARD findings: {hard_n}    SOFT findings: {soft_n}    (no aggregate score by design)")
    return "\n".join(L)
