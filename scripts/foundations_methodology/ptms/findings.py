"""findings.py — the Finding record and ordering helpers."""
from __future__ import annotations

from dataclasses import dataclass, field

HARD, SOFT = "hard", "soft"


@dataclass
class Finding:
    check_id: str
    severity: str          # "hard" | "soft"
    file: str              # "ANCHOR.md", "DRIFT.md", "registry", "" ...
    line: int              # 1-based; 0 if not line-located
    message: str
    evidence: str = ""

    def _key(self):
        return (0 if self.severity == HARD else 1, self.check_id, self.file, self.line, self.message)


def dedup(findings: list) -> list:
    seen, out = set(), []
    for f in findings:
        k = (f.check_id, f.severity, f.file, f.line, f.message)
        if k not in seen:
            seen.add(k)
            out.append(f)
    return out


def order(findings: list) -> list:
    return sorted(dedup(findings), key=lambda f: f._key())
