"""corpus.py — load the canon once, with a line index and structural finders.

Read-only. Provides: raw text + 1-based line access; top-level §-section windows
(for §14/§15 window-scoping); and the enclosing-block finder used by the sync check's
proximity window. Pure stdlib.
"""
from __future__ import annotations

import re
from pathlib import Path

CANON_FILES = ("ANCHOR.md", "DRIFT.md", "STATUS.md")

# top-level ANCHOR section header: '## 14. ...' or subsection '## 16.1 ...'
_SEC = re.compile(r"^##\s*(\d+)(?:\.(\d+))?")
# a line that ends a prose block (blank, any markdown header, or a table separator row)
_BLOCK_BREAK = re.compile(r"^\s*$|^#{1,6}\s|^\s*\|[\s:\-|]+\|\s*$")


class Corpus:
    def __init__(self, root: Path):
        self.root = Path(root)
        self._lines: dict[str, list[str]] = {}
        self._text: dict[str, str] = {}
        for f in CANON_FILES:
            p = self.root / f
            t = p.read_text(encoding="utf-8", errors="replace") if p.exists() else ""
            self._text[f] = t
            self._lines[f] = t.splitlines()

    # ---- basic access (1-based line numbers throughout the public API) ----
    def text(self, f: str) -> str:
        return self._text.get(f, "")

    def lines(self, f: str) -> list[str]:
        return self._lines.get(f, [])

    def n_lines(self, f: str) -> int:
        return len(self._lines.get(f, []))

    def line(self, f: str, n: int) -> str:
        ls = self._lines.get(f, [])
        return ls[n - 1] if 1 <= n <= len(ls) else ""

    def exists(self, f: str) -> bool:
        return bool(self._text.get(f))

    # ---- structure ----
    def headers(self, f: str):
        """List of (major:int, minor:int|None, line:int) for each '## N[.M]' header."""
        out = []
        for i, ln in enumerate(self._lines.get(f, []), start=1):
            m = _SEC.match(ln)
            if m:
                minor = int(m.group(2)) if m.group(2) else None
                out.append((int(m.group(1)), minor, i))
        return out

    def section_window(self, f: str, major: int):
        """1-based [start, end) line window for top-level section `major` (the '## major.'
        header up to the next top-level '## k.' header). Returns None if not found."""
        hs = self.headers(f)
        start = next((ln for mj, mn, ln in hs if mj == major and mn is None), None)
        if start is None:
            return None
        end = next((ln for mj, mn, ln in hs if mn is None and ln > start),
                   self.n_lines(f) + 1)
        return (start, end)

    def section_of(self, f: str, line: int):
        """Major number of the top-level '## N.' section containing `line` (None if before §0)."""
        major = None
        for mj, mn, ln in self.headers(f):
            if mn is None:
                if ln <= line:
                    major = mj
                else:
                    break
        return major

    def enclosing_block(self, f: str, line: int, cap: int = 8):
        """1-based [start, end] of the prose block containing `line`, bounded by blank
        lines / headers / table separators, capped at +/- `cap` lines. Used as the sync
        check's proximity window (paragraph-bounded, not a raw line count)."""
        ls = self._lines.get(f, [])
        n = len(ls)
        if not (1 <= line <= n):
            return (line, line)
        is_break = lambda i: bool(_BLOCK_BREAK.match(ls[i - 1]))
        start = line
        while start > 1 and start > line - cap and not is_break(start - 1):
            start -= 1
        end = line
        while end < n and end < line + cap and not is_break(end + 1):
            end += 1
        return (start, end)

    # ---- evidence (file existence only in v1; no execution) ----
    def script_exists(self, rel_path: str) -> bool:
        return (self.root / rel_path).exists()

    def script_path(self, rel_path: str) -> Path:
        """Absolute path for a repo-relative script (used by the Phase-2 evidence runner)."""
        return self.root / rel_path
