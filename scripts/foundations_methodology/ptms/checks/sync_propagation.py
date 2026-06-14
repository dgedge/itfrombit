"""checks/sync_propagation.py — the killer feature.

For each retired claim (status retired/superseded/withdrawn) with signature_strings, grep ANCHOR
for each signature. For every occurrence, decide severity:
  - locally flagged  (a retirement marker or the DRIFT code within +/-WINDOW lines)  -> no finding;
  - section-aware    (the enclosing top-level §-section references the DRIFT code, but this
                      restatement is not locally marked)                              -> SOFT (polish);
  - un-propagated    (the enclosing section has NO reference to the retraction at all) -> HARD.

This two-tier severity matches ANCHOR's structure: a long retired discussion notes the
retraction once at its head and restates the value in its body; only a survival in a section the
retraction never reached is a true propagation failure. DRIFT is never scanned (ANCHOR only).
"""
from __future__ import annotations

import re

from ..findings import Finding, HARD, SOFT

CHECK_ID = "sync_propagation"
WINDOW = 10                                   # flat +/- line window for a local marker
RETIRE_MARKERS = ("RETIRED", "WITHDRAWN", "SUPERSEDED", "FLAGGED", "REMOVED FROM LOCKED")
RETIRE_STATUSES = ("retired", "superseded", "withdrawn")


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def run(corpus, registry) -> list:
    A = "ANCHOR.md"
    alines = corpus.lines(A)
    norm = [_norm(l) for l in alines]
    n = len(alines)
    # cache section text (upper) per major for the awareness test
    sec_cache: dict = {}

    def section_aware(line: int, code: str) -> bool:
        # A retired value's enclosing §-section is "aware" of the retraction if it references
        # the DRIFT code OR carries ANY retirement marker. HARD is then reserved for a value
        # surviving in a section with ZERO retirement vocabulary (the true un-propagated case);
        # a restatement inside an aware section is SOFT (polish: not locally marked).
        major = corpus.section_of(A, line)
        if major is None:
            return False
        if major not in sec_cache:
            w = corpus.section_window(A, major)
            sec_cache[major] = " ".join(alines[w[0] - 1:w[1] - 1]).upper() if w else ""
        sect = sec_cache[major]
        return code in sect or any(m in sect for m in RETIRE_MARKERS)

    findings = []
    for c in registry:
        if c.kind != "drift-entry" or c.status not in RETIRE_STATUSES:
            continue
        code = c.id.split(":", 1)[1].upper()      # e.g. "M9"
        for sig in c.signature_strings:
            nsig = _norm(sig)
            if len(nsig) < 4:
                continue
            for i in range(1, n + 1):
                if nsig not in norm[i - 1]:
                    continue
                lo, hi = max(1, i - WINDOW), min(n, i + WINDOW)
                window = " ".join(alines[lo - 1:hi]).upper()
                if any(m in window for m in RETIRE_MARKERS) or code in window:
                    continue                        # locally flagged -> OK
                if section_aware(i, code):
                    findings.append(Finding(CHECK_ID, SOFT, A, i,
                        f"restatement of retired {c.id} signature '{sig}' not locally flagged "
                        f"(its §-section references {code}, but no marker within {WINDOW} lines)"))
                else:
                    findings.append(Finding(CHECK_ID, HARD, A, i,
                        f"un-propagated survival of retired {c.id} signature '{sig}'",
                        f"enclosing §-section has no reference to {code}"))
    return findings
