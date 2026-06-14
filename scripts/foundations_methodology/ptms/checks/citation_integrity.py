"""checks/citation_integrity.py — evidence/cross-ref integrity (no execution in v1).

HARD: a cited `python_code/X.py` that is absent on disk; a cross-ref (`§15 item N`, `DRIFT X`,
`§n`) resolving to no registry record; a duplicate §15 item number or DRIFT code. Section refs
resolve tolerantly to their top-level section (so `§6.6` is fine while `§99` is dangling).
"""
from __future__ import annotations

import re
from collections import Counter

from ..findings import Finding, HARD, SOFT

CHECK_ID = "citation_integrity"
ITEM = re.compile(r"^§15 item (\d+)$")
DRIFTREF = re.compile(r"^DRIFT ([A-Z]\d+)$")
SECREF = re.compile(r"^§(\d+)")


def run(corpus, registry) -> list:
    out = []
    item_nums = {int(c.id.rsplit(":", 1)[1]) for c in registry
                 if c.kind == "§15-item" and c.id.rsplit(":", 1)[1].isdigit()}
    drift_codes = {c.id.split(":", 1)[1] for c in registry if c.kind == "drift-entry"}
    sec_majors = {int(tail) for c in registry if c.kind == "anchor-section" and c.id.startswith("anchor:sec:")
                  for tail in [c.id.split(":")[-1].split(".")[0]] if tail.isdigit()}

    for c in registry:
        f, ln = c.location.get("file", ""), c.location.get("line", 0)
        for sc in c.cited_scripts:
            if not corpus.script_exists(sc):
                out.append(Finding(CHECK_ID, HARD, f, ln, f"cited script missing on disk: {sc}", c.id))
        for ref in c.cited_refs:
            # item/section dangling refs are SOFT: prose mentions (proposals, declined items)
            # are noisy and shouldn't gate. DRIFT-code refs and missing scripts are crisp -> HARD.
            m = ITEM.match(ref)
            if m:
                if int(m.group(1)) not in item_nums:
                    out.append(Finding(CHECK_ID, SOFT, f, ln, f"reference to non-existent §15 item: '{ref}'", c.id))
                continue
            m = DRIFTREF.match(ref)
            if m:
                if m.group(1) not in drift_codes:
                    out.append(Finding(CHECK_ID, HARD, f, ln, f"dangling cross-ref '{ref}' (no such DRIFT entry)", c.id))
                continue
            m = SECREF.match(ref)
            if m and int(m.group(1)) not in sec_majors:
                out.append(Finding(CHECK_ID, SOFT, f, ln, f"reference to non-existent section: '{ref}'", c.id))

    for num, k in Counter(c.id.rsplit(":", 1)[1] for c in registry if c.kind == "§15-item").items():
        if k > 1:
            out.append(Finding(CHECK_ID, HARD, "ANCHOR.md", 0, f"duplicate §15 item number {num} ({k}x)"))
    for code, k in Counter(c.id.split(":", 1)[1] for c in registry if c.kind == "drift-entry").items():
        if k > 1:
            out.append(Finding(CHECK_ID, HARD, "DRIFT.md", 0, f"duplicate DRIFT code {code} ({k}x)"))
    return out
