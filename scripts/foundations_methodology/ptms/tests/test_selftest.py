#!/usr/bin/env python3
"""Self-test — the linter eats its own dog food.

Runs the core checks against a fixture canon (tests/fixtures/ANCHOR.md) with KNOWN seeded
faults, and asserts the finding set equals the expected set EXACTLY — recall (it finds every
seeded fault) AND precision (it does NOT flag the correctly-propagated survivals). exit 0 iff
the checker's own logic is verified.

Seeded faults:
  - SIGALPHA: an un-propagated retired-claim survival (no marker)        -> sync HARD
  - SIGBETA / SIGGAMMA: retired survivals WITH a marker                  -> NO finding (precision)
  - drift:K9: multipart_count=2 but only 1 signature                     -> multipart HARD
  - drift:M9: multipart_count=2 with 2 signatures                        -> NO finding (precision)
  - anchor:locked: still [Tier: Locked]/derived, superseded by M9        -> tier_contradiction HARD
  - python_code/missing_stub.py cited but absent                         -> citation HARD
  - §15 item 999 referenced                                              -> citation SOFT
  - untagged §14 row                                                     -> tier_coverage SOFT
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))  # .../ai_methodology

from ptms.corpus import Corpus
from ptms.registry import Claim
from ptms.checks import (sync_propagation, multipart, tier_coverage,
                         tier_contradiction, citation_integrity)

FIX = pathlib.Path(__file__).resolve().parent / "fixtures"

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c:
        fails.append(m)


def line_of(corpus, sub):
    for i, l in enumerate(corpus.lines("ANCHOR.md"), 1):
        if sub in l:
            return i
    return 0


corpus = Corpus(FIX)
A = "ANCHOR.md"
registry = [
    Claim(id="anchor:§14:row:001", kind="§14-row", tier=None,
          location={"file": A, "line": line_of(corpus, "untagged fixture row"),
                    "anchor_text": "untagged fixture row"}),
    Claim(id="anchor:§15:1", kind="§15-item", status="live",
          location={"file": A, "line": line_of(corpus, "Fixture item one"),
                    "anchor_text": "Fixture item one"},
          cited_scripts=["python_code/missing_stub.py"], cited_refs=["§15 item 999"]),
    Claim(id="anchor:locked", kind="anchor-section", tier="Locked",
          location={"file": A, "line": line_of(corpus, "[Tier: Locked] and derived"),
                    "anchor_text": "This fixture result"}),
    Claim(id="drift:M9", kind="drift-entry", status="superseded",
          signature_strings=["SIGALPHA", "SIGBETA"], multipart_count=2,
          supersedes=["anchor:locked"],
          location={"file": "DRIFT.md", "line": 1, "anchor_text": ""}),
    Claim(id="drift:K9", kind="drift-entry", status="superseded",
          signature_strings=["SIGGAMMA"], multipart_count=2,
          location={"file": "DRIFT.md", "line": 2, "anchor_text": ""}),
]

found = []
for mod in (sync_propagation, multipart, tier_coverage, tier_contradiction, citation_integrity):
    found += mod.run(corpus, registry)


def hits(check_id, severity, needle=""):
    return [f for f in found if f.check_id == check_id and f.severity == severity and needle in f.message]


# --- recall: every seeded fault is found ---
check(len(hits("sync_propagation", "hard", "SIGALPHA")) == 1, "sync flags the un-propagated SIGALPHA survival")
check(len(hits("multipart", "hard", "drift:K9")) == 1, "multipart flags drift:K9 (K=2, 1 signature)")
check(len(hits("tier_coverage", "soft")) == 1, "tier_coverage flags the untagged §14 row")
check(len(hits("tier_contradiction", "hard", "anchor:locked")) == 1, "tier_contradiction flags the Locked survivor")
check(len(hits("citation_integrity", "hard", "missing_stub")) == 1, "citation flags the missing script")
check(len(hits("citation_integrity", "soft", "999")) == 1, "citation flags the dangling item 999 (soft)")

# --- precision: correctly-propagated things are NOT flagged ---
check(len(hits("sync_propagation", "hard", "SIGBETA")) == 0, "sync does NOT flag the correctly-retired SIGBETA")
check(len(hits("sync_propagation", "hard", "SIGGAMMA")) == 0, "sync does NOT flag the correctly-retired SIGGAMMA")
check(len(hits("multipart", "hard", "drift:M9")) == 0, "multipart does NOT flag drift:M9 (K matches signatures)")

# --- exact totals (recall + precision together) ---
hard = [f for f in found if f.severity == "hard"]
soft = [f for f in found if f.severity == "soft"]
check(len(hard) == 4, f"exactly 4 HARD findings (got {len(hard)})")
check(len(soft) == 2, f"exactly 2 SOFT findings (got {len(soft)})")

print()
if fails:
    print(f"SELFTEST: {len(fails)} FAILED")
    sys.exit(1)
print("SELFTEST: exit 0 — checker logic verified (recall + precision on the fixture canon).")
