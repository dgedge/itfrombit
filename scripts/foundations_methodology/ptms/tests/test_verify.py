#!/usr/bin/env python3
"""Phase-2 self-test — the evidence-runner eats its own dog food.

Builds a manifest over 4 fixture scripts with KNOWN outcomes and asserts the exact
evidence_regression finding set (recall + precision). exit 0 iff the runner's logic holds.

  pass_assert (assert, exit 0)  -> silent (precision)
  fail_assert (assert, exit 1)  -> 1 HARD (evidence regression)
  pure_print  (none,   exit 0)  -> SOFT (ran but asserts nothing)
  slow_stub   (timed out)       -> SOFT with no baseline; HARD against a seeded prior-passing baseline
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))  # .../ai_methodology

from ptms.corpus import Corpus
from ptms.registry import Claim
from ptms import manifest as M
from ptms.checks import evidence_regression as ev

FIX = pathlib.Path(__file__).resolve().parent / "fixtures"

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c:
        fails.append(m)


corpus = Corpus(FIX)
registry = [
    Claim(id="fixture:pass", kind="drift-entry", status="superseded",
          cited_scripts=["scripts/pass_assert.py"], location={"file": "DRIFT.md", "line": 1}),
    Claim(id="fixture:fail", kind="drift-entry", status="superseded",
          cited_scripts=["scripts/fail_assert.py"], location={"file": "DRIFT.md", "line": 2}),
    Claim(id="fixture:print", kind="drift-entry", status="superseded",
          cited_scripts=["scripts/pure_print.py"], location={"file": "DRIFT.md", "line": 3}),
    Claim(id="fixture:slow", kind="drift-entry", status="superseded",
          cited_scripts=["scripts/slow_stub.py"], location={"file": "DRIFT.md", "line": 4}),
]

data = M.build(corpus, registry, scope="cited", timeout=0.5, fast_threshold=0.5)


def has(findings, severity, needle):
    return [f for f in findings if f.severity == severity and needle in f.message]


# --- pass 1: no prior baseline ---
f1 = ev.run(corpus, registry, data, prior=None)
hard1 = [f for f in f1 if f.severity == "hard"]
check(len(has(f1, "hard", "fail_assert.py")) == 1, "fail_assert -> HARD evidence regression")
check(len(hard1) == 1, f"exactly 1 HARD with no baseline (got {len(hard1)})")
check(len(has(f1, "soft", "pure_print.py")) == 1, "pure_print -> SOFT (ran but asserts nothing)")
check(len(has(f1, "soft", "slow_stub.py")) == 1, "slow_stub (no baseline) -> SOFT (not verified)")
check(len(has(f1, "hard", "pass_assert.py")) == 0 and len(has(f1, "soft", "pass_assert.py")) == 0,
      "pass_assert -> silent (precision: a passing cited script yields no finding)")

# --- pass 2: seeded prior-passing baseline -> slow_stub timeout is now a regression ---
prior = {"scripts": {"scripts/slow_stub.py": {"exit_code": 0, "timed_out": False, "mtime": -1.0}}}
f2 = ev.run(corpus, registry, data, prior=prior)
check(len(has(f2, "hard", "slow_stub.py")) == 1,
      "slow_stub vs prior-passing baseline -> HARD (regression-to-timeout)")

print()
if fails:
    print(f"SELFTEST(verify): {len(fails)} FAILED")
    sys.exit(1)
print("SELFTEST(verify): exit 0 — evidence-runner logic verified (recall + precision).")
