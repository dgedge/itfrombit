"""checks/evidence_regression.py — Phase 2: bind cited claims to live exit-0 evidence.

Consumes the evidence manifest (already-measured exit codes) + the registry. The Phase-2
analogue of the sync killer feature: it catches the canon citing a proof script that has since
broken. A passing, assertion-bearing cited script is silent (precision). NOT run under `check`
(read-only); only under `ptms verify`.
"""
from __future__ import annotations

from ..findings import Finding, HARD, SOFT

CHECK_ID = "evidence_regression"


def run(corpus, registry, manifest, prior=None, run_all=False) -> list:
    scripts = manifest.get("scripts", {})
    prior_scripts = (prior or {}).get("scripts", {})
    findings = []
    cited = set()
    for c in registry:
        f, l = c.location.get("file", ""), c.location.get("line", 0)
        for sc in c.cited_scripts:
            cited.add(sc)
            rec = scripts.get(sc)
            if rec is None:
                findings.append(Finding(CHECK_ID, SOFT, f, l,
                    f"{sc} cited by {c.id} not in evidence manifest (run `ptms manifest --write`)", c.id))
                continue
            if rec.get("excluded"):
                findings.append(Finding(CHECK_ID, SOFT, f, l,
                    f"{sc} cited by {c.id} is hard-excluded ({rec.get('reason')}) — never auto-run", c.id))
                continue
            if rec.get("timed_out"):
                was_ok = (prior_scripts.get(sc, {}).get("exit_code") == 0
                          and not prior_scripts.get(sc, {}).get("timed_out"))
                if was_ok:
                    findings.append(Finding(CHECK_ID, HARD, f, l,
                        f"evidence regression: {sc} cited by {c.id} no longer completes within the timeout "
                        f"(was passing)", c.id))
                else:
                    findings.append(Finding(CHECK_ID, SOFT, f, l,
                        f"{sc} cited by {c.id} exceeds the timeout — not verified this pass", c.id))
                continue
            code = rec.get("exit_code")
            if code not in (0, None):
                findings.append(Finding(CHECK_ID, HARD, f, l,
                    f"evidence regression: {sc} cited by {c.id} now exits {code} (expected exit 0)",
                    (rec.get("stderr_tail") or "").strip()[-200:]))
                continue
            if rec.get("assertion_class") == "none":
                findings.append(Finding(CHECK_ID, SOFT, f, l,
                    f"{sc} cited by {c.id} ran but asserts nothing (no check/assert/sys.exit) — "
                    f"exit 0 is not evidence", c.id))
                continue
            # exit 0 + assertion-bearing (check/assert): the happy path — silent.

    if run_all:
        for sc, rec in scripts.items():
            if sc in cited or rec.get("excluded"):
                continue
            if rec.get("exit_code") == 0 and rec.get("assertion_class") in ("check", "assert"):
                findings.append(Finding(CHECK_ID, SOFT, sc, 0,
                    f"{sc} passes (exit 0) but is uncited — no claim binds to it"))
    return findings
