"""PTMS — a truth-maintenance canon-linter for the TCH ("Holographic Circlette") framework.

Inward, single-context, JTMS-style: it ENFORCES CONSISTENCY, NOT TRUTH. Every justification
edge is untrusted until bound to checkable evidence (cited script file resolves, cited
section/item resolves, a retired claim's signature is flagged at every survival site).

v1 (this package): mechanical consistency checks over the canon (ANCHOR/DRIFT/STATUS) plus a
sidecar claim registry (claims.jsonl). NO script execution yet (the evidence-runner is Phase 2).
On-demand, itemised report, no aggregate score. The journal screener is a separate tool.
"""
__version__ = "0.1.0"
