"""cli.py — `python3 -m ptms {check,extract,selftest}`.

`check`   : read-only over the canon; exit 1 iff >=1 HARD finding.
`extract` : refresh the registry's AUTO fields (the only sanctioned writer; --write to persist).
`selftest`: run the fixture self-test (recall + precision of the checker's own logic).
"""
from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path

from .corpus import Corpus
from .registry import load_registry, save_registry, index_by_id
from .findings import HARD
from .report import render
from . import checks as checks_pkg


def find_root(explicit: str | None) -> Path:
    if explicit:
        return Path(explicit).resolve()
    here = Path.cwd()
    for cand in (here, *here.parents):
        if (cand / "ANCHOR.md").exists():
            return cand
    # fall back: repo root is two levels up from this file (.../octahedrons/ai_methodology/ptms)
    return Path(__file__).resolve().parents[2]


def default_registry(root: Path) -> Path:
    return root / "ai_methodology" / "ptms" / "claims.jsonl"


def default_manifest(root: Path) -> Path:
    return root / "ai_methodology" / "ptms" / "evidence_manifest.json"


def cmd_check(args) -> int:
    root = find_root(args.root)
    corpus = Corpus(root)
    registry = load_registry(args.registry or default_registry(root))
    only = set(args.only.split(",")) if args.only else None
    findings, ran = [], []
    for m in checks_pkg.discover():
        if only and m.CHECK_ID not in only:
            continue
        ran.append(m.CHECK_ID)
        findings.extend(m.run(corpus, registry))
    meta = {"anchor_lines": corpus.n_lines("ANCHOR.md"),
            "drift_lines": corpus.n_lines("DRIFT.md"),
            "registry_n": len(registry), "checks_run": ran}
    print(render(findings, meta))
    return 1 if any(f.severity == HARD for f in findings) else 0


def cmd_extract(args) -> int:
    from .extract import extract, merge, summary
    root = find_root(args.root)
    corpus = Corpus(root)
    reg_path = Path(args.registry or default_registry(root))
    existing = load_registry(reg_path)
    fresh = extract(corpus)
    merged, stats = merge(existing, fresh)
    print(summary(corpus, merged))
    print(f"\n  merge: {stats['new']} new, {stats['kept']} updated (human fields kept), "
          f"{stats['stale']} stale (kept; flagged by registry_drift)")
    if args.write:
        save_registry(reg_path, merged)
        print(f"  wrote {len(merged)} records -> {reg_path}")
    else:
        print("  (dry run; pass --write to persist)")
    return 0


def cmd_manifest(args) -> int:
    from . import manifest as M
    root = find_root(args.root)
    corpus = Corpus(root)
    registry = load_registry(args.registry or default_registry(root))
    mpath = Path(args.manifest or default_manifest(root))
    prior = M.load_manifest(mpath)
    data = M.build(corpus, registry, scope=("all" if args.all else "cited"),
                   force=args.force, timeout=args.timeout, fast_threshold=args.fast_threshold,
                   prior=prior)
    print(M.summary(data))
    if args.write:
        M.save_manifest(mpath, data)
        print(f"  wrote {len(data['scripts'])} records -> {mpath}")
    else:
        print("  (dry run; pass --write to persist)")
    return 0


def cmd_verify(args) -> int:
    from . import manifest as M
    from .checks import evidence_regression as ev
    root = find_root(args.root)
    corpus = Corpus(root)
    registry = load_registry(args.registry or default_registry(root))
    mpath = Path(args.manifest or default_manifest(root))
    prior = M.load_manifest(mpath)
    # refresh the manifest in-memory (incremental) so verify reflects current disk state
    data = M.build(corpus, registry, scope=("all" if args.all else "cited"),
                   force=args.force, timeout=args.timeout, fast_threshold=args.fast_threshold,
                   prior=prior)
    if not args.no_write:
        M.save_manifest(mpath, data)
    findings = ev.run(corpus, registry, data, prior=prior, run_all=args.all)
    meta = {"anchor_lines": corpus.n_lines("ANCHOR.md"), "drift_lines": corpus.n_lines("DRIFT.md"),
            "registry_n": len(registry), "checks_run": [ev.CHECK_ID]}
    print(render(findings, meta))
    return 1 if any(f.severity == HARD for f in findings) else 0


def cmd_graph(args) -> int:
    from . import graph as G
    root = find_root(args.root)
    registry = load_registry(args.registry or default_registry(root))
    g = G.build(registry)
    out = G.to_dot(g, scope=args.scope, top=args.top) if args.format == "dot" else G.render_md(g, top=args.top)
    print(out)
    return 0                                       # analysis, not a gate — always exit 0


def cmd_selftest(args) -> int:
    tests_dir = Path(__file__).resolve().parent / "tests"
    rc = 0
    for name in ("test_selftest.py", "test_verify.py", "test_graph.py"):
        test = tests_dir / name
        if not test.exists():
            print(f"selftest not found at {test}"); rc = 1; continue
        print(f"--- {name} ---")
        try:
            runpy.run_path(str(test), run_name="__main__")
        except SystemExit as e:
            rc = rc or int(e.code or 0)
    return rc


def cmd_locate(args) -> int:
    """Retrieve via the index instead of grepping the 1.8MB monolith: filter claims.jsonl by
    id/keyword/--tier/--status/--kind and print id, status/tier, file:line, and dep/script counts."""
    root = find_root(args.root)
    registry = load_registry(args.registry or default_registry(root))
    q = (args.query or "").lower()
    hits = []
    for c in registry:
        if args.tier and (c.tier or "").lower() != args.tier.lower():
            continue
        if args.status and c.status.lower() != args.status.lower():
            continue
        if args.kind and c.kind != args.kind:
            continue
        hay = (c.id + " " + c.location.get("anchor_text", "")).lower()
        if q and q not in hay:
            continue
        hits.append(c)
    hits.sort(key=lambda c: c.id)
    for c in hits[:args.limit]:
        loc = c.location
        print(f"{c.id}\t[{c.status}/{c.tier}]\t{loc.get('file','')}:{loc.get('line',0)}\t"
              f"refs={len(c.cited_refs)} scripts={len(c.cited_scripts)}\t"
              f"{loc.get('anchor_text','')[:72]}")
    shown = min(len(hits), args.limit)
    tail = "" if shown == len(hits) else f"; showing {shown} (raise --limit)"
    print(f"# {len(hits)} match(es){tail}")
    return 0


def cmd_show(args) -> int:
    """Print one claim's metadata (status/tier/deps/scripts) + its exact ANCHOR/DRIFT block —
    so an LLM reads only the block it needs, never the whole file."""
    root = find_root(args.root)
    corpus = Corpus(root)
    registry = load_registry(args.registry or default_registry(root))
    c = index_by_id(registry).get(args.id)
    if c is None:
        near = [x.id for x in registry if args.id.lower() in x.id.lower()][:8]
        print(f"# no claim id '{args.id}'" + (f"; did you mean: {', '.join(near)}" if near else ""))
        return 1
    loc = c.location
    f, ln = loc.get("file", "ANCHOR.md"), loc.get("line", 0)
    print(f"# {c.id}  [status={c.status}, tier={c.tier}, kind={c.kind}]  {f}:{ln}")
    for label, vals in (("cited_refs", c.cited_refs), ("cited_scripts", c.cited_scripts),
                        ("supersedes", c.supersedes), ("superseded_by", c.superseded_by)):
        if vals:
            print(f"# {label}: {', '.join(vals)}")
    if ln and corpus.exists(f):
        s, e = corpus.enclosing_block(f, ln, cap=args.context)
        print(f"# block {f}:{s}-{e}\n")
        for i, line in enumerate(corpus.lines(f)[s - 1:e], s):
            print(f"{i}\t{line}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="ptms", description="canon-linter for the TCH framework")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("check", help="run consistency checks, emit report")
    c.add_argument("--root"); c.add_argument("--registry")
    c.add_argument("--only", help="comma-separated check ids")
    c.set_defaults(fn=cmd_check)

    e = sub.add_parser("extract", help="refresh the sidecar registry")
    e.add_argument("--root"); e.add_argument("--registry")
    e.add_argument("--write", action="store_true")
    e.set_defaults(fn=cmd_extract)

    def _runner_opts(parser):
        parser.add_argument("--root"); parser.add_argument("--registry"); parser.add_argument("--manifest")
        parser.add_argument("--cited", action="store_true", help="only scripts cited by claims (default)")
        parser.add_argument("--all", action="store_true", help="all python_code scripts, not just cited")
        parser.add_argument("--force", action="store_true", help="ignore the mtime cache; re-run all")
        parser.add_argument("--timeout", type=float, default=30.0)
        parser.add_argument("--fast-threshold", dest="fast_threshold", type=float, default=10.0)

    m = sub.add_parser("manifest", help="time + classify the cited scripts (evidence manifest)")
    _runner_opts(m); m.add_argument("--write", action="store_true")
    m.set_defaults(fn=cmd_manifest)

    v = sub.add_parser("verify", help="re-run cited self-asserting scripts; flag evidence regressions")
    _runner_opts(v); v.add_argument("--no-write", action="store_true", help="do not persist the refreshed manifest")
    v.set_defaults(fn=cmd_verify)

    gp = sub.add_parser("graph", help="claim dependency graph: contamination + frontier analyses, DOT export")
    gp.add_argument("--root"); gp.add_argument("--registry")
    gp.add_argument("--format", choices=["md", "dot"], default="md")
    gp.add_argument("--scope", choices=["full", "contamination", "frontier"], default="contamination")
    gp.add_argument("--top", type=int, default=15)
    gp.set_defaults(fn=cmd_graph)

    s = sub.add_parser("selftest", help="run the fixture self-test")
    s.set_defaults(fn=cmd_selftest)

    lo = sub.add_parser("locate", help="find claims via the index (id/keyword/--tier/--status/--kind)")
    lo.add_argument("query", nargs="?", default="", help="substring of id or anchor text")
    lo.add_argument("--root"); lo.add_argument("--registry")
    lo.add_argument("--tier"); lo.add_argument("--status"); lo.add_argument("--kind")
    lo.add_argument("--limit", type=int, default=40)
    lo.set_defaults(fn=cmd_locate)

    sh = sub.add_parser("show", help="print a claim's metadata + its exact ANCHOR/DRIFT block")
    sh.add_argument("id", help="registry id, e.g. 'anchor:§15:126' or 'drift:M9'")
    sh.add_argument("--root"); sh.add_argument("--registry")
    sh.add_argument("--context", type=int, default=8, help="enclosing-block cap in lines")
    sh.set_defaults(fn=cmd_show)

    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
