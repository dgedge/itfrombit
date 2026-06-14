"""extract.py — semi-automatic registry bootstrapper (read-only over the canon).

Builds the AUTO fields of claims.jsonl from ANCHOR/DRIFT with ZERO prose edits. Human-owned
fields (signature_strings, multipart_count, supersedes) are PRESERVED across re-extraction by
`id` (see registry.Claim.merge_auto). The extractor never guesses signature_strings; it leaves
them empty and lists the claims that need a human pass in a worklist.
"""
from __future__ import annotations

import re

from .corpus import Corpus
from .registry import Claim, SCHEMA_VERSION, index_by_id

SCRIPT_CITE = re.compile(r"python_code/([A-Za-z0-9_]+)\.py")          # path-prefixed ONLY
TIER        = re.compile(r"\[Tier:\s*([^\]]+)\]")
REF_ITEM    = re.compile(r"§\s*15\s*item\s*(\d+)", re.I)
REF_SEC     = re.compile(r"§\s*(\d+\.\d+(?:\.\d+)?)")                  # §6.6, §5.10.1 (subsections)
REF_DRIFT   = re.compile(r"\bDRIFT\s+([A-Z]\d+)\b")
ITEM15      = re.compile(r"^(\d+)\.\s+")                               # §15 item header (col 0)
DRIFT_HDR   = re.compile(r"^### ([A-Z])(\d+)\.\s")
STATUS_FLD  = re.compile(r"\*\*Status of older formulation\*\*", re.I)
ROW_SEP     = re.compile(r"^\s*\|[\s:\-|]+\|\s*$")
BOLD_STATUS = re.compile(r"\*\*\s*(superseded|flagged|withdrawn|retired)\s*\*\*", re.I)

RETIRE_WORDS = ("RETIRED", "WITHDRAWN", "SUPERSEDED", "FLAGGED")


def _scripts(text: str) -> list:
    return sorted({f"python_code/{m}.py" for m in SCRIPT_CITE.findall(text)})


def _refs(text: str) -> list:
    refs = [f"§15 item {n}" for n in REF_ITEM.findall(text)]
    refs += [f"§{s}" for s in REF_SEC.findall(text)]
    refs += [f"DRIFT {d}" for d in REF_DRIFT.findall(text)]
    return sorted(set(refs))


def _map_tier(block: str):
    m = TIER.search(block)
    if not m:
        return None
    r = m.group(1).lower()
    found = [t for t, w in (("Locked", "lock"), ("Proposition", "propos"), ("Open", "open")) if w in r]
    if "contested" in r and "Open" not in found:
        found.append("Open")
    if len(found) == 1:
        return found[0]
    if len(found) >= 2:
        return "compound"
    return None


def _anchor_status(block: str) -> str:
    up = block.upper()
    for w in RETIRE_WORDS:
        if w in up:
            return w.lower()
    return "live"


def _drift_status(block: str) -> str:
    m = STATUS_FLD.search(block)
    region = block[m.start():] if m else block
    up = region.upper()
    for w in ("WITHDRAWN", "RETIRED", "SUPERSEDED", "FLAGGED"):
        if w in up:
            return w.lower()
    if "OPEN" in up:
        return "open"
    if "PROPOSITION" in up:
        return "proposition"
    return "live"


def _drift_multipart_suggest(block: str) -> int:
    """AUTO-SUGGEST K = distinct bold status tokens in the Status field (human confirms)."""
    m = STATUS_FLD.search(block)
    region = block[m.start():m.start() + 800] if m else ""
    return len({t.lower() for t in BOLD_STATUS.findall(region)})


def extract(corpus: Corpus) -> list:
    A, D = "ANCHOR.md", "DRIFT.md"
    alines, dlines = corpus.lines(A), corpus.lines(D)
    claims: list = []

    # ANCHOR sections (cross-ref resolution targets)
    for mj, mn, ln in corpus.headers(A):
        sid = f"anchor:sec:{mj}" if mn is None else f"anchor:sec:{mj}.{mn}"
        claims.append(Claim(id=sid, kind="anchor-section",
                            location={"file": A, "line": ln, "anchor_text": corpus.line(A, ln)[:60]}))

    # §15 items (window-scoped)
    win = corpus.section_window(A, 15)
    if win:
        s, e = win
        item_lines = [(int(ITEM15.match(alines[i - 1]).group(1)), i)
                      for i in range(s, e) if ITEM15.match(alines[i - 1] or "")]
        for idx, (num, ln) in enumerate(item_lines):
            blk_end = item_lines[idx + 1][1] if idx + 1 < len(item_lines) else e
            block = "\n".join(alines[ln - 1:blk_end - 1])
            claims.append(Claim(id=f"anchor:§15:{num}", kind="§15-item",
                location={"file": A, "line": ln, "anchor_text": alines[ln - 1][:70]},
                tier=_map_tier(block), status=_anchor_status(block),
                cited_scripts=_scripts(block), cited_refs=_refs(block)))

    # §14 rows (window-scoped)
    win = corpus.section_window(A, 14)
    if win:
        s, e = win
        rown = 0
        for i in range(s, e):
            ln = alines[i - 1]
            if not ln.lstrip().startswith("|") or ROW_SEP.match(ln):
                continue
            if "Where to look" in ln and "Signature" in ln:   # the column-header row
                continue
            rown += 1
            claims.append(Claim(id=f"anchor:§14:row:{rown:03d}", kind="§14-row",
                location={"file": A, "line": i, "anchor_text": ln.strip()[:70]},
                tier=_map_tier(ln), status=_anchor_status(ln),
                cited_scripts=_scripts(ln), cited_refs=_refs(ln)))

    # DRIFT coded entries
    hdrs = [(i, DRIFT_HDR.match(dlines[i - 1])) for i in range(1, len(dlines) + 1)]
    hdrs = [(i, m) for i, m in hdrs if m]
    for idx, (ln, m) in enumerate(hdrs):
        end = hdrs[idx + 1][0] if idx + 1 < len(hdrs) else len(dlines) + 1
        block = "\n".join(dlines[ln - 1:end - 1])
        k_hint = _drift_multipart_suggest(block)
        claims.append(Claim(id=f"drift:{m.group(1)}{m.group(2)}", kind="drift-entry",
            location={"file": D, "line": ln, "anchor_text": dlines[ln - 1][:70]},
            status=_drift_status(block),
            cited_scripts=_scripts(block), cited_refs=_refs(block),
            auto_hint=(f"AUTO-SUGGEST multipart_count={k_hint}; " if k_hint else "")
                      + "needs human curation: signature_strings + supersedes"))
    return claims


def merge(existing: list, fresh: list):
    """Preserve human fields from `existing` by id; take auto fields from `fresh`. Existing
    claims absent from `fresh` are kept (stale) so annotations are never lost; the
    registry_drift check flags them."""
    ex = index_by_id(existing)
    fresh_ids = {c.id for c in fresh}
    out, new, kept = [], 0, 0
    for c in fresh:
        if c.id in ex:
            out.append(ex[c.id].merge_auto(c)); kept += 1
        else:
            out.append(c); new += 1
    stale = [c for c in existing if c.id not in fresh_ids]
    return out, {"new": new, "kept": kept, "stale": len(stale)}


def worklist(claims: list) -> list:
    """Drift entries that retire something but have no signatures registered yet. Entries
    deliberately curated as having no sweepable signature carry a `NO-SIG:` notes prefix
    (e.g. a value reinstated live elsewhere, or a pure principle) and are excluded."""
    return [c.id for c in claims
            if c.kind == "drift-entry" and c.status in ("retired", "superseded", "withdrawn")
            and not c.signature_strings
            and not (c.notes or "").startswith("NO-SIG")]


def summary(corpus: Corpus, claims: list) -> str:
    by_kind: dict = {}
    for c in claims:
        by_kind.setdefault(c.kind, []).append(c)
    items = sorted(int(c.id.rsplit(":", 1)[1]) for c in by_kind.get("§15-item", []))
    dups = sorted({n for n in items if items.count(n) > 1})
    gaps = [n for n in range(items[0], items[-1] + 1) if n not in items] if items else []
    cited = sorted({s for c in claims for s in c.cited_scripts})
    missing = [s for s in cited if not corpus.script_exists(s)]
    L = []
    L.append("=" * 70)
    L.append(" PTMS extractor — registry summary (computed, not asserted)")
    L.append("=" * 70)
    for k in ("anchor-section", "§15-item", "§14-row", "drift-entry"):
        L.append(f"  {k:16s}: {len(by_kind.get(k, []))}")
    L.append(f"  §15 items present : {len(items)}  range {items[0] if items else '-'}..{items[-1] if items else '-'}")
    L.append(f"  §15 duplicates    : {dups or 'none'}")
    L.append(f"  §15 gaps          : {gaps or 'none'}")
    tagged14 = sum(1 for c in by_kind.get("§14-row", []) if c.tier)
    L.append(f"  §14 rows tagged   : {tagged14}/{len(by_kind.get('§14-row', []))}")
    L.append(f"  distinct cited scripts: {len(cited)}   missing on disk: {missing or 'none'}")
    wl = worklist(claims)
    L.append(f"  drift retirements needing signature curation: {len(wl)}")
    if wl:
        L.append(f"    -> {', '.join(wl)}")
    return "\n".join(L)
