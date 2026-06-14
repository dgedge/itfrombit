#!/usr/bin/env python3
"""Baez crackpot-index lint over the canon (light-hearted; same discipline).
Letter-of-the-law greps for the index's mechanically checkable items; the
judgment items are scored in the accompanying discussion. exit 0 always —
this is a lint, not a gate (item 20: complaining about the index = 20 pts;
we shall not)."""
import re
from pathlib import Path

root = Path(__file__).parent.parent
corpus = ""
for f in ["ANCHOR.md", "DRIFT.md"]:
    corpus += (root / f).read_text(encoding="utf-8")
for f in (root / "technical_notes").glob("*.md"):
    corpus += f.read_text(encoding="utf-8")

print("[item 8] misspellings (5 pts each): Einstien / Feynmann / Hawkins:")
for w in ("Einstien", "Feynmann", "Hawkins"):
    n = corpus.count(w)
    print(f"    {w}: {n}  -> {5*n} pts")

print("\n[items 19/27/28/34/35/36/16/21] rhetoric red-flags (the ranty phrases):")
flags = ["paradigm shift", "hidebound", "self-appointed", "conspiracy",
         "Galileo", "Nobel", "only a theory", "establishment", "suppress",
         "revolutionary", "sham"]
total_hits = 0
for p in flags:
    n = len(re.findall(p, corpus, re.IGNORECASE))
    total_hits += n
    if n:
        print(f"    '{p}': {n} hit(s)")
print(f"    total rhetoric hits: {total_hits}")

print("\n[item 7, letter-of-the-law] ALL-CAPS words (5 pts each, 'defective")
print("keyboard' exemption pending): the canon's verdict typography:")
caps = re.findall(r"\b[A-Z]{4,}\b", corpus)
allowed = {"ANCHOR","DRIFT","HMRC","ITSA"}  # acronyms happen
caps = [c for c in caps if c not in allowed]
from collections import Counter
top = Counter(caps).most_common(8)
print(f"    all-caps tokens (4+ letters): {len(caps)}  -> {5*len(caps):,} pts (!!)")
print(f"    the leading offenders: {', '.join(f'{w}x{n}' for w,n in top)}")
print("\nexit 0 — lint complete; scorecard in the discussion.")
