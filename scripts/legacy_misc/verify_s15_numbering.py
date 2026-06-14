#!/usr/bin/env python3
import re, subprocess, sys
from collections import Counter
cur = open('ANCHOR.md', encoding='utf-8').read()
ref = subprocess.run(['git','show','dab77db:ANCHOR.md'],capture_output=True,text=True).stdout
hdr = lambda t: [int(n) for n in re.findall(r'(?m)^(1[0-9][0-9])\.\s+\*\*', t)]  # §15 100s headers
cur_n, ref_n = hdr(cur), hdr(ref)
dups    = sorted(n for n,c in Counter(cur_n).items() if c > 1)
missing = sorted(set(ref_n) - set(cur_n))
extra   = sorted(set(cur_n) - set(ref_n))   # new items since dab77db are OK — review each
print(f"100s headers: current {len(cur_n)} ({len(set(cur_n))} distinct); ref {len(set(ref_n))} distinct")
print(f"DUPLICATES (must be empty): {dups}")
print(f"MISSING vs ref (must be empty unless retired on purpose): {missing}")
print(f"EXTRA vs ref (new items — confirm intentional): {extra}")
sys.exit(0 if not dups and not missing else 1)
