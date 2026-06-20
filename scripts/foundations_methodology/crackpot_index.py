#!/usr/bin/env python3
"""
crackpot_index.py  --  an automated, physics-aware Crackpot Index scorer for PDFs.

Implements John Baez's "Crackpot Index" (https://math.ucr.edu/home/baez/crackpot.html)
as far as it can be honestly mechanised, and runs it over the *text of a paper* (PDF,
or .txt/.md). The intended use is triage: physicists get a torrent of self-published
"theories of everything" and have no time to read them. A low score does not mean a
theory is correct -- it means the document is free of the *rhetorical* tells that
correlate with cranks. A high score means: probably do not spend your afternoon on it.

With thanks to John C. Baez (UC Riverside) -- mathematical physicist, AMS Conant Prize
for expository writing, and author of the Crackpot Index. This little utility is an
homage to a body of work whose clarity it could only aspire to. All credit to him; any
clumsiness in the automation is mine.

WHAT IT DOES AND DOES NOT SCORE
-------------------------------
Baez's index has 37 items. Many require a physicist's judgement ("a statement widely
agreed to be false", "logically inconsistent", "no concrete testable predictions").
Those CANNOT be scored by a text scan without lying about it, so this tool:

  * AUTO-SCORES the mechanically detectable items (rhetoric phrases, ALL-CAPS shouting,
    misspelled famous names, self-comparison patterns, ...). These are the items that
    actually separate a crank manuscript from an ordinary technical paper.
  * LISTS the judgement-only items as a manual-review checklist, scored 0, never folded
    silently into the number.

Item 6 ("5 points for each word in all capital letters") is the one item that wrecks a
naive scan of any technical paper, because acronyms (QCD, PMNS, CSS) are all-caps. We
resolve it the honest way: an ALL-CAPS token scores ONLY if its lower-case form is a real
English word (i.e. someone is SHOUTING a word), checked against the system dictionary.
Acronyms are not words, so they do not score. Pass --acronyms to extend the ignore set.

Item 7 in Baez's list is the *misspelling* rule ("Einstien/Hawkins/Feynmann"), so a
correctly-spelled citation of Feynman or Einstein does NOT score -- as it should not.

USAGE
-----
  python3 crackpot_index.py paper.pdf
  python3 crackpot_index.py *.pdf                  # comparison table over many files
  python3 crackpot_index.py paper.pdf --verbose    # show every triggering match
  python3 crackpot_index.py paper.pdf --json
  python3 crackpot_index.py paper.pdf --acronyms ANCHOR DRIFT TCH

Exit code is 0 (it is a lint, not a gate). Requires `pdftotext` (poppler) for PDFs,
or one of pypdf / PyPDF2 / pdfminer.six as a fallback.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

BAEZ_URL = "https://math.ucr.edu/home/baez/crackpot.html"
STARTING_CREDIT = -5

# --------------------------------------------------------------------------------------
# PDF / text extraction
# --------------------------------------------------------------------------------------

def extract_text(path: Path) -> str:
    """Return the plain text of a .pdf / .txt / .md file."""
    suffix = path.suffix.lower()
    if suffix in (".txt", ".md", ".tex"):
        return path.read_text(encoding="utf-8", errors="replace")
    if suffix != ".pdf":
        raise ValueError(f"unsupported file type: {path.suffix} (use .pdf/.txt/.md)")

    # Primary: poppler's pdftotext (fast, reliable, preserves case).
    if shutil.which("pdftotext"):
        out = subprocess.run(
            ["pdftotext", "-q", str(path), "-"],
            capture_output=True, text=True,
        )
        if out.returncode == 0 and out.stdout.strip():
            return out.stdout

    # Fallbacks: any installed python PDF lib.
    try:
        import pypdf  # type: ignore
        return "\n".join(p.extract_text() or "" for p in pypdf.PdfReader(str(path)).pages)
    except Exception:
        pass
    try:
        import PyPDF2  # type: ignore
        return "\n".join(p.extract_text() or "" for p in PyPDF2.PdfReader(str(path)).pages)
    except Exception:
        pass
    try:
        from pdfminer.high_level import extract_text as _pm  # type: ignore
        return _pm(str(path))
    except Exception:
        pass
    raise RuntimeError(
        "No PDF backend available. Install poppler (`brew install poppler`) "
        "for `pdftotext`, or `pip install pypdf`."
    )


# --------------------------------------------------------------------------------------
# Dictionary (for the ALL-CAPS / shouting heuristic)
# --------------------------------------------------------------------------------------

def load_dictionary() -> set[str]:
    for p in ("/usr/share/dict/words", "/usr/dict/words"):
        fp = Path(p)
        if fp.exists():
            return {w.strip().lower() for w in fp.read_text(errors="replace").splitlines() if w.strip()}
    return set()  # heuristic degrades to "acronym allowlist only" if no dictionary

# A small built-in fallback of common shouted words, used only if no system dictionary.
_FALLBACK_WORDS = {
    "truth", "true", "false", "wrong", "right", "proof", "proven", "fact", "facts",
    "everything", "nothing", "impossible", "obvious", "nonsense", "real", "reality",
    "revolution", "revolutionary", "establishment", "conspiracy", "suppressed",
    "genius", "discovery", "breakthrough", "fundamental", "absolute", "ultimate",
    "must", "cannot", "never", "always", "the", "and", "but", "god", "universe",
}

# ALL-CAPS tokens that are legitimate physics/maths acronyms or labels (never score).
DEFAULT_ACRONYM_ALLOW = {
    "ANCHOR", "DRIFT", "TCH", "QCD", "QED", "QFT", "SM", "CKM", "PMNS", "CSS",
    "CPT", "CP", "GUT", "BCS", "OZI", "MOND", "BTFR", "CNOT", "SWAP", "XOR",
    "BCC", "FCC", "TCH", "HDR", "HQET", "DGG", "PDG", "CODATA", "BBN", "DESI",
    "LEGEND", "GERDA", "CUPID", "BICEP", "AME", "MUB", "RG", "IR", "UV", "EW",
    "SMG", "KMS", "AQUAL", "LNH", "ER", "EPR", "DNA", "RNA", "AI", "LLM", "GPU",
    "ISBN", "DOI", "URL", "PDF", "TEX", "HTML", "JSON", "API", "USA", "UK",
    "NESS", "FLDW", "TMS", "PTMS", "STATUS", "MEMORY",
    # acronyms whose lower-case form is, annoyingly, a real dictionary word:
    "EFT",   # effective field theory  (lc 'eft' = a juvenile newt)
    "ABE",   # ...
    "WIMP",  # weakly-interacting massive particle ('wimp')
    "MACHO", # massive compact halo object ('macho')
    "SUSY", "TOE", "BAO", "DOF", "VEV", "LSP", "GIM", "AdS",
}


# --------------------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------------------

class Hit:
    __slots__ = ("item", "points", "label", "matches")

    def __init__(self, item: str, points: int, label: str, matches: list[str]):
        self.item = item
        self.points = points
        self.label = label
        self.matches = matches


def _find(patterns, text, flags=re.IGNORECASE):
    out = []
    for pat in patterns:
        out += [m.group(0) for m in re.finditer(pat, text, flags)]
    return out


def score_text(text: str, acronyms: set[str], words: set[str]) -> tuple[list[Hit], list[str]]:
    """Return (auto-scored hits, manual-review checklist)."""
    hits: list[Hit] = []

    def add(item, per, label, matches):
        if matches:
            hits.append(Hit(item, per * len(matches), label, matches))

    # --- Item 6: ALL-CAPS shouting (5 pts/word). ---
    # A genuine "shout" is an ISOLATED caps word inside normal sentence-case text.
    # A RUN of consecutive caps words is a heading / title / emphasis block -- that is
    # typography, not ranting -- so those are reported separately and NOT scored.
    # Acronyms (QCD, PMNS) never count: only words whose lower-case form is in the
    # dictionary are treated as shoutable.
    tokens = re.findall(r"\S+", text)

    def _capsish(tok: str) -> bool:  # caps neighbour test for run detection (>=2 letters)
        core = re.sub(r"[^A-Za-z]", "", tok)
        return len(core) >= 2 and core.isupper()

    shouted, headingish = [], []
    for i, tok in enumerate(tokens):
        core = re.sub(r"[^A-Za-z]", "", tok)
        if len(core) < 3 or not core.isupper() or core in acronyms:
            continue
        low = core.lower()
        is_word = (low in words) if words else (low in _FALLBACK_WORDS)
        if not is_word:
            continue
        prev_caps = i > 0 and _capsish(tokens[i - 1])
        next_caps = i < len(tokens) - 1 and _capsish(tokens[i + 1])
        (headingish if (prev_caps or next_caps) else shouted).append(core)
    add("6  ALL-CAPS shouted words", 5, "shouting an isolated word in caps", shouted)
    if headingish:  # informational, 0 points, for transparency
        hits.append(Hit("6* caps in headings/titles", 0,
                        "stylistic caps runs ignored (typographic, not shouting)",
                        headingish))

    # --- Item 7: misspelled famous names (5 pts each). ---
    add("7  misspelled names", 5, "Einstien / Hawkins / Feynmann",
        _find([r"\bEinstien\b", r"\bHawkins\b", r"\bFeynmann\b", r"\bSchrodfinger\b",
               r"\bHeisenburg\b", r"\bGallileo\b"], text))

    # --- Item 8: quantum mechanics fundamentally misguided. ---
    add("8  QM called misguided", 10, "claims QM is fundamentally wrong/misguided",
        _find([r"quantum mechanics is (fundamentally )?(wrong|misguided|flawed|incorrect|broken)",
               r"(overthrow|refute|disprove)s? quantum mechanics"], text))

    # --- Item 9: schooling as evidence of sanity. ---
    add("9  schooling as proof of sanity", 10, "'I went to school, therefore...'",
        _find([r"I have (a (PhD|degree)|gone to (school|university)).{0,40}(sane|crazy|crank|credential)"], text))

    # --- Item 10: opening with how long you've worked on it. ---
    add("10 'worked on this for N years'", 10, "boasting of years spent",
        _find([r"(I have|I've) (been )?(work|spent|labou?red)(ing|ed)? (on (this|my theory) )?(for )?\d+ ?(years|decades)",
               r"after \d+ (years|decades) of (work|research) (on )?(my|this) (theory|idea)"], text))

    # --- Item 13: hard to mechanise (undefined neologisms) -> manual. ---

    # --- Item 14: "I'm not good at math but...". ---
    add("14 'not good at math, but...'", 10, "outsourcing the mathematics",
        _find([r"(I('?m| am) )?not (good|great) (at|with) (the )?math(s|ematics)?",
               r"all I need is (for )?someone to (express|write|put) (it|this) in.{0,20}(equations|math)"], text))

    # --- Item 15: "only a theory". ---
    add("15 'only a theory'", 10, "'just/only a theory' as a slur",
        _find([r"(just|only) a theory"], text))

    # --- Item 16: 'predicts but does not explain why / no mechanism' as an attack. ---
    add("16 'correct but no mechanism' jibe", 10, "deriding accepted theory for lacking a 'why'",
        _find([r"predicts.{0,40}but (does ?n'?t|fails to) (explain|provide).{0,20}(why|mechanism)",
               r"no (real )?(mechanism|explanation of why)"], text))

    # --- Item 17: favorable self-comparison to Einstein / relativity is wrong. ---
    add("17 self-compared to Einstein / SR-GR wrong", 10, "Einstein/relativity rhetoric",
        _find([r"(like|as great as|the new|another|a modern) Einstein",
               r"(special|general) relativity is (wrong|misguided|flawed|incorrect)",
               r"(overthrow|refute|disprove)s? (special |general )?relativity"], text))

    # --- Item 18: 'paradigm shift'. ---
    add("18 'paradigm shift'", 10, "cutting-edge paradigm-shift talk",
        _find([r"paradigm[- ]shift(ing)?"], text))

    # --- Item 20: deserves a Nobel. ---
    add("20 deserves a Nobel", 20, "Nobel-prize entitlement",
        _find([r"(deserve|worthy of|merit)s? (a |the )?Nobel",
               r"Nobel[- ]prize[- ]worthy"], text))

    # --- Item 21: self-compared to Newton / classical mechanics is wrong. ---
    add("21 self-compared to Newton", 20, "Newton rhetoric",
        _find([r"(like|as great as|the new|another) Newton",
               r"classical mechanics is (wrong|misguided|fundamentally flawed)"], text))

    # --- Item 22: science fiction / myth as fact. (light heuristic) ---
    add("22 sci-fi/myth as fact", 20, "fiction or myth cited as evidence",
        _find([r"as (shown|proven|predicted) (in|by) (Star ?Trek|Star ?Wars|the Bible|ancient (myths|texts))"], text))

    # --- Item 24: naming the theory after yourself. (manual; light heuristic only) ---

    # --- Item 26 / 27: the signature ranty phrases. ---
    add("26 'hidebound reactionary'", 20, "ranty phrase", _find([r"hidebound reactionar"], text))
    add("27 'self-appointed defender of the orthodoxy'", 20, "ranty phrase",
        _find([r"self[- ]appointed defender", r"defender of (the )?orthodoxy"], text))

    # --- Item 29: Einstein in his later years groping toward your ideas. ---
    add("29 'Einstein was groping toward this'", 30, "claiming Einstein anticipated you",
        _find([r"Einstein.{0,40}(later years|near the end).{0,40}(toward|towards|anticipat)"], text))

    # --- Item 30: theory from extraterrestrials. ---
    add("30 extraterrestrial origin", 30, "aliens gave you the theory",
        _find([r"(given|revealed|taught) (to me )?by (an? )?(extraterrestrial|alien)"], text))

    # --- Item 32: opponents = Nazis. ---
    add("32 opponents compared to Nazis", 40, "Nazi/stormtrooper comparison",
        _find([r"\b(Nazis?|stormtroopers?|brownshirts?|Gestapo)\b"], text))

    # --- Item 33: 'scientific establishment' 'conspiracy'. ---
    add("33 establishment 'conspiracy'", 40, "conspiracy against your work",
        _find([r"(scientific )?establishment.{0,40}conspiracy",
               r"conspiracy.{0,40}(to )?(suppress|silence|bury) (my|this) (work|theory)"], text))

    # --- Item 34: Galileo / modern Inquisition. ---
    add("34 Galileo / Inquisition", 40, "comparing yourself to Galileo",
        _find([r"(like|comparing myself to|another) Galileo",
               r"modern[- ]day Inquisition"], text))

    # --- Item 35: 'present-day science is a sham'. ---
    add("35 'science is a sham'", 40, "all of science is a sham",
        _find([r"(present[- ]day |modern )?science (will be seen|is) (for )?the sham",
               r"the sham (it|that) (it )?(truly|really) is"], text))

    # --- Item 36: revolutionary, no testable predictions: detect the boast only. ---
    add("36 'revolutionary' (check for predictions!)", 0,
        "FLAGGED not scored: confirm whether testable predictions exist",
        _find([r"\brevolutionar(y|ise|ize)s?\b"], text))

    # --- judgement-only items: never auto-scored ---
    manual = [
        "1  statements widely agreed to be false (+1 each)",
        "2  vacuous statements (+2 each)",
        "3  logically inconsistent statements (+3 each)",
        "4  errors adhered to despite correction (+5 each)",
        "5  thought experiment contradicting a real experiment (+5)",
        "13 new terms used without definition (+10 each)",
        "19 complaining about the crackpot index (+20)  [n/a to papers]",
        "23 citing past ridicule of your theories (+20)",
        "24 naming something after yourself (+20)",
        "25 praising the theory without ever explaining it (+20)",
        "28 famous figure secretly disbelieved a theory (+30)",
        "31 asylum-delay allusions (+30)",
        "36 revolutionary claim with NO testable predictions (+50)  <- decisive; judge manually",
    ]
    return hits, manual


# --------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------

def band(total: int) -> str:
    if total < 0:
        return "no crackpot markers detected (productive-physicist range)"
    if total <= 10:
        return "essentially clean; a stray stylistic flag or two"
    if total <= 30:
        return "minor rhetorical flags; read with mild caution"
    if total <= 75:
        return "several crank tells; triage carefully"
    if total <= 200:
        return "strong crank signature"
    return "textbook crackpot"


def report_one(path: Path, hits: list[Hit], manual: list[str], verbose: bool) -> int:
    total = STARTING_CREDIT + sum(h.points for h in hits)
    print(f"\n=== {path.name} ===")
    print(f"starting credit: {STARTING_CREDIT}")
    if hits:
        print("auto-scored items:")
        for h in sorted(hits, key=lambda x: -x.points):
            n = len(h.matches)
            print(f"  +{h.points:>4}  [{h.item}]  ({n} hit{'s' if n != 1 else ''}) -- {h.label}")
            if verbose:
                shown = ", ".join(sorted(set(h.matches))[:12])
                print(f"          e.g. {shown}")
    else:
        print("auto-scored items: none triggered")
    print(f"\nAUTO SCORE: {total}   -> {band(total)}")
    print("\nmanual-review checklist (NOT in the number above -- needs a physicist):")
    for m in manual:
        print(f"  [ ] item {m}")
    print(f"\n(Baez index: {BAEZ_URL}.  Low = free of rhetorical crank-tells; says")
    print(" NOTHING about whether the physics is correct. Judgement items above are yours.)")
    return total


def main(argv=None):
    ap = argparse.ArgumentParser(description="Automated Baez Crackpot Index scorer for PDFs.")
    ap.add_argument("files", nargs="+", help="PDF / txt / md files to score")
    ap.add_argument("--verbose", action="store_true", help="show example matches per item")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    ap.add_argument("--acronyms", nargs="*", default=[], help="extra ALL-CAPS tokens to ignore")
    args = ap.parse_args(argv)

    if not args.json:
        print("-" * 78)
        print("CRACKPOT INDEX (Baez) -- automated, light-hearted triage. NOT a verdict on")
        print("correctness: it scans for the *rhetorical* tells of cranks, nothing more. A")
        print("low score just means the prose is sane; the physics still has to be checked.")
        print("Better than rejecting on email domains or keywords -- it reads the actual text.")
        print("-" * 78)

    words = load_dictionary()
    acronyms = set(DEFAULT_ACRONYM_ALLOW) | {a.upper() for a in args.acronyms}

    results = {}
    for f in args.files:
        path = Path(f)
        try:
            text = extract_text(path)
        except Exception as e:
            print(f"!! {path.name}: {e}", file=sys.stderr)
            continue
        hits, manual = score_text(text, acronyms, words)
        total = STARTING_CREDIT + sum(h.points for h in hits)
        results[path.name] = {
            "auto_score": total,
            "band": band(total),
            "items": {h.item: {"points": h.points, "hits": len(h.matches)} for h in hits},
            "manual_checklist": manual,
        }
        if not args.json:
            report_one(path, hits, manual, args.verbose)

    if args.json:
        print(json.dumps(results, indent=2))
    elif len(results) > 1:
        print("\n=== comparison (auto score, lower is better) ===")
        for name, r in sorted(results.items(), key=lambda kv: kv[1]["auto_score"]):
            print(f"  {r['auto_score']:>5}  {name}   ({r['band']})")
        if not words:
            print("\nNOTE: no system dictionary found; ALL-CAPS heuristic ran in fallback mode.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
