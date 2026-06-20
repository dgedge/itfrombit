"""registry.py — the sidecar claim registry (claims.jsonl).

One JSON object per line (diff-friendly, append-friendly, stdlib only). The registry is a
SIDECAR: the human-facing ANCHOR/DRIFT/STATUS prose is never edited. Auto fields are produced
by extract.py; human-curated fields (signature_strings, multipart_count, supersedes) are
preserved across re-extraction by `id`.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path

SCHEMA_VERSION = "0.1"

# fields a re-extraction is allowed to overwrite; everything else is human-owned and preserved
AUTO_FIELDS = {"kind", "location", "tier", "status", "cited_scripts", "cited_refs",
               "auto_hint", "extractor_version"}
HUMAN_FIELDS = {"supersedes", "superseded_by", "signature_strings", "multipart_count", "notes", "status_override"}


@dataclass
class Claim:
    id: str                                   # anchor:§15:116 | anchor:§14:row:042 | drift:M9 | anchor:sec:6.6
    kind: str                                 # §15-item | §14-row | drift-entry | anchor-section
    location: dict = field(default_factory=dict)         # {file, line, anchor_text}
    tier: str | None = None                   # Locked | Proposition | Open | compound | None
    status: str = "live"                      # live|retired|superseded|withdrawn|flagged|proposition|open
    cited_scripts: list = field(default_factory=list)    # ["python_code/X.py", ...]
    cited_refs: list = field(default_factory=list)        # ["§15 item 86", "DRIFT M12", "§6.6"]
    supersedes: list = field(default_factory=list)        # registry ids this claim retires
    superseded_by: list = field(default_factory=list)
    signature_strings: list = field(default_factory=list)  # LOAD-BEARING: grep targets for sync
    multipart_count: int = 0                  # K = # distinct things a DRIFT entry retires (HUMAN-set)
    notes: str = ""                           # human note
    status_override: str = ""                 # HUMAN: pin status when auto-extraction mis-flags (live item w/ a retired sub-leg)
    auto_hint: str = ""                       # extractor suggestion (e.g. multipart K guess); not authoritative
    extractor_version: str = SCHEMA_VERSION

    def to_json(self) -> str:
        return json.dumps(asdict(self), ensure_ascii=False, sort_keys=True)

    @classmethod
    def from_dict(cls, d: dict) -> "Claim":
        known = set(cls.__dataclass_fields__)
        return cls(**{k: v for k, v in d.items() if k in known})

    def merge_auto(self, auto: "Claim") -> "Claim":
        """Return a copy with AUTO_FIELDS taken from `auto`, HUMAN_FIELDS kept from self."""
        merged = asdict(self)
        for f in AUTO_FIELDS:
            merged[f] = getattr(auto, f)
        if merged.get("status_override"):          # human pin overrides the auto-extracted status
            merged["status"] = merged["status_override"]
        return Claim.from_dict(merged)


def load_registry(path) -> list[Claim]:
    p = Path(path)
    if not p.exists():
        return []
    out = []
    for raw in p.read_text(encoding="utf-8").splitlines():
        s = raw.strip()
        if not s or s.startswith("//") or s.startswith("#"):
            continue
        out.append(Claim.from_dict(json.loads(s)))
    return out


def save_registry(path, claims: list[Claim]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(c.to_json() for c in claims) + "\n", encoding="utf-8")


def index_by_id(claims: list[Claim]) -> dict:
    return {c.id: c for c in claims}
