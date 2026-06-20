#!/usr/bin/env python3
"""Fetch live Zenodo records and match them to local publication inventories.

The script uses the `ZENODO_TOKEN` environment variable and writes:

- publication_refresh/zenodo/zenodo_depositions_live.json
- publication_refresh/inventory/zenodo_records.csv
- publication_refresh/zenodo/unmatched_zenodo_records.csv
- publication_refresh/zenodo/zenodo_match_report.txt

It deliberately does not read historical `papers/zenodo_state*.json` caches.
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "publication_refresh"
INV = OUT / "inventory"
ZEN = OUT / "zenodo"

API = "https://zenodo.org/api/deposit/depositions"
PAGE_SIZE = 100


def norm(value: str) -> str:
    value = value.lower()
    value = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^\]]*\])?\{([^{}]*)\}", r"\1", value)
    value = re.sub(r"\\[a-zA-Z]+", " ", value)
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def stem_words(path: str) -> str:
    return norm(Path(path).stem.replace("_", " ").replace("-", " "))


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


@dataclass
class LocalCandidate:
    source: str
    local_key: str
    title: str
    pdf: str
    norm_title: str
    file_name: str
    priority: int


def local_candidates() -> list[LocalCandidate]:
    candidates: list[LocalCandidate] = []

    for source, path, key in [
        ("paper", INV / "paper_status.csv", "folder"),
        ("legacy", INV / "legacy_paper_status.csv", "folder"),
    ]:
        for row in read_csv(path):
            pdf = row.get("pdf", "")
            candidates.append(LocalCandidate(
                source=source,
                local_key=row.get(key, "") if source == "paper" else row.get("main_tex", ""),
                title=row.get("title", ""),
                pdf=pdf,
                norm_title=norm(row.get("title", "")),
                file_name=Path(pdf).name if pdf else "",
                priority=3 if source == "paper" else 2,
            ))

    for row in read_csv(INV / "pdf_only_status.csv"):
        pdf = row.get("pdf", "")
        status = row.get("status", "")
        if status in {"ignored", "asset_colocated", "asset_source_unclear", "asset_nested"}:
            continue
        candidates.append(LocalCandidate(
            source="pdf_only",
            local_key=pdf,
            title=stem_words(pdf),
            pdf=pdf,
            norm_title=stem_words(pdf),
            file_name=Path(pdf).name if pdf else "",
            priority=1,
        ))

    return candidates


def request_json(url: str, token: str) -> Any:
    req = urllib.request.Request(url, headers={
        "Authorization": f"Bearer {token}",
        "Accept": "application/json",
        "User-Agent": "octahedrons-publication-refresh/1.0",
    })
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        detail = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"Zenodo API HTTP {e.code}: {detail[:500]}") from e


def fetch_depositions(token: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    page = 1
    while True:
        query = urllib.parse.urlencode({
            "page": page,
            "size": PAGE_SIZE,
            "all_versions": "true",
            "sort": "mostrecent",
        })
        data = request_json(f"{API}?{query}", token)
        if not isinstance(data, list):
            raise RuntimeError(f"Unexpected Zenodo API response: {type(data).__name__}")
        rows.extend(data)
        if len(data) < PAGE_SIZE:
            break
        page += 1
    return rows


def file_names(record: dict[str, Any]) -> list[str]:
    names = []
    for f in record.get("files", []) or []:
        if isinstance(f, dict):
            name = f.get("filename") or f.get("key") or ""
            if name:
                names.append(Path(name).name)
    return names


def record_publication_date(record: dict[str, Any]) -> str:
    md = record.get("metadata", {}) or {}
    return (
        md.get("publication_date")
        or record.get("submitted")
        or record.get("created", "")[:10]
        or record.get("modified", "")[:10]
    )


def record_row(record: dict[str, Any]) -> dict[str, Any]:
    md = record.get("metadata", {}) or {}
    rid = record.get("id", "")
    conceptrecid = record.get("conceptrecid") or record.get("concept_record_id") or md.get("conceptrecid", "")
    doi = record.get("doi") or md.get("doi", "")
    conceptdoi = record.get("conceptdoi") or md.get("conceptdoi", "")
    return {
        "title": md.get("title", ""),
        "doi": doi,
        "conceptdoi": conceptdoi,
        "record_id": rid,
        "conceptrecid": conceptrecid,
        "version": md.get("version", ""),
        "publication_date": record_publication_date(record),
        "local_folder": "",
        "notes": "",
        "files": ";".join(file_names(record)),
        "match_source": "",
        "match_score": "",
        "record_url": (record.get("links", {}) or {}).get("html")
            or (f"https://zenodo.org/record/{rid}" if rid else ""),
    }


def is_published(record: dict[str, Any]) -> bool:
    if record.get("submitted") is True:
        return True
    state = (record.get("state") or "").lower()
    return state == "done" or bool(record.get("doi") or (record.get("metadata", {}) or {}).get("doi"))


def score(record: dict[str, Any], candidate: LocalCandidate) -> tuple[float, str]:
    md = record.get("metadata", {}) or {}
    title = md.get("title", "")
    rt = norm(title)
    title_ratio = SequenceMatcher(None, rt, candidate.norm_title).ratio() if rt and candidate.norm_title else 0.0
    file_set = set(file_names(record))
    file_hit = candidate.file_name in file_set if candidate.file_name else False
    stem_ratio = SequenceMatcher(None, rt, stem_words(candidate.pdf)).ratio() if rt and candidate.pdf else 0.0

    best = title_ratio
    method = "title"
    if file_hit:
        # Exact uploaded filename is strong but not sufficient when old and new
        # papers reuse generic names.  Combine it with title similarity and
        # source priority rather than blindly accepting it.
        combined = 0.72 + 0.20 * max(title_ratio, stem_ratio) + 0.02 * candidate.priority
        if combined > best:
            best = combined
            method = "file+title"
    elif stem_ratio >= 0.82:
        best = max(best, stem_ratio - 0.05)
        method = "stem"

    return best, method


def match_records(records: list[dict[str, Any]], candidates: list[LocalCandidate]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    unmatched: list[dict[str, Any]] = []
    for record in records:
        row = record_row(record)
        if not is_published(record):
            row["notes"] = "unsubmitted draft or unpublished deposition"
            unmatched.append(row)
            rows.append(row)
            continue

        scored = []
        for cand in candidates:
            s, method = score(record, cand)
            scored.append((s, cand.priority, method, cand))
        scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
        best_s, _priority, method, best = scored[0] if scored else (0.0, 0, "none", None)

        if best is not None and best_s >= 0.78:
            row["local_folder"] = best.local_key
            row["match_source"] = best.source
            row["match_score"] = f"{best_s:.3f}"
            row["notes"] = f"matched by {method} to {best.pdf or best.local_key}"
        else:
            row["match_score"] = f"{best_s:.3f}" if scored else ""
            row["notes"] = "no confident local match"
            unmatched.append(row)
        rows.append(row)
    return rows, unmatched


def select_latest_per_concept(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    # Keep all unmatched rows for the Zenodo ledger, but if multiple records map
    # to the same local file/folder, put the latest one last so the inventory
    # generator's local_folder-keyed merge uses the newest live version.
    return sorted(rows, key=lambda r: (
        r.get("local_folder", ""),
        r.get("publication_date", ""),
        str(r.get("record_id", "")),
    ))


def write_report(rows: list[dict[str, Any]], unmatched: list[dict[str, Any]]) -> None:
    by_source = defaultdict(int)
    for row in rows:
        by_source[row.get("match_source", "") or "unmatched"] += 1
    text = [
        "Zenodo live inventory refresh",
        "=============================",
        "",
        f"published/draft records fetched: {len(rows)}",
        f"unmatched or unpublished records: {len(unmatched)}",
        "",
        "Match source counts:",
    ]
    for key in sorted(by_source):
        text.append(f"- {key}: {by_source[key]}")
    text.extend(["", "Unmatched records:"])
    for row in unmatched:
        text.append(f"- {row.get('record_id')} | {row.get('doi')} | {row.get('publication_date')} | {row.get('title')} | {row.get('notes')}")
    (ZEN / "zenodo_match_report.txt").write_text("\n".join(text) + "\n", encoding="utf-8")


def main() -> int:
    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        print("ZENODO_TOKEN is not set", file=sys.stderr)
        return 2

    ZEN.mkdir(parents=True, exist_ok=True)
    records = fetch_depositions(token)
    (ZEN / "zenodo_depositions_live.json").write_text(
        json.dumps(records, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    rows, unmatched = match_records(records, local_candidates())
    rows = select_latest_per_concept(rows)
    fields = [
        "title", "doi", "conceptdoi", "record_id", "conceptrecid", "version",
        "publication_date", "local_folder", "notes", "files", "match_source",
        "match_score", "record_url",
    ]
    write_csv(INV / "zenodo_records.csv", rows, fields)
    write_csv(ZEN / "unmatched_zenodo_records.csv", unmatched, fields)
    write_report(rows, unmatched)
    print(f"fetched {len(records)} Zenodo depositions")
    print(f"wrote {INV / 'zenodo_records.csv'}")
    print(f"wrote {ZEN / 'unmatched_zenodo_records.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
