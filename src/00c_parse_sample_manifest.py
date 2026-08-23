#!/usr/bin/env python3
"""Parse GSE174554_family.soft.gz into results/tables/sample_manifest.csv.

One row per GSM. This is a transcription of what GEO records, not an analysis:
no cohort rule is applied and no sample is filtered out. Selecting the discovery
cohort is a separate, written decision (CLAUDE.md 10.1) taken by Lanes 1 and 2
together; this script exists so that decision is made against the authoritative
record rather than against filename conventions.

Reads the gzipped SOFT file from data/raw/ without modifying it, and writes an
uncompressed copy to data/interim/ (CLAUDE.md 7.10).

Usage:  python3 src/00c_parse_sample_manifest.py
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _download_utils import read_id_column  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
SOFT_GZ = REPO / "data" / "raw" / "GSE174554" / "GSE174554_family.soft.gz"
INTERIM = REPO / "data" / "interim" / "GSE174554_family.soft"
OUT_CSV = REPO / "results" / "tables" / "sample_manifest.csv"

# Fields taken verbatim from the SOFT record.
SCALAR_FIELDS = {
    "!Sample_title": "title",
    "!Sample_platform_id": "platform_id",
    "!Sample_organism_ch1": "organism",
    "!Sample_library_strategy": "library_strategy",
    "!Sample_library_source": "library_source",
    "!Sample_type": "sample_type",
    "!Sample_source_name_ch1": "source_name",
    "!Sample_description": "description",
    "!Sample_instrument_model": "instrument_model",
}


def parse_soft(path: Path) -> list[dict]:
    """Split the SOFT family file into per-GSM records."""
    records: list[dict] = []
    cur: dict | None = None

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                if cur is not None:
                    records.append(cur)
                cur = {"gsm": line.split("=", 1)[1].strip(), "characteristics": []}
                continue
            if cur is None or "=" not in line:
                continue
            key, val = line.split("=", 1)
            key, val = key.strip(), val.strip()
            if key == "!Sample_characteristics_ch1":
                cur["characteristics"].append(val)
            elif key in SCALAR_FIELDS:
                # A few fields repeat; keep the first, which is the primary value.
                cur.setdefault(SCALAR_FIELDS[key], val)

    if cur is not None:
        records.append(cur)
    return records


def characteristic(chars: list[str], key: str) -> str:
    """Pull 'key: value' out of the characteristics list. '' if absent."""
    for c in chars:
        if ":" in c:
            k, v = c.split(":", 1)
            if k.strip().lower() == key.lower():
                return v.strip()
    return ""


def sample_id_from_title(title: str) -> str:
    """Titles look like 'GBM SF10099_snRNA' or 'GBM SF12704v2_snRNA'.

    Stop at the underscore: \w includes '_', so a naive \w+ swallows the assay
    suffix and yields 'SF12704v2_snRNA' instead of 'SF12704v2'.
    """
    m = re.search(r"\b(SF[^\s_]+)", title)
    return m.group(1) if m else ""


def assay_from(record: dict) -> str:
    """Coarse assay label from title/description/library strategy."""
    blob = " ".join([record.get("title", ""), record.get("description", ""),
                     record.get("library_strategy", "")]).lower()
    if "atac" in blob:
        return "snATAC-seq"
    if "snrna" in blob or "single nuc" in blob or "nuclei" in blob:
        return "snRNA-seq"
    if "rna-seq" in record.get("library_strategy", "").lower():
        return "RNA-seq (unspecified)"
    return record.get("library_strategy", "unknown")


def main() -> int:
    INTERIM.parent.mkdir(parents=True, exist_ok=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(SOFT_GZ, "rt", encoding="utf-8", errors="replace") as fin:
        INTERIM.write_text(fin.read(), encoding="utf-8")

    records = parse_soft(INTERIM)
    print(f"parsed {len(records)} GSM records from {SOFT_GZ.name}")

    rows = []
    for r in records:
        chars = r.get("characteristics", [])
        rows.append({
            "gsm": r["gsm"],
            "title": r.get("title", ""),
            "sample_id": sample_id_from_title(r.get("title", "")),
            "platform_id": r.get("platform_id", ""),
            "organism": r.get("organism", ""),
            "assay": assay_from(r),
            "library_strategy": r.get("library_strategy", ""),
            "sample_type": r.get("sample_type", ""),
            "source_name": r.get("source_name", ""),
            "patient_pair": read_id_column(characteristic(chars, "pair#")) or "",
            "patient_pair_raw": characteristic(chars, "pair#"),
            "progression": characteristic(chars, "progression"),
            "diagnosis": characteristic(chars, "diagnosis"),
            "age": characteristic(chars, "age"),
            "gender": characteristic(chars, "gender"),
            "description": r.get("description", ""),
            "characteristics_raw": " | ".join(chars),
        })

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT_CSV.relative_to(REPO)}  ({len(rows)} rows)")

    # ---------------- summary: report what the file says ----------------
    def show(title, counter):
        print(f"\n{title}")
        for k, n in counter.most_common():
            print(f"   {n:>4}  {k or '(empty)'}")

    show("organism", Counter(r["organism"] for r in rows))
    show("platform", Counter(r["platform_id"] for r in rows))
    show("assay", Counter(r["assay"] for r in rows))
    show("progression (timepoint)", Counter(r["progression"] for r in rows))
    show("diagnosis", Counter(r["diagnosis"] for r in rows))

    human_sn = [r for r in rows
                if r["organism"] == "Homo sapiens" and r["assay"] == "snRNA-seq"]
    print(f"\nhuman snRNA-seq GSMs: {len(human_sn)}")

    with_pair = [r for r in human_sn if r["patient_pair"]]
    print(f"  ...carrying a pair# : {len(with_pair)}")
    print(f"  ...distinct pair#   : {len({r['patient_pair'] for r in with_pair})}")

    by_pair: dict[str, set] = defaultdict(set)
    for r in with_pair:
        by_pair[r["patient_pair"]].add(r["progression"])

    both = {p for p, tp in by_pair.items()
            if any("primary" in t.lower() for t in tp)
            and any("recurren" in t.lower() for t in tp)}
    print(f"  ...pairs with BOTH primary and recurrent: {len(both)}")

    print("\n  timepoint sets observed per pair#:")
    for combo, n in Counter(tuple(sorted(tp)) for tp in by_pair.values()).most_common():
        print(f"     {n:>3}  {combo}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
