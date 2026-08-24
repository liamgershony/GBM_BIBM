#!/usr/bin/env python3
"""Extract cohort libraries from GSE174554_RAW.tar -> data/interim/GSE174554_RAW/.

data/raw/ is immutable (CLAUDE.md 7.4) and unpacked-but-untransformed content
lives in data/interim/ (7.10). Nothing here is transformed: members are copied out
of the archive and renamed to the standard 10x triplet so scanpy can read them.

Selection is BY EXPLICIT MEMBER NAME from results/tables/discovery_cohort.csv --
never by glob. Mouse, snATAC, proteomics and spatial content is excluded by
ASSERTION on the selected set, plus a per-library positive species check on the
gene symbols actually present.

Layout produced:
    data/interim/GSE174554_RAW/<batch_key>/{barcodes.tsv.gz,features.tsv.gz,matrix.mtx.gz}
"""

from __future__ import annotations

import csv
import gzip
import re
import shutil
import sys
import tarfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
COHORT = REPO / "results" / "tables" / "discovery_cohort.csv"
RAW_TAR = REPO / "data" / "raw" / "GSE174554" / "GSE174554_RAW.tar"
DEST = REPO / "data" / "interim" / "GSE174554_RAW"

# Content in the deposit that must never reach the pipeline.
FORBIDDEN = re.compile(r"Mouse|IR_|T5224|snATAC|Proteomics|Transcriptomics\.tar",
                       re.IGNORECASE)

# Positive species check, run per library on that library's own features file.
HUMAN_GENES = ("EGFR", "GFAP", "PTPRC")
MOUSE_CASED = ("Egfr", "Gfap", "Ptprc")

KIND = {"barcodes": "barcodes.tsv.gz",
        "features": "features.tsv.gz",
        "matrix": "matrix.mtx.gz"}


def kind_of(member: str) -> str | None:
    base = member.rsplit("/", 1)[-1]
    for k in KIND:
        if base.endswith(f"_{k}.tsv.gz") or base.endswith(f"_{k}.mtx.gz"):
            return k
    return None


def main() -> int:
    rows = list(csv.DictReader(open(COHORT)))
    print(f"cohort libraries: {len(rows)}")

    # ---- pre-extraction assertions ---------------------------------------
    selected = []
    for r in rows:
        members = r["tar_members"].split(";")
        assert len(members) == 3, \
            f"{r['batch_key']}: expected 3 tar members, got {len(members)}: {members}"
        kinds = sorted(filter(None, (kind_of(m) for m in members)))
        assert kinds == ["barcodes", "features", "matrix"], \
            f"{r['batch_key']}: members are not a 10x triplet -> {kinds} from {members}"
        selected.extend(members)

    bad = [m for m in selected if FORBIDDEN.search(m)]
    assert not bad, f"forbidden content in selection: {bad[:10]}"
    assert all(r["organism"] == "Homo sapiens" for r in rows)
    assert all(r["platform_id"] == "GPL24676" for r in rows)
    print(f"selected members: {len(selected)} (3 per library, verified per library)")
    print("no mouse / snATAC / proteomics / spatial members in selection")

    # ---- extract ----------------------------------------------------------
    DEST.mkdir(parents=True, exist_ok=True)
    wanted = {m: (r["batch_key"], kind_of(m)) for r in rows
              for m in r["tar_members"].split(";")}
    written = 0
    with tarfile.open(RAW_TAR, "r") as tf:
        for member in tf:
            if member.name not in wanted:
                continue
            batch_key, kind = wanted[member.name]
            outdir = DEST / batch_key
            outdir.mkdir(parents=True, exist_ok=True)
            target = outdir / KIND[kind]
            if target.exists() and target.stat().st_size > 0:
                written += 1
                continue
            src = tf.extractfile(member)
            assert src is not None, f"could not read member {member.name}"
            with open(target, "wb") as fh:
                shutil.copyfileobj(src, fh)
            written += 1
    print(f"extracted/verified {written} files into {DEST.relative_to(REPO)}")

    # ---- post-extraction: per-library species check ------------------------
    checked = 0
    for r in rows:
        feat = DEST / r["batch_key"] / "features.tsv.gz"
        assert feat.exists(), f"missing features for {r['batch_key']}"
        with gzip.open(feat, "rt", encoding="utf-8", errors="replace") as fh:
            syms = {line.split("\t")[1].strip() if "\t" in line else line.strip()
                    for line in fh}
        missing = [g for g in HUMAN_GENES if g not in syms]
        assert not missing, \
            f"{r['batch_key']}: human genes absent {missing} -- not a human library?"
        mouse = [g for g in MOUSE_CASED if g in syms]
        assert not mouse, \
            f"{r['batch_key']}: mouse-cased symbols present {mouse} -- wrong species"
        checked += 1
    print(f"per-library species check passed for all {checked} libraries")

    for r in rows:
        for fn in KIND.values():
            p = DEST / r["batch_key"] / fn
            assert p.exists() and p.stat().st_size > 0, f"empty/missing {p}"
    print(f"all {len(rows)} libraries have a complete non-empty 10x triplet")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
