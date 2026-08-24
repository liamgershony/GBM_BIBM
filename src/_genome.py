"""hg38 gene coordinates and p/q arms, from the UCSC files in data/raw/ucsc_hg38/.

Arms matter: disjoint_set_S names chr9**p**, an arm, not a chromosome. Genes on
9q remain eligible for Stage B selection (CLAUDE.md §3.4), and this is enforced by
tests/test_chr_disjoint.py rather than by care.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
CYTOBAND = REPO / "data" / "raw" / "ucsc_hg38" / "cytoBand.txt.gz"
REFFLAT = REPO / "data" / "raw" / "ucsc_hg38" / "refFlat.txt.gz"

MAIN_CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]


def centromere_boundaries() -> dict[str, int]:
    """chrom -> coordinate where the p arm ends (start of the q arm).

    UCSC marks the centromere with two `acen` bands, `pXX` then `qXX`; the p arm
    ends where the p-side acen band ends.
    """
    out: dict[str, int] = {}
    with gzip.open(CYTOBAND, "rt") as fh:
        for line in fh:
            chrom, start, end, name, stain = line.rstrip("\n").split("\t")
            if stain == "acen" and name.startswith("p"):
                out[chrom] = max(out.get(chrom, 0), int(end))
    return out


def gene_positions() -> pd.DataFrame:
    """gene_symbol -> chromosome, start, end, arm. One row per symbol.

    A symbol with several transcripts is collapsed to its widest span on its most
    common chromosome. Symbols on scaffolds or alt contigs are dropped.
    """
    rows = []
    with gzip.open(REFFLAT, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            sym, chrom, start, end = f[0], f[2], int(f[4]), int(f[5])
            if chrom in MAIN_CHROMS:
                rows.append((sym, chrom, start, end))
    df = pd.DataFrame(rows, columns=["gene", "chromosome", "start", "end"])

    # a symbol can appear on several contigs; keep its most frequent chromosome
    primary = (df.groupby(["gene", "chromosome"]).size()
                 .reset_index(name="n")
                 .sort_values(["gene", "n"], ascending=[True, False])
                 .drop_duplicates("gene")[["gene", "chromosome"]])
    df = df.merge(primary, on=["gene", "chromosome"], how="inner")
    df = (df.groupby(["gene", "chromosome"], as_index=False)
            .agg(start=("start", "min"), end=("end", "max")))

    cen = centromere_boundaries()
    mid = (df["start"] + df["end"]) / 2
    df["arm"] = [("p" if m < cen.get(c, 0) else "q")
                 for m, c in zip(mid, df["chromosome"])]
    return df.set_index("gene")


def in_region_set(chrom: str, arm: str, regions: list[dict]) -> bool:
    """Is this gene inside disjoint_set_S, respecting arm granularity?"""
    for r in regions:
        if chrom == r["chrom"] and (r["arm"] == "both" or arm == r["arm"]):
            return True
    return False


def annotate_var(var_names, regions: list[dict]) -> pd.DataFrame:
    """Build the .var annotation infercnvpy needs, plus in_disjoint_set_S."""
    pos = gene_positions()
    idx = pd.Index(var_names)
    ann = pos.reindex(idx)
    ann["in_disjoint_set_S"] = [
        (isinstance(c, str) and in_region_set(c, a, regions))
        for c, a in zip(ann["chromosome"], ann["arm"])
    ]
    return ann
