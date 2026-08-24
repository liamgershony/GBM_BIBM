#!/usr/bin/env python3
"""QC each cohort library -> data/processed/01_qc.h5ad.

Thresholds come from configs/pipeline_config.yaml (frozen) and
configs/runtime_thresholds.yaml (not frozen). No literal thresholds in this file.

ALL NUCLEI ARE RETAINED, malignant and non-malignant alike. inferCNV needs the
non-malignant nuclei as its reference on Day 2 (CLAUDE.md 3.2 / Day 2 plan), so
nothing is dropped on malignancy here. Only clause (d) counting, in
01c_clause_d_gate.py, is malignant-restricted.

Malignancy comes from the authors' own annotation (CLAUDE.md 4) -- never a
classifier we write. The authors annotated only the nuclei THEY retained, so
divergence from our QC is expected in both directions and is not an error:
nuclei we keep that they dropped are labelled `unknown`, never discarded, and are
excluded from both the malignant set and the inferCNV reference set.

Data is single-nucleus (snRNA-seq). The 5% mitochondrial ceiling is correct for
nuclear preps, which are mito-depleted by design (CLAUDE.md 7.5).
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import yaml

REPO = Path(__file__).resolve().parent.parent
COHORT = REPO / "results" / "tables" / "discovery_cohort.csv"
LIBDIR = REPO / "data" / "interim" / "GSE174554_RAW"
ANNOT = REPO / "data" / "interim" / "tumor_normal_metadata.txt"
CONF = REPO / "configs" / "pipeline_config.yaml"
RCONF = REPO / "configs" / "runtime_thresholds.yaml"
OUT_H5 = REPO / "data" / "processed" / "01_qc.h5ad"
OUT_CSV = REPO / "results" / "tables" / "qc_per_library.csv"

sc.settings.verbosity = 1

# Scrublet needs enough nuclei to build a simulated-doublet neighbourhood.
MIN_NUCLEI_FOR_SCRUBLET = 30


def read_library(libdir: Path) -> ad.AnnData:
    """Read one library's 10x triplet.

    GSE174554 does NOT ship the standard 3-column features.tsv (id, symbol, type):
    it ships a SINGLE column of gene symbols, so scanpy's read_10x_mtx fails with
    KeyError: 1. The matrix is also stored genes-x-cells and must be transposed to
    the cells-x-genes orientation AnnData expects.
    """
    with gzip.open(libdir / "features.tsv.gz", "rt", encoding="utf-8",
                   errors="replace") as fh:
        genes = [ln.split("\t")[0].strip() if "\t" in ln else ln.strip()
                 for ln in fh if ln.strip()]
    with gzip.open(libdir / "barcodes.tsv.gz", "rt", encoding="utf-8",
                   errors="replace") as fh:
        barcodes = [ln.strip() for ln in fh if ln.strip()]

    with gzip.open(libdir / "matrix.mtx.gz", "rb") as fh:
        m = scipy.io.mmread(fh)
    m = sp.csr_matrix(m)
    assert m.shape == (len(genes), len(barcodes)), (
        f"{libdir.name}: matrix {m.shape} does not match "
        f"{len(genes)} genes x {len(barcodes)} barcodes")

    a = ad.AnnData(X=sp.csr_matrix(m.T), dtype="float32")
    a.var_names = pd.Index(genes)
    a.obs_names = pd.Index(barcodes)
    return a


# 10x lane suffixes. The libraries use "-1"; the authors' annotation file uses
# "-1" for some specimens and ".1" for others (e.g. SF10108, SF7388, SF1343 --
# where barcode COUNTS match exactly but zero keys join). Both sides are
# normalised by stripping any trailing [-.]<digits> so the convention cannot
# silently split otherwise-identical nuclei.
LANE_SUFFIX = re.compile(r"[-.]\d+$")


def normalise_barcode(bc: str) -> str:
    return LANE_SUFFIX.sub("", bc)


def normalise_key(key: str) -> str:
    pre, _, bc = key.partition("_")
    return f"{pre}_{normalise_barcode(bc)}"


def load_annotation() -> tuple[dict[str, str], set[str], dict[str, int]]:
    """{SampleID}_{barcode} -> 'Tumor' | 'Normal'. Authors' own call.

    Returns (labels, sample_prefixes, rows_per_prefix).
    """
    out, prefixes, counts = {}, set(), {}
    with open(ANNOT, encoding="utf-8", errors="replace") as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) >= 2:
                key = normalise_key(parts[0].strip())
                out[key] = parts[1].strip()
                pre = key.split("_")[0]
                prefixes.add(pre)
                counts[pre] = counts.get(pre, 0) + 1
    return out, prefixes, counts


def annotation_prefix(sample_id: str, prefixes: set[str]) -> str | None:
    """Map a GEO sample_id to the prefix the annotation file actually uses.

    The authors' tumour/normal file does NOT use GEO's specimen naming
    consistently. For 9 cohort specimens GEO says `SF11916` while the annotation
    says `SF11916v2`; GEO has no `SF11916v2` GSM and the annotation has no
    `SF11916`, so the mapping is 1:1 and unambiguous -- verified per specimen
    below. `SF11981` appears under neither form and has no annotation at all.
    """
    if sample_id in prefixes:
        return sample_id
    alias = f"{sample_id}v2"
    if alias in prefixes:
        return alias
    return None


def main() -> int:
    conf = yaml.safe_load(open(CONF))
    rconf = yaml.safe_load(open(RCONF))
    min_genes = conf["qc"]["min_genes_per_nucleus"]
    max_mito = conf["qc"]["max_mito_fraction"]
    seed = conf["seed"]["master"]
    min_join = rconf["qc"]["annotation_join_min_rate"]
    print(f"thresholds -- min_genes={min_genes} max_mito={max_mito} "
          f"seed={seed} min_join_rate={min_join}")

    rows = list(csv.DictReader(open(COHORT)))
    annot, prefixes, ann_counts = load_annotation()
    print(f"libraries: {len(rows)}   annotation entries: {len(annot):,}")

    # ---- resolve the annotation prefix for each specimen, with assertions ----
    alias_of, unmapped = {}, []
    for sid in sorted({r["sample_id"] for r in rows}):
        pre = annotation_prefix(sid, prefixes)
        if pre is None:
            unmapped.append(sid)
        else:
            alias_of[sid] = pre
    used = [p for p in alias_of.values()]
    assert len(used) == len(set(used)), \
        f"two specimens resolved to the same annotation prefix: {used}"
    aliased = {k: v for k, v in alias_of.items() if k != v}
    print(f"annotation prefixes: {len(alias_of)} resolved, {len(unmapped)} unmapped")
    if aliased:
        print(f"  GEO id -> annotation id alias applied for {len(aliased)}: "
              f"{sorted(aliased.items())}")
    if unmapped:
        print(f"  NO annotation under any form (all nuclei -> unknown): {unmapped}")

    # ---- batch2 is NOT covered by the annotation -----------------------------
    # Row counts per specimen match the batch1 barcode count almost exactly, so the
    # authors annotated the primary library only. Worse, batch1 and batch2 can share
    # barcode SEQUENCES (SF7307: all 51 batch2 barcodes also occur in batch1), so a
    # {sample_id}_{barcode} key CANNOT distinguish them and any match is a
    # coincidence, not an identity. Joining batch2 would assign batch1 nuclei's
    # labels to different nuclei. batch2 is therefore always `unknown`.
    print("batch2 libraries are not covered by the annotation -> forced to unknown")

    stats, parts = [], []
    for i, r in enumerate(rows, 1):
        bk = r["batch_key"]
        a = read_library(LIBDIR / bk)
        a.var_names_make_unique()
        n_in = a.n_obs

        a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None,
                                   log1p=False, inplace=True)

        a = a[a.obs["n_genes_by_counts"] >= min_genes].copy()
        n_gene = a.n_obs
        a = a[a.obs["pct_counts_mt"] / 100.0 <= max_mito].copy()
        n_mito = a.n_obs

        # Doublet detection.
        #
        # This block previously caught `Exception`, which swallowed a
        # ModuleNotFoundError in all 60 libraries that reached it -- scrublet was
        # never installed -- and reported success while removing zero doublets.
        # A missing dependency is a fatal environment error, never a per-library
        # data condition, so ImportError is NOT caught here and will propagate.
        # Only genuine numerical failures on small/degenerate libraries are
        # tolerated, and the outcome is recorded in qc_per_library.csv so a
        # silent no-op is visible in the artifact rather than only in stdout.
        n_doub = 0
        if a.n_obs >= MIN_NUCLEI_FOR_SCRUBLET:
            try:
                sc.pp.scrublet(a, random_state=seed, verbose=False)
                n_doub = int(a.obs["predicted_doublet"].sum())
                a = a[~a.obs["predicted_doublet"]].copy()
                scrublet_status = "ran"
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as e:
                print(f"   [{bk}] scrublet numerical failure "
                      f"({type(e).__name__}: {e}); nuclei retained")
                a.obs["predicted_doublet"] = False
                scrublet_status = f"numerical_failure:{type(e).__name__}"
        else:
            a.obs["predicted_doublet"] = False
            scrublet_status = f"skipped_too_few(<{MIN_NUCLEI_FOR_SCRUBLET})"
        n_post = a.n_obs

        # annotation join AFTER our QC, on surviving nuclei
        bare = a.obs_names.str.replace(LANE_SUFFIX, "", regex=True)
        pre = alias_of.get(r["sample_id"])
        applicable = (r["library"] == "batch1") and (pre is not None)
        if applicable:
            keys = pre + "_" + bare
            labels = np.array([annot.get(k, "unknown") for k in keys])
        else:
            keys = pd.Series([""] * a.n_obs, index=a.obs_names)
            labels = np.array(["unknown"] * a.n_obs)
        a.obs["annotation_key"] = keys.values
        a.obs["tumor_normal_annotation"] = labels
        a.obs["is_malignant"] = labels == "Tumor"
        a.obs["is_reference_normal"] = labels == "Normal"
        n_annot = int((labels != "unknown").sum())
        rate = n_annot / n_post if n_post else 0.0

        for col, val in (("patient_id", r["patient_id"]), ("sample_id", r["sample_id"]),
                         ("gsm", r["gsm"]), ("timepoint", r["timepoint"]),
                         ("library", r["library"]), ("batch_key", bk)):
            a.obs[col] = val
        a.obs_names = [f"{bk}_{b}" for b in bare]

        stats.append({"batch_key": bk, "annotation_prefix": pre or "",
                      "annotation_applicable": applicable,
                      "patient_id": r["patient_id"],
                      "sample_id": r["sample_id"], "timepoint": r["timepoint"],
                      "library": r["library"], "n_input": n_in,
                      "n_post_gene_filter": n_gene, "n_post_mito_filter": n_mito,
                      "n_doublets_removed": n_doub,
                      "scrublet_status": scrublet_status, "n_post_qc": n_post,
                      "n_annotated": n_annot, "n_unknown": n_post - n_annot,
                      "annotation_rate": round(rate, 4),
                      "n_malignant": int((labels == "Tumor").sum()),
                      "n_normal": int((labels == "Normal").sum())})
        parts.append(a)
        print(f"[{i:>2}/{len(rows)}] {bk:<20} {n_in:>7,} -> {n_post:>7,}  "
              f"annot={rate:6.1%}  malignant={int((labels=='Tumor').sum()):>6,}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(stats[0].keys()))
        w.writeheader(); w.writerows(stats)
    print(f"\nwrote {OUT_CSV.relative_to(REPO)}")

    adata = ad.concat(parts, join="outer", label=None, index_unique=None)
    for c in ("patient_id", "sample_id", "timepoint", "library", "batch_key",
              "tumor_normal_annotation"):
        adata.obs[c] = adata.obs[c].astype("category")
    OUT_H5.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT_H5, compression="gzip")
    print(f"wrote {OUT_H5.relative_to(REPO)}  {adata.n_obs:,} nuclei x {adata.n_vars:,} genes")

    # ---- reporting + the bug-detector gate --------------------------------
    ran = [s for s in stats if s["scrublet_status"] == "ran"]
    print(f"\ndoublet detection ran in {len(ran)}/{len(stats)} libraries; "
          f"{sum(s['n_doublets_removed'] for s in stats):,} doublets removed")
    if not ran:
        print("FATAL: doublet detection ran in ZERO libraries. Refusing to report "
              "QC counts that silently skipped a pipeline stage.")
        return 1

    tot = sum(s["n_post_qc"] for s in stats)
    print(f"\ntotal post-QC nuclei : {tot:,}")
    print(f"  malignant          : {sum(s['n_malignant'] for s in stats):,}")
    print(f"  normal (reference) : {sum(s['n_normal'] for s in stats):,}")
    print(f"  unknown            : {sum(s['n_unknown'] for s in stats):,}")

    b2 = [s for s in stats if s["library"] == "batch2"]
    print("\nbatch2 libraries (annotation keying was unverified):")
    for s in b2:
        print(f"   {s['batch_key']:<20} post_qc={s['n_post_qc']:>7,} "
              f"annot={s['annotation_rate']:6.1%}")

    # the gate applies only where the annotation is applicable at all
    low = [s for s in stats
           if s["annotation_applicable"] and s["n_post_qc"] > 0
           and s["annotation_rate"] < min_join]
    empty = [s for s in stats if s["n_post_qc"] == 0]
    if empty:
        print(f"\nlibraries with ZERO nuclei surviving QC (join rate undefined): "
              f"{[e['batch_key'] for e in empty]}")
    exempt = [s for s in stats if not s["annotation_applicable"]]
    print(f"\nlibraries exempt from the join gate (no annotation coverage): "
          f"{len(exempt)}")
    for s in exempt:
        why = "batch2 not annotated" if s["library"] == "batch2" else "specimen absent from annotation"
        print(f"   {s['batch_key']:<20} {why}  ({s['n_post_qc']:,} nuclei -> unknown)")
    if low:
        print(f"\nANNOTATION JOIN BELOW {min_join:.0%} in {len(low)} librar(ies):")
        for s in low:
            print(f"   {s['batch_key']:<20} {s['annotation_rate']:6.1%} "
                  f"({s['n_annotated']:,}/{s['n_post_qc']:,})")
        print("This threshold is a bug detector, not a quality gate. Investigate the")
        print("barcode key format before trusting any downstream malignancy call.")
        return 1
    print(f"\nannotation join >= {min_join:.0%} in all {len(stats)} libraries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
