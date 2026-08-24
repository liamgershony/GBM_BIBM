#!/usr/bin/env python3
"""Harmony integration -> data/processed/02_integrated.h5ad, plus the Step 3 gate.

HARMONY INTEGRATES ON `patient_id`, NEVER ON `sample_id` OR `batch_key`.

This is the single most consequential line in the file. `sample_id` distinguishes
a patient's primary specimen from their recurrent specimen. Integrating on it
would instruct Harmony to remove the primary-vs-recurrent difference as though it
were technical noise -- and that difference IS RAS component T, the transcriptional
similarity between a primary cell and its patient's recurrent centroid. The
pipeline would run cleanly and return a null that looked methodologically sound.
See the DEVIATIONS.md entry of 2026-08-24 correcting the frozen config.

`batch_key` and `library` ride along in .obs as covariates for downstream use.

Step 3 gate. LISI is computed on `patient_id` (should be HIGH -- patients mixed,
correction worked) and on `timepoint` (should be LOW -- primary and recurrent
still separable, signal preserved). Raw LISI is bounded by category count -- 1..n
for patients but 1..2 for timepoint -- so raw values are NOT comparable and both
are normalised to (LISI - 1)/(k - 1) in [0, 1]. Integration FAILS if normalised
timepoint LISI >= lisi_timepoint_fail_ratio x normalised patient LISI.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scanpy.external as sce
import yaml
from harmonypy.lisi import compute_lisi

REPO = Path(__file__).resolve().parent.parent
QC_H5 = REPO / "data" / "processed" / "01_qc.h5ad"
COHORT_N = REPO / "results" / "tables" / "cohort_n.json"
CONF = REPO / "configs" / "pipeline_config.yaml"
RCONF = REPO / "configs" / "runtime_thresholds.yaml"
OUT_H5 = REPO / "data" / "processed" / "02_integrated.h5ad"
OUT_LISI = REPO / "results" / "tables" / "lisi_gate.csv"

sc.settings.verbosity = 1


def main() -> int:
    conf = yaml.safe_load(open(CONF))
    rconf = yaml.safe_load(open(RCONF))
    batch_key = conf["integration"]["batch_key"]
    seed = conf["seed"]["master"]
    n_hvg = conf["features"]["n_hvg"]
    fail_ratio = rconf["integration"]["lisi_timepoint_fail_ratio"]

    assert batch_key == "patient_id", (
        f"Harmony batch key must be patient_id, got {batch_key!r}. Integrating on "
        f"a specimen-level key would remove the primary-vs-recurrent difference "
        f"(RAS component T).")

    n_info = json.loads(COHORT_N.read_text())
    keep = set(n_info["patient_ids"])
    print(f"cohort_n.json: n_patients={n_info['n_patients']}")
    print(f"Harmony batch key: {batch_key}  (seed {seed}, {n_hvg} HVGs)")

    adata = ad.read_h5ad(QC_H5)
    print(f"01_qc.h5ad: {adata.n_obs:,} nuclei")
    adata = adata[adata.obs["patient_id"].isin(keep)].copy()
    print(f"after clause (d) subset: {adata.n_obs:,} nuclei, "
          f"{adata.obs['patient_id'].nunique()} patients")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=batch_key)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack", random_state=seed)

    print(f"running Harmony on {batch_key} ...")
    sce.pp.harmony_integrate(adata, key=batch_key, random_state=seed)
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", random_state=seed)

    # ---------------- Step 3 LISI gate ----------------
    emb = adata.obsm["X_pca_harmony"]
    meta = adata.obs[["patient_id", "timepoint"]].astype(str)
    lisi = compute_lisi(emb, meta, ["patient_id", "timepoint"])
    k_pat = meta["patient_id"].nunique()
    k_tp = meta["timepoint"].nunique()
    raw_pat, raw_tp = float(np.median(lisi[:, 0])), float(np.median(lisi[:, 1]))
    norm_pat = (raw_pat - 1) / (k_pat - 1)
    norm_tp = (raw_tp - 1) / (k_tp - 1)
    threshold = fail_ratio * norm_pat
    passed = norm_tp < threshold

    with open(OUT_LISI, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "n_categories", "median_lisi_raw",
                    "median_lisi_normalised"])
        w.writerow(["patient_id", k_pat, round(raw_pat, 4), round(norm_pat, 4)])
        w.writerow(["timepoint", k_tp, round(raw_tp, 4), round(norm_tp, 4)])
        w.writerow(["fail_ratio", "", "", fail_ratio])
        w.writerow(["fail_threshold_normalised", "", "", round(threshold, 4)])
        w.writerow(["gate", "", "", "PASS" if passed else "FAIL"])

    print("\n--- Step 3 LISI gate ---")
    print(f"  patient_id : raw {raw_pat:7.3f} / {k_pat:<3} -> normalised {norm_pat:.4f}"
          "   (want HIGH: patients mixed)")
    print(f"  timepoint  : raw {raw_tp:7.3f} / {k_tp:<3} -> normalised {norm_tp:.4f}"
          "   (want LOW: signal preserved)")
    print(f"  fail if normalised timepoint >= {fail_ratio} x {norm_pat:.4f} "
          f"= {threshold:.4f}")
    print(f"  -> {'PASS' if passed else 'FAIL'}")
    print(f"  wrote {OUT_LISI.relative_to(REPO)}")

    if not passed:
        print("\nGATE FAILED: the embedding has mixed the timepoints. RAS component T")
        print("is compromised. Do not proceed to RAS construction.")
        return 1

    adata.write_h5ad(OUT_H5, compression="gzip")
    print(f"\nwrote {OUT_H5.relative_to(REPO)}  "
          f"{adata.n_obs:,} nuclei x {adata.n_vars:,} HVGs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
