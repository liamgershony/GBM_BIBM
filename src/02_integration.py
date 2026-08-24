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

# Pre-specified in DEVIATIONS.md before any permutation value was computed.
N_PERMUTATIONS = 3
TIMEPOINT_BELOW_NULL = 0.95   # timepoint "below its null" iff ratio < this
PATIENT_NEAR_NULL = 0.80      # patient "near its null" iff ratio >= this


def main() -> int:
    conf = yaml.safe_load(open(CONF))
    rconf = yaml.safe_load(open(RCONF))
    batch_key = conf["integration"]["batch_key"]
    seed = conf["seed"]["master"]
    n_hvg = conf["features"]["n_hvg"]


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

    # ---------------- Step 3 LISI gate (permutation null) ----------------
    # Each label is compared against ITS OWN permutation null rather than against
    # the other label. The earlier (LISI-1)/(k-1) normalisation was invalid:
    # compute_lisi uses perplexity=30, so attainable LISI is bounded by
    # neighbourhood size as well as by category count, and dividing by (k-1)
    # inverted the comparison. See DEVIATIONS.md, 2026-08-24.
    emb = adata.obsm["X_pca_harmony"]
    rng = np.random.default_rng(seed)
    meta = adata.obs[["patient_id", "timepoint"]].astype(str).copy()

    labels = ["patient_id", "timepoint"]
    perm_cols = {l: [] for l in labels}
    for i in range(N_PERMUTATIONS):
        r = np.random.default_rng(seed + i)
        for l in labels:
            col = f"{l}__perm{i}"
            meta[col] = r.permutation(meta[l].values)
            perm_cols[l].append(col)

    # ONE compute_lisi call: the neighbourhood is built once and shared by the
    # observed and permuted columns, so the null differs only in the labels.
    all_cols = labels + [c for l in labels for c in perm_cols[l]]
    print(f"computing LISI for {len(all_cols)} label columns "
          f"({N_PERMUTATIONS} permutations per label) ...")
    lisi = compute_lisi(emb, meta, all_cols)
    med = {c: float(np.median(lisi[:, i])) for i, c in enumerate(all_cols)}

    results = {}
    for l in labels:
        obs = med[l]
        null = float(np.mean([med[c] for c in perm_cols[l]]))
        results[l] = {"observed": obs, "null": null,
                      "ratio": obs / null if null else float("nan"),
                      "k": int(meta[l].nunique())}

    r_tp = results["timepoint"]["ratio"]
    r_pat = results["patient_id"]["ratio"]

    if r_tp >= TIMEPOINT_BELOW_NULL:
        outcome, passed = "c", False
        verdict = ("REAL FAILURE -- timepoint LISI is at or above its null: the "
                   "embedding has mixed primary and recurrent. RAS component T is "
                   "compromised. Invoking STOP/GO gate 1.")
    elif r_pat >= PATIENT_NEAR_NULL:
        outcome, passed = "a", True
        verdict = "PASS -- timepoints separable, patients well mixed."
    else:
        outcome, passed = "b", True
        verdict = ("PROCEED WITH STATED LIMITATION -- timepoints are separable "
                   "(component T intact), but integration under-corrected across "
                   "patients. Report these LISI values as a limitation. Do NOT "
                   "retune Harmony theta to force mixing.")

    with open(OUT_LISI, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "k", "median_lisi_observed", "median_lisi_null",
                    "ratio_observed_over_null"])
        for l in labels:
            v = results[l]
            w.writerow([l, v["k"], round(v["observed"], 4), round(v["null"], 4),
                        round(v["ratio"], 4)])
        w.writerow(["n_permutations", "", "", "", N_PERMUTATIONS])
        w.writerow(["cutoff_timepoint_below_null", "", "", "", TIMEPOINT_BELOW_NULL])
        w.writerow(["cutoff_patient_near_null", "", "", "", PATIENT_NEAR_NULL])
        w.writerow(["outcome", "", "", "", outcome])
        w.writerow(["gate", "", "", "", "PASS" if passed else "FAIL"])

    print("\n--- Step 3 LISI gate (permutation null) ---")
    for l in labels:
        v = results[l]
        print(f"  {l:<11} k={v['k']:<3} observed {v['observed']:7.3f}   "
              f"null {v['null']:7.3f}   ratio {v['ratio']:.4f}")
    print(f"  cutoffs: timepoint below null if ratio < {TIMEPOINT_BELOW_NULL}; "
          f"patient near null if ratio >= {PATIENT_NEAR_NULL}")
    print(f"  -> pre-specified outcome ({outcome}): {verdict}")
    print(f"  wrote {OUT_LISI.relative_to(REPO)}")

    # Written regardless of outcome so a failure can be re-diagnosed without
    # re-running Harmony. Writing the embedding is not proceeding with RAS.
    adata.write_h5ad(OUT_H5, compression="gzip")
    print(f"  wrote {OUT_H5.relative_to(REPO)}  "
          f"{adata.n_obs:,} nuclei x {adata.n_vars:,} HVGs")

    if not passed:
        print("\nSTOP/GO gate 1 invoked. Do not proceed to RAS construction.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
