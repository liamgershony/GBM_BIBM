#!/usr/bin/env python3
"""Per-patient inferCNV on chr7/chr9p/chr10 -> clone_catalog.csv, 05_genotype.h5ad.

Clone identity G and clone abundance Ab(clone) for Tier C-disjoint are computed
from inferCNV signal restricted to disjoint_set_S (CLAUDE.md §3.4). Genes outside
S are excluded from the CNV input entirely, not merely from clustering, so no gene
eligible for Stage B selection contributes expression to the target.

chr9p is an ARM: chr9q genes stay eligible and are excluded from S here.

Nuclei per patient are capped at infercnv.max_nuclei_per_patient from the frozen
config (malignant only); ALL non-malignant nuclei of that patient serve as the
inferCNV reference. `unknown` nuclei are used for neither.

Clone calling has no cross-patient coupling, so patients run in parallel.

Per-patient convergence status and reference-nuclei count are written as COLUMNS
of clone_catalog.csv, the same way scrublet_status is recorded in
qc_per_library.csv: a silent failure must be visible in the artifact, not only in
stdout.
"""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _genome import annotate_var  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
QC_H5 = REPO / "data" / "processed" / "01_qc.h5ad"
COHORT_N = REPO / "results" / "tables" / "cohort_n.json"
CONF = REPO / "configs" / "pipeline_config.yaml"
IN_DIR = REPO / "data" / "interim" / "cnv_input"
OUT_DIR = REPO / "data" / "interim" / "cnv_out"
CATALOG = REPO / "results" / "tables" / "clone_catalog.csv"
ASSIGN = REPO / "results" / "tables" / "clone_assignments.csv"
OUT_H5 = REPO / "data" / "processed" / "05_genotype.h5ad"

# Implementation parameters, not in the frozen config. Recorded in the catalog so
# every reported clone count carries the settings that produced it.
WINDOW_SIZE = 100        # infercnvpy default
LEIDEN_RESOLUTION = 1.0  # scanpy default
MIN_REFERENCE_NUCLEI = 20

REGION_LABEL = {"chr7": "chr7", "chr9": "chr9p", "chr10": "chr10"}


def _leiden_on(matrix: np.ndarray, seed: int, tag: str) -> int:
    """Cluster one CNV block and return the number of clones it resolves."""
    if matrix.shape[1] == 0:
        return 0
    if matrix.shape[0] < 10:
        return 1
    tmp = ad.AnnData(np.asarray(matrix, dtype="float32"))
    n_comps = int(min(20, max(2, min(tmp.n_obs, tmp.n_vars) - 1)))
    sc.tl.pca(tmp, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(tmp, random_state=seed)
    sc.tl.leiden(tmp, resolution=LEIDEN_RESOLUTION, key_added="c",
                 random_state=seed, flavor="igraph", n_iterations=2,
                 directed=False)
    return int(tmp.obs["c"].nunique())


def run_patient(patient_id: str, seed: int) -> dict:
    """inferCNV + clone calling for one patient. Never raises: failures are
    returned as a status string so they land in the artifact."""
    import infercnvpy as icnv

    rec = {"patient_id": patient_id, "status": "not_run",
           "n_malignant_used": 0, "n_reference_nuclei": 0,
           "n_genes_in_S": 0, "n_cnv_windows": 0,
           "n_clones_combined": 0, "n_clones_chr7": 0,
           "n_clones_chr9p": 0, "n_clones_chr10": 0,
           "window_size": WINDOW_SIZE, "leiden_resolution": LEIDEN_RESOLUTION,
           "error": ""}
    try:
        a = ad.read_h5ad(IN_DIR / f"{patient_id}.h5ad")
        n_ref = int((a.obs["cnv_reference"] == "Normal").sum())
        n_mal = int((a.obs["cnv_reference"] == "Malignant").sum())
        rec.update(n_malignant_used=n_mal, n_reference_nuclei=n_ref,
                   n_genes_in_S=int(a.n_vars))

        if n_ref < MIN_REFERENCE_NUCLEI:
            rec["status"] = f"skipped_insufficient_reference(<{MIN_REFERENCE_NUCLEI})"
            return rec
        if n_mal < 10:
            rec["status"] = "skipped_insufficient_malignant(<10)"
            return rec

        icnv.tl.infercnv(a, reference_key="cnv_reference",
                         reference_cat=["Normal"], window_size=WINDOW_SIZE)
        X = np.asarray(a.obsm["X_cnv"].todense()
                       if hasattr(a.obsm["X_cnv"], "todense") else a.obsm["X_cnv"])
        rec["n_cnv_windows"] = int(X.shape[1])

        mal = (a.obs["cnv_reference"] == "Malignant").values
        Xm = X[mal]

        # clones from the whole of S (this is the operative clone_id)
        rec["n_clones_combined"] = _leiden_on(Xm, seed, "combined")

        # CLAUDE.md §10.2: each region checked INDEPENDENTLY, so a region that
        # contributes nothing cannot hide behind the combined result.
        chr_pos = dict(a.uns["cnv"]["chr_pos"])
        order = sorted(chr_pos.items(), key=lambda kv: kv[1])
        bounds = {}
        for i, (chrom, start) in enumerate(order):
            end = order[i + 1][1] if i + 1 < len(order) else X.shape[1]
            bounds[chrom] = (start, end)
        for chrom, (s, e) in bounds.items():
            label = REGION_LABEL.get(chrom)
            if label:
                rec[f"n_clones_{label}"] = _leiden_on(Xm[:, s:e], seed, label)
                rec[f"n_windows_{label}"] = int(e - s)

        # persist per-nucleus clone ids from the combined run
        tmp = ad.AnnData(np.asarray(Xm, dtype="float32"))
        n_comps = int(min(20, max(2, min(tmp.n_obs, tmp.n_vars) - 1)))
        sc.tl.pca(tmp, n_comps=n_comps, svd_solver="arpack", random_state=seed)
        sc.pp.neighbors(tmp, random_state=seed)
        sc.tl.leiden(tmp, resolution=LEIDEN_RESOLUTION, key_added="clone",
                     random_state=seed, flavor="igraph", n_iterations=2,
                     directed=False)
        pd.DataFrame({"nucleus_id": a.obs_names[mal],
                      "patient_id": patient_id,
                      "clone_id": [f"{patient_id}_c{c}" for c in tmp.obs["clone"]],
                      "timepoint": a.obs["timepoint"].values[mal]}
                     ).to_csv(OUT_DIR / f"{patient_id}_clones.csv", index=False)
        rec["status"] = "converged"
    except Exception as e:                                  # noqa: BLE001
        # Deliberately broad HERE ONLY: a worker must not take down the pool, and
        # the failure is recorded in the artifact rather than swallowed.
        rec["status"] = f"failed:{type(e).__name__}"
        rec["error"] = f"{e}"[:300]
        traceback.print_exc()
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--prepare-only", action="store_true")
    args = ap.parse_args()

    conf = yaml.safe_load(open(CONF))
    seed = conf["seed"]["master"]
    cap = conf["infercnv"]["max_nuclei_per_patient"]
    regions = conf["disjoint_set_S"]["regions"]
    patients = json.loads(COHORT_N.read_text())["patient_ids"]
    print(f"patients {len(patients)}  cap {cap} malignant/patient  seed {seed}")
    print(f"disjoint_set_S: {regions}")

    IN_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not all((IN_DIR / f"{p}.h5ad").exists() for p in patients):
        adata = ad.read_h5ad(QC_H5)
        ann = annotate_var(adata.var_names, regions)
        keep = (ann["in_disjoint_set_S"].values & ann["chromosome"].notna().values)
        print(f"genes: {adata.n_vars:,} -> {int(keep.sum()):,} inside disjoint_set_S")
        assert int(((ann["chromosome"] == "chr9") & (ann["arm"] == "q")
                    & ann["in_disjoint_set_S"]).sum()) == 0, \
            "chr9q gene marked in_disjoint_set_S -- arm exclusion is broken"

        adata = adata[:, keep].copy()
        for c in ("chromosome", "start", "end", "arm"):
            adata.var[c] = ann.loc[adata.var_names, c].values
        adata.var["start"] = adata.var["start"].astype(int)
        adata.var["end"] = adata.var["end"].astype(int)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        rng = np.random.default_rng(seed)
        for p in patients:
            sub = adata[adata.obs["patient_id"] == p]
            mal = np.where(sub.obs["is_malignant"].values)[0]
            ref = np.where(sub.obs["is_reference_normal"].values)[0]
            if len(mal) > cap:
                mal = rng.choice(mal, size=cap, replace=False)
            idx = np.sort(np.concatenate([mal, ref]))
            s = sub[idx].copy()
            s.obs["cnv_reference"] = np.where(
                s.obs["is_malignant"].values, "Malignant",
                np.where(s.obs["is_reference_normal"].values, "Normal", "unknown"))
            s = s[s.obs["cnv_reference"] != "unknown"].copy()
            s.obs["cnv_reference"] = s.obs["cnv_reference"].astype("category")
            s.write_h5ad(IN_DIR / f"{p}.h5ad")
            print(f"  {p:<5} malignant {int((s.obs['cnv_reference']=='Malignant').sum()):>5,}"
                  f"  reference {int((s.obs['cnv_reference']=='Normal').sum()):>5,}")
        del adata
    else:
        print("per-patient inputs already prepared")

    if args.prepare_only:
        return 0

    print(f"\nrunning inferCNV across {args.workers} workers ...")
    recs = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_patient, p, seed): p for p in patients}
        for fut in as_completed(futs):
            r = fut.result()
            recs.append(r)
            print(f"  [{len(recs):>2}/{len(patients)}] {r['patient_id']:<5} "
                  f"{r['status']:<40} ref={r['n_reference_nuclei']:>5,} "
                  f"mal={r['n_malignant_used']:>5,} clones={r['n_clones_combined']}")

    cat = pd.DataFrame(recs).sort_values(
        "patient_id", key=lambda s: s.astype(str).str.zfill(6))
    cat.to_csv(CATALOG, index=False)
    print(f"\nwrote {CATALOG.relative_to(REPO)}")

    parts = [pd.read_csv(f) for f in sorted(OUT_DIR.glob("*_clones.csv"))]
    if parts:
        pd.concat(parts).to_csv(ASSIGN, index=False)
        print(f"wrote {ASSIGN.relative_to(REPO)} "
              f"({sum(len(p) for p in parts):,} nuclei)")

    conv = (cat["status"] == "converged").sum()
    print(f"\nconverged {conv}/{len(cat)} patients")
    if conv == 0:
        print("FATAL: inferCNV converged for zero patients.")
        return 1
    bad = cat[cat["status"] != "converged"]
    if len(bad):
        print("non-converged patients (recorded in clone_catalog.csv):")
        for _, r in bad.iterrows():
            print(f"   {r['patient_id']:<5} {r['status']}  {r['error'][:80]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
