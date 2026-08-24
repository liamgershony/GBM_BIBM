#!/usr/bin/env python3
"""Is the clone structure real, or is Leiden partitioning noise?

02b reported 8-17 clones for EVERY region including chr9p, which carries a focal
~1 Mb event against a 23-window block. Near-identical clone counts from regions
with very different expected signal is the signature of a clustering artifact, not
of clonal structure. Two direct tests, neither of which depends on clustering:

1. SIGNAL MAGNITUDE. Mean inferCNV value per region, malignant vs reference. The
   canonical GBM events are chr7 GAIN (positive) and chr10 LOSS (negative). If
   those signs do not appear, the CNV signal itself is absent and no clustering of
   it can be meaningful.

2. PERMUTATION NULL for the clone count. Shuffle each CNV window independently
   across cells, destroying cell-cell covariance while preserving each window's
   marginal distribution, then cluster identically. If the null yields a similar
   clone count, the observed count carries no information.

Report only.
"""

from __future__ import annotations

import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
IN_DIR = REPO / "data" / "interim" / "cnv_input"
CONF = REPO / "configs" / "pipeline_config.yaml"
COHORT_N = REPO / "results" / "tables" / "cohort_n.json"
OUT = REPO / "results" / "tables" / "clone_validity_check.csv"

WINDOW_SIZE, LEIDEN_RESOLUTION = 100, 1.0
LABEL = {"chr7": "chr7", "chr9": "chr9p", "chr10": "chr10"}


def _n_clones(M, seed):
    if M.shape[0] < 10 or M.shape[1] == 0:
        return 1
    t = ad.AnnData(np.asarray(M, dtype="float32"))
    n = int(min(20, max(2, min(t.n_obs, t.n_vars) - 1)))
    sc.tl.pca(t, n_comps=n, svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(t, random_state=seed)
    sc.tl.leiden(t, resolution=LEIDEN_RESOLUTION, key_added="c", random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    return int(t.obs["c"].nunique())


def check(pid, seed):
    import infercnvpy as icnv
    a = ad.read_h5ad(IN_DIR / f"{pid}.h5ad")
    icnv.tl.infercnv(a, reference_key="cnv_reference", reference_cat=["Normal"],
                     window_size=WINDOW_SIZE)
    X = a.obsm["X_cnv"]
    X = np.asarray(X.todense() if hasattr(X, "todense") else X)
    mal = (a.obs["cnv_reference"] == "Malignant").values
    ref = ~mal
    cp = sorted(dict(a.uns["cnv"]["chr_pos"]).items(), key=lambda kv: kv[1])
    out = {"patient_id": pid}
    rng = np.random.default_rng(seed)
    for i, (chrom, s) in enumerate(cp):
        e = cp[i + 1][1] if i + 1 < len(cp) else X.shape[1]
        lab = LABEL.get(chrom)
        if not lab:
            continue
        out[f"{lab}_mean_malignant"] = float(X[mal, s:e].mean())
        out[f"{lab}_mean_reference"] = float(X[ref, s:e].mean())
        out[f"{lab}_delta"] = out[f"{lab}_mean_malignant"] - out[f"{lab}_mean_reference"]
        Mreal = X[mal, s:e]
        out[f"{lab}_clones_observed"] = _n_clones(Mreal, seed)
        Mnull = np.column_stack([rng.permutation(Mreal[:, j])
                                 for j in range(Mreal.shape[1])])
        out[f"{lab}_clones_null"] = _n_clones(Mnull, seed)
    return out


def main() -> int:
    conf = yaml.safe_load(open(CONF))
    seed = conf["seed"]["master"]
    import json
    pts = json.loads(COHORT_N.read_text())["patient_ids"]
    rows = []
    with ProcessPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(check, p, seed): p for p in pts}
        for f in as_completed(futs):
            rows.append(f.result())
            print(f"  {len(rows)}/{len(pts)} done", flush=True)
    df = pd.DataFrame(rows).sort_values("patient_id",
                                        key=lambda s: s.astype(str).str.zfill(6))
    df.to_csv(OUT, index=False)
    print(f"\nwrote {OUT.relative_to(REPO)}\n")

    print("1. SIGNAL MAGNITUDE  (expect chr7 delta > 0, chr10 delta < 0)")
    print(f"   {'region':<8}{'mean malignant':>16}{'mean reference':>16}"
          f"{'delta':>10}{'patients w/ expected sign':>28}")
    exp = {"chr7": 1, "chr9p": -1, "chr10": -1}
    for r in ("chr7", "chr9p", "chr10"):
        d = df[f"{r}_delta"]
        ok = (d > 0).sum() if exp[r] > 0 else (d < 0).sum()
        print(f"   {r:<8}{df[f'{r}_mean_malignant'].mean():>16.5f}"
              f"{df[f'{r}_mean_reference'].mean():>16.5f}{d.mean():>10.5f}"
              f"{ok:>20}/{len(df)}")

    print("\n2. CLONE COUNT vs PERMUTATION NULL")
    print(f"   {'region':<8}{'observed (med)':>16}{'null (med)':>14}{'ratio':>9}"
          f"   interpretation")
    for r in ("chr7", "chr9p", "chr10"):
        o, n = df[f"{r}_clones_observed"].median(), df[f"{r}_clones_null"].median()
        ratio = o / n if n else float("nan")
        verdict = ("clone count is INDISTINGUISHABLE from noise"
                   if ratio < 1.5 else "clone count exceeds the null")
        print(f"   {r:<8}{o:>16.1f}{n:>14.1f}{ratio:>9.2f}   {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
