#!/usr/bin/env python3
"""Stage B: ElasticNet discovery on metacells, with LOPO and stability selection.

    target_i = f(X_i) + eta_i

X is metacell-level log-normalised expression. HVGs are re-derived INSIDE each
training fold (CLAUDE.md §3.6) -- never once on the full matrix, which would leak
held-out patients into feature selection.

For the Tier C-disjoint arm, every gene on chr7, chr9p or chr10 is removed from X
(arm-aware: chr9q stays eligible). Enforced by assertion here and by
tests/test_chr_disjoint.py.

X NEVER contains PCA coordinates, transport mass, or abundance statistics -- those
are RAS's own inputs (§3.6).

Stability selection: resample patients with replacement, refit, keep genes chosen
at or above the threshold. Declared deviation (§9.1): the ElasticNet alpha is
tuned once per outer fold and reused across that fold's bootstraps, which slightly
understates selection variance.

Usage:
    python3 src/06_stage_b.py --target <col> --arm <name> [--exclude-disjoint]
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from sklearn.linear_model import ElasticNet, ElasticNetCV

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _genome import annotate_var  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
QC = REPO / "data" / "processed" / "01_qc.h5ad"
ASSIGN = REPO / "results" / "tables" / "metacell_assignments.csv"
RESID = REPO / "results" / "tables" / "stage_a_residuals.csv"
CONF = REPO / "configs" / "pipeline_config.yaml"
CACHE = REPO / "data" / "interim" / "metacell_expression.h5ad"


def build_metacell_expression(conf) -> ad.AnnData:
    """Aggregate raw counts per metacell, then log-normalise. Cached."""
    if CACHE.exists():
        return ad.read_h5ad(CACHE)
    a = ad.read_h5ad(QC)
    asg = pd.read_csv(ASSIGN)
    asg = asg[asg["nucleus_id"].isin(set(a.obs_names))]
    a = a[asg["nucleus_id"].values].copy()
    codes, uniq = pd.factorize(asg["metacell_id"].values)
    import scipy.sparse as sp
    M = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                      shape=(len(uniq), len(codes)))
    X = M @ a.X                                   # summed raw counts per metacell
    mc = ad.AnnData(X)
    mc.obs_names = uniq
    mc.var_names = a.var_names
    first = asg.drop_duplicates("metacell_id").set_index("metacell_id").loc[uniq]
    mc.obs["patient_id"] = first["patient_id"].astype(str).values
    mc.obs["timepoint"] = first["timepoint"].values
    mc.obs["n_nuclei"] = np.bincount(codes, minlength=len(uniq))
    sc.pp.normalize_total(mc, target_sum=1e4)
    sc.pp.log1p(mc)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    mc.write_h5ad(CACHE, compression="gzip")
    return mc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--exclude-disjoint", action="store_true")
    ap.add_argument("--bootstraps", type=int, default=None)
    args = ap.parse_args()

    conf = yaml.safe_load(open(CONF))
    seed = conf["seed"]["master"]
    n_hvg = conf["features"]["n_hvg"]
    n_boot = args.bootstraps or conf["stability_selection"]["n_bootstraps"]
    thr = [conf["stability_selection"]["threshold_primary"]] + \
          list(conf["stability_selection"]["thresholds_sensitivity"])
    regions = conf["disjoint_set_S"]["regions"]

    mc = build_metacell_expression(conf)
    resid = pd.read_csv(RESID)
    resid["patient_id"] = resid["patient_id"].astype(str)
    asg = pd.read_csv(ASSIGN)
    y_by_mc = (resid.merge(asg, on="nucleus_id", how="inner")
                    .groupby("metacell_id")[args.target].mean())     # §10.3: mean

    mc = mc[mc.obs_names.isin(y_by_mc.index)].copy()
    y = y_by_mc.loc[mc.obs_names].values
    print(f"arm={args.arm} target={args.target}")
    print(f"metacells with a target value: {mc.n_obs:,}  genes: {mc.n_vars:,}")
    print(f"patients: {mc.obs['patient_id'].nunique()}")
    print("per-patient metacell counts:")
    print(mc.obs["patient_id"].value_counts().sort_index().to_string())

    if args.exclude_disjoint:
        ann = annotate_var(mc.var_names, regions)
        drop = ann["in_disjoint_set_S"].fillna(False).values
        q9 = ((ann["chromosome"] == "chr9") & (ann["arm"] == "q")).fillna(False).values
        assert not (drop & q9).any(), "chr9q gene marked in disjoint set"
        mc = mc[:, ~drop].copy()
        print(f"disjoint exclusion: dropped {int(drop.sum()):,} genes on "
              f"chr7/chr9p/chr10 -> {mc.n_vars:,} eligible")

    pats = np.array(sorted(mc.obs["patient_id"].unique()))
    rng = np.random.default_rng(seed)
    counts = pd.Series(0, index=mc.var_names, dtype=int)
    n_fits = 0

    for fi, p in enumerate(pats):
        te = (mc.obs["patient_id"] == p).values
        tr = ~te
        sub = mc[tr].copy()
        sc.pp.highly_variable_genes(sub, n_top_genes=min(n_hvg, sub.n_vars - 1))
        hv = sub.var_names[sub.var["highly_variable"].values]
        Xtr = np.asarray(mc[tr, hv].X.todense() if hasattr(mc.X, "todense")
                         else mc[tr, hv].X)
        mu, sd = Xtr.mean(0), Xtr.std(0); sd[sd == 0] = 1
        Xtr = (Xtr - mu) / sd
        ytr = y[tr]
        cv = ElasticNetCV(l1_ratio=0.5, n_alphas=20, cv=3, random_state=seed,
                          max_iter=5000, n_jobs=-1).fit(Xtr, ytr)
        alpha = cv.alpha_                      # tuned once per outer fold (§9.1)
        tr_pat = mc.obs["patient_id"].values[tr]
        for b in range(n_boot):
            r = np.random.default_rng(seed + b)
            boot = r.choice(np.unique(tr_pat), size=len(np.unique(tr_pat)),
                            replace=True)
            idx = np.concatenate([np.where(tr_pat == q)[0] for q in boot])
            m = ElasticNet(alpha=alpha, l1_ratio=0.5, max_iter=5000,
                           warm_start=True, random_state=seed).fit(Xtr[idx], ytr[idx])
            counts[hv[np.abs(m.coef_) > 0]] += 1
            n_fits += 1
        print(f"  fold {fi+1}/{len(pats)} held out patient {p}: alpha={alpha:.5f}")

    freq = (counts / n_fits).sort_values(ascending=False)
    outdir = REPO / "results" / "gene_lists"
    outdir.mkdir(parents=True, exist_ok=True)
    freq.to_frame("selection_frequency").to_csv(outdir / f"{args.arm}_frequencies.csv")
    print(f"\nselection frequencies over {n_fits} fits")
    for t in thr:
        sel = freq[freq >= t].index.tolist()
        pd.Series(sel).to_csv(outdir / f"{args.arm}_genes_{int(t*100)}.csv",
                              index=False, header=["gene"])
        print(f"  threshold {int(t*100):>2}%: {len(sel):>5} genes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
