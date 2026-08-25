#!/usr/bin/env python3
"""Stage B for all six arms (3 H3 variants x 2 tiers), LOPO + stability selection.

X is metacell-level log-normalised expression. HVGs are re-derived INSIDE each
training fold (CLAUDE.md §3.6). For every Tier C arm, all chr7 / chr9p / chr10
genes are removed from X (arm-aware; chr9q stays eligible) and asserted.

X never contains PCA coordinates, transport mass or abundance statistics.

Declared deviation §9.1: ElasticNet alpha is tuned once per outer fold and reused
across that fold's bootstraps with warm starts.

Per-fold results are written individually -- §4.6 -- not collapsed to a mean.
"""
from __future__ import annotations
import argparse, json, sys, warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, yaml
from sklearn.linear_model import ElasticNet, ElasticNetCV
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _genome import annotate_var  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
QC = REPO / "data" / "processed" / "01_qc.h5ad"
ASSIGN = REPO / "results" / "tables" / "metacell_assignments.csv"
TARGETS = REPO / "results" / "tables" / "stage_b_targets.csv"
CONF = REPO / "configs" / "pipeline_config.yaml"
CACHE = REPO / "data" / "interim" / "metacell_expression.h5ad"
GL = REPO / "results" / "gene_lists"
FOLDS = REPO / "results" / "tables" / "stage_b_folds.csv"

ARMS = [("v1_tierA", False), ("v1_tierC", True),
        ("v2_tierA", False), ("v2_tierC", True),
        ("v3_tierA", False), ("v3_tierC", True)]


def build_cache(conf):
    if CACHE.exists():
        return
    a = ad.read_h5ad(QC)
    asg = pd.read_csv(ASSIGN)
    asg = asg[asg["nucleus_id"].isin(set(a.obs_names))]
    a = a[asg["nucleus_id"].values].copy()
    codes, uniq = pd.factorize(asg["metacell_id"].values)
    M = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                      shape=(len(uniq), len(codes)))
    mc = ad.AnnData(sp.csr_matrix(M @ a.X))
    mc.obs_names, mc.var_names = uniq, a.var_names
    f = asg.drop_duplicates("metacell_id").set_index("metacell_id").loc[uniq]
    mc.obs["patient_id"] = f["patient_id"].astype(str).values
    mc.obs["timepoint"] = f["timepoint"].values
    mc.obs["n_nuclei"] = np.bincount(codes, minlength=len(uniq))
    sc.pp.normalize_total(mc, target_sum=1e4); sc.pp.log1p(mc)
    mc.write_h5ad(CACHE, compression="gzip")


def run_fold(arm, exclude, p, seed, n_hvg, n_boot):
    mc = ad.read_h5ad(CACHE)
    tg = pd.read_csv(TARGETS)
    asg = pd.read_csv(ASSIGN)
    y_by = (tg.merge(asg, on="nucleus_id", how="inner")
              .groupby("metacell_id")[arm].mean())
    mc = mc[mc.obs_names.isin(y_by.index)].copy()
    y = y_by.loc[mc.obs_names].values
    if exclude:
        ann = annotate_var(mc.var_names, yaml.safe_load(open(CONF))["disjoint_set_S"]["regions"])
        drop = ann["in_disjoint_set_S"].fillna(False).values
        q9 = ((ann["chromosome"] == "chr9") & (ann["arm"] == "q")).fillna(False).values
        assert not (drop & q9).any(), "chr9q gene marked in disjoint set"
        mc = mc[:, ~drop].copy()
    te = (mc.obs["patient_id"] == p).values
    tr = ~te
    sub = mc[tr].copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=min(n_hvg, sub.n_vars - 1))
    hv = sub.var_names[sub.var["highly_variable"].values]
    X = mc[tr, hv].X
    X = np.asarray(X.todense() if sp.issparse(X) else X, dtype=np.float64)
    mu, sd = X.mean(0), X.std(0); sd[sd == 0] = 1
    X = (X - mu) / sd
    yt = y[tr]
    cv = ElasticNetCV(l1_ratio=0.5, n_alphas=10, cv=3, random_state=seed,
                      max_iter=2000, tol=1e-3, n_jobs=1).fit(X, yt)
    alpha = float(cv.alpha_)
    tp = mc.obs["patient_id"].values[tr]
    uniq_p = np.unique(tp)
    cnt = {}
    m = ElasticNet(alpha=alpha, l1_ratio=0.5, max_iter=2000, tol=1e-3,
                   warm_start=True, random_state=seed)
    for b in range(n_boot):
        r = np.random.default_rng(seed + b)
        boot = r.choice(uniq_p, size=len(uniq_p), replace=True)
        idx = np.concatenate([np.where(tp == q)[0] for q in boot])
        m.fit(X[idx], yt[idx])
        for g in hv[np.abs(m.coef_) > 0]:
            cnt[g] = cnt.get(g, 0) + 1
    return {"arm": arm, "heldout_patient": p, "alpha": alpha,
            "n_train_metacells": int(tr.sum()), "n_test_metacells": int(te.sum()),
            "n_features": int(len(hv)), "n_selected_any": len(cnt),
            "counts": cnt, "n_fits": n_boot}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--bootstraps", type=int, default=None)
    args = ap.parse_args()
    conf = yaml.safe_load(open(CONF))
    seed, n_hvg = conf["seed"]["master"], conf["features"]["n_hvg"]
    n_boot = args.bootstraps or conf["stability_selection"]["n_bootstraps"]
    thr = [conf["stability_selection"]["threshold_primary"]] + \
          list(conf["stability_selection"]["thresholds_sensitivity"])
    build_cache(conf)
    mc = ad.read_h5ad(CACHE)
    tg = pd.read_csv(TARGETS); asg = pd.read_csv(ASSIGN)
    y_by = (tg.merge(asg, on="nucleus_id", how="inner")
              .groupby("metacell_id")["v1_tierA"].mean())
    usable = mc.obs_names.isin(y_by.index)
    pats = sorted(mc.obs["patient_id"][usable].astype(str).unique())
    print(f"Stage B units: {int(usable.sum()):,} metacells with a target "
          f"(of {mc.n_obs:,} total), {len(pats)} patients, {n_boot} bootstraps")
    GL.mkdir(parents=True, exist_ok=True)

    jobs = [(a, e, p) for a, e in ARMS for p in pats]
    print(f"jobs: {len(jobs)} (6 arms x {len(pats)} folds)")
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_fold, a, e, p, seed, n_hvg, n_boot): (a, p)
                for a, e, p in jobs}
        for f in as_completed(futs):
            results.append(f.result())
            print(f"  {len(results)}/{len(jobs)} done", flush=True)

    rows = [{k: v for k, v in r.items() if k != "counts"} for r in results]
    pd.DataFrame(rows).sort_values(["arm", "heldout_patient"]).to_csv(FOLDS, index=False)
    print(f"\nwrote {FOLDS.relative_to(REPO)} (per-fold, not collapsed -- §4.6)")

    print(f"\n{'arm':<12}{'30%':>8}{'50%':>8}{'80%':>8}")
    summary = {}
    for arm, _ in ARMS:
        rs = [r for r in results if r["arm"] == arm]
        tot = sum(r["n_fits"] for r in rs)
        c = {}
        for r in rs:
            for g, n in r["counts"].items():
                c[g] = c.get(g, 0) + n
        freq = pd.Series(c, dtype=float).sort_values(ascending=False) / tot
        freq.to_frame("selection_frequency").to_csv(GL / f"{arm}_frequencies.csv")
        line = {}
        for t in thr:
            sel = freq[freq >= t].index.tolist()
            pd.Series(sel, dtype=object).to_csv(
                GL / f"{arm}_genes_{int(t*100)}.csv", index=False, header=["gene"])
            line[t] = len(sel)
        summary[arm] = line
        print(f"{arm:<12}{line[0.30]:>8}{line[0.50]:>8}{line[0.80]:>8}")
    json.dump(summary, open(REPO / "results" / "tables" / "stage_b_summary.json", "w"),
              indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
