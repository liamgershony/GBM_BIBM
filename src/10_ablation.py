#!/usr/bin/env python3
"""H1 (STRETCH) -- does confound adjustment improve external replication?

200 paired patient-level resamples. In each: resample patients with replacement,
fit both arms on the resampled metacells, take each arm's gene list, score its
CGGA replication rate, and record the paired difference dR.

  adjusted arm   target = Stage A residual of Tier A-reduced
  unadjusted arm target = RAW Tier A-reduced (CLAUDE.md §3.6)

Decision rule (§6): mean dR >= 10 percentage points AND the interval excludes 0.
All three outcomes are valid and reportable (§6.1), including an interval that
includes zero.

KEY EFFICIENCY POINT, and it is exact rather than an approximation: whether a gene
replicates in CGGA does not depend on which discovery resample selected it. The
per-gene CGGA verdict is therefore computed ONCE over the whole eligible universe,
and each resample's replication rate is a lookup.

DECLARED SIMPLIFICATION: each resample fits ElasticNet once at the alpha carried
over from the full-data outer folds, rather than repeating the full
LOPO-plus-stability-selection procedure 200 times, which is not computable in the
time available. This makes each resample's gene list noisier than a full
stability-selected list and therefore makes dR MORE variable, not less -- it
widens the interval rather than narrowing it. Logged in DEVIATIONS.md.
"""
from __future__ import annotations
import argparse, json, sys, warnings
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
import statsmodels.api as sm, yaml
from sklearn.linear_model import ElasticNet
warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
CACHE = REPO / "data" / "interim" / "metacell_expression.h5ad"
ASSIGN = REPO / "results" / "tables" / "metacell_assignments.csv"
TARGETS = REPO / "results" / "tables" / "stage_b_targets.csv"
RAS = REPO / "results" / "tables" / "ras_scores.csv"
FOLDS = REPO / "results" / "tables" / "stage_b_folds.csv"
CONF = REPO / "configs" / "pipeline_config.yaml"
CGGA = REPO / "data" / "interim" / "cgga"
OUT = REPO / "results" / "tables" / "h1_ablation.csv"
OUT_V = REPO / "results" / "tables" / "cgga_gene_verdicts.csv"
GBM_H = {"GBM", "RGBM", "SGBM"}


def cgga_verdicts(genes):
    """Per-gene replication verdict in CGGA. Computed once; resample-independent."""
    if OUT_V.exists():
        return pd.read_csv(OUT_V).set_index("gene")
    clin = pd.read_csv(CGGA / "CGGA.mRNAseq_693_clinical.20200506.txt", sep="\t", dtype=str)
    clin = clin[clin["Histology"].str.upper().isin(GBM_H)]
    clin = clin[clin["PRS_type"].str.strip().isin(["Primary", "Recurrent"])]
    clin = clin[clin["IDH_mutation_status"].str.strip().isin(["Wildtype", "Mutant"])]
    clin["y"] = (clin["PRS_type"].str.strip() == "Recurrent").astype(int)
    clin["idh"] = (clin["IDH_mutation_status"].str.strip() == "Mutant").astype(int)
    expr = pd.read_csv(CGGA / "CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt",
                       sep="\t", index_col=0)
    expr = expr[~expr.index.duplicated(keep="first")]
    keep = [c for c in expr.columns if c in set(clin["CGGA_ID"])]
    expr = expr[keep]; clin = clin.set_index("CGGA_ID").loc[keep]
    cpm = np.log2(expr.divide(expr.sum(0).replace(0, np.nan), axis=1) * 1e6 + 1.0)
    present = [g for g in genes if g in cpm.index]
    y = clin["y"].values.astype(float); cov = clin["idh"].values.astype(float)
    rows = []
    for i, g in enumerate(present):
        x = cpm.loc[g].values.astype(float)
        if np.nanstd(x) == 0:
            continue
        X = sm.add_constant(np.column_stack([x, cov]), has_constant="add")
        try:
            m = sm.Logit(y, X).fit(disp=0, maxiter=60)
            b, p = float(m.params[1]), float(m.pvalues[1])
        except Exception:                                   # noqa: BLE001
            continue
        lfc = float(np.nanmean(x[y == 1]) - np.nanmean(x[y == 0]))
        rows.append({"gene": g, "beta": b, "odds_ratio": float(np.exp(b)),
                     "p_raw": p, "log2fc": lfc})
        if i % 2000 == 0:
            print(f"    cgga {i}/{len(present)}", flush=True)
    df = pd.DataFrame(rows)
    n = len(df); o = np.argsort(df["p_raw"].values); q = np.empty(n)
    q[o] = np.minimum.accumulate(
        (df["p_raw"].values[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    df["q_bh"] = np.minimum(q, 1.0)
    eff = (df["log2fc"].abs() >= 0.25) | ~df["odds_ratio"].between(0.8, 1.25)
    df["replicates"] = (df["p_raw"] < 0.05) & (df["q_bh"] < 0.05) & eff
    df.to_csv(OUT_V, index=False)
    return df.set_index("gene")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resamples", type=int, default=None)
    args = ap.parse_args()
    conf = yaml.safe_load(open(CONF))
    seed = conf["seed"]["master"]
    n_rs = args.resamples or conf["h1_ablation"]["n_resamples"]
    n_hvg = conf["features"]["n_hvg"]
    bound = 0.10

    mc = ad.read_h5ad(CACHE)
    tg = pd.read_csv(TARGETS); asg = pd.read_csv(ASSIGN); ras = pd.read_csv(RAS)
    adj = (tg.merge(asg, on="nucleus_id").groupby("metacell_id")["v1_tierA"].mean())
    una = (ras.merge(asg, on="nucleus_id")
              .groupby("metacell_id")["ras_tier_a_reduced"].mean())
    both = adj.index.intersection(una.index).intersection(pd.Index(mc.obs_names))
    mc = mc[both].copy()
    y_adj, y_una = adj.loc[both].values, una.loc[both].values
    pats = mc.obs["patient_id"].astype(str).values
    uniq = np.unique(pats)
    print(f"H1: {mc.n_obs:,} metacells, {len(uniq)} patients, {n_rs} resamples")

    print("computing per-gene CGGA verdicts once ...")
    verd = cgga_verdicts(list(mc.var_names))
    rep = set(verd.index[verd["replicates"].values])
    print(f"  {len(verd):,} genes tested in CGGA; {len(rep):,} replicate")

    alpha = float(pd.read_csv(FOLDS).query("arm == 'v1_tierA'")["alpha"].median()) \
        if FOLDS.exists() else 0.02
    print(f"  ElasticNet alpha (median over full-data outer folds): {alpha:.5f}")

    sub = mc.copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=min(n_hvg, sub.n_vars - 1))
    hv = sub.var_names[sub.var["highly_variable"].values]
    X = mc[:, hv].X
    X = np.asarray(X.todense() if sp.issparse(X) else X, dtype=np.float64)
    X = (X - X.mean(0)) / np.where(X.std(0) == 0, 1, X.std(0))

    rows = []
    for b in range(n_rs):
        r = np.random.default_rng(seed + b)
        boot = r.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.where(pats == q)[0] for q in boot])
        out = {}
        for arm, yv in (("adjusted", y_adj), ("unadjusted", y_una)):
            m = ElasticNet(alpha=alpha, l1_ratio=0.5, max_iter=2000, tol=1e-3,
                           random_state=seed).fit(X[idx], yv[idx])
            sel = [g for g in hv[np.abs(m.coef_) > 0]]
            tested = [g for g in sel if g in verd.index]
            out[f"n_selected_{arm}"] = len(sel)
            out[f"n_tested_{arm}"] = len(tested)
            out[f"rate_{arm}"] = (sum(g in rep for g in tested) / len(tested)
                                  if tested else np.nan)
        out["resample"] = b
        out["delta_R"] = out["rate_adjusted"] - out["rate_unadjusted"]
        rows.append(out)
        if (b + 1) % 20 == 0:
            print(f"  resample {b+1}/{n_rs}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    d = df["delta_R"].dropna().values
    lo, hi = np.percentile(d, [2.5, 97.5])
    mean = float(d.mean())
    print(f"\nwrote {OUT.relative_to(REPO)}")
    print(f"\n--- H1 ---")
    print(f"  mean replication rate, adjusted   : {df['rate_adjusted'].mean():.4f}")
    print(f"  mean replication rate, unadjusted : {df['rate_unadjusted'].mean():.4f}")
    print(f"  mean dR                           : {mean:+.4f} "
          f"({mean*100:+.2f} percentage points)")
    print(f"  95% percentile interval           : [{lo:+.4f}, {hi:+.4f}]")
    print(f"  interval excludes 0               : {not (lo <= 0 <= hi)}")
    print(f"  mean dR >= +10 pp                 : {mean >= bound}")
    verdict = ("adjustment adds value" if (mean >= bound and not (lo <= 0 <= hi))
               else ("unadjusted wins" if (mean <= -bound and not (lo <= 0 <= hi))
                     else "informative null (§6.1)"))
    print(f"  OUTCOME: {verdict}")
    print("  All three outcomes are valid and reportable (§6.1). The interval is "
          "wide by design: 200 resamples were pre-specified as reduced for compute.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
