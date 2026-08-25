#!/usr/bin/env python3
"""H3 -- the circularity check. Jaccard overlap vs a 1,000-permutation null.

Run in three variants, reported together. Variant 1 is PRIMARY by
pre-registration; variants 2 and 3 are DECLARED SENSITIVITY ANALYSES. No winner
is designated.

  v1 PRIMARY      Stage A residuals exactly as pre-registered (held-out random
                  effect predicted as 0, no centring).
  v2 SENSITIVITY  residuals centred within patient, realising the stated intent of
                  §3.5's patient intercept that cross-fitting with RE=0 does not.
  v3 SENSITIVITY  the shared z(G) term removed from both tiers and renormalised,
                  so the two targets share no component.

The null samples gene sets of the OBSERVED sizes from each arm's OWN eligible
universe -- Tier C's universe excludes chr7/chr9p/chr10, so a shared universe
would misstate the chance overlap.
"""
from __future__ import annotations
import json, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd, yaml
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _genome import annotate_var  # noqa: E402
import anndata as ad

REPO = Path(__file__).resolve().parent.parent
GL = REPO / "results" / "gene_lists"
CACHE = REPO / "data" / "interim" / "metacell_expression.h5ad"
CONF = REPO / "configs" / "pipeline_config.yaml"
TARGETS = REPO / "results" / "tables" / "stage_b_targets.csv"
OUT = REPO / "results" / "tables" / "circularity_check.csv"

VARIANTS = [("v1", "PRIMARY (pre-registered)"),
            ("v2", "SENSITIVITY (patient-centred)"),
            ("v3", "SENSITIVITY (G-removed)")]


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a or b) else 0.0


def main() -> int:
    conf = yaml.safe_load(open(CONF))
    seed = conf["seed"]["master"]
    n_perm = conf["h3_circularity"]["n_permutations"]
    alpha = conf["h3_circularity"]["alpha"]
    thr = int(conf["stability_selection"]["threshold_primary"] * 100)

    mc = ad.read_h5ad(CACHE)
    ann = annotate_var(mc.var_names, conf["disjoint_set_S"]["regions"])
    in_S = ann["in_disjoint_set_S"].fillna(False).values
    universe_A = list(mc.var_names)
    universe_C = list(mc.var_names[~in_S])
    print(f"eligible universes: Tier A {len(universe_A):,} genes, "
          f"Tier C {len(universe_C):,} (chr7/chr9p/chr10 excluded)")

    tg = pd.read_csv(TARGETS)
    rng = np.random.default_rng(seed)
    rows = []
    for v, label in VARIANTS:
        fa, fc = GL / f"{v}_tierA_genes_{thr}.csv", GL / f"{v}_tierC_genes_{thr}.csv"
        if not (fa.exists() and fc.exists()):
            print(f"  {v}: gene lists missing, skipping")
            continue
        A = pd.read_csv(fa)["gene"].dropna().tolist()
        C = pd.read_csv(fc)["gene"].dropna().tolist()
        obs = jaccard(A, C)
        null = np.empty(n_perm)
        for i in range(n_perm):
            ra = rng.choice(universe_A, size=min(len(A), len(universe_A)), replace=False)
            rc = rng.choice(universe_C, size=min(len(C), len(universe_C)), replace=False)
            null[i] = jaccard(ra, rc)
        p = (1 + int((null >= obs).sum())) / (1 + n_perm)
        corr = tg[f"{v}_tierA"].corr(tg[f"{v}_tierC"])
        rows.append({"variant": v, "role": label,
                     "target_correlation": round(float(corr), 4),
                     "n_genes_tierA": len(A), "n_genes_tierC": len(C),
                     "n_overlap": len(set(A) & set(C)),
                     "jaccard_observed": round(float(obs), 6),
                     "jaccard_null_mean": round(float(null.mean()), 6),
                     "jaccard_null_sd": round(float(null.std()), 6),
                     "jaccard_null_p95": round(float(np.percentile(null, 95)), 6),
                     "p_value": round(float(p), 6),
                     "alpha": alpha,
                     "significant": bool(p < alpha)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(f"\nwrote {OUT.relative_to(REPO)}\n")
    if df.empty:
        return 1
    print(f"{'variant':<6}{'role':<32}{'corr':>7}{'|A|':>7}{'|C|':>7}{'ovl':>6}"
          f"{'Jobs':>9}{'Jnull':>9}{'p':>9}  sig")
    print("-" * 104)
    for _, r in df.iterrows():
        print(f"{r['variant']:<6}{r['role']:<32}{r['target_correlation']:>7.3f}"
              f"{r['n_genes_tierA']:>7}{r['n_genes_tierC']:>7}{r['n_overlap']:>6}"
              f"{r['jaccard_observed']:>9.4f}{r['jaccard_null_mean']:>9.4f}"
              f"{r['p_value']:>9.4f}  {r['significant']}")
    print(f"\nalpha = {alpha} (Bonferroni over the H1/H3 family)")
    print("Variant 1 is primary by pre-registration; 2 and 3 are declared "
          "sensitivity analyses. No winner is designated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
