#!/usr/bin/env python3
"""CGGA bulk replication of discovery candidate genes (CLAUDE.md §4.3).

MANDATORY CAVEAT, restated wherever this output is discussed: CGGA replication is
BULK-LEVEL SUPPORTING EVIDENCE ONLY. It is never single-cell validation of the
specific persister population a gene was discovered in.

Test, per candidate gene g and cohort c:
    logit(recurrent) ~ expression(g) + covariates
Direction must match discovery, with the sign taken from mRNAseq_693.
Significance: raw p < 0.05 in mRNAseq_693, then Benjamini-Hochberg at 5% across
all tested genes JOINTLY (not per cohort).
Effect size: |log2 FC| >= 0.25 OR odds ratio outside [0.8, 1.25].
mRNAseq_325 is supportive: it must not contradict, and its significance is not
independently required. Genes where the two batches significantly disagree in
direction are EXCLUDED and reported separately as discordant, never silently
dropped.

COVARIATES. §4.3 names sequencing platform, IDH status and tumour purity.
  * IDH status  -- present, used.
  * tumour purity -- ABSENT from the CGGA download. Dropped and declared per
    CLAUDE.md §10.4. No proxy is substituted.
  * sequencing platform -- ABSENT as a column, and constant within each cohort
    (platform differs BETWEEN mRNAseq_693 and _325, which are analysed
    separately), so it cannot enter a within-cohort model. Dropped and declared.
Whichever covariates were actually used are written to the output table.

Cohort definition: Histology in {GBM, rGBM, sGBM}. Recurrent tumours carry an "r"
prefix; filtering on Histology == "GBM" alone yields 140 primaries and ZERO
recurrent samples, i.e. a constant outcome.

Usage:
    python3 src/09_cgga_replication.py --genes <file>     # one symbol per line
    python3 src/09_cgga_replication.py --self-test        # dummy list, end to end
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
CGGA = REPO / "data" / "interim" / "cgga"
OUT = REPO / "results" / "tables" / "cgga_replication.csv"
OUT_DISC = REPO / "results" / "tables" / "cgga_discordant.csv"

COHORTS = {
    "mRNAseq_693": ("CGGA.mRNAseq_693_clinical.20200506.txt",
                    "CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt"),
    "mRNAseq_325": ("CGGA.mRNAseq_325_clinical.20200506.txt",
                    "CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt"),
}
PRIMARY_COHORT = "mRNAseq_693"
GBM_HISTOLOGY = {"GBM", "RGBM", "SGBM"}
RAW_P, FDR_Q = 0.05, 0.05
MIN_ABS_LOG2FC, OR_BAND = 0.25, (0.8, 1.25)


def load_cohort(name: str, genes: list[str]):
    clin_f, expr_f = COHORTS[name]
    clin = pd.read_csv(CGGA / clin_f, sep="\t", dtype=str)
    clin = clin[clin["Histology"].str.upper().isin(GBM_HISTOLOGY)].copy()
    clin["recurrent"] = (clin["PRS_type"].str.strip() == "Recurrent").astype(int)
    # Secondary GBM is neither primary-untreated nor a recurrence of a treated
    # GBM; excluded so the outcome is a clean primary-vs-recurrent contrast.
    clin = clin[clin["PRS_type"].str.strip().isin(["Primary", "Recurrent"])].copy()
    clin["idh"] = clin["IDH_mutation_status"].str.strip()
    clin = clin[clin["idh"].isin(["Wildtype", "Mutant"])].copy()
    clin["idh_mut"] = (clin["idh"] == "Mutant").astype(int)

    expr = pd.read_csv(CGGA / expr_f, sep="\t", index_col=0)
    expr = expr[~expr.index.duplicated(keep="first")]
    keep = [c for c in expr.columns if c in set(clin["CGGA_ID"])]
    expr = expr[keep]
    clin = clin[clin["CGGA_ID"].isin(keep)].set_index("CGGA_ID").loc[keep]

    # counts -> log2 CPM
    lib = expr.sum(axis=0).replace(0, np.nan)
    cpm = np.log2(expr.divide(lib, axis=1) * 1e6 + 1.0)
    present = [g for g in genes if g in cpm.index]
    return clin, cpm.loc[present], present


def test_cohort(name, genes, log):
    clin, cpm, present = load_cohort(name, genes)
    n_r = int(clin["recurrent"].sum())
    log(f"  {name}: {len(clin)} GBM samples ({len(clin)-n_r} primary, {n_r} recurrent), "
        f"{len(present)}/{len(genes)} genes present")
    assert 0 < n_r < len(clin), f"{name}: outcome is constant -- check the cohort filter"

    covars, used = [], []
    if clin["idh_mut"].nunique() > 1:
        covars.append(clin["idh_mut"].astype(float).values); used.append("idh_status")
    y = clin["recurrent"].values.astype(float)

    rows = []
    for g in present:
        x = cpm.loc[g].values.astype(float)
        if np.nanstd(x) == 0:
            continue
        X = np.column_stack([x] + covars) if covars else x.reshape(-1, 1)
        X = sm.add_constant(X, has_constant="add")
        try:
            m = sm.Logit(y, X).fit(disp=0, maxiter=100)
            beta, p = float(m.params[1]), float(m.pvalues[1])
        except Exception:                                    # noqa: BLE001
            continue
        lfc = float(np.nanmean(x[y == 1]) - np.nanmean(x[y == 0]))  # log2 CPM diff
        rows.append({"gene": g, "cohort": name, "beta": beta,
                     "odds_ratio": float(np.exp(beta)), "p_raw": p,
                     "log2fc_recurrent_minus_primary": lfc,
                     "direction": int(np.sign(beta)),
                     "covariates_used": ";".join(used) or "(none)",
                     "n_samples": len(clin), "n_recurrent": n_r})
    return pd.DataFrame(rows)


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes")
    ap.add_argument("--discovery-direction",
                    help="TSV: gene<TAB>direction(+1/-1) from discovery")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    log = print

    if args.self_test:
        genes = ["EGFR", "CHI3L1", "CD44", "GFAP", "PDGFRA", "OLIG1", "AQP4",
                 "SOX4", "DLL3", "VIM", "MKI67", "NOTAREALGENE123"]
        disc = {g: 1 for g in genes}
        log(f"SELF-TEST with {len(genes)} dummy genes (one deliberately absent)\n")
    else:
        assert args.genes, "--genes or --self-test required"
        genes = [l.strip() for l in open(args.genes) if l.strip()]
        disc = {}
        if args.discovery_direction:
            for l in open(args.discovery_direction):
                a = l.split()
                if len(a) >= 2:
                    disc[a[0]] = int(float(a[1]))
        log(f"{len(genes)} candidate genes\n")

    res = {c: test_cohort(c, genes, log) for c in COHORTS}
    prim = res[PRIMARY_COHORT]
    if prim.empty:
        log("no genes testable in the primary cohort")
        return 1

    # BH across all tested genes JOINTLY across cohorts (CLAUDE.md §4.3)
    allr = pd.concat(res.values(), ignore_index=True)
    allr["q_bh_joint"] = bh(allr["p_raw"].values)
    log(f"\nBH-FDR applied jointly across {len(allr)} gene x cohort tests")

    p = allr[allr.cohort == PRIMARY_COHORT].set_index("gene")
    s = allr[allr.cohort != PRIMARY_COHORT].set_index("gene")

    out, discord = [], []
    for g, r in p.iterrows():
        eff_ok = (abs(r["log2fc_recurrent_minus_primary"]) >= MIN_ABS_LOG2FC
                  or not (OR_BAND[0] <= r["odds_ratio"] <= OR_BAND[1]))
        dir_ok = (disc.get(g, r["direction"]) == r["direction"]) if disc else True
        sup = s.loc[g] if g in s.index else None
        contradicts = False
        if sup is not None:
            contradicts = (sup["direction"] != r["direction"]
                           and sup["p_raw"] < RAW_P and r["p_raw"] < RAW_P)
        rec = {"gene": g, "beta_693": r["beta"], "or_693": r["odds_ratio"],
               "p_raw_693": r["p_raw"], "q_bh_joint_693": r["q_bh_joint"],
               "log2fc_693": r["log2fc_recurrent_minus_primary"],
               "direction_693": r["direction"],
               "direction_325": (int(sup["direction"]) if sup is not None else np.nan),
               "p_raw_325": (float(sup["p_raw"]) if sup is not None else np.nan),
               "covariates_used": r["covariates_used"],
               "passes_raw_p": bool(r["p_raw"] < RAW_P),
               "passes_fdr": bool(r["q_bh_joint"] < FDR_Q),
               "passes_effect_size": bool(eff_ok),
               "direction_matches_discovery": bool(dir_ok),
               "supportive_contradicts": bool(contradicts)}
        rec["REPLICATES"] = bool(rec["passes_raw_p"] and rec["passes_fdr"]
                                 and eff_ok and dir_ok and not contradicts)
        (discord if contradicts else out).append(rec)

    df = pd.DataFrame(out)
    df.to_csv(OUT, index=False)
    pd.DataFrame(discord).to_csv(OUT_DISC, index=False)
    log(f"\nwrote {OUT.relative_to(REPO)} ({len(df)} genes) and "
        f"{OUT_DISC.relative_to(REPO)} ({len(discord)} discordant)")
    if not df.empty:
        log(f"\n  raw p<0.05 in 693      : {int(df.passes_raw_p.sum())}/{len(df)}")
        log(f"  BH-FDR q<0.05 (joint)  : {int(df.passes_fdr.sum())}/{len(df)}")
        log(f"  effect size threshold  : {int(df.passes_effect_size.sum())}/{len(df)}")
        log(f"  REPLICATE (all criteria): {int(df.REPLICATES.sum())}/{len(df)}")
        log(f"  covariates used        : {df.covariates_used.iloc[0]}")
        log("\n  purity: ABSENT from the CGGA download, dropped and declared "
            "(CLAUDE.md §10.4)")
        log("  platform: no column, and constant within each cohort, so it cannot "
            "enter a within-cohort model")
        log("\n  CAVEAT: bulk-level supporting evidence only; never single-cell "
            "validation of the persister population.")
        log("\n" + df[["gene", "or_693", "p_raw_693", "q_bh_joint_693", "log2fc_693",
                       "REPLICATES"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
