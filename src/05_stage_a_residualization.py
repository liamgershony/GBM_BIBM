#!/usr/bin/env python3
"""Stage A residualization -> results/tables/stage_a_residuals.csv.

    RAS_i = b0 + sum_k b_k 1[state_i = k] + g*genotype_i + u_p(i) + e_i

Fixed effects: cell state (4 levels, all evaluable at n=21) and genotype
(chr7xchr10 class). Random effect: patient intercept.

CROSS-FITTED per CLAUDE.md §3.5: for held-out patient p the fixed-effect
coefficients are estimated on the other 20 patients only, and p's random effect is
predicted as 0 (standard mixed-model out-of-sample practice). The residual is
RAS - fixed-effect prediction.

Also reports the IN-SAMPLE R2 per tier, which is the STOP/GO check: if state +
genotype + patient explained essentially all of a tier's variance, the residual
would be noise and there would be nothing for Stage B to find.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
RAS = REPO / "results" / "tables" / "ras_scores.csv"
OUT = REPO / "results" / "tables" / "stage_a_residuals.csv"
OUT_R2 = REPO / "results" / "tables" / "stage_a_r2.csv"

TIERS = ["ras_tier_a_reduced", "ras_tier_c_disjoint"]
FORMULA = "{y} ~ C(cell_state) + C(genotype_class)"


def main() -> int:
    d = pd.read_csv(RAS)
    d["patient_id"] = d["patient_id"].astype(str)
    patients = sorted(d["patient_id"].unique())
    print(f"{len(d):,} nuclei, {len(patients)} patients")
    print(f"states: {sorted(d['cell_state'].unique())}")
    print(f"genotype classes: {sorted(d['genotype_class'].unique())}\n")

    r2_rows = []
    for tier in TIERS:
        # ---- in-sample fit: the STOP/GO check ----
        full = smf.mixedlm(FORMULA.format(y=tier), d, groups=d["patient_id"]).fit(reml=False)
        fe_only = np.asarray(full.fittedvalues)     # already includes RE in statsmodels
        ols_fe = smf.ols(FORMULA.format(y=tier), d).fit()
        ols_fe_pat = smf.ols(FORMULA.format(y=tier) + " + C(patient_id)", d).fit()
        r2_rows.append({
            "tier": tier,
            "r2_state_genotype_only": round(float(ols_fe.rsquared), 4),
            "r2_state_genotype_plus_patient": round(float(ols_fe_pat.rsquared), 4),
            "variance_remaining": round(float(1 - ols_fe_pat.rsquared), 4),
            "ras_sd": round(float(d[tier].std()), 4),
        })

        # ---- cross-fitted residuals ----
        resid = np.full(len(d), np.nan)
        pred = np.full(len(d), np.nan)
        fold = np.full(len(d), -1, dtype=int)
        for i, p in enumerate(patients):
            te = (d["patient_id"] == p).values
            tr = ~te
            m = smf.mixedlm(FORMULA.format(y=tier), d[tr],
                            groups=d.loc[tr, "patient_id"]).fit(reml=False)
            # held-out random effect predicted as 0: fixed effects only
            X = np.asarray(smf.ols(FORMULA.format(y=tier), d).exog)[te]
            names = smf.ols(FORMULA.format(y=tier), d).exog_names
            beta = np.array([m.fe_params.get(n, 0.0) for n in names])
            yhat = X @ beta
            pred[te] = yhat
            resid[te] = d.loc[te, tier].values - yhat
            fold[te] = i
        d[f"{tier}__pred"] = pred
        d[f"{tier}__residual"] = resid
        d[f"{tier}__fold"] = fold
        print(f"{tier}")
        print(f"   in-sample R2, state+genotype        : {ols_fe.rsquared:.4f}")
        print(f"   in-sample R2, +patient              : {ols_fe_pat.rsquared:.4f}")
        print(f"   variance remaining for Stage B      : {1-ols_fe_pat.rsquared:.2%}")
        print(f"   cross-fitted residual sd            : {np.nanstd(resid):.4f} "
              f"(RAS sd {d[tier].std():.4f})")
        wp = pd.Series(resid).groupby(d["patient_id"].values).std(ddof=0)
        print(f"   median within-patient residual sd   : {wp.median():.4f}, "
              f"patients with zero: {int((wp == 0).sum())}/{len(patients)}\n")

    keep = ["nucleus_id", "patient_id", "cell_state", "genotype_class"] + TIERS + \
           [f"{t}__{s}" for t in TIERS for s in ("pred", "residual", "fold")]
    d[keep].to_csv(OUT, index=False)
    pd.DataFrame(r2_rows).to_csv(OUT_R2, index=False)
    print(f"wrote {OUT.relative_to(REPO)} and {OUT_R2.relative_to(REPO)}")

    print("\n--- STOP/GO ---")
    for r in r2_rows:
        trig = r["variance_remaining"] < 0.05
        print(f"  {r['tier']:<22} {r['variance_remaining']:.2%} of variance remains "
              f"-> {'TRIGGERED (nothing left to find)' if trig else 'NOT triggered'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
