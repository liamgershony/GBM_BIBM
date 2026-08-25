#!/usr/bin/env python3
"""Build the six Stage B targets: 3 H3 variants x 2 tiers.

Variant 1 PRIMARY (pre-registered): Stage A residual, held-out random effect = 0.
Variant 2 patient-centred SENSITIVITY: variant-1 residual centred within patient.
Variant 3 G-removed SENSITIVITY: scores rebuilt without the shared z(G) term,
          renormalised, then residualised with the same Stage A model.

Variants 2 and 3 are declared sensitivity analyses, not replacements.
"""
from __future__ import annotations
import warnings
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.formula.api as smf
warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
RAS = REPO / "results" / "tables" / "ras_scores.csv"
RESID = REPO / "results" / "tables" / "stage_a_residuals.csv"
OUT = REPO / "results" / "tables" / "stage_b_targets.csv"
FORM = "{y} ~ C(cell_state) + C(genotype_class)"


def crossfit(d, ycol):
    pats = sorted(d["patient_id"].unique())
    res = np.full(len(d), np.nan)
    base = smf.ols(FORM.format(y=ycol), d)
    names, exog = base.exog_names, np.asarray(base.exog)
    for p in pats:
        te = (d["patient_id"] == p).values
        m = smf.mixedlm(FORM.format(y=ycol), d[~te],
                        groups=d.loc[~te, "patient_id"]).fit(reml=False)
        beta = np.array([m.fe_params.get(n, 0.0) for n in names])
        res[te] = d.loc[te, ycol].values - exog[te] @ beta
    return res


def main() -> int:
    r = pd.read_csv(RAS); r["patient_id"] = r["patient_id"].astype(str)
    d = pd.read_csv(RESID); d["patient_id"] = d["patient_id"].astype(str)
    assert (r["nucleus_id"].values == d["nucleus_id"].values).all()
    out = d[["nucleus_id", "patient_id", "cell_state", "genotype_class"]].copy()

    # --- variant 1: pre-registered ---
    out["v1_tierA"] = d["ras_tier_a_reduced__residual"].values
    out["v1_tierC"] = d["ras_tier_c_disjoint__residual"].values

    # --- variant 2: patient-centred ---
    for t in ("tierA", "tierC"):
        v = out[f"v1_{t}"]
        out[f"v2_{t}"] = (v - v.groupby(out["patient_id"]).transform("mean")).values

    # --- variant 3: G removed, renormalised, then Stage A ---
    z = lambda v: (v - v.mean()) / v.std()
    r["noG_tierA"] = 0.5 * r["z_T"] + 0.5 * r["z_Ab_state"]
    r["noG_tierC"] = r["z_Ab_clone"]
    tmp = r[["patient_id", "cell_state", "genotype_class",
             "noG_tierA", "noG_tierC"]].copy()
    out["v3_tierA"] = crossfit(tmp, "noG_tierA")
    out["v3_tierC"] = crossfit(tmp, "noG_tierC")

    out.to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(REPO)} ({len(out):,} nuclei)")
    print("\ntarget correlations (nucleus level):")
    for v, lab in (("v1", "PRIMARY pre-registered"),
                   ("v2", "SENSITIVITY patient-centred"),
                   ("v3", "SENSITIVITY G-removed")):
        c = out[f"{v}_tierA"].corr(out[f"{v}_tierC"])
        print(f"  {v} {lab:<30} corr(TierA, TierC) = {c:+.4f}")
    print("\nraw-score correlations for reference:")
    print(f"  with G    : {r['ras_tier_a_reduced'].corr(r['ras_tier_c_disjoint']):+.4f}")
    print(f"  without G : {r['noG_tierA'].corr(r['noG_tierC']):+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
