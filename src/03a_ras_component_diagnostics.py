#!/usr/bin/env python3
"""Pre-RAS diagnostics: self-residualisation, z-scope, and Tier A variance.

Builds NO RAS score and writes nothing into the pipeline. Three questions that
must be answered before any RAS value exists:

1. SELF-RESIDUALISATION. RAS_C is now essentially Ab(clone), a deterministic
   function of (patient, chr7xchr10 class). If Stage A's genotype fixed effect is
   that same class, it can explain nearly all of RAS_C and the residual goes to
   zero -- the STOP/GO "nothing left to find" gate.
2. Z-SCOPE. The frozen config does not say whether z() is within-patient or
   pooled. With G constant in 19/21 patients, within-patient z is a divide by
   zero.
3. TIER A COST. How much of RAS_A's variance does a constant G actually carry?

NOT COMPUTABLE YET, and reported as such rather than guessed:
  * Ab(state) and Stage A's cell-state fixed effect -- blocked on the Neftel
    metamodules (see src/04_states.py).
  * O (optimal transport) -- not built; Day 3.
Omitting the state term can only UNDERSTATE the R2 of the full Stage A model, so
the self-residualisation result below is a lower bound.
"""

from __future__ import annotations

import json
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import yaml

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
IN_DIR = REPO / "data" / "interim" / "cnv_input"
INTEG = REPO / "data" / "processed" / "02_integrated.h5ad"
CONF = REPO / "configs" / "pipeline_config.yaml"
COHORT_N = REPO / "results" / "tables" / "cohort_n.json"
OUT_CELLS = REPO / "results" / "tables" / "ras_component_cells.csv"
OUT_Z = REPO / "results" / "tables" / "zscope_comparison.csv"
OUT_VAR = REPO / "results" / "tables" / "tier_a_variance.csv"

WINDOW_SIZE, K_SD = 100, 2.0


def components_for_patient(pid: str) -> pd.DataFrame:
    """Per-nucleus G, Ab(clone) and genotype class for one patient's PRIMARY cells."""
    import infercnvpy as icnv
    a = ad.read_h5ad(IN_DIR / f"{pid}.h5ad")
    icnv.tl.infercnv(a, reference_key="cnv_reference", reference_cat=["Normal"],
                     window_size=WINDOW_SIZE)
    X = a.obsm["X_cnv"]
    X = np.asarray(X.todense() if hasattr(X, "todense") else X)
    mal = (a.obs["cnv_reference"] == "Malignant").values
    ref = ~mal
    cp = sorted(dict(a.uns["cnv"]["chr_pos"]).items(), key=lambda kv: kv[1])
    b = {}
    for i, (c, s) in enumerate(cp):
        b[c] = (s, cp[i + 1][1] if i + 1 < len(cp) else X.shape[1])
    rm = lambda rows, c: X[rows, b[c][0]:b[c][1]].mean(axis=1)
    thr7 = rm(ref, "chr7").mean() + K_SD * rm(ref, "chr7").std()
    thr10 = rm(ref, "chr10").mean() - K_SD * rm(ref, "chr10").std()
    g7, l10 = rm(mal, "chr7") > thr7, rm(mal, "chr10") < thr10
    cls = np.array([f"{'+7' if g else 'no7'}_{'-10' if l else 'no10'}"
                    for g, l in zip(g7, l10)])
    tp = a.obs["timepoint"].values[mal]
    ids = np.asarray(a.obs_names)[mal]

    prim_m, rec_m = tp == "Primary", tp == "Recurrent"
    prim, rec = cls[prim_m], cls[rec_m]
    rec_set = set(rec.tolist())
    pc, rc = pd.Series(prim).value_counts(), pd.Series(rec).value_counts()
    eps = 1.0
    ab = {c: float(np.log2(((rc.get(c, 0) + eps) / (len(rec) + eps)) /
                           ((pc.get(c, 0) + eps) / (len(prim) + eps))))
          for c in sorted(set(pc.index) | set(rc.index))}
    return pd.DataFrame({
        "nucleus_id": ids[prim_m], "patient_id": pid,
        "genotype_class": prim,
        "chr7_gain": g7[prim_m].astype(int),
        "chr10_loss": l10[prim_m].astype(int),
        "G": [1 if c in rec_set else 0 for c in prim],
        "Ab_clone": [ab[c] for c in prim],
    })


def zscore(v, ddof=0):
    v = np.asarray(v, dtype=float)
    s = v.std(ddof=ddof)
    return np.full_like(v, np.nan) if s == 0 else (v - v.mean()) / s


def main() -> int:
    conf = yaml.safe_load(open(CONF))
    pts = json.loads(COHORT_N.read_text())["patient_ids"]

    parts = []
    with ProcessPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(components_for_patient, p): p for p in pts}
        for f in as_completed(futs):
            parts.append(f.result())
            print(f"  cnv {len(parts)}/{len(pts)}", flush=True)
    df = pd.concat(parts, ignore_index=True)

    # ---- T: cosine similarity to the patient's recurrent centroid -----------
    ad_i = ad.read_h5ad(INTEG)
    emb = ad_i.obsm["X_pca_harmony"]
    o = ad_i.obs
    T = {}
    for p in pts:
        rec = (o["patient_id"].astype(str) == p) & (o["timepoint"] == "Recurrent") \
              & o["is_malignant"].values
        if rec.sum() == 0:
            continue
        cen = emb[rec.values].mean(axis=0)
        cn = np.linalg.norm(cen)
        sel = (o["patient_id"].astype(str) == p) & (o["timepoint"] == "Primary") \
              & o["is_malignant"].values
        V = emb[sel.values]
        sims = (V @ cen) / (np.linalg.norm(V, axis=1) * cn + 1e-12)
        for nid, s in zip(o.index[sel.values], sims):
            T[nid] = float(s)
    df["T"] = df["nucleus_id"].map(T)
    n_missing = int(df["T"].isna().sum())
    df = df.dropna(subset=["T"]).copy()
    print(f"\nprimary malignant cells with all available components: {len(df):,} "
          f"({n_missing:,} dropped for missing embedding)")
    df.to_csv(OUT_CELLS, index=False)

    # ================= 2. Z-SCOPE =================
    print("\n=== 2. z-scoring scope ===")
    comps = ["T", "G", "Ab_clone"]
    rows = []
    for c in comps:
        pooled = zscore(df[c].values)
        pooled_bad = int(np.isnan(pooled).sum())
        within_nan_pts, within_all = 0, []
        for p, g in df.groupby("patient_id"):
            z = zscore(g[c].values)
            if np.isnan(z).all():
                within_nan_pts += 1
            within_all.append(z)
        within = np.concatenate(within_all)
        rows.append({
            "component": c,
            "pooled_sd": round(float(df[c].std()), 5),
            "pooled_all_nan": pooled_bad > 0,
            "n_patients_zero_variance": int(sum(
                g[c].std(ddof=0) == 0 for _, g in df.groupby("patient_id"))),
            "n_patients_total": df["patient_id"].nunique(),
            "within_patient_nan_cells": int(np.isnan(within).sum()),
            "within_patient_nan_frac": round(float(np.isnan(within).mean()), 4),
        })
    zdf = pd.DataFrame(rows)
    zdf.to_csv(OUT_Z, index=False)
    print(zdf.to_string(index=False))

    # ================= 3. TIER A VARIANCE =================
    print("\n=== 3. Tier A variance contribution (T and G only; O and Ab_state "
          "not yet built) ===")
    vrows = []
    for c in ["T", "G"]:
        pv = df.groupby("patient_id")[c].var(ddof=0)
        vrows.append({"component": c,
                      "pooled_variance": round(float(df[c].var(ddof=0)), 6),
                      "median_within_patient_variance": round(float(pv.median()), 6),
                      "n_patients_zero_within_variance": int((pv == 0).sum()),
                      "n_patients": len(pv)})
    vdf = pd.DataFrame(vrows)
    vdf.to_csv(OUT_VAR, index=False)
    print(vdf.to_string(index=False))
    zt, zg = zscore(df["T"].values), zscore(df["G"].values)
    if not np.isnan(zg).all():
        partial = 0.25 * zt + 0.25 * zg
        print(f"\n  var(0.25*z(T))           = {np.var(0.25*zt):.6f}")
        print(f"  var(0.25*z(G))           = {np.var(0.25*zg):.6f}")
        print(f"  G share of the two-term partial RAS_A variance = "
              f"{np.var(0.25*zg)/np.var(partial):.1%}")
        print("  (a four-component RAS_A would divide these shares further; O and "
              "Ab_state are not yet built)")

    # ================= 1. SELF-RESIDUALISATION =================
    print("\n=== 1. Stage A self-residualisation on Tier C ===")
    d = df.copy()
    zg_ = zscore(d["G"].values)
    if np.isnan(zg_).all():
        print("  z(G) is undefined pooled; using raw G for RAS_C construction")
        zg_ = d["G"].values.astype(float)
    d["RAS_C"] = 0.5 * zg_ + 0.5 * zscore(d["Ab_clone"].values)
    d["patient_id"] = d["patient_id"].astype(str)

    def pseudo_r2(formula, label):
        m = smf.mixedlm(formula, d, groups=d["patient_id"]).fit(reml=False)
        fitted_fe = m.fittedvalues
        re = m.random_effects
        fitted_full = fitted_fe + d["patient_id"].map(
            {k: float(np.asarray(v).ravel()[0]) for k, v in re.items()}).values
        r2_fe = 1 - np.var(d["RAS_C"] - fitted_fe) / np.var(d["RAS_C"])
        r2_full = 1 - np.var(d["RAS_C"] - fitted_full) / np.var(d["RAS_C"])
        print(f"  {label}")
        print(f"     pseudo-R2 fixed effects only      : {r2_fe:.4f}")
        print(f"     pseudo-R2 fixed + patient random  : {r2_full:.4f}")
        print(f"     residual SD                       : "
              f"{np.std(d['RAS_C'] - fitted_full):.5f}  "
              f"(RAS_C SD {np.std(d['RAS_C']):.5f})")
        return r2_full

    a1 = pseudo_r2("RAS_C ~ C(genotype_class)",
                   "WITH genotype = chr7xchr10 class (4 levels), + (1|patient)")
    a2 = pseudo_r2("RAS_C ~ chr7_gain",
                   "WITH genotype = chr7 gain only (EGFR-amp proxy, §3.5), + (1|patient)")
    a3 = pseudo_r2("RAS_C ~ 1",
                   "OMITTING genotype entirely, + (1|patient)")
    print(f"\n  Both reported, no choice made between them.")
    print(f"  variance of RAS_C left after the class model: {1-a1:.4%}")
    print(f"  variance of RAS_C left after no-genotype model: {1-a3:.4%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
