#!/usr/bin/env python3
"""Would a chr7-gain x chr10-loss genotype class survive as a clone definition?

02c showed the CNV signal is real (chr7 gain and chr10 loss in 21/21 patients)
but that Leiden clone calls are indistinguishable from noise. The natural
replacement is to threshold the two canonical events into up to 4 classes:

    (+7,-10)   (+7, .)   ( . ,-10)   ( . , . )

But chr7 gain and chr10 loss are near-universal in GBM. If essentially every
malignant cell lands in (+7,-10), then:

  * G is constant  -- every primary cell's class is also seen at recurrence, so
                      G == 1 for all cells and contributes no variance;
  * Ab(clone) is constant per patient -- one class means one abundance ratio.

RAS Tier C-disjoint = 0.5*z(G) + 0.5*z(Ab_clone) would then have ZERO variance and
H3 would be untestable. This is the CLAUDE.md §10.2 degeneracy arriving from the
other direction, and it must be measured before RAS is built on it.

Thresholds are REFERENCE-CALIBRATED, not absolute: for each patient and region,
gain/loss is called against that patient's own non-malignant nuclei
(mean +/- K*sd), so no global cutoff is invented.

REPORT ONLY. Builds no RAS component and writes no genotype into the pipeline.
"""

from __future__ import annotations

import json
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
IN_DIR = REPO / "data" / "interim" / "cnv_input"
CONF = REPO / "configs" / "pipeline_config.yaml"
COHORT_N = REPO / "results" / "tables" / "cohort_n.json"
OUT = REPO / "results" / "tables" / "genotype_class_degeneracy.csv"
OUT_CELLS = REPO / "results" / "tables" / "genotype_class_per_patient.csv"

WINDOW_SIZE = 100
K_SD = 2.0          # gain/loss called at mean +/- K_SD of the patient's own reference
MIN_CLASS_NUCLEI = 20


def classify_patient(pid: str) -> dict:
    import infercnvpy as icnv
    a = ad.read_h5ad(IN_DIR / f"{pid}.h5ad")
    icnv.tl.infercnv(a, reference_key="cnv_reference", reference_cat=["Normal"],
                     window_size=WINDOW_SIZE)
    X = a.obsm["X_cnv"]
    X = np.asarray(X.todense() if hasattr(X, "todense") else X)
    mal = (a.obs["cnv_reference"] == "Malignant").values
    ref = ~mal
    cp = sorted(dict(a.uns["cnv"]["chr_pos"]).items(), key=lambda kv: kv[1])
    bounds = {}
    for i, (chrom, s) in enumerate(cp):
        e = cp[i + 1][1] if i + 1 < len(cp) else X.shape[1]
        bounds[chrom] = (s, e)

    def region_mean(rows, chrom):
        s, e = bounds[chrom]
        return X[rows, s:e].mean(axis=1)

    g7_ref, g10_ref = region_mean(ref, "chr7"), region_mean(ref, "chr10")
    g7_mal, g10_mal = region_mean(mal, "chr7"), region_mean(mal, "chr10")
    thr7 = g7_ref.mean() + K_SD * g7_ref.std()
    thr10 = g10_ref.mean() - K_SD * g10_ref.std()

    gain7 = g7_mal > thr7
    loss10 = g10_mal < thr10
    cls = np.array([f"{'+7' if g else '..'}/{'-10' if l else '...'}"
                    for g, l in zip(gain7, loss10)])
    tp = a.obs["timepoint"].values[mal]

    prim, rec = cls[tp == "Primary"], cls[tp == "Recurrent"]
    rec_classes = set(rec.tolist())
    # G per CLAUDE.md §3.2: 1 if the primary cell's class is also seen at recurrence
    G = np.array([1 if c in rec_classes else 0 for c in prim])

    # Ab(clone) per CLAUDE.md §3.2: log2 fold-change in the abundance of the
    # cell's class, primary -> recurrent, within the patient. If G is constant,
    # this is the ONLY remaining source of variance in RAS Tier C-disjoint.
    pc, rc = pd.Series(prim).value_counts(), pd.Series(rec).value_counts()
    all_cls = sorted(set(pc.index) | set(rc.index))
    eps = 1.0
    ab = {c: float(np.log2(((rc.get(c, 0) + eps) / (len(rec) + eps)) /
                           ((pc.get(c, 0) + eps) / (len(prim) + eps))))
          for c in all_cls}
    ab_prim = np.array([ab[c] for c in prim]) if len(prim) else np.array([])

    counts = pd.Series(cls).value_counts()
    big = counts[counts >= MIN_CLASS_NUCLEI]
    return {
        "patient_id": pid,
        "n_malignant": int(mal.sum()),
        "n_primary": int((tp == "Primary").sum()),
        "n_recurrent": int((tp == "Recurrent").sum()),
        "n_classes_observed": int(len(counts)),
        "n_classes_ge20": int(len(big)),
        "dominant_class": str(counts.index[0]),
        "n_dominant": int(counts.iloc[0]),
        "frac_dominant": round(float(counts.iloc[0] / len(cls)), 4),
        "classes_at_recurrence": ";".join(sorted(rec_classes)),
        "G_mean": round(float(G.mean()), 4) if len(G) else float("nan"),
        "G_varies_within_patient": bool(len(G) and 0 < G.mean() < 1),
        "n_primary_G0": int((G == 0).sum()),
        "n_primary_G1": int((G == 1).sum()),
        "class_distribution": ";".join(f"{k}={v}" for k, v in counts.items()),
        "n_classes_primary": int(len(pc)),
        "n_classes_recurrent": int(len(rc)),
        "classes_missing_at_recurrence": ";".join(sorted(set(pc.index) - set(rc.index))) or "(none)",
        "ab_clone_distinct_values": int(len(set(np.round(ab_prim, 6)))) if len(ab_prim) else 0,
        "ab_clone_sd_over_primary": round(float(ab_prim.std()), 4) if len(ab_prim) else 0.0,
        "ab_clone_range": round(float(ab_prim.max() - ab_prim.min()), 4) if len(ab_prim) else 0.0,
        "thr_chr7_gain": round(float(thr7), 5),
        "thr_chr10_loss": round(float(thr10), 5),
    }


def main() -> int:
    conf = yaml.safe_load(open(CONF))
    pts = json.loads(COHORT_N.read_text())["patient_ids"]
    rows = []
    with ProcessPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(classify_patient, p): p for p in pts}
        for f in as_completed(futs):
            rows.append(f.result())
            print(f"  {len(rows)}/{len(pts)}", flush=True)
    df = pd.DataFrame(rows).sort_values(
        "patient_id", key=lambda s: s.astype(str).str.zfill(6))
    df.to_csv(OUT, index=False)
    df[["patient_id", "class_distribution", "classes_at_recurrence",
        "G_mean", "G_varies_within_patient"]].to_csv(OUT_CELLS, index=False)
    print(f"\nwrote {OUT.relative_to(REPO)}\n")

    print(f"{'pt':<5}{'malig':>7}{'prim':>6}{'rec':>6}{'cls':>5}{'cls>=20':>9}"
          f"{'dominant':>12}{'frac':>7}{'G mean':>8}{'G varies':>10}")
    print("-" * 82)
    for _, r in df.iterrows():
        print(f"{r['patient_id']:<5}{r['n_malignant']:>7,}{r['n_primary']:>6,}"
              f"{r['n_recurrent']:>6,}{r['n_classes_observed']:>5}"
              f"{r['n_classes_ge20']:>9}{r['dominant_class']:>12}"
              f"{r['frac_dominant']:>7.0%}{r['G_mean']:>8.2f}"
              f"{str(r['G_varies_within_patient']):>10}")

    n = len(df)
    multi = int((df["n_classes_ge20"] > 1).sum())
    gvar = int(df["G_varies_within_patient"].sum())
    print(f"\n--- verdict ---")
    print(f"  patients with >1 class holding >={MIN_CLASS_NUCLEI} nuclei : {multi}/{n}")
    print(f"  median fraction in the dominant class               : "
          f"{df['frac_dominant'].median():.0%}")
    print(f"  patients where G VARIES within the patient          : {gvar}/{n}")
    print(f"  patients where G is constant (no variance)          : {n - gvar}/{n}")
    print(f"\n  --- Ab(clone): the other half of RAS Tier C-disjoint ---")
    print(f"  median distinct Ab(clone) values per patient : "
          f"{df['ab_clone_distinct_values'].median():.0f}")
    print(f"  median SD of Ab(clone) across primary cells  : "
          f"{df['ab_clone_sd_over_primary'].median():.3f}")
    print(f"  patients with zero Ab(clone) variance        : "
          f"{int((df['ab_clone_sd_over_primary'] == 0).sum())}/{n}")
    print(f"  patients missing >=1 primary class at recurrence: "
          f"{int((df['classes_missing_at_recurrence'] != '(none)').sum())}/{n}")

    if gvar == 0:
        print("\n  G IS CONSTANT IN EVERY PATIENT. RAS Tier C-disjoint would have zero")
        print("  variance from the G term and H3 could not be tested on it.")
    elif gvar < n / 2:
        print(f"\n  G varies in a MINORITY of patients ({gvar}/{n}). Tier C-disjoint")
        print("  would carry information for those patients only.")
    else:
        print(f"\n  G varies in {gvar}/{n} patients -- the class definition carries "
              "within-patient variance.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
