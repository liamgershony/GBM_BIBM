#!/usr/bin/env python3
"""Build RAS Tier A-reduced and Tier C-disjoint -> results/tables/ras_scores.csv.

Tier A-reduced = (1/3)z(T) + (1/3)z(G) + (1/3)z(Ab_state)
Tier C-disjoint = 0.5 z(G) + 0.5 z(Ab_clone)

O is DROPPED per the §9.1 contingency: balanced EMD conserves mass, so O as
specified is identically the metacell's own size (corr 1.000000). The score is
named `ras_tier_a_reduced`, never `ras_tier_a`. See DEVIATIONS.md.

z() is POOLED across the cohort (DEVIATIONS.md): within-patient z is a
divide-by-zero for G in 19/21 patients.

Weights are the frozen equal values. G is constant within 19/21 patients, so after
Stage A removes between-patient variance Tier A-reduced behaves as a two-component
score -- reported, not reweighted.

Scope: primary malignant nuclei that entered the inferCNV input, i.e. capped at
infercnv.max_nuclei_per_patient. G and Ab(clone) are undefined outside that set.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
COMPONENTS = REPO / "results" / "tables" / "ras_component_cells.csv"
STATES = REPO / "results" / "tables" / "state_assignments.csv"
OUT = REPO / "results" / "tables" / "ras_scores.csv"

W_A = {"T": 1/3, "G": 1/3, "Ab_state": 1/3}
W_C = {"G": 0.5, "Ab_clone": 0.5}


def z_pooled(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    sd = v.std()
    if sd == 0:
        raise ValueError("pooled SD is zero -- component carries no variance at all")
    return (v - v.mean()) / sd


def main() -> int:
    comp = pd.read_csv(COMPONENTS)
    st = pd.read_csv(STATES)
    st = st.rename(columns={st.columns[0]: "nucleus_id"})
    st["patient_id"] = st["patient_id"].astype(str)
    comp["patient_id"] = comp["patient_id"].astype(str)

    # Ab(state): log2 FC in the abundance of the cell's state, primary -> recurrent
    eps = 1.0
    ab_state = {}
    for p, g in st.groupby("patient_id"):
        pr = g[g["timepoint"] == "Primary"]["cell_state"].value_counts()
        rc = g[g["timepoint"] == "Recurrent"]["cell_state"].value_counts()
        npr, nrc = max(pr.sum(), 1), max(rc.sum(), 1)
        for s_ in set(pr.index) | set(rc.index):
            ab_state[(p, s_)] = float(np.log2(((rc.get(s_, 0) + eps) / (nrc + eps)) /
                                              ((pr.get(s_, 0) + eps) / (npr + eps))))

    state_of = dict(zip(st["nucleus_id"], st["cell_state"]))
    comp["cell_state"] = comp["nucleus_id"].map(state_of)
    missing = comp["cell_state"].isna().sum()
    print(f"component rows {len(comp):,}; without a state call: {missing:,}")
    comp = comp.dropna(subset=["cell_state"]).copy()
    comp["Ab_state"] = [ab_state[(p, s_)] for p, s_ in
                        zip(comp["patient_id"], comp["cell_state"])]

    for c in ("T", "G", "Ab_state", "Ab_clone"):
        comp[f"z_{c}"] = z_pooled(comp[c].values)

    comp["ras_tier_a_reduced"] = sum(w * comp[f"z_{c}"] for c, w in W_A.items())
    comp["ras_tier_c_disjoint"] = sum(w * comp[f"z_{c}"] for c, w in W_C.items())

    cols = ["nucleus_id", "patient_id", "cell_state", "genotype_class",
            "chr7_gain", "chr10_loss", "T", "G", "Ab_state", "Ab_clone",
            "z_T", "z_G", "z_Ab_state", "z_Ab_clone",
            "ras_tier_a_reduced", "ras_tier_c_disjoint"]
    comp[cols].to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(REPO)} ({len(comp):,} primary malignant nuclei, "
          f"{comp['patient_id'].nunique()} patients)")

    print("\ncomponent summary (pooled z):")
    for c in ("T", "G", "Ab_state", "Ab_clone"):
        wp = comp.groupby("patient_id")[c].std(ddof=0)
        print(f"  {c:<9} pooled sd {comp[c].std():8.4f}   median within-patient sd "
              f"{wp.median():8.4f}   patients with zero within-patient sd "
              f"{int((wp == 0).sum())}/{comp['patient_id'].nunique()}")
    print("\nscores:")
    for s_ in ("ras_tier_a_reduced", "ras_tier_c_disjoint"):
        wp = comp.groupby("patient_id")[s_].std(ddof=0)
        print(f"  {s_:<22} mean {comp[s_].mean():7.4f}  sd {comp[s_].std():7.4f}  "
              f"median within-patient sd {wp.median():7.4f}")
    print(f"\n  corr(Tier A-reduced, Tier C-disjoint) = "
          f"{comp['ras_tier_a_reduced'].corr(comp['ras_tier_c_disjoint']):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
