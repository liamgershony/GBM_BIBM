#!/usr/bin/env python3
"""Clause (d) gate -- the step where n stops being 29.

docs/COHORT_RULE.md clause (d): a patient enters the discovery cohort only if it
has >=100 usable nuclei at BOTH timepoints after QC. This is deliberately a
standalone, logged gate rather than a filter inside a preprocessing loop, because
every frozen formula downstream reads the resulting n:

    n_folds        = n
    nb_correction  = 1/n + 1/(n-1)          (Nadeau-Bengio, CLAUDE.md 3.7)
    evaluability   = floor(n/2) + 1         (CLAUDE.md 6.2)

"Usable" means MALIGNANT nuclei surviving QC. Non-malignant nuclei are retained in
01_qc.h5ad as the inferCNV reference (Day 2) but do not count toward clause (d);
`unknown` nuclei count toward neither, since neither status can be asserted.

Libraries are POOLED per (patient_id, timepoint) per the agreed rule.

Every one of the 29 candidate patients appears in the output table, passing or
failing, with its counts and reason. Nothing is silently dropped.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import anndata as ad
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
QC_H5 = REPO / "data" / "processed" / "01_qc.h5ad"
COHORT = REPO / "results" / "tables" / "discovery_cohort.csv"
MANIFEST = REPO / "results" / "tables" / "sample_manifest.csv"
OUT_COUNTS = REPO / "results" / "tables" / "clause_d_counts.csv"
OUT_FLOW = REPO / "results" / "tables" / "cohort_flow.csv"
OUT_N = REPO / "results" / "tables" / "cohort_n.json"

MIN_NUCLEI = 100   # docs/COHORT_RULE.md clause (d)


def main() -> int:
    adata = ad.read_h5ad(QC_H5)
    obs = adata.obs
    print(f"01_qc.h5ad: {adata.n_obs:,} nuclei x {adata.n_vars:,} genes")
    print(f"  malignant {int(obs['is_malignant'].sum()):,} | "
          f"normal {int(obs['is_reference_normal'].sum()):,} | "
          f"unknown {int((obs['tumor_normal_annotation'] == 'unknown').sum()):,}")

    cohort = list(csv.DictReader(open(COHORT)))
    patients = sorted({r["patient_id"] for r in cohort},
                      key=lambda x: int(x) if x.isdigit() else 10**6)
    print(f"candidate patients (clauses a-c): {len(patients)}")

    mal = obs[obs["is_malignant"]]
    counts = (mal.groupby(["patient_id", "timepoint"], observed=True)
                 .size().unstack(fill_value=0))
    for tp in ("Primary", "Recurrent"):
        if tp not in counts.columns:
            counts[tp] = 0

    rows = []
    for p in patients:
        n_p = int(counts.loc[p, "Primary"]) if p in counts.index else 0
        n_r = int(counts.loc[p, "Recurrent"]) if p in counts.index else 0
        ok_p, ok_r = n_p >= MIN_NUCLEI, n_r >= MIN_NUCLEI
        if ok_p and ok_r:
            reason = "pass"
        else:
            failed = [f"{t} {n} < {MIN_NUCLEI}" for t, n, ok in
                      (("Primary", n_p, ok_p), ("Recurrent", n_r, ok_r)) if not ok]
            reason = "FAIL: " + "; ".join(failed)
        specs = sorted({r["sample_id"] for r in cohort if r["patient_id"] == p})
        rows.append({"patient_id": p, "n_malignant_primary": n_p,
                     "n_malignant_recurrent": n_r,
                     "passes_clause_d": ok_p and ok_r, "reason": reason,
                     "specimens": ";".join(specs)})

    OUT_COUNTS.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_COUNTS, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    passed = [r for r in rows if r["passes_clause_d"]]
    failed = [r for r in rows if not r["passes_clause_d"]]
    n = len(passed)

    print(f"\nwrote {OUT_COUNTS.relative_to(REPO)} ({len(rows)} rows: "
          f"{n} pass, {len(failed)} fail)")
    if failed:
        print("\nFAILED clause (d) -- reported, not silently dropped:")
        for r in failed:
            print(f"   patient {r['patient_id']:<5} P={r['n_malignant_primary']:>6,} "
                  f"R={r['n_malignant_recurrent']:>6,}   {r['reason']}")

    # ---- cohort flow (Figure 1) -------------------------------------------
    gsms = list(csv.DictReader(open(MANIFEST)))
    flow = [
        ("GSMs in GSE174554", len(gsms)),
        ("human snRNA-seq GSMs", sum(1 for g in gsms
                                     if g["organism"] == "Homo sapiens"
                                     and g["assay"] == "snRNA-seq")),
        ("matched pairs (clauses a+b)", 30),
        ("IDH-wildtype pairs (clause c)", len(patients)),
        ("libraries loaded", len({r["batch_key"] for r in cohort})),
        ("nuclei after QC", int(adata.n_obs)),
        ("malignant nuclei", int(obs["is_malignant"].sum())),
        (f"patients passing clause (d) (>={MIN_NUCLEI} malignant both timepoints)", n),
    ]
    with open(OUT_FLOW, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["stage", "count"]); w.writerows(flow)
    print(f"\ncohort flow -> {OUT_FLOW.relative_to(REPO)}")
    for label, v in flow:
        print(f"   {label:<62} {v:>9,}")

    # ---- the single source of n -------------------------------------------
    payload = {
        "n_patients": n,
        "patient_ids": [r["patient_id"] for r in passed],
        "min_nuclei_per_timepoint": MIN_NUCLEI,
        "n_folds": n,
        "nadeau_bengio_correction": (1 / n + 1 / (n - 1)) if n > 1 else None,
        "evaluability_floor_patients": math.floor(n / 2) + 1,
        "source": "src/01c_clause_d_gate.py",
    }
    OUT_N.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {OUT_N.relative_to(REPO)} -- THE single source of n downstream")
    print(f"   n_patients          = {n}")
    print(f"   n_folds             = {n}")
    print(f"   Nadeau-Bengio       = 1/{n} + 1/{n-1} = "
          f"{payload['nadeau_bengio_correction']:.6f}")
    print(f"   evaluability floor  = floor({n}/2)+1 = "
          f"{payload['evaluability_floor_patients']} patients per state")
    print(f"   H1 admissible (n>=16)? {'yes' if n >= 16 else 'NO -- H1 is dropped'}")

    assert len(rows) == len(patients), "clause (d) table lost a candidate patient"
    assert len(passed) + len(failed) == len(patients)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
