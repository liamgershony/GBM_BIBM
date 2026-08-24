#!/usr/bin/env python3
"""Diagnose the clause (d) failures -> results/tables/clause_d_failure_diagnosis.csv.

REPORT ONLY. This script does not change the cohort, does not rewrite
cohort_n.json, and does not alter any threshold. It cross-tabulates each failing
patient's libraries against qc_per_library.csv so the failure can be classified as
genuine biological/technical sparsity or as an artifact of our own pipeline.

Patient 8 has 1 malignant primary nucleus against 1,969 recurrent, and patient 28
has 3 against 2,748. Those ratios are implausible as biology and warrant a
mechanism, not an assumption.

Classification is by explicit rule, applied per failing timepoint:

  ANNOTATION_ARTIFACT   the timepoint retained nuclei through QC, but few or none
                        carry a Tumor/Normal label -- the shortfall is missing
                        annotation coverage, not missing cells.
  QC_ATTRITION          nuclei were present in the input but were removed by our
                        own QC gates (gene count, mito fraction, doublets).
  GENUINE_SPARSITY      few nuclei were deposited for that timepoint at all.
  ANNOTATED_NON_TUMOUR  nuclei survived and were annotated, but the authors called
                        them Normal -- a real biological/compositional result.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
QC = REPO / "results" / "tables" / "qc_per_library.csv"
CLAUSE_D = REPO / "results" / "tables" / "clause_d_counts.csv"
COHORT = REPO / "results" / "tables" / "discovery_cohort.csv"
OUT = REPO / "results" / "tables" / "clause_d_failure_diagnosis.csv"

MIN_NUCLEI = 100


def classify(inp, post_qc, annotated, malignant, normal) -> tuple[str, str]:
    if inp < MIN_NUCLEI:
        return ("GENUINE_SPARSITY",
                f"only {inp} nuclei deposited for this timepoint")
    if post_qc == 0:
        return ("QC_ATTRITION",
                f"all {inp} input nuclei removed by our QC")
    if annotated == 0:
        return ("ANNOTATION_ARTIFACT",
                f"{post_qc} nuclei survived QC but 0 carry an author annotation")
    if annotated / post_qc < 0.5:
        return ("ANNOTATION_ARTIFACT",
                f"{annotated}/{post_qc} annotated ({annotated/post_qc:.0%}) -- "
                f"malignant count limited by annotation coverage")
    if post_qc < inp * 0.25:
        return ("QC_ATTRITION",
                f"{inp} -> {post_qc} nuclei ({post_qc/inp:.0%} retained) through QC")
    if malignant < MIN_NUCLEI <= normal + malignant:
        return ("ANNOTATED_NON_TUMOUR",
                f"{post_qc} nuclei annotated, but only {malignant} called Tumor "
                f"({normal} Normal) -- authors' own call")
    return ("GENUINE_SPARSITY",
            f"{inp} in -> {post_qc} post-QC -> {malignant} malignant")


def main() -> int:
    qc = list(csv.DictReader(open(QC)))
    cohort = list(csv.DictReader(open(COHORT)))
    failures = [r for r in csv.DictReader(open(CLAUSE_D))
                if r["passes_clause_d"] in ("False", "false")]
    alias_needed = {r["batch_key"]: (r["annotation_prefix"] not in ("", r["sample_id"]))
                    for r in qc}

    by_pt = defaultdict(list)
    for r in qc:
        by_pt[(r["patient_id"], r["timepoint"])].append(r)

    rows = []
    for f in failures:
        pid = f["patient_id"]
        for tp, n_mal_col in (("Primary", "n_malignant_primary"),
                              ("Recurrent", "n_malignant_recurrent")):
            libs = by_pt.get((pid, tp), [])
            i = lambda k: sum(int(x[k]) for x in libs)
            inp, post_gene = i("n_input"), i("n_post_gene_filter")
            post_mito, post_qc = i("n_post_mito_filter"), i("n_post_qc")
            annotated, unknown = i("n_annotated"), i("n_unknown")
            mal, nor = i("n_malignant"), i("n_normal")
            doub = i("n_doublets_removed")
            failed_here = int(f[n_mal_col]) < MIN_NUCLEI
            kind, why = classify(inp, post_qc, annotated, mal, nor) if failed_here \
                else ("PASSES", f"{f[n_mal_col]} malignant >= {MIN_NUCLEI}")
            rows.append({
                "patient_id": pid, "timepoint": tp,
                "is_failing_timepoint": failed_here,
                "specimens": ";".join(sorted({x["sample_id"] for x in libs})),
                "n_libraries": len(libs),
                "n_input": inp,
                "n_post_gene_filter": post_gene,
                "lost_to_gene_filter": inp - post_gene,
                "n_post_mito_filter": post_mito,
                "lost_to_mito_filter": post_gene - post_mito,
                "n_doublets_removed": doub,
                "n_post_qc": post_qc,
                "qc_retention": round(post_qc / inp, 4) if inp else 0.0,
                "n_annotated": annotated,
                "n_unknown": unknown,
                "annotation_coverage": round(annotated / post_qc, 4) if post_qc else 0.0,
                "n_malignant": mal, "n_normal": nor,
                "alias_mapping_used": any(alias_needed.get(x["batch_key"], False)
                                          for x in libs),
                "diagnosis": kind, "explanation": why,
            })

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"wrote {OUT.relative_to(REPO)} ({len(rows)} rows, "
          f"{len(failures)} failing patients x 2 timepoints)\n")

    bad = [r for r in rows if r["is_failing_timepoint"]]
    print(f"{'pt':<5}{'timepoint':<11}{'input':>7}{'postQC':>8}{'ret':>6}"
          f"{'annot':>7}{'cov':>7}{'malig':>7}{'norm':>7}  {'alias':<6}diagnosis")
    print("-" * 108)
    for r in sorted(bad, key=lambda x: int(x["patient_id"])):
        print(f"{r['patient_id']:<5}{r['timepoint']:<11}{r['n_input']:>7,}"
              f"{r['n_post_qc']:>8,}{r['qc_retention']:>6.0%}"
              f"{r['n_annotated']:>7,}{r['annotation_coverage']:>7.0%}"
              f"{r['n_malignant']:>7,}{r['n_normal']:>7,}  "
              f"{str(r['alias_mapping_used']):<6}{r['diagnosis']}")
    print()
    for r in sorted(bad, key=lambda x: int(x["patient_id"])):
        print(f"  patient {r['patient_id']:<4} {r['timepoint']:<10} "
              f"[{r['specimens']}] {r['explanation']}")

    print("\nsummary by diagnosis:")
    tally = defaultdict(int)
    for r in bad:
        tally[r["diagnosis"]] += 1
    for k, v in sorted(tally.items(), key=lambda x: -x[1]):
        print(f"   {v:>2}  {k}")
    print("\nREPORT ONLY -- cohort unchanged, cohort_n.json untouched (n = 21).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
