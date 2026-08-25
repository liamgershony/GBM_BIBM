"""Compare a fresh pipeline run against the baseline snapshot.

Reports every number RESULTS_SUMMARY.md cites, baseline vs re-run, and flags any
that differ. Exit code is non-zero if any headline number moved.
"""
import json, sys
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
BASE = REPO / "results" / "_repro_baseline"
NEW = REPO / "results" / "tables"
GL_B, GL_N = BASE / "gene_lists", REPO / "results" / "gene_lists"


def load(d, name):
    f = d / name
    return pd.read_csv(f) if f.exists() else None


def cmp_scalar(label, a, b, out, tol=0.0):
    same = (a == b) if tol == 0 else (a is not None and b is not None
                                      and abs(float(a) - float(b)) <= tol)
    out.append({"quantity": label, "baseline": a, "rerun": b,
                "match": bool(same)})


def main() -> int:
    rows = []

    # cohort
    for f, key in (("cohort_n.json", "n_patients"),):
        b = json.loads((BASE / f).read_text()) if (BASE / f).exists() else {}
        n = json.loads((NEW / f).read_text()) if (NEW / f).exists() else {}
        cmp_scalar(f"cohort_n.{key}", b.get(key), n.get(key), rows)

    for name, col, label in (
        ("cohort_flow.csv", None, "cohort_flow"),
        ("discovery_cohort.csv", None, "discovery_cohort rows"),
        ("ras_scores.csv", None, "ras_scores rows"),
        ("stage_a_r2.csv", None, None),
    ):
        b, n = load(BASE, name), load(NEW, name)
        if b is None or n is None:
            rows.append({"quantity": name, "baseline": "present" if b is not None
                         else "MISSING", "rerun": "present" if n is not None
                         else "MISSING", "match": False})
            continue
        if label:
            cmp_scalar(label, len(b), len(n), rows)

    # stage A R2
    b, n = load(BASE, "stage_a_r2.csv"), load(NEW, "stage_a_r2.csv")
    if b is not None and n is not None:
        for _, rb in b.iterrows():
            rn = n[n.tier == rb.tier]
            if len(rn):
                cmp_scalar(f"stage_a R2 +patient [{rb.tier}]",
                           rb.r2_state_genotype_plus_patient,
                           float(rn.iloc[0].r2_state_genotype_plus_patient),
                           rows, tol=1e-4)

    # stage B gene counts
    bs = json.loads((BASE / "stage_b_summary.json").read_text()) \
        if (BASE / "stage_b_summary.json").exists() else {}
    ns = json.loads((NEW / "stage_b_summary.json").read_text()) \
        if (NEW / "stage_b_summary.json").exists() else {}
    for arm in sorted(set(bs) | set(ns)):
        for t in ("0.3", "0.5", "0.8"):
            cmp_scalar(f"stage_b {arm} @{t}",
                       (bs.get(arm) or {}).get(t), (ns.get(arm) or {}).get(t), rows)

    # H3
    b, n = load(BASE, "circularity_check.csv"), load(NEW, "circularity_check.csv")
    if b is not None and n is not None:
        for _, rb in b.iterrows():
            rn = n[n.variant == rb.variant]
            if len(rn):
                rn = rn.iloc[0]
                for c, tol in (("n_genes_tierA", 0), ("n_genes_tierC", 0),
                               ("n_overlap", 0), ("jaccard_observed", 1e-6),
                               ("p_value", 1e-3), ("target_correlation", 1e-3)):
                    cmp_scalar(f"H3 {rb.variant}.{c}", rb[c], rn[c], rows, tol=tol)

    # H1
    b, n = load(BASE, "h1_ablation.csv"), load(NEW, "h1_ablation.csv")
    if b is not None and n is not None:
        for c in ("rate_adjusted", "rate_unadjusted", "delta_R"):
            cmp_scalar(f"H1 mean {c}", round(float(b[c].mean()), 4),
                       round(float(n[c].mean()), 4), rows, tol=5e-3)

    # gene lists
    for arm in ("v1_tierA", "v1_tierC", "v3_tierA", "v3_tierC"):
        fb, fn = GL_B / f"{arm}_genes_30.csv", GL_N / f"{arm}_genes_30.csv"
        if fb.exists() and fn.exists():
            gb = set(pd.read_csv(fb)["gene"].dropna())
            gn = set(pd.read_csv(fn)["gene"].dropna())
            rows.append({"quantity": f"gene list {arm}@30 identical",
                         "baseline": len(gb), "rerun": len(gn),
                         "match": gb == gn})

    df = pd.DataFrame(rows)
    df.to_csv(REPO / "results" / "tables" / "reproducibility_check.csv", index=False)
    ok = int(df["match"].sum()); tot = len(df)
    print(f"{'quantity':<42}{'baseline':>14}{'rerun':>14}  match")
    print("-" * 78)
    for _, r in df.iterrows():
        flag = "ok" if r["match"] else "**DIFFERS**"
        print(f"{r['quantity']:<42}{str(r['baseline']):>14}{str(r['rerun']):>14}  {flag}")
    print(f"\n{ok}/{tot} quantities reproduce exactly")
    bad = df[~df["match"]]
    if len(bad):
        print("\nDIFFERENCES:")
        for _, r in bad.iterrows():
            print(f"  {r['quantity']}: {r['baseline']} -> {r['rerun']}")
    return 0 if ok == tot else 1


if __name__ == "__main__":
    raise SystemExit(main())
