#!/usr/bin/env python3
"""CLAUDE.md §10.2 degeneracy check -- each region of S, independently.

"If clone assignment is degenerate (e.g. one clone per patient), Tier C-disjoint
carries no information and H3 cannot be tested. Check this on Day 2, not Day 4."
and: "Check each region independently, not S as a whole ... If 9p is
uninformative, say so in the paper rather than implying all three regions
contributed."

chr7 gain and chr10 loss are whole-chromosome events and should resolve from
windowed expression. CDKN2A deletion at 9p21 is focal (~1 Mb) against a ~100-gene
window and may vanish entirely.

Report only. Does not alter the clone catalog.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
CATALOG = REPO / "results" / "tables" / "clone_catalog.csv"
OUT = REPO / "results" / "tables" / "region_degeneracy_check.csv"

REGIONS = ["chr7", "chr9p", "chr10"]


def main() -> int:
    cat = pd.read_csv(CATALOG)
    conv = cat[cat["status"] == "converged"].copy()
    print(f"clone_catalog.csv: {len(cat)} patients, {len(conv)} converged")
    print(f"window_size={cat['window_size'].iloc[0]}  "
          f"leiden_resolution={cat['leiden_resolution'].iloc[0]}\n")

    rows = []
    for r in REGIONS:
        col, wcol = f"n_clones_{r}", f"n_windows_{r}"
        n_clones = conv[col]
        n_win = conv[wcol] if wcol in conv.columns else pd.Series([0] * len(conv))
        rows.append({
            "region": r,
            "n_patients_converged": len(conv),
            "median_windows": float(n_win.median()),
            "median_clones": float(n_clones.median()),
            "min_clones": int(n_clones.min()),
            "max_clones": int(n_clones.max()),
            "n_patients_single_clone": int((n_clones <= 1).sum()),
            "n_patients_multi_clone": int((n_clones > 1).sum()),
            "frac_patients_degenerate": round(float((n_clones <= 1).mean()), 4),
            "informative": bool((n_clones > 1).mean() >= 0.5),
        })
    combined = conv["n_clones_combined"]
    rows.append({
        "region": "ALL_S_combined", "n_patients_converged": len(conv),
        "median_windows": float(conv["n_cnv_windows"].median()),
        "median_clones": float(combined.median()),
        "min_clones": int(combined.min()), "max_clones": int(combined.max()),
        "n_patients_single_clone": int((combined <= 1).sum()),
        "n_patients_multi_clone": int((combined > 1).sum()),
        "frac_patients_degenerate": round(float((combined <= 1).mean()), 4),
        "informative": bool((combined > 1).mean() >= 0.5),
    })

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(f"{'region':<16}{'windows':>9}{'clones med':>12}{'min':>6}{'max':>6}"
          f"{'1-clone pts':>13}{'degenerate':>12}  informative")
    print("-" * 88)
    for _, r in df.iterrows():
        print(f"{r['region']:<16}{r['median_windows']:>9.0f}{r['median_clones']:>12.1f}"
              f"{r['min_clones']:>6}{r['max_clones']:>6}"
              f"{r['n_patients_single_clone']:>7}/{r['n_patients_converged']:<5}"
              f"{r['frac_patients_degenerate']:>12.0%}  {r['informative']}")
    print(f"\nwrote {OUT.relative_to(REPO)}")

    print("\nper-patient clone counts by region:")
    cols = ["patient_id", "n_cnv_windows"] + [f"n_windows_{r}" for r in REGIONS] \
           + [f"n_clones_{r}" for r in REGIONS] + ["n_clones_combined"]
    cols = [c for c in cols if c in conv.columns]
    print(conv[cols].to_string(index=False))

    print("\n--- verdict ---")
    for _, r in df[df["region"] != "ALL_S_combined"].iterrows():
        if r["n_patients_multi_clone"] == 0:
            print(f"  {r['region']}: CONTRIBUTES NOTHING -- every patient resolves a "
                  f"single clone. Must be stated in the paper; the combined result "
                  f"must not imply this region contributed.")
        elif not r["informative"]:
            print(f"  {r['region']}: LARGELY UNINFORMATIVE -- "
                  f"{r['n_patients_single_clone']}/{r['n_patients_converged']} "
                  f"patients degenerate.")
        else:
            print(f"  {r['region']}: informative -- multi-clone in "
                  f"{r['n_patients_multi_clone']}/{r['n_patients_converged']} patients, "
                  f"median {r['median_clones']:.0f} clones.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
