#!/usr/bin/env python3
"""Join GEO sample manifest to Wang et al. Supplementary Table 1 for IDH status.

Per-patient IDH status is not in the GEO deposit (CLAUDE.md §4). This script
joins results/tables/sample_manifest.csv to Supplementary Table 1 of
Wang et al. 2022 (PMC9767870, CC BY) on the specimen ID, and writes
results/tables/patient_idh_status.csv.

This is a transcription and a join. NO COHORT RULE IS APPLIED and no sample is
excluded — selecting the discovery cohort is a separate written decision under
CLAUDE.md §10.1.

Usage:  python3 src/00e_join_idh_status.py
"""

from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _download_utils import read_id_column  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
MANIFEST = REPO / "results" / "tables" / "sample_manifest.csv"
SUPP = REPO / "data" / "interim" / "wang2022" / "43018_2022_475_MOESM2_ESM.xlsx"
OUT = REPO / "results" / "tables" / "patient_idh_status.csv"


def read_supp_table1() -> dict[str, dict]:
    wb = openpyxl.load_workbook(SUPP, read_only=True, data_only=True)
    ws = wb["Table 1"]
    rows = [[("" if c is None else str(c).strip()) for c in r]
            for r in ws.iter_rows(values_only=True)]
    wb.close()
    rows = [r for r in rows if any(r)]
    header, body = rows[0], rows[1:]
    idx = {name: i for i, name in enumerate(header)}
    out = {}
    for r in body:
        def get(col):
            i = idx.get(col)
            return r[i] if i is not None and i < len(r) else ""
        sid = read_id_column(get("ID"))
        if sid:
            out[sid] = {
                "supp_stage": get("Stage"),
                "supp_pair": read_id_column(get("Pair#")) or "",
                "supp_diagnosis": get("Diagnosis"),
                "idh": get("IDH"),
                "age": get("Age"),
                "sex": get("Sex"),
                "tumor_site": get("Tumor site"),
                "snRNAseq_flag": get("snRNA-seq"),
            }
    return out


def main() -> int:
    supp = read_supp_table1()
    print(f"Supplementary Table 1: {len(supp)} specimen rows")

    rows = list(csv.DictReader(open(MANIFEST)))
    hs = [r for r in rows
          if r["organism"] == "Homo sapiens" and r["assay"] == "snRNA-seq"]
    print(f"GEO manifest: {len(rows)} GSMs, {len(hs)} human snRNA-seq")

    # GEO pair# -> {timepoint: [sample_id]}
    by_pair: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    for r in hs:
        if r["patient_pair"]:
            by_pair[r["patient_pair"]][r["progression"]].append(r["sample_id"])
    both = {p: v for p, v in by_pair.items()
            if "Primary" in v and "Recurrent" in v}
    print(f"GEO pairs with both timepoints: {len(both)}")

    out_rows, joined_pairs, idh_by_pair = [], set(), {}
    for pair, tps in sorted(both.items()):
        sids = tps["Primary"] + tps["Recurrent"]
        hits = {s: supp.get(s) for s in sids}
        matched = {s: v for s, v in hits.items() if v}
        idh_vals = {v["idh"] for v in matched.values() if v["idh"]}
        if matched:
            joined_pairs.add(pair)
        idh_by_pair[pair] = idh_vals
        for s in sids:
            v = hits[s] or {}
            out_rows.append({
                "geo_pair": pair,
                "sample_id": s,
                "geo_timepoint": next(t for t in tps if s in tps[t]),
                "joined_to_supp": bool(hits[s]),
                "supp_stage": v.get("supp_stage", ""),
                "supp_pair": v.get("supp_pair", ""),
                "idh": v.get("idh", ""),
                "supp_diagnosis": v.get("supp_diagnosis", ""),
                "age": v.get("age", ""),
                "sex": v.get("sex", ""),
                "tumor_site": v.get("tumor_site", ""),
            })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {OUT.relative_to(REPO)} ({len(out_rows)} rows)\n")

    # ---------------- report ----------------
    n_specimen = len(out_rows)
    n_joined = sum(1 for r in out_rows if r["joined_to_supp"])
    print(f"specimens inside the {len(both)} pairs : {n_specimen}")
    print(f"  joined to Supplementary Table 1     : {n_joined}")
    print(f"  NOT found in Supplementary Table 1  : {n_specimen - n_joined}")
    print(f"\npairs with >=1 specimen joined        : {len(joined_pairs)}")
    fully = sum(1 for p, tps in both.items()
                if all(supp.get(s) for s in tps["Primary"] + tps["Recurrent"]))
    print(f"pairs with ALL specimens joined       : {fully}")

    print("\nIDH values observed (per specimen):")
    for k, n in Counter(r["idh"] for r in out_rows).most_common():
        print(f"   {n:>4}  {k or '(no join / blank)'}")

    wt = {p for p, v in idh_by_pair.items() if v and v <= {"IDH wildtype"}}
    print(f"\npairs where every joined specimen is IDH wildtype: {len(wt)}")
    other = {p: v for p, v in idh_by_pair.items() if v and not (v <= {"IDH wildtype"})}
    print(f"pairs with a non-wildtype or mixed IDH value     : {len(other)}")
    for p, v in sorted(other.items()):
        print(f"   pair {p}: {sorted(v)}")
    none = {p for p, v in idh_by_pair.items() if not v}
    print(f"pairs with NO IDH value available               : {len(none)}"
          + (f"  -> {sorted(none)}" if none else ""))

    print("\nGEO pair# vs Supplementary pair# (they are different numbering schemes):")
    shown = 0
    for r in out_rows:
        if r["joined_to_supp"] and shown < 5:
            print(f"   {r['sample_id']:<10} GEO pair {r['geo_pair']:<6} supp pair {r['supp_pair']}")
            shown += 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
