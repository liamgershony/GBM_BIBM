#!/usr/bin/env python3
"""Compare two pairing schemes for GSE174554 and write pairing_comparison.csv.

Two independent records claim to pair primary and recurrent specimens:

  GEO `pair#`     -- a characteristic on each GSM, i.e. the submission form.
  Supplementary Table 1 `Pair#` -- the authors' own record in the published paper.

GEO's series summary claims 40 matched IDH-wildtype pairs; grouping human
snRNA-seq GSMs by GEO pair# yields only 30, and 12 human snRNA-seq GSMs carry no
usable pair# at all (literal "NA" or absent). This script re-derives pairing from
Supplementary Table 1, joined to GEO by specimen ID, and reports both schemes
side by side.

Transcription and join only. NO COHORT RULE IS APPLIED and nothing is filtered.

Usage:  python3 src/00f_compare_pairing.py
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import openpyxl

REPO = Path(__file__).resolve().parent.parent
MANIFEST = REPO / "results" / "tables" / "sample_manifest.csv"
SUPP = REPO / "data" / "interim" / "wang2022" / "43018_2022_475_MOESM2_ESM.xlsx"
OUT = REPO / "results" / "tables" / "pairing_comparison.csv"


def read_supp_table1() -> dict[str, dict]:
    wb = openpyxl.load_workbook(SUPP, read_only=True, data_only=True)
    ws = wb["Table 1"]
    rows = [[("" if c is None else str(c).strip()) for c in r]
            for r in ws.iter_rows(values_only=True)]
    wb.close()
    rows = [r for r in rows if any(r)]
    header, body = rows[0], rows[1:]
    idx = {n: i for i, n in enumerate(header)}
    get = lambda r, c: (r[idx[c]] if c in idx and idx[c] < len(r) else "")
    return {get(r, "ID"): {"stage": get(r, "Stage"),
                           "pair": clean_pair(get(r, "Pair#")),
                           "pair_raw": get(r, "Pair#"),
                           "idh": get(r, "IDH"), "diagnosis": get(r, "Diagnosis")}
            for r in body if get(r, "ID")}


def clean_pair(value: str) -> str:
    """Supplementary Table 1 uses a literal 'NA' in Pair#, exactly as GEO does.

    Left as a string it fuses every unpaired specimen into one fabricated patient
    -- here 3 primaries and 7 recurrents -- which then looks like a legitimate
    matched pair and inflates the count. Same trap as GEO's pair#; both records
    must be normalised or neither.
    """
    return "" if value.strip().upper() in {"NA", "N/A", "NONE", ""} else value.strip()


def norm_stage(s: str) -> str:
    s = s.lower()
    if s.startswith("primary"):
        return "Primary"
    if s.startswith("recurren"):
        return "Recurrent"
    return s or "(blank)"


def main() -> int:
    supp = read_supp_table1()
    geo_rows = list(csv.DictReader(open(MANIFEST)))
    hs = [r for r in geo_rows
          if r["organism"] == "Homo sapiens" and r["assay"] == "snRNA-seq"]
    geo_by_sid = {r["sample_id"]: r for r in hs}

    print(f"Supplementary Table 1 specimens : {len(supp)}")
    print(f"GEO human snRNA-seq GSMs        : {len(hs)}")
    print(f"  with a usable GEO pair#       : {sum(1 for r in hs if r['patient_pair'])}")
    print(f"  without                       : {sum(1 for r in hs if not r['patient_pair'])}")

    # ---- scheme A: GEO pair# ------------------------------------------------
    geo_groups = defaultdict(lambda: defaultdict(list))
    for r in hs:
        if r["patient_pair"]:
            geo_groups[r["patient_pair"]][r["progression"]].append(r["sample_id"])
    geo_matched = {p for p, t in geo_groups.items()
                   if "Primary" in t and "Recurrent" in t}

    # ---- scheme B: Supplementary Pair#, restricted to GEO human snRNA-seq ----
    supp_groups = defaultdict(lambda: defaultdict(list))
    for sid, v in supp.items():
        if sid in geo_by_sid and v["pair"]:
            supp_groups[v["pair"]][norm_stage(v["stage"])].append(sid)
    supp_matched = {p for p, t in supp_groups.items()
                    if "Primary" in t and "Recurrent" in t}

    print(f"\nmatched pairs -- GEO pair# scheme          : {len(geo_matched)}")
    print(f"matched pairs -- Supplementary Pair# scheme: {len(supp_matched)}")

    # ---- per-patient comparison rows ---------------------------------------
    sid_geo_pair = {r["sample_id"]: r["patient_pair"] for r in hs}
    sid_supp_pair = {s: v["pair"] for s, v in supp.items()}

    rows = []
    for sp in sorted(supp_matched, key=lambda x: (len(x), x)):
        t = supp_groups[sp]
        sids = t.get("Primary", []) + t.get("Recurrent", [])
        geo_pairs = sorted({sid_geo_pair.get(s, "") for s in sids} - {""})
        idh = sorted({supp[s]["idh"] for s in sids if supp[s]["idh"]})
        rows.append({
            "supp_pair": sp,
            "geo_pair": ";".join(geo_pairs) or "(none)",
            "in_geo_matched": any(g in geo_matched for g in geo_pairs),
            "supp_primary": ";".join(sorted(t.get("Primary", []))),
            "supp_recurrent": ";".join(sorted(t.get("Recurrent", []))),
            "geo_primary": ";".join(sorted(
                s for g in geo_pairs for s in geo_groups[g].get("Primary", []))),
            "geo_recurrent": ";".join(sorted(
                s for g in geo_pairs for s in geo_groups[g].get("Recurrent", []))),
            "idh": ";".join(idh) or "(none)",
            "n_specimens": len(sids),
        })

    # GEO-matched pairs that the supplementary scheme does not reproduce
    supp_covered_geo = {g for r in rows for g in r["geo_pair"].split(";") if g}
    geo_only = sorted(geo_matched - supp_covered_geo)
    for gp in geo_only:
        t = geo_groups[gp]
        sids = t.get("Primary", []) + t.get("Recurrent", [])
        rows.append({
            "supp_pair": "(none)",
            "geo_pair": gp,
            "in_geo_matched": True,
            "supp_primary": "", "supp_recurrent": "",
            "geo_primary": ";".join(sorted(t.get("Primary", []))),
            "geo_recurrent": ";".join(sorted(t.get("Recurrent", []))),
            "idh": ";".join(sorted({supp[s]["idh"] for s in sids
                                    if s in supp and supp[s]["idh"]})) or "(none)",
            "n_specimens": len(sids),
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT.relative_to(REPO)} ({len(rows)} rows)")

    # ---- report -------------------------------------------------------------
    recovered = [r for r in rows if r["supp_pair"] != "(none)" and not r["in_geo_matched"]]
    print(f"\npairs the supplementary scheme RECOVERS that GEO's loses : {len(recovered)}")
    for r in recovered:
        print(f"   supp pair {r['supp_pair']:<6} geo_pair={r['geo_pair']:<8} "
              f"P={r['supp_primary']:<22} R={r['supp_recurrent']:<22} {r['idh']}")

    print(f"\npairs matched in GEO but NOT reproduced by the supplementary table: {len(geo_only)}")
    for gp in geo_only:
        t = geo_groups[gp]
        print(f"   geo pair {gp}: P={t.get('Primary')} R={t.get('Recurrent')}")

    # ---- the pair#-NA GSMs --------------------------------------------------
    na = [r for r in hs if not r["patient_pair"]]
    print(f"\n=== the {len(na)} human snRNA-seq GSMs with no usable GEO pair# ===")
    rec_na = 0
    for r in sorted(na, key=lambda x: x["sample_id"]):
        sid = r["sample_id"]
        in_supp = sid in supp
        sp = supp[sid]["pair"] if in_supp else ""
        sp_raw = supp[sid]["pair_raw"] if in_supp else ""
        in_matched = bool(sp) and sp in supp_matched
        rec_na += in_matched
        if not in_supp:
            why = "ABSENT from Supp Table 1"
        elif not sp:
            why = f"in Supp Table 1 but Pair# is literal '{sp_raw}'"
        elif not in_matched:
            why = f"Supp Pair# {sp}, but that pair has no partner timepoint in GEO snRNA-seq"
        else:
            why = "RECOVERED into a matched pair"
        print(f"   {sid:<12} {r['progression']:<10} geo_pair={r['patient_pair_raw'] or '(absent)':<4} "
              f"| {why}")
    print(f"\n   recovered into a matched pair by the supplementary scheme: {rec_na} of {len(na)}")

    idh_wt = [r for r in rows if r["idh"] == "IDH wildtype"]
    print(f"\nsupplementary-scheme matched pairs that are IDH wildtype: "
          f"{sum(1 for r in rows if r['supp_pair'] != '(none)' and r['idh'] == 'IDH wildtype')}"
          f" of {len(supp_matched)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
