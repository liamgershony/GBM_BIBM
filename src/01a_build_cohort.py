#!/usr/bin/env python3
"""Materialise the discovery cohort -> results/tables/discovery_cohort.csv.

Applies clauses (a)-(c) of docs/COHORT_RULE.md, agreed by Lane 1 and Lane 2
(DEVIATIONS.md, 23 Aug 2026):

  (a) Supplementary Table 1 `Pair#` links a Primary and a Recurrent specimen.
      GEO's `pair#` is a cross-check only, never the source of pairing.
  (b) Both specimens present in GSE174554 as human snRNA-seq GSMs.
  (c) IDH-wildtype per Supplementary Table 1.

Clause (d) (>=100 usable nuclei per timepoint after QC) is NOT applied here --
it lives in 01c_clause_d_gate.py, which is where n stops being 29.

Two agreed amendments are implemented here:
  * a GSM maps to the patient/timepoint of its specimen ID with a trailing
    re-sample suffix (vN) stripped, if the stripped id is in Supp Table 1;
  * batch_key = {sample_id}__{library}, library in {batch1, batch2}.

Output is ONE ROW PER LIBRARY. Expected: 29 patients / 61 specimens / 68 libraries.
"""

from __future__ import annotations

import csv
import re
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _download_utils import read_id_column  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
MANIFEST = REPO / "results" / "tables" / "sample_manifest.csv"
SUPP = REPO / "data" / "interim" / "wang2022" / "43018_2022_475_MOESM2_ESM.xlsx"
RAW_TAR = REPO / "data" / "raw" / "GSE174554" / "GSE174554_RAW.tar"
OUT = REPO / "results" / "tables" / "discovery_cohort.csv"

EXPECT_PATIENTS, EXPECT_SPECIMENS, EXPECT_LIBRARIES = 29, 61, 68


def strip_resample_suffix(sample_id: str) -> str:
    """SF6118v2 -> SF6118.  SF6118 -> SF6118.  Agreed amendment, 23 Aug 2026."""
    return re.sub(r"v\d+$", "", sample_id)


def read_supp_table1() -> dict[str, dict]:
    wb = openpyxl.load_workbook(SUPP, read_only=True, data_only=True)
    ws = wb["Table 1"]
    rows = [[("" if c is None else str(c).strip()) for c in r]
            for r in ws.iter_rows(values_only=True)]
    wb.close()
    rows = [r for r in rows if any(r)]
    hdr, body = rows[0], rows[1:]
    idx = {n: i for i, n in enumerate(hdr)}
    get = lambda r, c: (r[idx[c]] if c in idx and idx[c] < len(r) else "")
    out = {}
    for r in body:
        sid = read_id_column(get(r, "ID"))
        if not sid:
            continue
        stage = get(r, "Stage").split()[0] if get(r, "Stage") else ""
        out[sid] = {"pair": read_id_column(get(r, "Pair#")),
                    "stage": stage, "idh": get(r, "IDH")}
    return out


def tar_index() -> dict[str, list[str]]:
    """GSM -> its 10x member paths inside GSE174554_RAW.tar."""
    idx = defaultdict(list)
    with tarfile.open(RAW_TAR, "r") as tf:
        for name in tf.getnames():
            m = re.match(r"^(GSM\d+)_", name)
            if m:
                idx[m.group(1)].append(name)
    return idx


def library_of(member: str) -> str | None:
    """batch2 members carry '_batch2_'; the primary library carries nothing."""
    base = member.rsplit("/", 1)[-1]
    if not base.endswith((".mtx.gz", ".tsv.gz")):
        return None
    return "batch2" if "_batch2_" in base else "batch1"


def main() -> int:
    supp = read_supp_table1()
    geo = [r for r in csv.DictReader(open(MANIFEST))
           if r["organism"] == "Homo sapiens" and r["assay"] == "snRNA-seq"]
    tars = tar_index()
    print(f"Supp Table 1 specimens : {len(supp)}")
    print(f"GEO human snRNA-seq    : {len(geo)}")

    # ---- clause (a): supp pairs with both timepoints, restricted to (b) ----
    geo_by_stripped = defaultdict(list)
    for r in geo:
        geo_by_stripped[strip_resample_suffix(r["sample_id"])].append(r)

    pair_stages = defaultdict(set)
    for sid, v in supp.items():
        if v["pair"] and sid in geo_by_stripped:
            pair_stages[v["pair"]].add(v["stage"])
    matched = {p for p, st in pair_stages.items()
               if "Primary" in st and "Recurrent" in st}
    print(f"clause (a)+(b) matched pairs : {len(matched)}")

    # ---- clause (c): IDH wildtype ----------------------------------------
    def pair_idh(pair):
        return {v["idh"] for s, v in supp.items()
                if v["pair"] == pair and s in geo_by_stripped and v["idh"]}
    cohort_pairs, excluded = set(), []
    for p in matched:
        idh = pair_idh(p)
        if idh == {"IDH wildtype"}:
            cohort_pairs.add(p)
        else:
            excluded.append((p, sorted(idh)))
    print(f"clause (c) IDH-wildtype pairs: {len(cohort_pairs)}")
    for p, idh in excluded:
        sids = sorted(s for s, v in supp.items()
                      if v["pair"] == p and s in geo_by_stripped)
        print(f"   EXCLUDED pair {p}: {idh}  specimens={sids}")

    # ---- emit one row per library ----------------------------------------
    rows = []
    for sid, v in sorted(supp.items()):
        if v["pair"] not in cohort_pairs:
            continue
        for g in geo_by_stripped.get(sid, []):
            members = tars.get(g["gsm"], [])
            libs = defaultdict(list)
            for m in members:
                lib = library_of(m)
                if lib:
                    libs[lib].append(m)
            for lib, mem in sorted(libs.items()):
                rows.append({
                    "patient_id": v["pair"],
                    "sample_id": g["sample_id"],
                    "supp_specimen_id": sid,
                    "gsm": g["gsm"],
                    "timepoint": v["stage"],
                    "library": lib,
                    "batch_key": f"{g['sample_id']}__{lib}",
                    "idh": v["idh"],
                    "geo_pair": g["patient_pair"] or "",      # cross-check only
                    "platform_id": g["platform_id"],
                    "organism": g["organism"],
                    "tar_members": ";".join(sorted(mem)),
                    "n_tar_members": len(mem),
                })

    rows.sort(key=lambda r: (int(r["patient_id"]) if r["patient_id"].isdigit()
                             else 10**6, r["timepoint"], r["sample_id"], r["library"]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    n_pat = len({r["patient_id"] for r in rows})
    n_spec = len({r["sample_id"] for r in rows})
    n_lib = len(rows)
    resampled = sorted({r["sample_id"] for r in rows
                        if strip_resample_suffix(r["sample_id"]) != r["sample_id"]})
    b2 = sorted({r["sample_id"] for r in rows if r["library"] == "batch2"})
    print(f"\nwrote {OUT.relative_to(REPO)}")
    print(f"  patients  : {n_pat}")
    print(f"  specimens : {n_spec}")
    print(f"  libraries : {n_lib}")
    print(f"  recovered by vN suffix-stripping : {resampled}")
    print(f"  specimens with a batch2 library  : {len(b2)} {b2}")

    # ---- exit assertions --------------------------------------------------
    assert n_pat == EXPECT_PATIENTS, f"expected {EXPECT_PATIENTS} patients, got {n_pat}"
    assert n_spec == EXPECT_SPECIMENS, f"expected {EXPECT_SPECIMENS} specimens, got {n_spec}"
    assert n_lib == EXPECT_LIBRARIES, f"expected {EXPECT_LIBRARIES} libraries, got {n_lib}"
    by_pat = defaultdict(set)
    for r in rows:
        by_pat[r["patient_id"]].add(r["timepoint"])
    bad = {p: sorted(t) for p, t in by_pat.items()
           if not {"Primary", "Recurrent"} <= t}
    assert not bad, f"patients missing a timepoint: {bad}"
    assert all(r["idh"] == "IDH wildtype" for r in rows), "non-wildtype specimen present"
    assert all(r["organism"] == "Homo sapiens" for r in rows), "non-human specimen present"
    assert all(r["platform_id"] == "GPL24676" for r in rows), "unexpected platform"
    assert all(r["n_tar_members"] == 3 for r in rows), \
        f"library without exactly 3 tar members: " \
        f"{[(r['batch_key'], r['n_tar_members']) for r in rows if r['n_tar_members'] != 3]}"
    print("\nall exit assertions passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
