#!/usr/bin/env python3
"""Neftel four-state scoring -> data/processed/04_states.h5ad.

BLOCKED ON AN INPUT FILE. See §"Required input" below.

Assigns each malignant nucleus one of the four Neftel states (MES, AC, OPC, NPC)
by scoring the published metamodules with sc.tl.score_genes and taking the argmax.
Scoring is per nucleus and involves no cross-patient information, so it is safe to
run before any fold structure exists.

Required input
--------------
    data/raw/neftel_signatures/neftel_metamodules.tsv
    two columns, tab-separated, with a header:  state<TAB>gene
    state in {MES, AC, OPC, NPC}

CLAUDE.md §4 lists the Neftel signatures as "Open — published supplementary
tables". That is not true for automated retrieval, and this was verified on
2026-08-24:
  * Europe PMC: "Article with id PMC6703186 is not open access one"
  * PMC OA service: "identifier 'PMC6703186' is not Open Access"
  * Neftel's own GEO deposit GSE131928 carries only GSE131928_RAW.tar and a
    sample-name spreadsheet — no metamodule gene lists.

The gene lists are therefore NOT auto-downloadable and must be placed at the path
above by hand, from Neftel et al. 2019 Cell, Table S. This script refuses to run
without them and will never substitute a guessed gene list: an invented signature
would silently define the cell-state covariate that Stage A residualises on.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

REPO = Path(__file__).resolve().parent.parent
IN_H5 = REPO / "data" / "processed" / "02_integrated.h5ad"
SIG = REPO / "data" / "raw" / "neftel_signatures" / "neftel_metamodules.tsv"
CONF = REPO / "configs" / "pipeline_config.yaml"
OUT_H5 = REPO / "data" / "processed" / "04_states.h5ad"
OUT_CSV = REPO / "results" / "tables" / "state_assignments.csv"

EXPECTED_STATES = {"MES", "AC", "OPC", "NPC"}


def load_signatures() -> dict[str, list[str]]:
    if not SIG.exists():
        raise SystemExit(
            f"MISSING REQUIRED INPUT: {SIG}\n\n"
            "The Neftel metamodules are not auto-downloadable — Neftel et al. 2019\n"
            "is not open access on PMC and GSE131928 does not deposit the gene\n"
            "lists. Place a two-column TSV (state<TAB>gene, header included) at the\n"
            "path above, taken from the paper's supplementary tables.\n\n"
            "This script will not proceed with a guessed signature: the state call\n"
            "becomes Stage A's fixed effect, so an invented gene list would silently\n"
            "define the covariate the whole adjustment rests on.")
    df = pd.read_csv(SIG, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    assert "state" in cols and "gene" in cols, \
        f"{SIG} must have columns 'state' and 'gene', found {list(df.columns)}"
    df = df.rename(columns={cols["state"]: "state", cols["gene"]: "gene"})
    df["state"] = df["state"].str.strip().str.upper()
    sigs = {s: sorted(set(g.strip() for g in d["gene"] if isinstance(g, str)))
            for s, d in df.groupby("state")}
    unknown = set(sigs) - EXPECTED_STATES
    missing = EXPECTED_STATES - set(sigs)
    assert not unknown, f"unexpected states in signature file: {sorted(unknown)}"
    assert not missing, f"signature file is missing states: {sorted(missing)}"
    return sigs


def main() -> int:
    conf = yaml.safe_load(open(CONF))
    seed = conf["seed"]["master"]
    sigs = load_signatures()
    for s, g in sigs.items():
        print(f"  {s:<4} {len(g):>4} genes")

    adata = ad.read_h5ad(IN_H5)
    print(f"{IN_H5.name}: {adata.n_obs:,} nuclei")
    # score on the full log-normalised gene space, not the HVG subset
    full = adata.raw.to_adata() if adata.raw is not None else adata
    full.obs = adata.obs.copy()

    mal = full[full.obs["is_malignant"].values].copy()
    print(f"malignant nuclei scored: {mal.n_obs:,}")

    cols = []
    for state, genes in sorted(sigs.items()):
        present = [g for g in genes if g in mal.var_names]
        cov = len(present) / len(genes) if genes else 0.0
        print(f"  {state}: {len(present)}/{len(genes)} signature genes present "
              f"({cov:.0%})")
        assert cov >= 0.5, (
            f"only {cov:.0%} of the {state} signature is present in the data — "
            f"check the gene symbol namespace before trusting any state call")
        sc.tl.score_genes(mal, present, score_name=f"score_{state}",
                          random_state=seed)
        cols.append(f"score_{state}")

    scores = mal.obs[cols].to_numpy()
    winner = np.array([c.replace("score_", "") for c in cols])[scores.argmax(axis=1)]
    mal.obs["cell_state"] = pd.Categorical(winner, categories=sorted(EXPECTED_STATES))
    # margin between best and runner-up: low margin = ambiguous assignment
    srt = np.sort(scores, axis=1)
    mal.obs["state_margin"] = srt[:, -1] - srt[:, -2]

    out = mal.obs[["patient_id", "sample_id", "timepoint", "cell_state",
                   "state_margin"] + cols].copy()
    out.index.name = "nucleus_id"
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV)
    mal.write_h5ad(OUT_H5, compression="gzip")
    print(f"\nwrote {OUT_CSV.relative_to(REPO)} and {OUT_H5.relative_to(REPO)}")

    print("\nstate composition overall:")
    print(mal.obs["cell_state"].value_counts().to_string())
    print("\nper-patient state counts (CLAUDE.md §6.2 evaluability input):")
    tab = pd.crosstab(mal.obs["patient_id"], mal.obs["cell_state"])
    print(tab.to_string())
    floor = len(tab) // 2 + 1
    print(f"\nevaluability floor at n={len(tab)}: floor(n/2)+1 = {floor} patients "
          f"with >=20 nuclei in a state")
    for st in tab.columns:
        n_ok = int((tab[st] >= 20).sum())
        print(f"  {st:<4} {n_ok:>3}/{len(tab)} patients >=20 nuclei  "
              f"-> {'EVALUABLE' if n_ok >= floor else 'EVALUABILITY-FAILED'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
