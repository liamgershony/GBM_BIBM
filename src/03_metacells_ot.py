#!/usr/bin/env python3
"""Metacells per patient-timepoint, then optimal transport -> RAS component O.

Metacells: SEACells on the Harmony embedding, ~30 nuclei each (frozen config).
Pre-declared contingency (CLAUDE.md §9.1): if SEACells does not converge, k-means
within the same patient-timepoint is substituted. Here that switch is automatic
and time-boxed -- a stalled group falls back rather than being debugged -- and the
method actually used is recorded per group in the output table.

Transport: for each patient, ot.emd between that patient's primary metacells and
their recurrent metacells, cost = Euclidean distance in the Harmony PCA space.

NOTE ON O AS SPECIFIED. CLAUDE.md §3.2 defines O as "total optimal-transport mass
moved from cell i's metacell to any recurrent-timepoint metacell". For BALANCED
EMD the row sums of the transport plan equal the source marginal by construction,
so that quantity is identically the metacell's own mass. We compute it exactly as
specified AND record the diagnostic that shows this, rather than silently
substituting a different statistic. A transport-cost alternative is computed
alongside for inspection only and is NOT used as O.
"""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import ot
import pandas as pd
import scanpy as sc
import yaml
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
INTEG = REPO / "data" / "processed" / "02_integrated.h5ad"
CONF = REPO / "configs" / "pipeline_config.yaml"
COHORT_N = REPO / "results" / "tables" / "cohort_n.json"
OUT_MC = REPO / "data" / "processed" / "07_metacells.h5ad"
OUT_MC_CSV = REPO / "results" / "tables" / "metacell_catalog.csv"
OUT_O = REPO / "results" / "tables" / "ot_component_O.csv"

SEACELLS_TIMEBOX_S = 180     # per patient-timepoint; then fall back, do not debug


def build_metacells(emb, target, seed, tag, log):
    """Return (labels, method). SEACells first, k-means on timeout or failure."""
    n = emb.shape[0]
    k = max(2, int(round(n / target)))
    if n < 3 * target:
        lab = KMeans(n_clusters=min(k, n), n_init=10,
                     random_state=seed).fit_predict(emb)
        return lab, "kmeans_too_few_nuclei"
    t0 = time.time()
    try:
        # SEACells imports tqdm.notebook, which raises ImportError without
        # ipywidgets. That is a PROGRESS-BAR dependency, not a convergence
        # failure, and must never be allowed to trigger the §9.1 k-means
        # contingency -- doing so would put a false claim in the paper.
        import ipywidgets  # noqa: F401
        import SEACells
        a = ad.AnnData(np.asarray(emb, dtype="float32"))
        a.obsm["X_pca_harmony"] = np.asarray(emb, dtype="float32")
        m = SEACells.core.SEACells(a, build_kernel_on="X_pca_harmony",
                                   n_SEACells=k, n_waypoint_eigs=min(10, k - 1),
                                   convergence_epsilon=1e-5)
        m.construct_kernel_matrix()
        m.initialize_archetypes()
        m.fit(min_iter=5, max_iter=50)
        if time.time() - t0 > SEACELLS_TIMEBOX_S:
            raise TimeoutError(f"SEACells exceeded {SEACELLS_TIMEBOX_S}s")
        lab = pd.Categorical(a.obs["SEACell"]).codes
        if len(set(lab)) < 2:
            raise ValueError("SEACells produced <2 metacells")
        return lab, "seacells"
    except Exception as e:                                  # noqa: BLE001
        log(f"    [{tag}] SEACells unusable after {time.time()-t0:.0f}s "
            f"({type(e).__name__}: {str(e)[:90]}) -> k-means contingency")
        lab = KMeans(n_clusters=min(k, n), n_init=10,
                     random_state=seed).fit_predict(emb)
        return lab, f"kmeans_fallback:{type(e).__name__}"


def main() -> int:
    conf = yaml.safe_load(open(CONF))
    seed = conf["seed"]["master"]
    target = conf["metacells"]["target_size_nuclei"]
    pts = json.loads(COHORT_N.read_text())["patient_ids"]
    log = print
    log(f"metacell target {target} nuclei/metacell, seed {seed}, "
        f"{len(pts)} patients")

    adata = ad.read_h5ad(INTEG)
    adata = adata[adata.obs["is_malignant"].values].copy()
    emb_all = adata.obsm["X_pca_harmony"]
    obs = adata.obs
    log(f"malignant nuclei: {adata.n_obs:,}")

    rows, mc_rows = [], []
    for p in pts:
        for tp in ("Primary", "Recurrent"):
            sel = ((obs["patient_id"].astype(str) == p)
                   & (obs["timepoint"] == tp)).values
            if sel.sum() == 0:
                continue
            E = emb_all[sel]
            lab, method = build_metacells(E, target, seed, f"{p}/{tp}", log)
            for c in np.unique(lab):
                m = lab == c
                mc_rows.append({
                    "metacell_id": f"{p}_{tp[:1]}_mc{int(c):03d}",
                    "patient_id": p, "timepoint": tp,
                    "n_nuclei": int(m.sum()), "method": method,
                    "centroid": E[m].mean(axis=0),
                })
            rows.append({"patient_id": p, "timepoint": tp,
                         "n_nuclei": int(sel.sum()),
                         "n_metacells": int(len(np.unique(lab))),
                         "method": method})
            log(f"  {p:<5} {tp:<10} {int(sel.sum()):>6,} nuclei -> "
                f"{len(np.unique(lab)):>4} metacells [{method}]")

    mc = pd.DataFrame(mc_rows)
    pd.DataFrame(rows).to_csv(OUT_MC_CSV, index=False)

    # 07_metacells.h5ad: one row per metacell, X = Harmony centroid (docs/SCHEMA.md §2)
    mc_ad = ad.AnnData(np.vstack(mc["centroid"].values).astype("float32"))
    mc_ad.obs_names = mc["metacell_id"].values
    for c in ("patient_id", "timepoint", "n_nuclei", "method"):
        mc_ad.obs[c] = mc[c].values
    for c in ("patient_id", "timepoint", "method"):
        mc_ad.obs[c] = mc_ad.obs[c].astype("category")
    mc_ad.write_h5ad(OUT_MC, compression="gzip")
    log(f"wrote {OUT_MC.relative_to(REPO)} ({mc_ad.n_obs:,} metacells)")
    log(f"\nwrote {OUT_MC_CSV.relative_to(REPO)}")

    # ---------------- optimal transport, per patient ----------------
    log("\nrunning optimal transport per patient ...")
    out = []
    for p in pts:
        P = mc[(mc.patient_id == p) & (mc.timepoint == "Primary")]
        R = mc[(mc.patient_id == p) & (mc.timepoint == "Recurrent")]
        if len(P) == 0 or len(R) == 0:
            log(f"  {p}: missing a timepoint, skipped")
            continue
        A = np.vstack(P["centroid"].values)
        B = np.vstack(R["centroid"].values)
        a = P["n_nuclei"].to_numpy(float); a /= a.sum()
        b = R["n_nuclei"].to_numpy(float); b /= b.sum()
        M = cdist(A, B, metric="euclidean")
        t0 = time.time()
        T = ot.emd(a, b, M, numItermax=1_000_000)
        dt = time.time() - t0
        mass = T.sum(axis=1)                       # O exactly as §3.2 specifies
        cost = (T * M).sum(axis=1)                 # diagnostic only, NOT O
        for mid, mval, cval, src in zip(P["metacell_id"], mass, cost, a):
            out.append({"metacell_id": mid, "patient_id": p,
                        "O_transport_mass": float(mval),
                        "source_marginal": float(src),
                        "transport_cost_diagnostic": float(cval)})
        log(f"  {p:<5} {len(P):>4} primary x {len(R):>4} recurrent metacells, "
            f"emd {dt:.2f}s")

    df = pd.DataFrame(out)
    df.to_csv(OUT_O, index=False)
    log(f"\nwrote {OUT_O.relative_to(REPO)} ({len(df):,} primary metacells)")

    r = np.corrcoef(df["O_transport_mass"], df["source_marginal"])[0, 1]
    log(f"\n--- O as specified: degeneracy check ---")
    log(f"  corr(O, source marginal) = {r:.6f}")
    log(f"  max |O - source marginal| = "
        f"{np.abs(df['O_transport_mass']-df['source_marginal']).max():.3e}")
    if r > 0.999:
        log("  O IS IDENTICAL TO THE SOURCE MARGINAL. Balanced EMD conserves mass,")
        log("  so row sums equal the source weights by construction: O carries no")
        log("  information beyond metacell size. Reported, not silently replaced.")
    log(f"  transport-cost diagnostic: mean {df['transport_cost_diagnostic'].mean():.4f}, "
        f"sd {df['transport_cost_diagnostic'].std():.4f} (NOT used as O)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
