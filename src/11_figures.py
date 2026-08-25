#!/usr/bin/env python3
"""Figures 1 and 2, sized and styled for IEEE two-column print in greyscale.

Greyscale constraint: identity is never carried by tone alone. Every series also
carries a hatch texture and a direct label, so the figures survive photocopying
and colour-blind readers alike.

Figure 1  cohort flow + RAS distribution by cell state
Figure 2  the three-variant H3 panel -- observed Jaccard against its permutation
          null, with the target correlation annotated, so the reader sees not just
          that overlap exists but what it is made of.
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, yaml
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

REPO = Path(__file__).resolve().parent.parent
# RUN PROVENANCE. The paper reports RUN 1, snapshotted in results/_repro_baseline/.
# results/tables/ currently holds RUN 2 from the reproducibility check, whose
# Stage B numbers differ (see RESULTS_SUMMARY.md §7b). Both figures must read from
# ONE run; mixing them would put run-1 gene counts beside run-2 Jaccards.
RUN_LABEL = "run 1 (results/_repro_baseline/)"
SRC = REPO / "results" / "_repro_baseline"
if not (SRC / "cohort_flow.csv").exists():
    SRC = REPO / "results" / "tables"
    RUN_LABEL = "results/tables/ (baseline snapshot absent)"
FIG = REPO / "results" / "figures"

IEEE_W = 7.16          # full two-column width, inches
INK, MUTED, GRID = "#111111", "#555555", "#CCCCCC"
# Four greys with wide, even lightness separation; each also gets a hatch.
TONES = ["#1a1a1a", "#5c5c5c", "#9a9a9a", "#d0d0d0"]
HATCH = ["", "///", "...", "xxx"]

plt.rcParams.update({
    "font.size": 7.5, "axes.labelsize": 7.5, "axes.titlesize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.labelcolor": INK, "figure.dpi": 400, "savefig.dpi": 400,
    "axes.grid": False, "font.family": "DejaVu Sans",
})


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=2.5, width=0.6)


# ----------------------------------------------------------------- Figure 1
def figure1():
    flow = pd.read_csv(SRC / "cohort_flow.csv")
    ras = pd.read_csv(SRC / "ras_scores.csv")

    # read from cohort_flow.csv rather than hardcoding, so the figure cannot
    # drift from the table it claims to depict
    want = {"GSMs in GSE174554": "GSMs in GSE174554",
            "human snRNA-seq GSMs": "Human snRNA-seq",
            "matched pairs (clauses a+b)": "Matched pairs",
            "IDH-wildtype pairs (clause c)": "IDH-wildtype pairs"}
    lut = dict(zip(flow["stage"], flow["count"]))
    steps = [(lbl, int(lut[k])) for k, lbl in want.items()]
    dkey = [k for k in lut if k.startswith("patients passing clause (d)")][0]
    steps.append(("Pass clause (d)", int(lut[dkey])))
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(IEEE_W, 2.55), gridspec_kw={"width_ratios": [1.05, 1]})

    # (a) cohort funnel -- magnitude by length, one axis, direct labels
    y = np.arange(len(steps))[::-1]
    vals = [v for _, v in steps]
    ax1.barh(y, vals, height=0.6, color=TONES[1], edgecolor=INK, linewidth=0.5)
    for yi, (lab, v) in zip(y, steps):
        ax1.text(v + 2.5, yi, str(v), va="center", ha="left", fontsize=7.5,
                 color=INK, fontweight="bold")
    ax1.set_yticks(y); ax1.set_yticklabels([s for s, _ in steps])
    ax1.set_xlim(0, 128); ax1.set_xlabel("samples / patients")
    ax1.set_title("(a) Cohort selection", loc="left", fontweight="bold")
    style(ax1)
    # annotate the two exclusions that matter
    ax1.annotate("1 IDH-mutant pair", xy=(29, y[3]), xytext=(52, y[3] - 0.42),
                 fontsize=6.2, color=MUTED,
                 arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED))
    ax1.annotate("8 fail $\\geq$100 nuclei", xy=(21, y[4]), xytext=(52, y[4] - 0.42),
                 fontsize=6.2, color=MUTED,
                 arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED))

    # (b) RAS distribution by state -- box + strip, tone AND hatch AND label
    order = ["AC", "MES", "NPC", "OPC"]
    data = [ras.loc[ras.cell_state == s, "ras_tier_a_reduced"].values for s in order]
    bp = ax2.boxplot(data, positions=np.arange(len(order)), widths=0.55,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color=INK, lw=1.1),
                     whiskerprops=dict(color=MUTED, lw=0.6),
                     capprops=dict(color=MUTED, lw=0.6),
                     boxprops=dict(edgecolor=INK, lw=0.5))
    for patch, tone, h in zip(bp["boxes"], TONES, HATCH):
        patch.set_facecolor(tone); patch.set_hatch(h)
    for i, s in enumerate(order):
        n = int((ras.cell_state == s).sum())
        ax2.text(i, ax2.get_ylim()[1], f"n={n:,}", ha="center", va="bottom",
                 fontsize=6.2, color=MUTED)
    ax2.axhline(0, color=GRID, lw=0.6, zorder=0)
    ax2.set_xticks(np.arange(len(order))); ax2.set_xticklabels(order)
    ax2.set_ylabel("RAS Tier A-reduced")
    ax2.set_xlabel("Neftel cell state")
    ax2.set_title("(b) RAS by cell state", loc="left", fontweight="bold")
    style(ax2)

    fig.tight_layout(pad=0.5, w_pad=1.6)
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"figure1_cohort_and_ras.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"  figure1: {[f'{s}={v}' for s, v in steps]}")
    print(f"  RAS by state medians: "
          f"{ {s: round(float(np.median(d)), 3) for s, d in zip(order, data)} }")


# ----------------------------------------------------------------- Figure 2
def figure2():
    cc = pd.read_csv(SRC / "circularity_check.csv")
    conf = yaml.safe_load(open(REPO / "configs" / "pipeline_config.yaml"))
    seed = conf["seed"]["master"]; n_perm = conf["h3_circularity"]["n_permutations"]
    alpha = conf["h3_circularity"]["alpha"]
    # eligible universes derived, not hardcoded; both are run-invariant
    import anndata as ad_
    from _genome import annotate_var as _av
    _mc = ad_.read_h5ad(REPO / "data" / "interim" / "metacell_expression.h5ad")
    _ann = _av(_mc.var_names, conf["disjoint_set_S"]["regions"])
    nA = int(_mc.n_vars)
    nC = int((~_ann["in_disjoint_set_S"].fillna(False).values).sum())
    assert (nA, nC) == (33694, 31744), f"universes changed: {(nA, nC)}"

    def jac(a, b):
        a, b = set(a), set(b)
        return len(a & b) / len(a | b) if (a or b) else 0.0

    fig, axes = plt.subplots(1, 3, figsize=(IEEE_W, 2.5), sharey=True, sharex=True)
    rng = np.random.default_rng(seed)
    xmax = float(cc["jaccard_observed"].max()) * 1.45
    for k, (ax, (_, r)) in enumerate(zip(axes, cc.iterrows())):
        null = np.empty(n_perm)
        for i in range(n_perm):
            ra = rng.choice(nA, size=int(r.n_genes_tierA), replace=False)
            rc = rng.choice(nC, size=int(r.n_genes_tierC), replace=False)
            null[i] = jac(ra, rc)
        # The null is highly DISCRETE -- with lists of 3 and 11 genes only a few
        # Jaccard values are attainable at all. A histogram bins a point mass into
        # a wide bar and reads as a broad distribution, which is false. Stems at
        # the values actually attained tell the truth.
        vals, cnts = np.unique(np.round(null, 6), return_counts=True)
        ax.vlines(vals, 0, cnts, color=TONES[2], lw=2.2,
                  label="permutation null" if k == 0 else None)
        ax.plot(vals, cnts, "o", ms=2.6, color=TONES[1], zorder=3)
        ax.axvline(r.jaccard_observed, color=INK, lw=1.7, zorder=4,
                   label="observed" if k == 0 else None)
        ax.set_yscale("symlog", linthresh=1)
        ax.set_xlabel("Jaccard index")
        role = "PRIMARY" if r.variant == "v1" else "sensitivity"
        ax.set_title(f"({'abc'[k]}) {r.variant} — {role}", loc="left",
                     fontweight="bold", fontsize=8)
        sig = "*" if r.p_value < alpha else "n.s."
        txt = (f"corr(targets) {r.target_correlation:.3f}\n"
               f"|A| {int(r.n_genes_tierA)}   |C| {int(r.n_genes_tierC)}   "
               f"overlap {int(r.n_overlap)}\n"
               f"J = {r.jaccard_observed:.4f}\n"
               f"p = {r.p_value:.4f}  {sig}")
        ax.text(0.96, 0.97, txt, transform=ax.transAxes, ha="right", va="top",
                fontsize=6.1, color=INK, linespacing=1.4, zorder=6,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=GRID, lw=0.5,
                          alpha=1.0))
        if r.n_overlap == 0:
            ax.annotate("observed lies\ninside the null",
                        xy=(float(r.jaccard_observed), 1.6), xytext=(0.30, 0.30),
                        textcoords="axes fraction", fontsize=6.0, color=MUTED,
                        ha="left", arrowprops=dict(arrowstyle="->", lw=0.5,
                                                   color=MUTED))
        style(ax)
        ax.set_xlim(-0.0025, xmax)
    axes[0].set_ylabel("permutations (symlog)")
    axes[0].legend(frameon=False, loc="lower left", fontsize=6.3,
                   handlelength=1.5, borderaxespad=0.3,
                   bbox_to_anchor=(0.02, 0.03))
    fig.suptitle("H3: gene-list overlap against a 1,000-permutation null "
                 f"($\\alpha$ = {alpha}); v1 is primary by pre-registration",
                 fontsize=8, y=1.03, x=0.01, ha="left")
    fig.tight_layout(pad=0.5, w_pad=0.9)
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"figure2_h3_permutation.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("  figure2 rows:")
    print(cc[["variant", "target_correlation", "n_genes_tierA", "n_genes_tierC",
              "n_overlap", "jaccard_observed", "p_value"]].to_string(index=False))


def main() -> int:
    FIG.mkdir(parents=True, exist_ok=True)
    print(f"RUN PROVENANCE: every value in both figures comes from {RUN_LABEL}")
    print(f"source tables : {SRC.relative_to(REPO)}")
    figure1(); figure2()
    print(f"\nwrote {FIG.relative_to(REPO)}/figure1_cohort_and_ras.[pdf|png]")
    print(f"wrote {FIG.relative_to(REPO)}/figure2_h3_permutation.[pdf|png]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
