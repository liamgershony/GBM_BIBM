"""CLAUDE.md §7.3: chromosome disjointness is enforced by TEST, not by care.

No gene on chr7, chr9p or chr10 may appear in any Tier C-disjoint feature matrix
or selected gene list. chr9 is ARM-aware: chr9q genes remain eligible.
"""
import sys
from pathlib import Path
import pandas as pd, yaml, anndata as ad

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from _genome import annotate_var  # noqa: E402

CONF = yaml.safe_load(open(REPO / "configs" / "pipeline_config.yaml"))
REGIONS = CONF["disjoint_set_S"]["regions"]
GL = REPO / "results" / "gene_lists"
TIER_C_ARMS = ["v1_tierC", "v2_tierC", "v3_tierC"]


def _ann(genes):
    return annotate_var(pd.Index(genes), REGIONS)


def test_no_disjoint_genes_in_tier_c_lists():
    for arm in TIER_C_ARMS:
        for t in (30, 50, 80):
            f = GL / f"{arm}_genes_{t}.csv"
            if not f.exists():
                continue
            genes = pd.read_csv(f)["gene"].dropna().tolist()
            if not genes:
                continue
            a = _ann(genes)
            bad = a.index[a["in_disjoint_set_S"].fillna(False).values].tolist()
            assert not bad, f"{arm}@{t}%: chr7/chr9p/chr10 genes present: {bad}"


def test_chr9q_remains_eligible():
    """The exclusion must be arm-level; excluding all of chr9 would be wrong."""
    mc = ad.read_h5ad(REPO / "data" / "interim" / "metacell_expression.h5ad")
    a = _ann(mc.var_names)
    q9 = (a["chromosome"] == "chr9") & (a["arm"] == "q")
    assert q9.sum() > 0, "no chr9q genes found at all -- annotation is broken"
    assert not a.loc[q9, "in_disjoint_set_S"].fillna(False).any(), \
        "chr9q genes marked as inside the disjoint set"


def test_disjoint_set_covers_expected_regions():
    mc = ad.read_h5ad(REPO / "data" / "interim" / "metacell_expression.h5ad")
    a = _ann(mc.var_names)
    ins = a["in_disjoint_set_S"].fillna(False)
    for chrom, arm in (("chr7", None), ("chr9", "p"), ("chr10", None)):
        m = (a["chromosome"] == chrom)
        if arm:
            m &= (a["arm"] == arm)
        assert ins[m].all(), f"{chrom}{arm or ''} genes missing from the disjoint set"
