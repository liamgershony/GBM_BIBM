"""CLAUDE.md §7.2: no cell from a held-out patient may influence its own prediction."""
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
FOLDS = REPO / "results" / "tables" / "stage_b_folds.csv"
ASSIGN = REPO / "results" / "tables" / "metacell_assignments.csv"


def test_every_patient_held_out_exactly_once_per_arm():
    f = pd.read_csv(FOLDS)
    for arm, g in f.groupby("arm"):
        assert g["heldout_patient"].nunique() == len(g), f"{arm}: duplicate folds"


def test_train_and_test_metacell_counts_are_disjoint_and_complete():
    f = pd.read_csv(FOLDS)
    for arm, g in f.groupby("arm"):
        tot = (g["n_train_metacells"] + g["n_test_metacells"]).unique()
        assert len(tot) == 1, f"{arm}: fold sizes inconsistent -> {tot}"


def test_metacells_never_span_patients():
    """A metacell drawn from two patients would leak across every fold."""
    a = pd.read_csv(ASSIGN)
    n = a.groupby("metacell_id")["patient_id"].nunique()
    assert (n == 1).all(), f"metacells spanning patients: {n[n > 1].index.tolist()}"


def test_metacells_never_span_timepoints():
    a = pd.read_csv(ASSIGN)
    n = a.groupby("metacell_id")["timepoint"].nunique()
    assert (n == 1).all(), f"metacells spanning timepoints: {n[n > 1].index.tolist()}"
