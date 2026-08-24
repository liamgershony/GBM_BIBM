# Negative controls on clone identity: what the disjoint control does and does not support

**Candidate Results subsection. Day 2, 24 August 2026.**
Every number below is produced by a script in `src/` and written to a table in
`results/tables/`; provenance is given per claim.

---

## Summary

Tier C-disjoint requires a clone identity derived only from copy-number signal on
S = {chr7, chr9p, chr10}. We constructed that identity, then subjected it to three
controls before using it. **All three failed**, and the failures are informative
rather than fatal: the copy-number signal itself is sound, but the *clone
definition* built on it does not survive a null, and the genomic-compatibility
term **G** carries no within-patient information at any granularity we could
construct. We report Tier C-disjoint as effectively single-component.

This section exists because a negative control that is run and reported is worth
more than a gene list obtained without one.

---

## 1. Unsupervised clone calls are indistinguishable from noise

inferCNV converged for **21/21 patients** (`clone_catalog.csv`). Clustering the
resulting copy-number profiles with Leiden at resolution 1.0 resolved 8–17 clones
per patient in every region of S — including chr9p, which spans only 23 smoothing
windows.

Near-identical clone counts from regions carrying very different expected signal is
the signature of the clustering procedure, not of clonal structure. We tested this
against a permutation null: each copy-number window was shuffled independently
across cells, destroying cell–cell covariance while preserving every window's
marginal distribution, and the identical clustering was re-applied.

| Region | Observed clones (median) | Permutation null (median) | Ratio |
|---|---|---|---|
| chr7 | 12 | 14 | 0.86 |
| chr9p | 10 | 20 | 0.50 |
| chr10 | 13 | 14 | 0.93 |

**Every region's observed clone count falls at or below its own null.** Leiden
recovers as many clusters from shuffled copy-number values as from real ones. The
clone counts carry no information.

*Source: `src/02c_clone_validity_check.py` → `results/tables/clone_validity_check.csv`.*

---

## 2. chr9p contributes no detectable copy-number signal

The same analysis measured signal magnitude directly, independently of any
clustering: the mean inferCNV value per region in malignant nuclei against that
patient's own non-malignant reference nuclei.

| Region | Δ (malignant − reference) | Patients with the expected sign |
|---|---|---|
| chr7 (expect gain) | **+0.0622** | **21/21** |
| chr10 (expect loss) | **−0.0340** | **21/21** |
| chr9p (expect loss) | −0.0005 | **9/21** — indistinguishable from chance |

The canonical whole-chromosome events are recovered cleanly and unanimously. The
chr9p signal is absent, at chance direction. This is the expected consequence of
window smoothing: CDKN2A deletion at 9p21 is focal, on the order of a megabase,
against a ~100-gene window, and only 228 of the genes in our matrix map to 9p.

**Clone identity in this study is therefore defined by chr7 and chr10 alone.** The
paper must not describe S as three contributing regions. It is two, with chr9p
retained in the gene-exclusion rule — where it still does useful work, since Stage B
must exclude 9p genes regardless of whether 9p informed the target.

*Source: `src/02c_clone_validity_check.py`; gene counts from `src/_genome.py`.*

---

## 3. Coarsening the clone definition makes G constant

Replacing unsupervised clustering with a thresholded genotype — chr7 gain × chr10
loss, calling gain and loss against each patient's own reference nuclei at ±2 SD —
yields a well-spread partition. All four classes appear in all 21 patients, all four
hold ≥20 malignant nuclei in every patient, and the median dominant-class fraction
is 48% (range 32–68%). **The classes do not collapse.**

G nevertheless degenerates. Per protocol, G is a **set-membership test**: 1 if the
primary cell's clone is also observed among that patient's recurrent clones. With
four coarse classes and hundreds to thousands of recurrent nuclei per patient, every
class is trivially observed at recurrence.

> **G is constant (= 1) in 19 of 21 patients.**

This is structural, not empirical. Fine-grained clone definitions fail control 1;
coarse definitions make G certain. **Coarsening the clone definition makes G more
certain, not less** — the two failure modes are opposite ends of the same axis, and
no granularity between them was available from two informative regions.

*Source: `src/02d_genotype_class_degeneracy.py` → `results/tables/genotype_class_degeneracy.csv`.*

---

## 4. The two exceptions are confounded with sequencing depth

G varies in exactly two patients. In both, a primary class is *absent* at recurrence
(patient 19: `+7/−10` missing; patient 5: `../−10` missing), and both sit in the
shallow tail of recurrent-timepoint depth (113 and 150 recurrent malignant nuclei).

Depth alone does not determine it — patients 29 (95 recurrent nuclei) and 14 (113)
retain all four classes with G ≡ 1 — but the only two patients contributing
within-patient variance to G are drawn from the shallowest end of the cohort.

> **Where G varies, it is measuring how deeply the recurrent timepoint was
> sequenced, not whether a clone persisted.**

A component that fires when fewer recurrent nuclei were captured is not measuring
persistence. We do not treat these two patients as evidence that G is informative.

*Source: `results/tables/genotype_class_degeneracy.csv`.*

---

## 5. G survives z-scoring as variance, but not as information

RAS components are z-scored before weighting. Pooled across the cohort (the scope
fixed in `DEVIATIONS.md` on numerical grounds — within-patient z is a divide-by-zero
for G in 19/21 patients, undefined for 94.99% of nuclei), each component contributes
equal variance **by construction**.

| Component | Pooled variance | Median within-patient variance | Patients with zero within-patient variance |
|---|---|---|---|
| T | 0.031655 | 0.018768 | 0/21 |
| **G** | 0.017690 | **0.000000** | **19/21** |

After z-scoring, **G supplies 41.8% of the two-term partial Tier A variance while
being constant within 19 of 21 patients.** All of G's variance is *between* patients.

Stage A then regresses RAS on cell state, genotype and a **patient random
intercept** — which removes exactly the between-patient variance that is all G has.
G therefore looks like a full quarter of Tier A in the raw score and contributes
essentially nothing to the residual that Stage B predicts.

> **Tier A is four equally *weighted* components, not four equally *informative*
> ones.** The equal weighting is real and pre-registered; the equal information is
> not, and reporting the former without the latter would misdescribe the score.

*Source: `src/03a_ras_component_diagnostics.py` → `results/tables/tier_a_variance.csv`,
`results/tables/zscope_comparison.csv`.*

---

## What this means for the analysis

- **Tier C-disjoint is reported as effectively single-component.** Ab(clone) retains
  variance in every patient (median 4 distinct values, median SD 1.336); G does not.
  `RAS_C = 0.5·z(G) + 0.5·z(Ab_clone)` reduces in practice to `0.5·z(Ab_clone)` plus a
  per-patient constant in 19 of 21 patients.
- **The weights are not changed.** Reweighting after observing which component
  degenerates would be a data-dependent choice, which is exactly what fixing the
  weights a priori exists to prevent. We report the degeneracy instead.
- **H3 remains testable.** Tier C-disjoint retains within-patient variance through
  Ab(clone), so the disjoint gene list is still constructed and the permutation test
  still runs. What narrows is the *claim*: the control tests a clone-abundance
  signal derived from chr7 and chr10, not a two-component genomic score derived from
  three regions.
- **Stage A does not residualize Tier C to nothing.** Fitting RAS_C on genotype with
  a patient term leaves 44.3% of its variance (R² = 0.5569); omitting genotype leaves
  61.6% (R² = 0.3843). The "nothing left to find" stop condition is not triggered.

*Stage A figures from `src/03a_ras_component_diagnostics.py`, patient as a fixed
effect. Cell state is not yet available, so the state term is omitted; including it
could only raise these R² values, making the reported residual an upper bound.*
