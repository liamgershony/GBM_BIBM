# SCHEMA.md — Lane 1 → Lane 2 handoff contract

**Status: DRAFT, 23 Aug 2026 (Day 0).** Arnav is coding Stage A/B against simulated
data until Day 3. This file is what the simulator must match, so that modelling
code is already debugged when the real object lands.

Every field is tagged:

- **CONFIRMED** — read directly from a downloaded file. Will not change.
- **PROVISIONAL** — inferred from the protocol or from filename conventions, and
  **not yet verified against GSE174554's authoritative sample metadata.** May change.
  Code defensively against these; do not hardcode their values.

Changes to a CONFIRMED field require agreement from both lanes. Changes to a
PROVISIONAL field are expected and will be announced.

---

## 0. What has actually been read so far

Downloaded and checksummed (see `data/raw/PROVENANCE.md`):

| Fact | Value | Source |
|---|---|---|
| Tumour/normal annotation columns | `Sample#_Barcode`, `Tumor_Normal_annotation` | CONFIRMED — file header |
| Annotation values | `Tumor`, `Normal` | CONFIRMED |
| Annotated nuclei | 254,288 (126,773 Tumor / 127,515 Normal) | CONFIRMED — row counts |
| Distinct sample IDs in that file | 100 | CONFIRMED |
| Barcode field format | `<SampleID>_<16bp barcode>`, e.g. `SF10022_CTATCTAAGCAAGCCA` | CONFIRMED |
| `GSE174554_RAW.tar` members | 329 | CONFIRMED — archive index |

**Not yet read:** `GSE174554_family.soft.gz` and the three `*_series_matrix.txt.gz`
files, which carry per-GSM sample characteristics. **These are the authoritative
source for timepoint and patient mapping.** Everything marked PROVISIONAL below
becomes CONFIRMED or gets corrected once they are parsed.

> **Warning — the deposit is not all human GBM.** GEO lists this series under three
> platforms including **GPL24247 (*Mus musculus*)**, and the archive contains mouse
> irradiation samples (`Mouse_IR_Treated`, `T5224_3days_treated`), snATAC, spatial
> transcriptomics and proteomics. The discovery cohort is a **subset**. Any code that
> globs the archive and assumes human snRNA-seq will silently ingest mouse data.

---

## 1. Identifier conventions

| Field | Dtype | Status | Definition |
|---|---|---|---|
| `sample_id` | `str` / `category` | **CONFIRMED** | Text before the underscore in `Sample#_Barcode`, e.g. `SF10022`, `SF10099v2`. |
| `nucleus_id` | `str` | **CONFIRMED** | Full `Sample#_Barcode` string. Unique across the object. Use as `adata.obs_names`. |
| `patient_id` | `str` / `category` | **PROVISIONAL** | Intended: `sample_id` with the timepoint suffix stripped. **Unverified.** |
| `timepoint` | `category` | **PROVISIONAL** | Intended: `{"primary", "recurrent"}`. **Encoding unverified — see §1.1.** |
| `is_malignant` | `bool` | **CONFIRMED** | `Tumor_Normal_annotation == "Tumor"`. Authors' own call; we do not reclassify (CLAUDE.md §4). |

### 1.1 PROVISIONAL — the timepoint encoding is not settled

Observed suffix patterns across the 100 sample IDs:

- 77 bare (`SF10022`), 23 carrying a **`v2`** suffix (`SF10099v2`).
- Applying "strip `v2` → patient, `v2` → recurrent" yields **13 patients with both
  timepoints present**, plus **10 `v2` samples whose base ID is absent** from the
  tumour/normal file.
- At least one different convention exists: **`SF9259R` / `SF9259S`**, which does not
  fit the `v2` scheme at all.

**`v2` has not been shown to mean "recurrent."** It could equally denote a second
library, a re-prep, or a re-sequencing of the same timepoint. Do not build this
assumption into anything that is hard to unwind.

### 1.2 PROVISIONAL — patient count, and why it matters

`configs/pipeline_config.yaml` records `expected_n_patients: 19` with
`n_folds: "equals_n_patients"`, and CLAUDE.md §10.1 explicitly forbids assuming 19
is achievable before the metadata is read.

The naive `v2` rule gives **13**, not 19. That is a *provisional* number from a
*provisional* rule against a file that is not the authoritative sample manifest —
it is **not** a finding yet. But it is close enough to the line that Lane 2 should
treat `n_patients` as a **runtime value**, never a literal:

```python
n_patients = adata.obs["patient_id"].nunique()
n_folds    = n_patients
nb_correction = 1 / n_folds + 1 / (n_folds - 1)   # Nadeau-Bengio
```

Whatever the realised count turns out to be, it goes in `DEVIATIONS.md`, not into
the frozen config.

---

## 2. `data/processed/07_metacells.h5ad` — the Stage B unit

Metacells, ~30 nuclei each (CLAUDE.md §3.8). Built **within** `(patient_id, timepoint)`,
never across patients.

### `.obs`

| Column | Dtype | Status | Notes |
|---|---|---|---|
| `metacell_id` | `str` | CONFIRMED (convention) | `obs_names`. Format `{sample_id}_mc{NNN}`. |
| `patient_id` | `category` | PROVISIONAL | **The grouping key for every fold operation.** |
| `timepoint` | `category` | PROVISIONAL | `primary` / `recurrent`. |
| `sample_id` | `category` | CONFIRMED | |
| `n_nuclei` | `int32` | CONFIRMED | Nuclei aggregated into this metacell. |
| `cell_state` | `category` | PROVISIONAL | Neftel: `MES`, `AC`, `OPC`, `NPC`. 4 levels; §3.5 allows 4–6. |
| `genotype_egfr_amp` | `int8` | PROVISIONAL | 0/1. Stage A fixed effect. |
| `clone_id` | `category` | PROVISIONAL | From inferCNV on chr7/9p/10 only. May be degenerate — CLAUDE.md §10.2. |
| `stage_a_residual` | `float64` | CONFIRMED (definition) | **Mean** of per-nucleus residuals (§10.3). Stage B's adjusted target. |

### `.var`

| Column | Dtype | Status | Notes |
|---|---|---|---|
| `gene_symbol` | `str` | CONFIRMED | `var_names`. |
| `chrom` | `category` | CONFIRMED (definition) | `chr1`…`chrX`, cytoband hg38. |
| `arm` | `category` | CONFIRMED (definition) | `p` / `q`. **Required** — chr9**p** is an arm, not a chromosome. |
| `in_disjoint_set_S` | `bool` | CONFIRMED (definition) | True iff chr7, chr9p, or chr10. **Excluded from X in the Tier C-disjoint run.** 9q must be False. |

---

## 3. `results/tables/ras_scores.csv`

One row per **nucleus**. Components z-scored before weighting; weights fixed a
priori and never fit (§3.3).

| Column | Dtype | Notes |
|---|---|---|
| `nucleus_id` | `str` | Join key. |
| `patient_id` | `str` | PROVISIONAL. |
| `T`, `G`, `O`, `Ab_state` | `float64` | Tier A components. `G` is 0/1. |
| `Ab_clone` | `float64` | Tier C-disjoint only. |
| `ras_tier_a` | `float64` | `0.25·z(T)+0.25·z(G)+0.25·z(O)+0.25·z(Ab_state)`. |
| `ras_tier_c_disjoint` | `float64` | `0.5·z(G)+0.5·z(Ab_clone)`, both from chr7/9p/10 signal. |

If optimal transport misses its Day 3 hard stop, `O` is dropped and the column is
named **`ras_tier_a_reduced`** — never silently reweighted under the old name.

## 4. `results/tables/stage_a_residuals.csv`

Per **nucleus**. Aggregated to metacells by mean before Stage B.

| Column | Dtype | Notes |
|---|---|---|
| `nucleus_id` | `str` | |
| `patient_id` | `str` | PROVISIONAL. |
| `ras_tier_a` | `float64` | Model input. |
| `fixed_effect_pred` | `float64` | Cross-fitted: coefficients from the other patients only; held-out random effect = 0. |
| `residual` | `float64` | `ras_tier_a − fixed_effect_pred`. |
| `heldout_fold` | `int16` | Which LOPO fold treated this nucleus as test. |

## 5. Stage B feature matrix `X`

Rows = metacells, aligned to `07_metacells.h5ad.obs_names`.

Blocks: `hvg_expression` (2,000, **re-derived inside each training fold**),
`hallmark_pathway_scores`, `dorothea_regulon_scores`.

**Forbidden — enforced by test, not by care:**

- PCA coordinates, transport mass, abundance statistics. These are RAS's own inputs;
  including them lets the model reconstruct the target (§3.6).
- For the Tier C-disjoint run: every gene with `in_disjoint_set_S == True`.
  `tests/test_chr_disjoint.py` asserts this at **arm** granularity — genes on **9q
  remain eligible**. Read the region list from
  `configs/pipeline_config.yaml: disjoint_set_S.regions`; do not hardcode it.

---

## 6. Contract rules

1. **`patient_id` is the grouping key everywhere** — HVG selection, scaling, tuning,
   stability resampling. No cell from a held-out patient may influence its own
   prediction (CLAUDE.md §7.2).
2. **`n_patients` is read at runtime**, never assumed (§1.2).
3. **Categories are `category` dtype**, not `object`, so level order is stable across folds.
4. **Seeds come from `configs/pipeline_config.yaml`** — master 42, bootstraps at
   `master + resample_index`. No hardcoded seeds.
5. **A renamed column is a breaking change.** Announce it in both lanes' sessions.
