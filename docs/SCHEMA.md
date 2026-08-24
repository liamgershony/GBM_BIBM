# SCHEMA.md — Lane 1 → Lane 2 handoff contract

**Status: UPDATED 24 Aug 2026 (Day 1 complete). Supersedes the Day 0 draft.**

> **Lane 2 — read this first.** Four things changed materially since the Day 0 draft:
> 1. **`n = 21`, provisional** (not 19, not 29). Read it from
>    `results/tables/cohort_n.json`; never hardcode. Provisional because doublet
>    detection did not run — see §0.2.
> 2. **`patient_id` and `timepoint` are now CONFIRMED**, from the SOFT file's
>    `pair#` and `progression` characteristics. The `v2`-suffix guess was wrong.
> 3. **New `.obs` columns**: `batch_key`, `library`, `is_malignant`,
>    `is_reference_normal`, `tumor_normal_annotation`, `annotation_key`, `gsm`.
> 4. **`batch2` libraries contribute ZERO usable malignant nuclei** — they are
>    unannotated by construction. If your simulator gives them malignant cells,
>    it does not match the real object.

**Original status: DRAFT, 23 Aug 2026 (Day 0).** Arnav is coding Stage A/B against simulated
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

## 0. Realised cohort and pipeline facts (Day 1)

### 0.1 Confirmed counts

| Fact | Value | Source |
|---|---|---|
| GSMs in GSE174554 | 113 | `sample_manifest.csv` |
| Human snRNA-seq GSMs | 81 | `sample_manifest.csv` |
| Matched pairs (clauses a+b) | 30 | `pairing_comparison.csv` |
| IDH-wildtype pairs (clause c) | **29** | `discovery_cohort.csv` |
| Specimens / libraries in cohort | **61 / 68** | `discovery_cohort.csv` |
| Nuclei after QC | 142,726 | `01_qc.h5ad` |
| — malignant | 80,201 | |
| — normal (inferCNV reference) | 54,986 | |
| — unknown | 7,539 | |
| **Patients passing clause (d)** | **21** | `cohort_n.json` |
| Nuclei in the 21-patient cohort | 121,268 | |

At `n = 21`: `n_folds = 21`, Nadeau-Bengio `1/21 + 1/20 = 0.097619`, evaluability
floor `floor(21/2)+1 = 11` patients per state, **H1 remains admissible**.

### 0.2 Why `n = 21` is PROVISIONAL

**Doublet detection never ran.** `sc.pp.scrublet` raised `ModuleNotFoundError` in
all 60 libraries that reached it — scanpy ships the wrapper, not the
implementation, and the `scrublet` package was not installed. The exception was
caught and every nucleus retained. All QC counts, and therefore clause (d) and
`n = 21`, were computed **without doublet removal**. Expect `n` to move once this
is fixed. **Read `n` from `cohort_n.json` at runtime — never hardcode 21.**

### 0.3 File-format facts that cost a day, so you don't repeat them

- **`features.tsv.gz` is a SINGLE column of gene symbols**, not the standard 3-column
  10x file, and the matrix is stored **genes × cells**. `sc.read_10x_mtx` fails with
  `KeyError: 1`. See `read_library()` in `src/01_qc_integration.py`.
- **The authors' annotation renames 9 cohort specimens**: GEO `SF11916` is annotation
  `SF11916v2`. Verified 1:1 and unambiguous per specimen.
- **Barcode lane suffixes are inconsistent** — libraries use `-1`; the annotation
  uses `-1` for some specimens and `.1` for others. Barcode counts match exactly
  while zero keys join. Both sides are normalised by stripping `[-.]\d+$`.
- **`batch2` is not annotated at all.** Annotation row counts match `batch1` barcode
  counts, and batch1/batch2 share barcode *sequences* (`SF7307`: all 51 batch2
  barcodes occur in batch1), so `{sample_id}_{barcode}` cannot distinguish them.
  batch2 is forced to `unknown` — a match would assign one nucleus's label to a
  different nucleus.
- **The mitochondrial filter is a safety net, not an active filter.** Observed
  `pct_counts_mt` runs median 0.0–0.14%, max 4.8% — nuclear preps are mito-depleted,
  exactly as CLAUDE.md §7.5 says. 0 nuclei removed at the 5% ceiling is correct.

> **Warning — the deposit is not all human GBM.** GEO lists three platforms including
> **GPL24247 (*Mus musculus*)**; the archive holds mouse irradiation samples, snATAC,
> spatial and proteomics. Code that globs the archive will ingest mouse data.

## 1. Identifier conventions

| Field | Dtype | Status | Definition |
|---|---|---|---|
| `sample_id` | `str` / `category` | **CONFIRMED** | Text before the underscore in `Sample#_Barcode`, e.g. `SF10022`, `SF10099v2`. |
| `nucleus_id` | `str` | **CONFIRMED** | Full `Sample#_Barcode` string. Unique across the object. Use as `adata.obs_names`. |
| `patient_id` | `category` | **CONFIRMED** | Wang et al. Supplementary Table 1 `Pair#`. **Not** GEO's `pair#` — the two use different numbering (GEO #1 = supp 28). |
| `timepoint` | `category` | **CONFIRMED** | `Primary` / `Recurrent`, from the SOFT `progression` characteristic. Capitalised exactly so. |
| `gsm` | `category` | **CONFIRMED** | GEO accession of the source sample. |
| `library` | `category` | **CONFIRMED** | `batch1` / `batch2`. 7 specimens carry a second library. |
| `batch_key` | `category` | **CONFIRMED** | `{sample_id}__{library}`. Covariate only — **never** Harmony's key. |
| `is_malignant` | `bool` | **CONFIRMED** | `Tumor_Normal_annotation == "Tumor"`. Authors' own call; we do not reclassify (CLAUDE.md §4). |

### 1.1 RESOLVED — timepoint encoding

The `v2`-suffix guess in the Day 0 draft was **wrong**. The SOFT file carries
`progression: Primary|Recurrent` per GSM directly. `v2` marks a **second sample at
the same timepoint**, not recurrence: all three `v2` samples are Recurrent and sit
alongside a non-`v2` recurrent sample from the same patient.

A GSM maps to the patient/timepoint of its specimen ID **with a trailing `vN`
stripped**, if the stripped id is in Supplementary Table 1 (agreed amendment,
DEVIATIONS.md). This recovers `SF6118v2` and `SF9715v2`.

### 1.2 `n` is read at runtime, always

```python
import json
info = json.load(open("results/tables/cohort_n.json"))
n_patients    = info["n_patients"]          # 21 today; will move
n_folds       = info["n_folds"]
nb_correction = info["nadeau_bengio_correction"]   # 1/n + 1/(n-1)
floor_        = info["evaluability_floor_patients"]  # floor(n/2)+1
```

Never hardcode 19, 21 or 29. A realised `n` other than the pre-specified 19 is a
DEVIATIONS entry, never an edit to the frozen config.

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
| `library` | `category` | **CONFIRMED** | `batch1`/`batch2`. Covariate only. |
| `batch_key` | `category` | **CONFIRMED** | `{sample_id}__{library}`. Covariate only. |
| `gsm` | `category` | **CONFIRMED** | Source GEO accession. |

### 2.1 Nucleus-level `.obs` in `01_qc.h5ad` (upstream of metacells)

`01_qc.h5ad` retains **all** nuclei, malignant and non-malignant — inferCNV needs
the non-malignant ones as reference on Day 2. Only clause (d) counting is
malignant-restricted.

| Column | Dtype | Notes |
|---|---|---|
| `tumor_normal_annotation` | `category` | `Tumor` / `Normal` / `unknown`. Authors' own call (CLAUDE.md §4). |
| `is_malignant` | `bool` | `annotation == "Tumor"`. **This is the Stage A/B analysis set.** |
| `is_reference_normal` | `bool` | `annotation == "Normal"`. **inferCNV reference set.** |
| `annotation_key` | `str` | `{annotation_prefix}_{normalised_barcode}`; empty where not applicable. |

`unknown` nuclei belong to **neither** set — neither status can be asserted for
them. All `batch2` nuclei are `unknown`, as is every nucleus of `SF11981`, which is
absent from the annotation file under any name.

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
1b. **Harmony integrates on `patient_id`, never on `sample_id` or `batch_key`.**
   `sample_id` separates a patient's primary specimen from their recurrent one, so
   integrating on it would regress out the primary-vs-recurrent difference — which
   **is RAS component T**. The frozen config said `sample_id`; that was an error,
   corrected in DEVIATIONS.md. `02_integration.py` asserts the key at runtime.
2. **`n_patients` is read at runtime**, never assumed (§1.2).
3. **Categories are `category` dtype**, not `object`, so level order is stable across folds.
4. **Seeds come from `configs/pipeline_config.yaml`** — master 42, bootstraps at
   `master + resample_index`. No hardcoded seeds.
5. **A renamed column is a breaking change.** Announce it in both lanes' sessions.
