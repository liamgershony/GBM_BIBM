# DEVIATIONS.md

Append-only log of every departure from the protocol frozen in `CLAUDE.md` and
`configs/pipeline_config.yaml`.

**Rules.** Append, never edit or delete a prior entry. Every entry needs a UTC
timestamp, the parameter or procedure affected, the reason, and who authorised it.
Declared deviations are fine; undeclared ones invalidate the paper (CLAUDE.md §7.9).

Format:

```
## YYYY-MM-DD HH:MM UTC — <short title>
- **Affects:** <config key / script / procedure>
- **From → To:** <old> → <new>
- **Reason:** <why>
- **Authorised by:** <name(s)>
- **Paper:** <where this is disclosed, e.g. "Methods §III-C" or "Limitations">
```

---

## Deviations already declared in the protocol (Revision 2, pre-registered)

These were fixed before any data was processed and are recorded here so the paper's
deviation list has a single source. They are NOT post-hoc changes.

- Stability bootstraps reduced 300 → 100, for compute.
- ElasticNet alpha tuned once per outer fold and reused across that fold's
  bootstraps with warm starts, rather than a full inner CV inside every bootstrap.
  Slightly understates selection variance.
- XGBoost dropped; linear-only discovery.
- H1 ablation at 200 resamples rather than 1,000, pre-specified as reduced for compute.
- Tier C-strict (exome-derived clone identity) replaced by Tier C-disjoint, because
  patient-level exome is not in the public GSE174554 deposit.
- Tier B dropped for scope.
- GLASS replaced by CGGA mRNAseq_693 / _325, because GLASS requires a Synapse
  account and data-use agreement.
- Max mitochondrial fraction 20% → 5%, correcting for single-nucleus chemistry.
- H2, H4, H5 removed from scope.

## Contingencies pre-declared (trigger, then log here if they fire)

- If SEACells does not converge by end of Day 3, substitute k-means metacells
  within each patient-timepoint.
- If optimal transport does not converge by end of Day 3, drop component **O** and
  use **Tier A-reduced** (T, G, Ab_state), named as such. Never silently drop a component.
- If CGGA tumour purity is unavailable, drop it from the covariate set and declare it.
  Do not substitute a proxy without saying so.

---

## Log

<!-- Append new entries below this line. -->

## 2026-08-24 00:25 UTC — CORRECTION to frozen config: Harmony batch key was wrong

- **Affects:** `configs/pipeline_config.yaml` -> `integration.batch_key`.
- **From -> To:** `"sample_id"` -> `"patient_id"`.
- **Authorised by:** Liam Gershony (Lane 1), Arnav Vishwakarma (Lane 2).
- **Paper:** Methods (integration).

**This is a correction of an error, not a change of plan.** CLAUDE.md never
specifies Harmony's batch key. The line `batch_key: "sample_id"` was introduced by
Claude when the config was first written on 23 Aug 2026 and had no basis in the
protocol.

**Why it was wrong, and why it mattered.** `sample_id` distinguishes a patient's
primary specimen from their recurrent specimen. Integrating on it would have
instructed Harmony to remove the primary-versus-recurrent difference as though it
were technical batch. That difference **is RAS component T** -- the transcriptional
similarity between a primary cell and the patient's recurrent centroid -- and it is
the signal the whole study is built on. The pipeline would have run cleanly and
produced a null that looked methodologically sound.

Caught during Day 1 planning review, **before any preprocessing was executed.** No
result was produced under the incorrect key.

**Correct behaviour.** Harmony integrates on `patient_id`. `batch_key`
(`{sample_id}__{library}`) and `library` are retained in `.obs` as covariates
for downstream use and are **never** passed to Harmony.

**Gate added (protocol Step 3).** After integration, LISI is computed on
`patient_id` and on `timepoint`. Because LISI is bounded by category count -- 1..n
for patients but 1..2 for timepoint -- raw values are not comparable; both are
normalised to `(LISI - 1)/(k - 1)`. Integration **fails** if normalised timepoint
LISI >= 0.90 x normalised patient LISI, i.e. if the timepoints have been mixed.
Reported in `results/tables/lisi_gate.csv`.

## 2026-08-24 00:25 UTC — Two amendments to the agreed cohort rule

- **Affects:** `docs/COHORT_RULE.md` clause (b) and the pooling clause.
- **Authorised by:** Liam Gershony (Lane 1), Arnav Vishwakarma (Lane 2).
- **Status:** agreed before any preprocessing was executed.

**1. Re-sample suffix handling.** `SF6118v2` and `SF9715v2` are human snRNA-seq
GSMs in GSE174554 carrying real nuclei from cohort patients at the Recurrent
timepoint, but Supplementary Table 1 does not index them. A literal reading of
clause (b) discarded them on a bookkeeping technicality.

> **Sub-rule.** A GSM is assigned to the patient and timepoint of its specimen ID
> **with a trailing re-sample suffix (`vN`) stripped**, provided the stripped ID
> appears in Supplementary Table 1 for that timepoint.

Recovers both. The cohort is **29 patients / 61 specimens**. `SF12704v2` is *not*
recovered -- its stripped ID `SF12704` carries a literal `NA` pair and is not in a
matched pair.

**2. Batch key composition.** 7 cohort specimens (`SF10099`, `SF10433`, `SF12243`,
`SF4209`, `SF4449`, `SF4810`, `SF7307`) carry a **second 10x library** (`batch2`)
under the same GSM -- a technical batch *inside* one specimen. The agreed wording
"sample ID retained as batch key" cannot distinguish them, so a real technical
batch would have been invisible.

> **Amendment.** The batch key is **`{sample_id}__{library}`** where `library`
> is `batch1` or `batch2`. `sample_id` remains its own `.obs` column.

The cohort spans **68 libraries** across 61 specimens. Per the correction above,
this key is a covariate only and is never Harmony's integration key.

## 2026-08-24 00:25 UTC — Non-frozen runtime thresholds file introduced

- **Affects:** new `configs/runtime_thresholds.yaml`.
- **Reason:** implementation surfaced two operational thresholds that are not part
  of the pre-registered analysis and must not sit in the frozen config.
- **Authorised by:** Liam Gershony (Lane 1).

- `qc.annotation_join_min_rate: 0.80` — a **bug detector**, not a quality gate.
  The authors annotated only the nuclei they retained, so divergence from our QC is
  expected and lands in the high 90s; a barcode key-format mismatch yields near 0%.
  0.80 separates those two regimes. 0.95 was considered and rejected because it
  would fire on ordinary divergence.
- `integration.lisi_timepoint_fail_ratio: 0.90` — the Step 3 gate described above.

**Constraint on this file:** every value in it must be a bug detector or an
engineering guard. Anything capable of influencing a reported statistic belongs in
the frozen config with its own DEVIATIONS entry.


## 2026-08-23 23:58 UTC — Discovery cohort rule agreed (CLAUDE.md §10.1 requirement satisfied)

- **Affects:** definition of the discovery cohort; `n_patients` and every quantity
  derived from it.
- **Decision:** `docs/COHORT_RULE.md` is **agreed by Liam Gershony (Lane 1) and
  Arnav Vishwakarma (Lane 2)**, 23 Aug 2026. CLAUDE.md §10.1 requires the subset to
  be defined by a written, reproducible rule agreed by both lanes *before anything
  downstream is computed*. That condition is now met. No preprocessing had been run
  at the time of this entry.
- **Authorised by:** Liam Gershony and Arnav Vishwakarma.
- **Paper:** Methods (cohort definition), and Results (cohort flow, Figure 1).

**The rule.** A patient enters the discovery cohort if all four hold:

- **(a)** Wang et al. **Supplementary Table 1 `Pair#`** links the patient to at
  least one `Primary` and at least one `Recurrent` specimen. GEO's `pair#`
  characteristic is recorded as a cross-check, **not** as the source of pairing.
- **(b)** Both specimens are present in **GSE174554 as human snRNA-seq** GSMs.
- **(c)** **IDH-wildtype** per Supplementary Table 1.
- **(d)** **≥100 usable nuclei at both timepoints after QC.**

Multiple specimens at one timepoint are **pooled**, with `sample_id` retained as a
batch key.

**Pre-QC ceiling: n = 29.** Clauses (a)–(c) are evaluable from metadata alone and
yield 29 matched pairs. Clause (d) requires Day 1 QC and **can only reduce** this
number. The realised n is therefore not known at the time of this entry, and every
formula that depends on it — fold count, the Nadeau-Bengio correction
`1/n + 1/(n-1)`, and the §6.2 evaluability floor `floor(n/2) + 1` — is evaluated at
runtime against the realised value. **The frozen config is not edited.**

**Explicit exclusion — one IDH-mutant pair.** Supplementary Table 1 records
**`SF8963` (Primary) and `SF12165` (Recurrent)**, GEO `pair#` #33, as **IDH mutant**
despite a diagnosis of Glioblastoma. This pair satisfies clauses (a) and (b) and is
excluded solely by clause (c). It is the only difference between the 30 matched
pairs recoverable from the deposit and the 29 entering the cohort. Recorded here so
the 30 → 29 step is traceable to a stated criterion and not to an unexplained count.

**Note on the deposit.** GEO's series summary claims 40 matched IDH-wildtype pairs.
That figure is not reproducible from either pairing record. Supplementary Table 1
supports 36 matched pairs flagged `snRNA-seq = Y`, but only 30 have both specimens
present in GEO as human snRNA-seq GSMs — 86 specimens are flagged snRNA-seq in the
paper while 78 appear as such in the deposit. This is a deposit-versus-publication
discrepancy, not an analysis choice, and belongs in the paper's cohort description.


## 2026-08-23 23:25 UTC — Pre-specified cascade if the discovery cohort has n != 19 patients

- **Affects:** `cross_validation.n_folds`, `variance_correction`, the analysis
  evaluability floor, and hypothesis H1's admissibility.
- **Status:** Written **before the patient count is known.** At the time of this
  entry the authoritative sample manifest (`GSE174554_family.soft.gz`) has not been
  parsed. A provisional reading of sample-ID suffixes suggests the count may be
  below 19, which is why this is being fixed now rather than after the number lands.
- **Authorised by:** Liam Gershony (Lane 1), 23 Aug 2026.

**1. Fold count and variance correction — already formula-encoded, no change needed.**
`configs/pipeline_config.yaml` records `expected_n_patients: 19` alongside
`n_folds: "equals_n_patients"` and
`variance_correction.formula: "1/n_folds + 1/(n_folds - 1)"`. Both therefore track
the realised count automatically. The literal `1/19 + 1/18` in CLAUDE.md 5 is the
pre-specified *expectation*, not the operative definition. **The frozen config is
not edited when n differs**; the realised n is reported in the paper alongside the
expectation.

**2. Evaluability floor.** The floor is generalised from **"10 of 19"** to
**`floor(n/2) + 1`**, evaluated against the realised n. This is a restatement of the
same rule, not a change of stringency: at n = 19 it gives `floor(19/2) + 1 = 10`,
reproducing the pre-specified expectation exactly. At n = 30 it gives 16.

> **Source.** This rule comes from **Step 12 of the original protocol**, not from
> `CLAUDE.md` — the committed CLAUDE.md carried no evaluability floor at the time of
> this entry (its §12 is "How to work with us"). The rule has since been added to
> CLAUDE.md as **§6.2** so that the operative document states it. Cite the protocol
> Step 12 as the origin and CLAUDE.md §6.2 as the current statement.

> **Correction, same session:** this entry originally recorded the formula as
> `ceil(n/2) + 1`, which yields 11 at n = 19 and therefore contradicts the
> pre-specified "10 of 19". Corrected to `floor(n/2) + 1` before any evaluability
> result existed.

**3. H1 is formally dropped if n < 16.** H1 (CLAUDE.md 6) requires a paired
patient-level resample interval to support the decision rule
"mean dR >= 10 pp AND interval excludes 0". Below 16 patients the interval is too
wide for that rule to discriminate, so H1 would be reported as an uninformative
null caused by cohort size rather than by biology. **H1 is therefore abandoned
outright at n < 16 and not reported at all** — including not reported as a null.
H3 remains the core deliverable and is unaffected: it is a permutation test on gene
list overlap and does not depend on patient count in the same way.

This threshold is set **before** the count is known so that it cannot be chosen to
include or exclude H1 on the basis of a result.

## 2026-08-23 23:25 UTC — CGGA expression matrix: Read_Counts primary, RSEM sensitivity only

- **Affects:** CLAUDE.md 4.3 replication test input.
- **Reason:** CGGA publishes two expression matrices per cohort — `Read_Counts`
  (2022-06-20) and `RSEM` (2020-05-06). CLAUDE.md 4 refers to read counts but does
  not explicitly exclude RSEM. Both are archived in `data/raw/CGGA/` with checksums.
- **Decision:** The **`Read_Counts` tables are the protocol matrix** for the
  replication test. The **RSEM tables are a pre-specified sensitivity analysis**.
- **Hard constraint:** RSEM **must never be substituted for Read_Counts after
  Read_Counts results have been seen.** Reporting the RSEM result in place of the
  Read_Counts result, or selecting between them on the basis of which replicates
  better, is outcome-dependent analysis selection and would invalidate the
  replication claim. If both are reported, Read_Counts is identified as primary and
  RSEM as sensitivity, in that order, regardless of which is more favourable.
- **Authorised by:** Liam Gershony (Lane 1), 23 Aug 2026, **before any replication
  result existed** — no CGGA file has been decompressed or parsed at this entry.
- **Paper:** Methods (replication), and Limitations if the two disagree.

