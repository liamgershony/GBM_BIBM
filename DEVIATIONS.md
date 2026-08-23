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

**2. Evaluability floor.** The floor changes from **"10 of 19"** to
**`ceil(n/2) + 1`**, evaluated against the realised n. "10 of 19" is recorded here
as the pre-specified expectation. For n = 19 the new rule gives
`ceil(19/2) + 1 = 11`, which is *stricter* than the original 10 — the rule is not a
relaxation designed to keep a marginal cohort evaluable, and that is the point of
fixing it before the count is known.

> **Citation note.** The "10 of 19" floor does **not appear in the committed
> `CLAUDE.md`** as of this entry — §12 is "How to work with us" and no evaluability
> floor is stated anywhere in the file. It is recorded here as a decision taken on
> 23 Aug 2026 rather than as a quotation from a section that a reader cannot find.
> If the rule exists in the original protocol draft, CLAUDE.md should be amended to
> carry it and this note updated to cite the correct section.

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

