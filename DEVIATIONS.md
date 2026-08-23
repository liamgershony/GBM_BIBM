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
