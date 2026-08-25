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

## 2026-08-25 17:54 UTC — Stage A is predominantly patient-level centring, and the two tiers share G

- **Affects:** interpretation of H1 and, more seriously, of **H3**. No parameter changes.
- **Authorised by:** Liam Gershony (Lane 1).
- **Belongs in the paper's interpretation, not only in this log.**

### (a) What Stage A actually removes

| Tier | R², state + genotype | R², + patient | Variance remaining |
|---|---|---|---|
| A-reduced | **0.1112** | 0.5252 | 47.48% |
| C-disjoint | **0.1772** | 0.5572 | 44.28% |

State and genotype together explain only 11–18% of RAS. Nearly all of Stage A's
explanatory power comes from the patient term. **The confound adjustment is
predominantly patient-level centring rather than removal of state and genotype.**
H1 asks whether confound adjustment improves external replication; if H1 shows a
benefit, the honest reading is that most of it comes from patient centring, not
from removing the cell-state and genotype confounds the protocol names.

### (b) CORRECTION: Stage A does NOT remove between-patient variance

An earlier entry in this log stated that Stage A's patient random intercept removes
G's purely between-patient variance, and therefore that Tier A-reduced behaves as
`0.5·z(T) + 0.5·z(Ab_state)` after Stage A. **That was wrong.**

§3.5 specifies that the held-out patient's random effect is predicted as **0**.
Only the state and genotype *fixed* effects are subtracted, so the patient offset
`u_p` remains in the residual. Measured: Tier A-reduced is 46.4% between-patient
variance raw and **47.7% after Stage A** — the patient structure is not reduced.

Consequently **G does contribute to the Stage A residual that Stage B fits**, and
the earlier "effectively two-component" description applies to the raw score's
within-patient behaviour, not to the residual. The corrected statement: G is
constant within 19/21 patients, so it contributes no *within-patient* variance at
any stage, but its *between-patient* variance survives Stage A intact.

### (c) The two tiers are NOT independent by construction

`Tier A-reduced = (1/3)z(T) + (1/3)z(G) + (1/3)z(Ab_state)`
`Tier C-disjoint = 0.5·z(G) + 0.5·z(Ab_clone)`

**Both contain z(G).** Measured contribution:

| quantity | value |
|---|---|
| corr(Tier A, Tier C) | **0.5634** |
| corr with the shared G removed from both | **0.1366** |
| z(G) share of Tier A variance | 23.0% |
| z(G) share of Tier C variance | 38.2% |
| corr(Stage A residuals) | **0.6521** |
| corr(Stage A residuals), patient-centred | 0.4119 |

**Roughly three quarters of the correlation between the two tiers is the shared G
term.** H3 tests whether genes selected under Tier A overlap those selected under
Tier C-disjoint above a permutation baseline. Because the two targets share a
component, some overlap is expected **by construction**, and a permutation null
that shuffles labels does not account for correlation between the targets
themselves.

This is a property of the protocol as frozen in §3.3, not a change introduced here,
and it is recorded **before H3 is run** so the result cannot be reinterpreted
afterwards. **H3's Jaccard index must be reported alongside this shared-component
correlation**, and the claim narrowed accordingly: a significant overlap does not
by itself demonstrate a shared biological programme, because part of the overlap is
attributable to the shared G term.

One mitigation is real and worth stating: the *cis* route is still closed for Tier
C. Stage B's feature matrix for Tier C-disjoint excludes every gene on chr7, chr9p
and chr10, so genes whose expression tracks the copy number that defines G cannot
be selected in that arm. What remains shared is the *trans* route and the
patient-level offset.


## 2026-08-25 17:36 UTC — Metacells: SEACells in 33/42 groups, k-means contingency in 9

- **Affects:** `data/processed/07_metacells.h5ad`, `results/tables/metacell_catalog.csv`.
- **Authorised by:** Liam Gershony (Lane 1).

SEACells converged for **33 of 42** patient-timepoint groups. The §9.1 k-means
contingency fired for **9**: 7 exceeded the 180 s per-group time box and 2 raised a
numerical `RuntimeWarning`. The method used is recorded per group in the `method`
column of `metacell_catalog.csv`, so no group's provenance is ambiguous.
2,426 metacells over 72,857 malignant nuclei; median 30.0 nuclei per metacell
against a target of 30.

**An earlier run must not be cited as evidence about SEACells.** It fell back in
all 42 groups on `ImportError: IProgress not found` — a missing `ipywidgets`
progress-bar dependency, not non-convergence. `ipywidgets==8.1.5` is now pinned and
imported explicitly. The 33/42 figure above is from the first run in which SEACells
was actually exercised.


## 2026-08-25 17:36 UTC — CONTINGENCY INVOKED: component O dropped, Tier A-reduced adopted

- **Affects:** CLAUDE.md §3.2 component **O**; §3.3 Tier A weights; the name of the
  Tier A score throughout the paper.
- **Class:** **specification error**, not a data-driven choice. Pre-declared
  contingency in CLAUDE.md §9.1 invoked.
- **Authorised by:** Liam Gershony (Lane 1).

**O is degenerate by construction, and no data could have made it otherwise.**
§3.2 defines O as "total optimal-transport mass moved from cell *i*'s metacell to
any recurrent-timepoint metacell". For **balanced** optimal transport the transport
plan's marginals are fixed: the row sums of the plan equal the source weights
exactly. O is therefore identically the metacell's own mass — its size — and
carries no information about transport at all.

Measured (`results/tables/ot_component_O.csv`, `src/03_metacells_ot.py`):

| quantity | value |
|---|---|
| corr(O, source marginal) | **1.000000** |
| max abs difference | **2.25e-16** (floating-point noise) |

This is an algebraic identity of `ot.emd`, not a property of GBM, this cohort, or
our preprocessing. It would hold for any dataset. **The error is in the protocol's
definition of O, not in the data.**

**Action, per §9.1:** component **O is dropped** and the score is
**Tier A-reduced = (T, G, Ab_state)**, named as such in the paper. Never silently
reweighted under the Tier A name. The RAS column is `ras_tier_a_reduced`.

**A transport-cost alternative was computed and NOT substituted.** The cost
actually incurred by each metacell, `sum_j T[i,j] * M[i,j]`, is non-degenerate
(mean 0.1991, sd 0.5188) and is recorded in `ot_component_O.csv` as
`transport_cost_diagnostic`. It is **not** used as O and does not enter any score.
Replacing a pre-registered component with a different statistic after discovering
the original was degenerate would be a post-hoc redefinition, which is exactly what
freezing the components a priori exists to prevent. **O is dropped, not redesigned.**

The optimal transport machinery itself ran correctly and is retained for the
record: SEACells metacells per patient-timepoint and `ot.emd` on a PCA-Euclidean
cost matrix, over all 21 patients. Nothing failed to converge.

## 2026-08-25 17:36 UTC — Tier A-reduced is in practice a two-component score

- **Affects:** how Tier A-reduced is described in the paper. **No weight is changed.**
- **Authorised by:** Liam Gershony (Lane 1).

With O dropped, Tier A-reduced is
`(1/3)·z(T) + (1/3)·z(G) + (1/3)·z(Ab_state)` at the frozen equal weighting.

**G is constant within 19 of 21 patients** (`results/tables/genotype_class_degeneracy.csv`;
see `results/tables/negative_control_failures.md` §3). Its variance is almost
entirely *between* patients: pooled variance 0.017690, median within-patient
variance **0.000000**. Stage A regresses RAS on state, genotype and a **patient
random intercept**, which removes precisely that between-patient variance.

> **After Stage A, Tier A-reduced behaves as `0.5·z(T) + 0.5·z(Ab_state)` in
> effect**, because G contributes essentially nothing to the residual that Stage B
> predicts.

**The weights are NOT changed.** They stay at the frozen equal thirds. Reweighting
after observing which component degenerates would be a data-dependent choice, and
the equal-weighting decision in CLAUDE.md §3.3 exists specifically to remove that
route. We report the effective structure and leave the specification alone.

The paper must state the weighting *and* the effective structure. Reporting three
equally weighted components without noting that one is within-patient constant
would misdescribe the score.


## 2026-08-24 07:57 UTC — MANUAL INPUT: Neftel Table S2 cannot be fetched automatically

- **Affects:** `data/raw/neftel_signatures/NIHMS1532254-supplement-9.xlsx`; the
  cell-state calls in `src/04_states.py`; Stage A's state fixed effect; Ab(state).
- **Class:** documented departure from CLAUDE.md §4's requirement that every input
  be direct-download. Access is open; **automated retrieval is blocked.**
- **Authorised by:** Liam Gershony (Lane 1).

**This is the ONE pipeline input that cannot be regenerated by running `src/`.**
Everything else in `data/raw/` is reproducible from `src/00*.py`.

| | |
|---|---|
| File | `NIHMS1532254-supplement-9.xlsx` |
| URL | `https://pmc.ncbi.nlm.nih.gov/articles/instance/6703186/bin/NIHMS1532254-supplement-9.xlsx` |
| Bytes | 44,805 |
| SHA256 | `208e73ab3d22c494caf85c867d69dc6be38df3fc62ab1f043d7fcc5441066277` |
| Retrieved | 2026-08-24, manually, in a browser |
| Source | Neftel et al. 2019 *Cell*; PMID 31327527; PMC6703186; Table S2 |

**Why it is manual.** That path returns an HTML proof-of-work interstitial
(`cloudpmc-viewer-pow`, difficulty 4) to non-browser clients rather than the file —
a 1,817-byte page titled "Preparing to download ...". Every alternative was tried
and refused: the PMC OA service and Europe PMC both report the article as not open
access; `www.ncbi.nlm.nih.gov/pmc/articles/`, `pmc.ncbi.nlm.nih.gov/articles/`,
`europepmc.org/articles/` and the FTP tree all return 404; and the NIHMS manuscript
PDF contains only ~5 marker genes per module in prose, stating "see Table S2 for a
full list".

**No proof-of-work solver was implemented, by deliberate choice.** The challenge is
an anti-automation control. The content is open and we are entitled to it, but
building a solver into the pipeline to defeat a bot-detection measure is not an
appropriate way to satisfy a reproducibility requirement. A browser download is one
step, and the SHA256 above pins exactly what was used.

**Verification instead of regeneration.** `src/00h_build_neftel_signatures.py`
asserts this SHA256 before reading the file, so a substituted or corrupted table
fails loudly. Anyone reproducing this work must fetch the file from the URL above
and check the hash matches.

**Related hardening.** The first automated attempt saved the 1,817-byte HTML
interstitial under the `.xlsx` name and recorded its checksum in `PROVENANCE.md` —
a failure written down as a success. `_download_utils.verify_file_type()` now
checks magic bytes against the extension and quarantines mismatches as
`<name>.rejected`.


## 2026-08-24 04:13 UTC — SPECIFICATION GAP RESOLVED: z() scope is POOLED ACROSS THE COHORT

- **Affects:** every `z()` in CLAUDE.md §3.3 — Tier A `0.25·z(T)+0.25·z(G)+0.25·z(O)+0.25·z(Ab_state)` and Tier C-disjoint `0.5·z(G)+0.5·z(Ab_clone)`.
- **Class:** specification gap in the frozen config, resolved. Not a parameter change.
- **Authorised by:** Liam Gershony (Lane 1).
- **Decision:** `z()` is computed **pooled across the whole discovery cohort**, once
  per component, not within patient.

**The gap.** `configs/pipeline_config.yaml` records `ras.standardize: "zscore"` and
CLAUDE.md §3.3 writes `z(...)` without stating the scope. Pooled and within-patient
give materially different answers, so the pipeline cannot proceed without fixing it.

**Numerical grounds, measured before any RAS value was examined**
(`results/tables/zscope_comparison.csv`, from `src/03a_ras_component_diagnostics.py`,
over 14,710 primary malignant nuclei in 21 patients):

| component | pooled SD | patients with zero within-patient variance | cells that would be NaN under within-patient z |
|---|---|---|---|
| T | 0.17793 | 0/21 | 0 (0.00%) |
| G | 0.13301 | **19/21** | **13,973 (94.99%)** |
| Ab_clone | 1.96764 | 0/21 | 0 (0.00%) |

**Within-patient z-scoring is a divide-by-zero for G in 19 of 21 patients**, which
would render 95% of all cells' G undefined and make Tier C-disjoint uncomputable
for almost the whole cohort. Pooled z is well defined for every component in every
patient. The decision is therefore forced numerically, not chosen on preference.

**This was decided on the definedness of the statistic, not on any RAS value or any
downstream result.** No RAS score existed when it was made, and no H3 or H1 output
existed. `O` and `Ab(state)` are not yet built and inherit the same pooled rule.

**Consequence that must be stated in the paper.** Pooled z-scoring means each
component's variance includes a between-patient part. Stage A's patient random
intercept (§3.5) then absorbs that part. For a component whose variance is almost
entirely between-patient — which G now is — the contribution to the Stage A
**residual** is close to zero even though its contribution to raw RAS is large.
See the Tier A variance entry in `results/tables/tier_a_variance.csv`.


## 2026-08-24 02:39 UTC — Operationalising the LISI permutation gate (recorded BEFORE computing it)

- **Affects:** the Step 3 gate in `src/02_integration.py`.
- **Status:** written while **no permutation-null LISI value exists**. The QC re-run
  has completed (n = 21) but `02_integration.py` has not been re-run.

The decision rule recorded earlier is qualitative ("below its null", "near its
null"). Those phrases need numeric form, and choosing that form after seeing the
values would defeat the pre-specification. The cutoffs are therefore fixed here.

**Statistic.** For each label L, `ratio(L) = median observed LISI(L) / median null
LISI(L)`, where the null is LISI recomputed on the same embedding and the same
neighbourhood with L randomly permuted. **3 permutations**, seeded
`master_seed + i`, null taken as their mean. Both labels' observed and permuted
values come from a single `compute_lisi` call so the neighbourhood is identical
throughout, and the perplexity ceiling cancels within each ratio.

**Cutoffs.**
- `timepoint` is **"below its null"** iff `ratio(timepoint) < 0.95`.
- `patient_id` is **"near its null"** iff `ratio(patient_id) >= 0.80`;
  **"far below its null"** iff `ratio(patient_id) < 0.80`.

**Mapping to the pre-specified outcomes.**

| Condition | Outcome | Response |
|---|---|---|
| `ratio(tp) < 0.95` and `ratio(pat) >= 0.80` | **(a)** | PASS. Proceed to RAS. |
| `ratio(tp) < 0.95` and `ratio(pat) < 0.80` | **(b)** | PROCEED, report LISI as a stated limitation. **Do not retune Harmony theta.** |
| `ratio(tp) >= 0.95` | **(c)** | REAL FAILURE. Invoke STOP/GO gate 1 and stop. |

Condition (c) is evaluated first and dominates: if the timepoints are mixed, the
patient-side result is irrelevant.

The embedding is now written to disk **regardless of gate outcome**, so a failure
can be re-diagnosed without re-running Harmony. Writing the embedding is not
proceeding with the analysis; outcome (c) still halts RAS construction.


## 2026-08-24 02:29 UTC — BUG FIX: doublet detection never ran (scrublet was not installed)

- **Affects:** `envs/environment.yml`, `envs/environment.lock.yml`,
  `src/01_qc_integration.py`; all QC counts; clause (d); `n`.
- **Class:** **bug fix / correction of record.** Not a protocol change, not a
  cohort expansion.
- **Authorised by:** Liam Gershony (Lane 1).

**What happened.** `envs/environment.yml` omitted a `scrublet` pin, with the
comment *"scanpy 1.10 ships sc.pp.scrublet -- no separate scrublet pin."* That
comment was **false**: scanpy ships the wrapper, not the implementation. The
`scrublet` package was never installed, so `sc.pp.scrublet` raised
`ModuleNotFoundError` in **all 60 libraries that reached it**. The call sat inside
`except Exception`, which swallowed the error, retained every nucleus, and let the
pipeline report success. **0 doublets were removed from 170,266 input nuclei**, and
`n = 21` was computed without doublet removal.

**Why this is a correction and not a cohort expansion.** Doublet removal can only
**remove** nuclei. Re-running it can therefore only lower per-timepoint counts,
never raise them, so it **cannot cause any patient to newly satisfy clause (d)'s
>=100-nucleus threshold**. The fix is incapable of recovering a failed patient,
which is precisely what makes it safe to apply after the failures were inspected.
Any change in `n` will be downward or nil.

**Fix.**
1. `scrublet==0.2.3` pinned in the pip block; the false comment removed.
2. `envs/environment.lock.yml` regenerated from the repaired environment.
3. `src/01_qc_integration.py`: `ImportError` is **no longer caught** -- a missing
   dependency is a fatal environment error, never a per-library data condition.
   Only `ValueError`/`ArithmeticError`/`LinAlgError` are tolerated, and the outcome
   is written to a new `scrublet_status` column in `qc_per_library.csv` so a
   silent no-op is visible in the artifact. The script now **exits non-zero if
   doublet detection ran in zero libraries.**

**Related hardening, same session.** Every `except Exception` in `src/` was audited
and narrowed, because this bug shares its shape with two earlier ones -- the
`Content-Length: 0` HEAD response that made an absent `.part` look complete, and
the `batch2` barcode collisions that produced 51 false annotation joins. All three
turned a failure into a plausible-looking success. `remote_size` and every download
handler now catch only network errors; `00_download.py` treats an unreachable
`suppl/` index as fatal rather than silently downloading a shorter file list.
No `except Exception` remains in `src/`.

**Consequence.** `n = 21` is superseded by the value in
`results/tables/cohort_n.json` after the re-run.


## 2026-08-24 02:08 UTC — Step 3 LISI gate: redesign, and PRE-SPECIFIED decision rule

- **Affects:** the Step 3 integration gate in `src/02_integration.py`;
  `configs/runtime_thresholds.yaml -> integration.lisi_timepoint_fail_ratio`.
- **Authorised by:** Liam Gershony (Lane 1).
- **Status:** the decision rule below is recorded **before the permutation-null
  LISI values have been computed.** No redesigned-gate output existed when this
  entry was written.

### The original gate, and its failing values

The gate as first implemented computed LISI on `patient_id` and on `timepoint` in
the Harmony embedding, normalised each to `(LISI - 1)/(k - 1)`, and failed if
normalised timepoint LISI >= `0.90 x` normalised patient LISI.

It **fired on the first real run**, with these values (recorded verbatim, from
`results/tables/lisi_gate.csv`, commit b79dbad):

| label | k | median LISI (raw) | normalised `(LISI-1)/(k-1)` |
|---|---|---|---|
| `patient_id` | 21 | **2.747** | **0.0874** |
| `timepoint` | 2 | **1.263** | **0.2626** |

Fail threshold `0.90 x 0.0874 = 0.0786`; observed `0.2626` -> **FAIL**. The run
stopped fail-closed and `02_integrated.h5ad` was not written.

### Why the `(k-1)` normalisation was mathematically invalid

`(LISI - 1)/(k - 1)` presumes LISI can approach its category count `k`. It cannot.
`harmonypy.lisi.compute_lisi` uses **`perplexity = 30`**, so the effective
neighbourhood is roughly 30 nuclei and the attainable LISI is bounded by
**neighbourhood size as well as by `k`**. With `k = 21` patients the patient-LISI
ceiling lies far below 21, whereas with `k = 2` a timepoint LISI approaching 2 is
readily attainable. Dividing by `(k - 1)` therefore deflates the high-`k` variable
and inflates the low-`k` one, **inverting the comparison the gate was built to
make**. The failure was a property of the statistic, not of the embedding.

### Disclosure: the redesign happened AFTER the gate fired

This redesign was undertaken **after** the original gate failed on real data, by
the same author who specified the original gate. That ordering is disclosed
because it is exactly the circumstance in which a threshold change is
untrustworthy. Two constraints follow and are binding:

1. The replacement is justified **only** by the mathematical argument above, which
   is independent of the observed values and would hold had the gate passed.
2. **The decision rule for every outcome is fixed below, before the replacement
   statistic is computed.** Whatever the permutation nulls show, the response is
   already written.

### The replacement gate

For each label, compare observed LISI against **its own permutation null**: LISI
recomputed on the same embedding and neighbourhood structure with that label
randomly shuffled (seeded from `seed.master`). The null is the fully-mixed
reference *for that label under the actual neighbourhood size*, so no cross-`k`
comparison arises and the perplexity ceiling cancels.

- `patient_id`: observed should approach its null (patients mixed -> correction worked).
- `timepoint`: observed should sit well below its null (primary and recurrent still
  separable -> RAS component T preserved).

### PRE-SPECIFIED decision rule

| Outcome | Response |
|---|---|
| **(a)** timepoint **below** its null AND patient **near** its null | **PASS.** Proceed to RAS construction. |
| **(b)** timepoint **below** its null BUT patient **far below** its null | **PROCEED, with a stated limitation.** Integration under-corrected across patients. Report the LISI values in the paper as a limitation. **Do NOT retune Harmony `theta` to force mixing.** |
| **(c)** timepoint **at or above** its null | **REAL FAILURE.** Invoke **STOP/GO gate 1** and stop. Component T is compromised; no RAS work proceeds on this embedding. |

Outcome (b) is explicitly a *proceed* condition. Under-correction across patients
is a weaker embedding, not a circular one: it does not manufacture the
primary-vs-recurrent similarity that RAS component T measures. Tuning `theta`
upward until patients mix would be fitting a preprocessing parameter to make a
diagnostic look better, with no outcome-independent justification — precisely the
data-dependent choice the equal-weighting decision in CLAUDE.md §3.3 exists to
avoid. The honest report of a weaker embedding is preferred.

`lisi_timepoint_fail_ratio` is retired from `configs/runtime_thresholds.yaml`; the
permutation gate has no tunable ratio, which is part of its appeal.


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

