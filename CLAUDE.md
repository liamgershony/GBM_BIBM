# CLAUDE.md — GBM Persister Project

**Read this file completely before doing anything. It is the source of truth for scope, parameters, and rules.**

Last updated: 23 August 2026 (Revision 2)

---

## 0. Quick orientation

We are three students writing a 5-page computational-biology methods paper for the IEEE BIBM 2026 Undergraduate & High School Symposium (UGHS). The submission deadline is **31 August 2026, 11:59 PM AoE**. We started on 23 August. Nothing had been executed before that date.

The paper asks whether glioblastoma cells that survive treatment are already molecularly distinguishable inside the untreated primary tumor — and, more precisely, whether a confound-adjustment step produces candidate genes that survive a circularity control.

**The core deliverable is H3 (the circularity check). H1 (the ablation) is a stretch goal. H2, H4 and H5 are out of scope.**

If you are ever unsure whether something is in scope, the answer is almost certainly no. We are eight days from a deadline with a 5-page limit.

---

## 1. Team and lanes

| Person | Lane | Owns |
|---|---|---|
| **Oren** | **Lane 1 — Data & preprocessing** | Everything from raw download to labelled, genotyped, metacell-aggregated AnnData. Then the CGGA replication harness and Figure 1. **This is the person you are usually talking to.** |
| **Arnav** | Lane 2 — Modelling & statistics | RAS construction, Stage A residualization, Stage B, stability selection, the permutation test, leakage/disjointness tests, Figure 2. |
| **Devarsh** | Lane 3 — Paper, figures, venue | IEEE template, all prose, references, page-count discipline, release forms and parent documentation. Also floats to help Lanes 1 and 2 when they are blocked. |

Each of us runs Claude Code in a separate clone of the same repo. This file is committed and shared, so all three sessions inherit the same rules.

**Lane 2 works against simulated data until Day 3.** Arnav is writing and testing Stage A/B against a synthetic matrix with the same schema as the real object, so the modelling code is already debugged when Lane 1 hands over. If you are helping Lane 2, do not wait for real data. If you are helping Lane 1, the schema you produce must match what Lane 2 simulated — column names and dtypes agreed in writing, in `docs/SCHEMA.md`.

**Open admin item:** the author list on the original protocol draft reads "Devarsh Aswin, Arnav Vishwakarma, and Liam Gershony." That does not match the current team. Devarsh must reconcile the author list before submission, and the affiliation must state that the first author is a high school student.

---

## 2. The scientific question

### 2.1 Clinical background

Glioblastoma (GBM) is the most aggressive common brain tumour. Standard care is maximal safe surgical resection, then radiation with concurrent temozolomide. Recurrence is essentially universal, usually within six to ten months. A second surgery yields a **matched pair**: tissue from the same patient before treatment and after recurrence.

Prior longitudinal work found that primary and recurrent tumours from one patient are genetically similar — recurrence is not mostly a new mutant clone taking over. It is a shift in the *mix* of transcriptional states, most consistently toward a mesenchymal-like phenotype.

### 2.2 The question that follows

If recurrence is a shift in state mixture rather than genetic selection, the interesting question moves upstream: **are the cells that persist through treatment already distinguishable, molecularly, inside the untreated primary tumour?**

Answering that requires separating a genuine persistence-associated programme from two confounders that are easy to mistake for it:

1. **Cell state** — the transcriptional "job role" a cell has adopted.
2. **Genomic background** — which copy-number clone the cell belongs to.

A naive comparison of persisting versus non-persisting cells recovers those two long before it recovers anything specific to persistence.

### 2.3 What kind of paper this actually is

**This is a methods paper wearing a biology paper's clothes.** The headline is not "we found the recurrence genes." It is: *analyses that construct a target from expression data and then search expression data for correlates of that target are structurally circular; we built a control that removes this route by construction, and report what survives it.*

Keep this straight in all prose. It is what makes a null result publishable rather than embarrassing, and it is the reason the paper can be five pages.

---

## 3. Methodology

### 3.1 The proxy problem

We cannot follow a cell forward in time — sequencing destroys it. So we build a proxy: the **Recurrence Association Score (RAS)**, one number per primary-tumour cell, meaning roughly *how much does this cell resemble what this same patient's tumour looked like when it grew back?*

### 3.2 RAS components

| Component | Definition | Source |
|---|---|---|
| **T** — transcriptional similarity | Mean cosine similarity between cell *i*'s Harmony-corrected PCA embedding and the centroid of patient *p(i)*'s recurrent-timepoint cells in the same embedding. | Expression |
| **G** — genomic compatibility | 1 if cell *i*'s clone ID is also observed among patient *p(i)*'s recurrent-timepoint clones, else 0. | CNV |
| **O** — transport mass | Total optimal-transport mass moved from cell *i*'s metacell to any recurrent-timepoint metacell of patient *p(i)*, using a PCA-Euclidean cost matrix. | Expression |
| **Ab(state)** — state abundance shift | log2 fold-change in abundance of cell *i*'s assigned Neftel state, primary to recurrent, for patient *p(i)*. | Expression |
| **Ab(clone)** — clone abundance shift | log2 fold-change in abundance of cell *i*'s clone ID, primary to recurrent, for patient *p(i)*. | CNV |

**Optimal transport**, for reference: given piles of sand at some locations and holes at others, what is the cheapest way to move sand into holes? Primary metacells are the piles, recurrent metacells the holes, distance in expression space is the cost.

### 3.3 The two tiers we actually build

```
Tier A (full):        RAS_A = 0.25·z(T) + 0.25·z(G) + 0.25·z(O) + 0.25·z(Ab_state)
Tier C-disjoint:      RAS_C = 0.5·z(G) + 0.5·z(Ab_clone)
                              where G and Ab_clone are computed from inferCNV signal
                              restricted to chr7, chr9p, chr10 ONLY
```

Weights are **fixed a priori and never fit against any outcome**. This is deliberate: fitting them would open an avoidable data-dependent leakage route. Equal weighting sacrifices some predictive optimality to remove that risk entirely.

Tier B from the original protocol is dropped for scope. Say nothing about it in the paper.

### 3.4 The circularity problem, and why Tier C-disjoint exists

Three of Tier A's four components are built from expression data. Stage B then asks a model to *predict RAS from expression data*. The target is constructed from the same material used to predict it, so of course genes will be found. That is a tautology with p-values attached, not a discovery.

The original protocol's defence was "Tier C-strict": clone identity from patient-matched **exome** sequencing, containing no expression-derived quantity. **That is not buildable from open data.** Wang et al. used exome-seq only as internal validation of CNVs called from the nuclei; patient-level exome is not in the public deposit.

Our replacement is stronger in one specific way and weaker in another, and both must be stated honestly:

> **Tier C-disjoint.** Let S = {chr7, chr9p, chr10} — the canonical GBM copy-number regions (chr7 gain, chr10 loss, CDKN2A deletion at 9p21). Clone identity **G** and clone abundance change **Ab(clone)** are computed from inferCNV signal restricted to S. Stage B's feature matrix X **excludes every gene located in any region in S, where S = {chr7 (whole), chr9p (arm only — 9q genes remain eligible), chr10 (whole)}, resolved against cytoband hg38**.

- **Stronger than Tier C-strict:** no gene that Stage B can possibly select contributed any expression signal to constructing the target. This closes the *cis*-transcriptional route that the original protocol explicitly disclosed it could not close.
- **Weaker than Tier C-strict:** clone identity is still expression-inferred, just from expression on chromosomes whose genes are excluded from selection. The claim narrows from "no expression-derived quantity, direct or indirect" to **"no expression from any gene eligible for selection."** *Trans* effects remain unresolved.

Write the weaker sentence. Never the stronger one.

### 3.5 Stage A — residualization

Compare high-RAS to low-RAS cells directly and you mostly rediscover that mesenchymal cells look mesenchymal and EGFR-amplified cells look EGFR-amplified. Residualization removes that.

Intuition: predict height from age. A ten-year-old measures 150 cm, the model predicts 140, so the **residual** is +10 cm — meaning "tall *for their age*."

```
RAS_i = β0 + Σ_k β_k·1[state_i = k] + γ·genotype_i + u_p(i) + ε_i

u_p(i) ~ N(0, σ²_u)     ε_i ~ N(0, σ²_ε)
```

- **Fixed effects:** cell state (4–6 levels), genotype/EGFR-amplification status.
- **Random effect:** patient intercept, absorbing patient-level baseline shifts including unmeasured technical batch.
- **Cross-fitting:** for held-out patient *p*, fixed-effect coefficients are estimated on the other 18 patients only. The held-out patient's random effect is predicted as 0 (standard mixed-model out-of-sample practice). Residual = RAS − fixed-effect prediction.

The residual means: *more recurrence-like than you would expect knowing only this cell's type, its mutations, and which patient it came from.*

### 3.6 Stage B — discovery

```
target_i = f(X_i) + η_i
```

- **Adjusted arm target:** the Stage A residual.
- **Unadjusted arm target:** raw RAS_A (H1 only).
- **X:** gene expression (2,000 HVGs re-derived per training fold), Hallmark pathway scores, DoRothEA regulon scores. **Excludes** all genes on chr7/chr9p/chr10 for the Tier C-disjoint run.
- **f:** ElasticNet only. XGBoost dropped for compute — this is a declared deviation.

**X must never contain PCA coordinates, transport mass, or abundance statistics.** Those are RAS's own inputs; including them would let the model trivially reconstruct the target.

### 3.7 Guard rails

- **LOPO (leave-one-patient-out).** 19 patients, 19 folds, train on 18, test on 1. The unit is the *patient*, not the cell, because two cells from one patient are not independent observations.
- **Nadeau-Bengio variance correction.** The 19 folds share 17 of 18 training patients, so they are heavily correlated. Correction factor `1/19 + 1/18`. Applied to every LOPO-derived confidence interval.
- **Stability selection.** With 2,000 genes and 19 patients, ElasticNet's exact gene list is unstable. Resample patients 100 times, refit, keep genes chosen at least 30% of the time. The 30% comes from published small-sample benchmarking, not from taste.
- **Frozen thresholds.** Every cutoff was fixed before any result existed. See §5.

### 3.8 Metacell aggregation

Stage B fits **metacells (~30 nuclei each, ~3,000 units), not individual nuclei (~90,000).** This is a ~30x compute reduction that *also* directly reduces the within-patient pseudoreplication the original protocol disclosed as unresolved. Rare case where the cheap option is the more defensible one. Declare it in the paper as a deliberate design choice, never as a shortcut.

---

## 4. Datasets — all open, none gated

**Constraint: every dataset must be direct-download with no account, no data-use agreement, and no approval queue.** This was verified on 23 August 2026.

| Source | Access | Notes |
|---|---|---|
| **GSE174554** (GEO) | Open, anonymous FTP/HTTPS | Wang et al. 2022, *Nature Cancer*. 86 primary-recurrent matched specimens, **snRNA-seq**. 76 IDHwt; 52 carry matched-pair identifiers. Our discovery cohort. |
| **`GSE174554_Tumor_normal_metadata.txt`** | Open, same series | The authors' own malignant vs non-malignant annotation. **Use this. Do not roll our own classifier.** |
| **CGGA mRNAseq_693** | Open, direct from cgga.org.cn | Read counts opened to public access June 2022. GBM: **140 primary, 109 recurrent**. Mandatory replication cohort. |
| **CGGA mRNAseq_325** | Open, same portal | GBM: **85 primary, 24 recurrent, 30 secondary**. Supportive replication cohort. |
| Neftel four-state signatures | Open | Published supplementary tables. |
| Hallmark / DoRothEA | Open | MSigDB; omnipath via decoupler. |

### 4.1 Datasets explicitly NOT used

- **GLASS.** Hosted on Synapse; requires an account and a data-use agreement. **Removed.** Do not reference it in code or prose.
- **GSE174554 exome data.** Not in the public deposit. **Removed.** This is why Tier C-disjoint replaced Tier C-strict.
- **Spatial subcohort.** Out of scope (was H4).

### 4.2 The honest cost of the GLASS→CGGA swap

GLASS and CGGA would have been two genuinely independent consortia. The two CGGA batches are two sequencing batches from **one** consortium and one population. They are only semi-independent.

This exact sentence, or one like it, goes in the limitations section:

> *Replication was assessed in two batches of a single cohort, which controls for batch but not for population or ascertainment; independent replication in a second consortium remains untested.*

### 4.3 CGGA replication test definition

For candidate gene *g* in cohort *c*:

- **Test:** logistic regression of recurrence status on expression(*g*) plus covariates (sequencing platform, IDH status, tumour purity).
- **Direction:** sign in discovery must match sign in mRNAseq_693. mRNAseq_325 must not contradict; its significance is not independently required.
- **Significance:** raw p < 0.05 in mRNAseq_693, then Benjamini-Hochberg FDR at 5% applied across all tested genes jointly, not per-cohort.
- **Effect size:** |log2 FC| ≥ 0.25, or odds ratio outside [0.8, 1.25].
- **Discordance:** if the two batches significantly disagree in direction, the gene is excluded and reported separately as discordant — never silently dropped.

**Mandatory caveat, every time this is discussed:** CGGA replication is bulk-level supporting evidence only. It is never single-cell validation of the specific persister population a gene was discovered in.

---

## 5. Frozen parameters

**All of these live in `configs/pipeline_config.yaml`. That file is frozen. It was committed before any data was processed, and the git timestamp is the proof.**

| Parameter | Value | Changed in Rev 2? |
|---|---|---|
| Master seed | 42 | no |
| Min genes per nucleus | 500 | no |
| **Max mitochondrial fraction** | **5%** | **yes** — was 20%; snRNA correction |
| Highly variable genes | 2,000, re-derived per training fold | no |
| Metacell size | ~30 nuclei | no |
| Tier A weights | 0.25 each on T, G, O, Ab(state) | no |
| **Tier C-disjoint chromosome set S** | **chr7, chr9p, chr10** | **new** |
| Stability threshold (primary) | 30% | no |
| Stability thresholds (sensitivity) | 50%, 80% | no |
| **Stability bootstraps** | **100** | **yes** — was 300 |
| Cross-validation | LOPO, 19 folds | no |
| Variance correction | Nadeau-Bengio, 1/19 + 1/18 | no |
| Permutation iterations (H3) | 1,000 | no |
| **Ablation resamples (H1, stretch)** | **200** | **yes** — was 1,000 |
| **Confirmatory alpha** | **0.05 / 2 = 0.025** | **yes** — family is now H1, H3 only |
| H1 decision bound | mean ΔR ≥ 10 pp AND interval excludes 0 | no |
| **inferCNV nuclei cap per patient** | **1,500** | **new** |

**On the alpha change:** Bonferroni divides by the size of the confirmatory family. With H2 out of scope the family is H1 and H3, so per-hypothesis alpha is 0.025, not the 0.0167 in the original protocol. Getting this right and explaining why it changed signals that we understand the correction rather than having copied it.

---

## 6. Hypotheses and scope

| ID | Status | Statement |
|---|---|---|
| **H3** | **Confirmatory — CORE DELIVERABLE** | Genes selected under Tier A overlap genes selected under Tier C-disjoint at a rate significantly above a permutation-derived chance baseline. Jaccard index vs 1,000-permutation null, alpha 0.025. |
| **H1** | Confirmatory — **STRETCH** | The confound-adjusted model's external replication rate exceeds the unadjusted baseline's by at least 10 percentage points, via 200 paired patient-level resamples, alpha 0.025. |
| H2 | **OUT OF SCOPE** | Dropped. Its equivalence margin (d = ±0.2) was justified by "half the pilot's observed MES effect, d = 0.4" — a pilot that does not exist. Dropping H2 removes our single worst credibility risk. |
| H4 | **OUT OF SCOPE** | Spatial immune-context interaction. No time. |
| H5 | **OUT OF SCOPE** | Persister Fraction Index survival model. No time. |

### 6.1 All three H1 outcomes are valid and reportable

- Mean ΔR ≥ 10 pp and interval excludes zero → confound-adjustment demonstrably adds value.
- Interval includes zero, or mean below 10 pp → the persister signal is likely dominated by state and genotype. **A legitimate, informative null.**
- Mean ΔR ≤ −10 pp and interval excludes zero → the unadjusted model wins. Also a finding, also reported.

With 19 patients and 200 resamples the interval will be wide and will very likely include zero. Say so, and say the resample count was pre-specified as reduced for compute. **A wide honest interval is a result. A narrow dishonest one is misconduct.**

---

## 7. Hard rules

These are not stylistic preferences. Violating any of them invalidates the paper.

1. **`configs/pipeline_config.yaml` is frozen.** Never edit it to make a result look better. If a parameter genuinely must change, append to `DEVIATIONS.md` with the reason and a timestamp. **Never edit silently.**
2. **Patient-level grouping is mandatory everywhere.** No cell from a held-out patient may influence its own prediction — including HVG selection, scaling, hyperparameter tuning, and stability-selection resampling. `tests/test_no_patient_leakage.py` asserts zero patient overlap between train and test in every fold at every step.
3. **Chromosome disjointness is enforced by test, not by care.** `tests/test_chr_disjoint.py` asserts that no gene on chr7, chr9p or chr10 appears in the Tier C-disjoint feature matrix.
4. **`data/raw/` is immutable.** Never write to it, never edit in place. Every downstream file must be regenerable from `data/raw/` + `src/` + `configs/`.
5. **Data is snRNA-seq (single-nucleus), not scRNA-seq.** Say "single-nucleus" in all prose. Nuclear preps are depleted of mitochondrial transcripts by design, which is why the QC ceiling is 5% and not 20%.
6. **Notebooks are never the source of a reported number.** Every figure and statistic in the paper comes from a script in `src/` that runs start to finish. Notebooks are exploratory only.
7. **Every stochastic step reads its seed from the config.** Master seed 42. Bootstraps increment the seed by resample index for reproducible-but-distinct draws.
8. **No causal language, ever.** No clinical-utility claims. No "we found the recurrence genes." Association only.
9. **Declared deviations are fine; undeclared ones are not.** Reviewers punish silence, not honesty.

---

## 8. Repository structure

```
gbm-persister/
  CLAUDE.md                      <- this file
  DEVIATIONS.md                  <- running log, append-only
  README.md                      <- how to reproduce every result from raw
  docs/
    SCHEMA.md                    <- agreed column names/dtypes for the Lane1->Lane2 handoff
  data/
    raw/                         <- IMMUTABLE
      GSE174554/
      CGGA/
      PROVENANCE.md              <- accession, download date, SHA256 per file
    processed/                   <- generated, regenerable
      01_qc.h5ad
      02_integrated.h5ad
      03_labeled.h5ad
      04_states.h5ad
      05_genotype.h5ad
      07_metacells.h5ad
  configs/
    pipeline_config.yaml         <- FROZEN
  envs/
    environment.yml              <- pinned versions
  src/
    00_download.py
    01_qc_integration.py
    02_cnv_genotype.py
    03_metacells_ot.py
    04_ras_construction.py
    05_stage_a_residualization.py
    06_stage_b.py
    07_stability_selection.py
    08_circularity_check.py
    09_cgga_replication.py
    10_ablation.py               <- H1 stretch goal only
  tests/
    test_no_patient_leakage.py
    test_chr_disjoint.py
  notebooks/                     <- exploratory only, never source of truth
  results/
    tables/
    figures/
    gene_lists/
```

`.gitignore` must exclude `data/raw/`, `data/processed/`, and `*.h5ad`. Do not commit 40 GB of AnnData.

---

## 9. Timeline and current position

**Today is 23 August 2026 (Day 0). Nothing has been executed yet.**

| Day | Date | Goal | Artefact |
|---|---|---|---|
| **0** | Sun 23 Aug | Freeze config; start GEO + CGGA downloads; define the 19-patient rule; repo and environment | Config committed; all three datasets on disk |
| 1 | Mon 24 Aug | QC (500 genes, 5% mito), Scrublet, tumour/normal merge + join-drop audit, Harmony, Neftel scoring; verify CGGA parses | `01_qc.h5ad`, `02_integrated.h5ad`, `04_states.h5ad` |
| 2 | Tue 25 Aug | infercnvpy capped at 1,500 nuclei/patient on chr7/9p/10; Leiden clones; metacells; scCODA abundance | `05_genotype.h5ad`, `clone_catalog.csv`, `07_metacells.h5ad` |
| 3 | Wed 26 Aug | Optimal transport (**hard stop 6pm**); RAS Tier A + Tier C-disjoint; Stage A; **run both tests today** | `ras_scores.csv`, `stage_a_residuals.csv`, tests passing |
| 4 | Thu 27 Aug | Stage B, stability selection (100 boots), **circularity check** | `circularity_check.csv`, Figures 1 and 2 — **core deliverable done** |
| 5 | Fri 28 Aug | CGGA replication harness; **launch H1 ablation before sleeping** (~2 hrs). Not running by 10pm → abandon H1 | H1 running or formally abandoned |
| 6 | Sat 29 Aug | Full 5-page draft in IEEE template; limitations section | Complete draft |
| 7 | Sun 30 Aug | Three read-throughs; verify every number against its artefact; trim to 5 pages; **submit** | Submitted PDF + submission ID |
| — | Mon 31 Aug | **Buffer only.** Do not plan to use it. | — |

### 9.1 Declared deviations to write into the paper

- Stability bootstraps reduced 300 → 100, for compute.
- ElasticNet alpha tuned once per outer fold and reused across that fold's bootstraps with warm starts, rather than a full inner CV inside every bootstrap. Slightly understates selection variance.
- XGBoost dropped; linear-only discovery.
- H1 at 200 resamples rather than 1,000, pre-specified as reduced for compute.
- If SEACells does not converge by end of Day 3, k-means metacells within each patient-timepoint are substituted.
- If optimal transport does not converge by end of Day 3, component **O** is dropped and **Tier A-reduced** (T, G, Ab(state)) is used and named as such. Never silently drop a component.

---

## 10. Known risks and open decisions

1. **The 19-patient discovery subset does not yet exist as a defined set.** The original protocol inherited it from a pilot that was never run. It must be derived by a **written, reproducible rule applied to the downloaded metadata**, decided by Oren and Arnav together *before* anything downstream is computed, and stated in the paper. Do not assume 19 is achievable until the metadata is read.
2. **inferCNV clone quality on chr7/9p/10 only** is unverified. If clone assignment is degenerate (e.g. one clone per patient), Tier C-disjoint carries no information and H3 cannot be tested. Check this on Day 2, not Day 4.

   **Check each region independently, not S as a whole.** Chr7 gain and chr10 loss are whole-chromosome events and will resolve fine from windowed expression. **CDKN2A deletion at 9p21 is often focal — on the order of a megabase — and inferCNV's window smoothing across ~100 genes can erase it entirely.** So 9p may contribute nothing, leaving clones defined off chr7 and chr10 alone. That is still workable — those are the canonical events — but it must be a Day 2 finding, not a Day 4 discovery. Before accepting `clone_catalog.csv`, verify that **each** of the three regions independently produces non-degenerate clone structure, and record per-region results. If 9p is uninformative, say so in the paper rather than implying all three regions contributed.
3. **Metacell aggregation changes the Stage A/Stage B interface.** Stage A residuals are per-nucleus; Stage B fits per-metacell. The aggregation rule (mean residual per metacell) must be written down and agreed, not improvised.
4. **CGGA covariate availability.** Tumour purity may not ship with the CGGA download. If it is absent, drop it from the covariate set and declare the change — do not substitute a proxy without saying so.
5. **In-person attendance is mandatory.** The symposium is in Dallas, TX. At least one student author must attend and present or the paper does not enter the IEEE proceedings. High school authors additionally need a parent present, a signed release form, and a signed information document. Devarsh owns this; it must be settled with families this week, not after the 22 September notification.

---

## 11. The paper

5 pages maximum, including figures, tables and references (high school first author; undergraduate would be 6). IEEE Computer Society Proceedings template. Single-blind review — no anonymisation needed.

| Section | Budget | Content |
|---|---|---|
| I. Introduction | 0.6 pp | Recurrence is universal; prior work describes the recurrent tumour after regrowth; the prior question is whether persisters are already distinguishable. Frame as methodological from the first paragraph. |
| II. Background | 0.4 pp | Neftel states, matched primary/recurrent snRNA design, the confounding problem. |
| III. Methods | 1.7 pp | RAS and components. **Tier A vs Tier C-disjoint with the chromosome-exclusion rule stated explicitly — this is our novelty, give it room.** Stage A with the equation. Metacell aggregation and why. LOPO and patient-level grouping. Stability selection. Permutation test. |
| IV. Results | 1.0 pp | Figure 1: cohort flow + RAS distribution by state. Figure 2: observed Jaccard vs permutation null. Gene counts at all three thresholds. Both tests reported as passing. H1 if it completed. |
| V. Discussion & Limitations | 0.8 pp | What the disjoint control does and does not establish. *Trans* effects unresolved. Single-consortium replication. Every declared deviation. |
| References | 0.5 pp | 7–10 entries. Add CGGA. **Drop GLASS.** |

### 11.1 The sentence the paper is built around

> *Analyses that construct a target from expression data and then search expression data for correlates of that target are structurally circular; we built a chromosome-disjoint negative control that removes this route by construction, and report what survives it.*

True whether the overlap is significant or not. Defensible with data we can actually download. More interesting than a gene list would have been.

---

## 12. How to work with us

- **Explain before writing.** Say what a script will do, and flag anything that could leak test data into training, before producing code.
- **Prefer small readable scripts** over one large one.
- **When a number lands, say which file produced it.** Provenance per claim.
- **Use plan mode for anything touching preprocessing or fold structure.** A wrong preprocessing decision silently poisons everything downstream.
- **Push back.** If an instruction here conflicts with something you find in the data, say so rather than working around it. If a result looks too good, say that too — on a project whose entire credibility rests on pre-registration, a suspiciously clean number is a bug report.
- **Do not invent numbers.** If a value is unknown, say it is unknown. Never fill a gap with a plausible-looking figure.
- **We are the authors.** The students must be able to explain every method in this document to a reviewer, in person, at a poster in Dallas. Write code we can read and understand, not code that merely runs.
