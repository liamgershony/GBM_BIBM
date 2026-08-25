# RESULTS SUMMARY — every number the paper will cite

**Generated 25 August 2026.** Written for someone who was not in the analysis
session. Every figure names the script that produced it and the table it lives in.
Nothing here is quoted from memory.

> **Read §7 before writing the Results section.** Several results are statistically
> significant but rest on very small gene lists, and the paper must say so.

---

## 1. Cohort

| Quantity | Value | Source |
|---|---|---|
| GSMs in GSE174554 | 113 | `sample_manifest.csv` ← `00c` |
| Human snRNA-seq GSMs | 81 | same |
| Matched pairs (clauses a+b) | 30 | `pairing_comparison.csv` ← `00f` |
| IDH-wildtype pairs (clause c) | **29** | `discovery_cohort.csv` ← `01a` |
| Specimens / sequencing libraries | 61 / 68 | same |
| Nuclei after QC | 141,113 | `01_qc.h5ad` ← `01_qc_integration.py` |
| — malignant | 79,757 | same |
| — non-malignant (inferCNV reference) | 53,931 | same |
| — unknown | 7,425 | same |
| **Patients passing clause (d)** | **21** | `cohort_n.json` ← `01c` |
| Nuclei in the 21-patient cohort | 121,268 | `02_integrated.h5ad` |

**Derived at n = 21:** folds = 21; Nadeau-Bengio `1/21 + 1/20 = 0.097619`;
evaluability floor `floor(21/2)+1 = 11`. H1 remains admissible (threshold n < 16).

**Cohort discrepancy to disclose.** GEO's series summary claims 40 matched
IDH-wildtype pairs. That is not reproducible: Supplementary Table 1 supports 36
pairs flagged snRNA-seq, but only 30 have both specimens deposited as human
snRNA-seq GSMs — 86 specimens are flagged snRNA-seq in the paper against 78 in the
deposit. A deposit-versus-publication discrepancy, not an analysis choice.

**Clause (d) exclusions:** 8 of 29 patients. Diagnosed individually in
`clause_d_failure_diagnosis.csv` (`01d`): 5 QC attrition (all from the <500-gene
filter, on genuinely low-complexity libraries), 2 annotated non-tumour (the authors
called ~all nuclei Normal — patient 28: 3,595 of 3,598), 1 genuine sparsity, 1
annotation gap (`SF11981` is absent from the authors' annotation under any name).

---

## 2. Cell states

Neftel six meta-modules collapsed to four **before** the argmax
(`00h_build_neftel_signatures.py`): MES ← MES1+MES2 (95 unique genes), NPC ←
NPC1+NPC2 (89), AC (39), OPC (50). Cell-cycle programs G1/S and G2/M excluded.

Signature-gene coverage: AC 38/39 (97%), MES 94/95 (99%), NPC 87/89 (98%),
OPC 48/50 (96%).

**Evaluability against the §6.2 floor of 11 patients** (`state_assignments.csv` ← `04_states.py`):

| State | Patients with ≥20 nuclei | Verdict |
|---|---|---|
| AC | 20/21 | EVALUABLE |
| MES | 17/21 | EVALUABLE |
| NPC | 21/21 | EVALUABLE |
| OPC | 21/21 | EVALUABLE |

**No state was excluded**, so Stage A's fixed effect carries all four levels.
23.1% of nuclei have a top-to-runner-up score margin below 0.05 — consistent with
Neftel's hybrid fraction; those cells receive a hard label that Stage A treats as
certain.

---

## 3. Negative controls — three failures, all reportable

Full write-up in `negative_control_failures.md`. Summary:

**(a) Unsupervised clone calls are indistinguishable from noise** (`02c`). Leiden
at resolution 1.0 versus a permutation null that shuffles each CNV window across
cells:

| Region | Observed clones (median) | Null (median) | Ratio |
|---|---|---|---|
| chr7 | 12 | 14 | 0.86 |
| chr9p | 10 | 20 | 0.50 |
| chr10 | 13 | 14 | 0.93 |

Every region's observed count is at or **below** its null.

**(b) chr9p carries no detectable CNV signal** (`02c`):

| Region | Δ (malignant − reference) | Patients with expected sign |
|---|---|---|
| chr7 (gain) | **+0.0622** | **21/21** |
| chr10 (loss) | **−0.0340** | **21/21** |
| chr9p (loss) | −0.0005 | 9/21 — chance |

**Clone identity rests on chr7 and chr10 only.** The paper must not describe S as
three contributing regions. chr9p is retained in the *gene-exclusion* rule, where
it still does work.

**(c) G is constant in 19/21 patients** (`02d`). Thresholded chr7×chr10 classes do
not collapse — all 4 classes appear in all 21 patients, median dominant-class
fraction 48% — but G is a *set-membership* test, so with coarse classes every class
is trivially observed at recurrence. The two exceptions (patients 5 and 19) both
sit in the shallow tail of recurrent depth, so where G varies it partly tracks
sequencing depth rather than persistence.

**(d) Component O is degenerate by construction** (`03_metacells_ot.py`).
`corr(O, source marginal) = 1.000000`, max difference 2.25e-16. Balanced EMD
conserves mass, so the transport-plan row sums equal the source weights as an
algebraic identity. **O was dropped** per the §9.1 contingency and the score is
named **Tier A-reduced**. A transport-cost alternative was computed and explicitly
**not** substituted.

**Metacells:** SEACells converged in 33/42 patient-timepoint groups; the k-means
contingency fired in 9 (7 time-box, 2 numerical), recorded per group in
`metacell_catalog.csv`. 2,426 metacells, median 30.0 nuclei each.

---

## 4. RAS and Stage A

`ras_tier_a_reduced = (1/3)z(T) + (1/3)z(G) + (1/3)z(Ab_state)`
`ras_tier_c_disjoint = 0.5·z(G) + 0.5·z(Ab_clone)`
z() pooled across the cohort (within-patient z is a divide-by-zero for G in 19/21).
14,710 primary malignant nuclei, 21 patients (`ras_scores.csv` ← `04`).

**Stage A** (`stage_a_r2.csv` ← `05`), cross-fitted LOPO, held-out random effect = 0:

| Tier | R², state+genotype | R², +patient | Variance remaining |
|---|---|---|---|
| A-reduced | **0.1112** | 0.5252 | 47.48% |
| C-disjoint | **0.1772** | 0.5572 | 44.28% |

**STOP/GO not triggered on either tier.**

**Interpretation the paper must carry:** state and genotype explain only 11–18%;
nearly all of Stage A's power is the patient term. **The adjustment is
predominantly patient-level centring, not removal of the named confounds.**

**The two tiers are not independent** — both contain z(G):

| | |
|---|---|
| corr(Tier A, Tier C) raw | +0.5634 |
| corr with shared G removed | **+0.1366** |
| corr of Stage A residuals | +0.6521 |

Roughly three quarters of the between-tier correlation is the shared G term.

---

## 5. Stage B

1,413 metacells carry a target (of 2,426); 2,000 features per fold — **p > n**.
Chromosome disjointness enforced by `tests/test_chr_disjoint.py`; **7/7 tests
pass**, including that chr9q remains eligible.

**Selected gene counts** (`stage_b_summary.json`, `results/gene_lists/`):

| Arm | 30% | 50% | 80% |
|---|---|---|---|
| v1_tierA (primary) | **46** | 4 | 0 |
| v1_tierC (primary) | **3** | 0 | 0 |
| v2_tierA (centred) | 3 | 0 | 0 |
| v2_tierC (centred) | 1 | 0 | 0 |
| v3_tierA (G-removed) | **43** | 5 | 0 |
| v3_tierC (G-removed) | **11** | 1 | 0 |

**No arm selects a single gene at the 80% threshold.**

**Per-fold results are in `stage_b_folds.csv` (§4.6), not collapsed.** The fold
imbalance is visible directly in the tuned alpha: folds holding out large patients
(9 with 206 test metacells, 14 with 189) tune to alpha ≈ 0.005–0.02 and select
1,400–1,900 genes, while folds holding out small patients (10 with 4, 36 with 5)
tune to alpha ≈ 0.20 and select ~240. **Nadeau-Bengio assumes exchangeable folds;
these are not exchangeable, so every LOPO interval is narrower than the truth.** No
capping or rebalancing was applied — that would be a post-hoc change to fold
structure.

---

## 6. Hypothesis results

### H3 — circularity check (CORE DELIVERABLE)

`circularity_check.csv` ← `08`. α = 0.025 (Bonferroni over the H1/H3 family).
1,000 permutations, sampling from **each arm's own** eligible universe (Tier A
33,694 genes; Tier C 31,744 after chr7/chr9p/chr10 exclusion).

| Variant | Role | Target corr | \|A\| | \|C\| | Overlap | Jaccard obs | Null mean | p | Significant |
|---|---|---|---|---|---|---|---|---|---|
| **v1** | **PRIMARY (pre-registered)** | 0.652 | 46 | 3 | **2** | **0.0426** | 0.0001 | **0.0010** | **YES** |
| v2 | SENSITIVITY (patient-centred) | 0.412 | 3 | 1 | 0 | 0.0000 | 0.0000 | 1.0000 | no |
| v3 | SENSITIVITY (G-removed) | 0.276 | 43 | 11 | **2** | 0.0385 | 0.0005 | **0.0020** | YES |

**v1 is the registered result and stands.** v2 and v3 are declared sensitivity
analyses; no winner is designated.

Overlapping genes — **note they are not the same set across variants**:
- v1: `EXTL1`, `RP11-290O12.2`
- v3: `CH17-80A12.1`, `P2RY14`

### H1 — ablation (STRETCH)

`h1_ablation.csv` ← `10`. 200 paired patient-level resamples.

| Quantity | Value |
|---|---|
| Mean replication rate, **adjusted** | **0.2801** |
| Mean replication rate, **unadjusted** | **0.2461** |
| Mean ΔR | **+0.0341 (+3.41 pp)** |
| 95% percentile interval | **[−0.0376, +0.1122]** = [−3.76, +11.22] pp |
| Interval excludes 0 | **No** |
| Mean ΔR ≥ +10 pp | **No** |
| Resamples with adjusted > unadjusted | 165/200 (82.5%) |
| Median genes selected per resample | 91 adjusted, 88 unadjusted |

**OUTCOME: informative null (§6.1, second case).** The interval includes zero and
the mean is below the 10 pp bound. This is a valid, pre-specified, reportable
result.

**Essential context:** **24.1% of all genes tested in CGGA (7,182 of 29,796) meet
the replication criteria.** The unadjusted arm's 24.6% is therefore
indistinguishable from the background rate, and the adjusted arm's 28.0% is 3.9 pp
above it. A difference between two rates that both sit near a permissive background
means far less than the same difference between two selective ones.

### CGGA replication harness

`09_cgga_replication.py`, self-tested end to end. Cohorts reproduce §4 exactly:
693 = 140 primary / 109 recurrent; 325 = 85 / 24 / 30 secondary. Filtering on
`Histology == "GBM"` alone yields 140 primaries and **zero** recurrent, because
recurrent tumours carry an `r` prefix — the filter is `{GBM, rGBM, sGBM}`.

**Covariates:** IDH status used. **Tumour purity is absent from the CGGA download**
— dropped and declared per §10.4, no proxy substituted. **Sequencing platform has
no column and is constant within each cohort**, so it cannot enter a within-cohort
model. Both recorded in the output table.

---

## 7. What the Results section must say

1. **H3 v1 is significant (p = 0.0010) but the effect is two genes.** Tier C's 30%
   list has **three** genes; the overlap is two of them. A p-value computed against
   a near-zero null is easy to clear with a tiny list. Report the Jaccard, the list
   sizes, and the overlapping gene names together — the p-value alone overstates it.
2. **The overlapping genes differ between v1 and v3.** The overlap is not a stable
   gene set, so it does not support a claim about specific genes.
3. **The two tiers share z(G) by construction**, contributing ~3/4 of their
   correlation. A permutation null over gene labels does not account for correlated
   targets, so some overlap is expected structurally. The *cis* route is still
   closed for Tier C (its features exclude all chr7/chr9p/chr10 genes); what remains
   shared is the *trans* route and the patient offset.
4. **v2 found essentially nothing** (3 and 1 genes, zero overlap). When patient
   structure is removed — the thing §3.5 says the patient intercept is for — the
   signal does not survive.
5. **H1 is an informative null**, and both arms sit near CGGA's 24.1% background.
6. **No arm selects any gene at the 80% stability threshold**, and only 4–5 genes
   at 50%. The discovery is fragile at every threshold above the primary one.
7. **Three RAS components degenerated**: O by construction, G within-patient, and
   chr9p contributed no CNV signal. Tier A-reduced is three equally *weighted*
   components of which one is within-patient constant. Weights were never changed.

---

## 8. Reproducing this

Every input is regenerable by `src/00*.py` **except one**: Neftel Table S2
(`NIHMS1532254-supplement-9.xlsx`, SHA256
`208e73ab3d22c494caf85c867d69dc6be38df3fc62ab1f043d7fcc5441066277`) must be
downloaded in a browser — PMC gates the path with a proof-of-work anti-bot
challenge. `00h` asserts that hash before reading.

Order: `00_download` → `00b` → `00d` → `00g` → `00h` → `00c` → `00e` → `00f` →
`01a` → `01b` → `01_qc_integration` → `01c` → `02_cnv_genotype` → `02b`/`02c`/`02d`
→ `03_metacells_ot` → `03a` → `04_states` → `04_ras_construction` → `05` → `06a` →
`06` → `08` → `09` → `10`.

All deviations, corrections and specification gaps are in `DEVIATIONS.md`,
append-only and dated.
