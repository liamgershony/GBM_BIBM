# Discovery cohort rule — DRAFT, awaiting Lane 2 sign-off

CLAUDE.md §10.1 requires the discovery subset to be defined by a **written,
reproducible rule applied to the downloaded metadata**, agreed by Lane 1 and Lane 2
**before anything downstream is computed**, and stated in the paper.

**Status:** drafted by Liam (Lane 1), 23 Aug 2026. **Not yet agreed with Arnav.**
No downstream computation may use it until it is.

---

## The rule

A patient enters the discovery cohort if **all four** hold:

- **(a) Pairing.** Wang et al. **Supplementary Table 1 `Pair#`** links the patient to
  at least one `Primary` and at least one `Recurrent` specimen.
  *GEO's `pair#` characteristic is recorded as a cross-check, not as the source.*
- **(b) Assay.** Both specimens are present in **GSE174554 as human snRNA-seq** GSMs.
- **(c) Genotype.** **IDH-wildtype** per Supplementary Table 1.
- **(d) Depth.** **≥100 usable nuclei at both timepoints after QC**
  (≥500 genes/nucleus, ≤5% mitochondrial, Scrublet doublets removed).

**Multiple specimens at one timepoint are pooled**, with `sample_id` retained as a
batch key so pooling never hides a technical batch.

### Why pairing comes from the supplementary table

Both records were compared (`results/tables/pairing_comparison.csv`, produced by
`src/00f_compare_pairing.py`). They agree on 30 matched pairs and neither recovers a
pair the other loses. The supplementary table is preferred because it is the
authors' curated record rather than the submission form, and because it carries IDH
status, which GEO does not. GEO `pair#` is retained per patient for cross-check.

**Both records use a literal `"NA"` in their pair column.** Treated as an identifier
it fuses every unpaired specimen into one fabricated patient — 12 specimens under
GEO's scheme, 10 under the supplementary scheme. Any implementation must normalise
`"NA"` to missing in **both** files. This has already produced two wrong counts
during exploration and is the single easiest error to make here.

---

## What the rule yields on today's metadata

Clauses (a)–(c) are evaluable now. Clause (d) is not — it needs Day 1 QC.

| Filter | Matched pairs |
|---|---|
| Supp Table 1, any assay | 45 |
| Supp Table 1, flagged `snRNA-seq = Y` | 36 |
| …and IDH-wildtype | 34 |
| …**and both specimens present in GEO as human snRNA-seq** | 30 |
| **…and IDH-wildtype → clauses (a)+(b)+(c)** | **29** |

**Pre-QC ceiling: n = 29.** Clause (d) can only reduce it.

### Consequences at n = 29

- `n_folds = 29`; Nadeau-Bengio correction `1/29 + 1/28` (computed at runtime — the
  frozen config encodes the formula, not the literal).
- Evaluability floor `floor(29/2) + 1 = 15` patients per cell state (CLAUDE.md §6.2).
- **H1 remains admissible** — the DEVIATIONS threshold drops H1 only at n < 16.
- inferCNV budget ≈ **1.53×** the 19-patient plan.

---

## Open item for Lane 2

The pre-specified n was 19; the realised pre-QC ceiling is 29. Confirm the
DEVIATIONS.md cascade is the agreed handling, and that no analysis parameter other
than fold count and the evaluability floor scales with n.
