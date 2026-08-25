#!/usr/bin/env bash
# Full pipeline re-run from data/raw/. Every derived artefact is purged first so
# nothing can be silently reused. Stage order follows the real dependency graph.
set -uo pipefail
cd "$(dirname "$0")/.."
PY="$1"
echo "REPRO RUN started $(date -u +%Y-%m-%dT%H:%M:%SZ)"

echo "--- step 0: verify data/raw checksums ---"
"$PY" -u scripts/verify_raw.py

echo "--- purging all derived artefacts ---"
rm -rf data/interim/GSE174554_RAW data/interim/cnv_input data/interim/cnv_out \
       data/interim/metacell_expression.h5ad data/processed/*.h5ad \
       results/tables/*.csv results/tables/*.json results/gene_lists/*.csv \
       data/raw/neftel_signatures/neftel_metamodules.tsv 2>/dev/null
echo "  purged"

run () { echo "--- $1 ---"; "$PY" -u "src/$1" ${2:-} || echo "  !! $1 exited non-zero"; }

run 00h_build_neftel_signatures.py
run 00c_parse_sample_manifest.py
run 00e_join_idh_status.py
run 00f_compare_pairing.py
run 01a_build_cohort.py
run 01b_unpack_raw.py
run 01_qc_integration.py
run 01c_clause_d_gate.py
run 01d_diagnose_clause_d_failures.py
run 02_integration.py
run 02_cnv_genotype.py "--workers 5"
run 02b_clone_degeneracy_check.py
run 02c_clone_validity_check.py
run 02d_genotype_class_degeneracy.py
run 03_metacells_ot.py
run 03a_ras_component_diagnostics.py
run 04_states.py
run 04_ras_construction.py
run 05_stage_a_residualization.py
run 06a_build_stage_b_targets.py
run 06_stage_b.py "--workers 8"
run 08_circularity_check.py
run 10_ablation.py
echo "REPRO RUN complete $(date -u +%Y-%m-%dT%H:%M:%SZ)"
