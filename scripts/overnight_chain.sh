#!/usr/bin/env bash
# Wait for Stage B, then run H3, then H1. One overnight chain.
# Emits a heartbeat immediately and while waiting, so an idle wait is
# distinguishable from a dead process.
set -uo pipefail
cd "$(dirname "$0")/.."
PY="$1"
echo "chain started $(date -u +%Y-%m-%dT%H:%M:%SZ); waiting for Stage B (pid $(cat logs/06_stage_b.pid))"
while kill -0 "$(cat logs/06_stage_b.pid)" 2>/dev/null; do
  echo "  waiting: $(grep -c done logs/06_stage_b.log)/126 Stage B jobs done ($(date -u +%H:%M:%SZ))"
  sleep 120
done
echo "=== Stage B finished $(date -u +%H:%M:%SZ), running H3 ==="
"$PY" -u src/08_circularity_check.py
echo "=== H3 done $(date -u +%H:%M:%SZ), running H1 ==="
"$PY" -u src/10_ablation.py
echo "=== chain complete $(date -u +%H:%M:%SZ) ==="
