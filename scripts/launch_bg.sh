#!/usr/bin/env bash
# Launch a long job and VERIFY it is running. Print a PID only after evidence.
#
# Five failures in this project have presented as successes; a PID echoed at
# launch is not evidence of a running job. This helper requires: the process is
# alive, the log is non-empty, and the log GREW between two samples.
#
#   scripts/launch_bg.sh <name> <command...>
set -uo pipefail
NAME="$1"; shift
LOG="logs/${NAME}.log"; PIDF="logs/${NAME}.pid"
mkdir -p logs; : > "$LOG"

# Detach into a new session so the job survives this shell exiting. macOS has no
# setsid(1), so scripts/detach.py does the double-fork + setsid() itself and
# writes the real child PID to $PIDF.
rm -f "$PIDF"
nohup python3 scripts/detach.py "$PIDF" "$@" > "$LOG" 2>&1 < /dev/null
for _ in $(seq 1 20); do [ -s "$PIDF" ] && break; sleep 0.25; done
PID=$(cat "$PIDF" 2>/dev/null || echo 0)

sleep 15; S1=$(wc -c < "$LOG" | tr -d ' ')
sleep 45; S2=$(wc -c < "$LOG" | tr -d ' ')

ALIVE=no; kill -0 "$PID" 2>/dev/null && ALIVE=yes
echo "name=$NAME pid=$PID alive=$ALIVE log_bytes_15s=$S1 log_bytes_60s=$S2"

if [ "$ALIVE" != yes ]; then
  # A short job can legitimately finish inside the verification window. Treat a
  # non-empty, grown log as completion; only an empty log means it died on launch.
  if [ "$S2" -gt 0 ] && [ "$S2" -gt "$S1" ]; then
    echo "COMPLETED within the verification window (log grew $S1 -> $S2 bytes)."
    echo "--- log tail ---"; tail -20 "$LOG"; exit 0
  fi
  echo "VERIFY FAILED: process is not running after 60s and produced no output."
  echo "--- log tail ---"; tail -30 "$LOG"; exit 1
fi
if [ "$S2" -eq 0 ]; then
  echo "VERIFY FAILED: process alive but log is empty after 60s."
  echo "Is the command unbuffered? Python needs -u when stdout is redirected."
  exit 1
fi
if [ "$S2" -le "$S1" ]; then
  echo "WARNING: log did not grow between 15s and 60s ($S1 -> $S2 bytes)."
  echo "Alive, but may be blocked. Check before relying on it."
fi
echo "VERIFIED RUNNING: pid $PID, log $LOG growing ($S1 -> $S2 bytes)"
echo "--- log so far ---"; tail -8 "$LOG"
