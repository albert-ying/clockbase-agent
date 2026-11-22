#!/usr/bin/env bash
# Sequential full run: DNAm first (fast), then RNA (long).
# Each modality gets its own output / status dir so they don't collide.
set -uo pipefail
cd "$(dirname "$0")"

mkdir -p logs
LOG_DNAM="logs/run_hs_dnam_$(date +%Y%m%d_%H%M%S).log"
LOG_RNA="logs/run_hs_rna_$(date +%Y%m%d_%H%M%S).log"

# ---- DNAm (307 rows) ----
echo "=== $(date)  starting DNAm run ===" | tee -a "$LOG_DNAM"
CSV_MODALITY=hs_dnam NUM_WORKERS=4 \
  .venv/bin/python src/run.py >> "$LOG_DNAM" 2>&1
echo "=== $(date)  DNAm finished (exit $?) ===" | tee -a "$LOG_DNAM"
# Preserve DNAm outputs so RNA doesn't overwrite them
mv output        output_hs_dnam        2>/dev/null
mv status.json   status_hs_dnam.json   2>/dev/null

# ---- RNA (3,494 rows) ----
echo "=== $(date)  starting RNA run ===" | tee -a "$LOG_RNA"
CSV_MODALITY=hs_rna NUM_WORKERS=12 \
  .venv/bin/python src/run.py >> "$LOG_RNA" 2>&1
echo "=== $(date)  RNA finished (exit $?) ===" | tee -a "$LOG_RNA"
mv output        output_hs_rna         2>/dev/null
mv status.json   status_hs_rna.json    2>/dev/null

echo "=== $(date)  ALL DONE ===" | tee -a logs/summary.log
