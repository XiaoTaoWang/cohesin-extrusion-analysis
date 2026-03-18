#!/bin/bash

# 1. ChIA-PET RAD21 degron

# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(dirname "$(dirname "$SCRIPT_DIR")")}"
DATA_ROOT="${DATA_ROOT:-$(dirname "$PROJECT_ROOT")}"
GENOME_ROOT="${GENOME_ROOT:-$(dirname "$DATA_ROOT")/genome}"
# ──────────────────────────────────────────────────────────────────

python ${PROJECT_ROOT}/ChIA-PET/1.example/RAD21_degron/RAD21_6h_9h_1kb.py &

# 2. ChIA-PET RAD21 degron enhance
python ${PROJECT_ROOT}/ChIA-PET/1.example/RAD21_degron_enhance/RAD21_6h_9h_1kb.py &

# 3. Micro-C RAD21 degron enhance
python ${PROJECT_ROOT}/micro-c/1.example/RAD21_degron_enhance/RAD21_6h_9h_1kb.py &

# 4. Micro-C RAD21 degron
python ${PROJECT_ROOT}/micro-c/1.example/RAD21_degron/RAD21_6h_9h.py &

# 5. Micro-C CTCF degron enhance
python ${PROJECT_ROOT}/micro-c/1.example/CTCF_degron_enhance/CTCF_6h_9h_1kb.py &

# 6. Micro-C CTCF degron
python ${PROJECT_ROOT}/micro-c/1.example/CTCF_degron/CTCF_6h_9h_1kb.py &

# Wait for all background tasks to finish
wait
echo "All plotting tasks finished."
