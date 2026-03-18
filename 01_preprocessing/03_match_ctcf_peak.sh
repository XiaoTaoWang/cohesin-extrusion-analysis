
# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(dirname "$(dirname "$SCRIPT_DIR")")}"
DATA_ROOT="${DATA_ROOT:-$(dirname "$PROJECT_ROOT")}"
GENOME_ROOT="${GENOME_ROOT:-$(dirname "$DATA_ROOT")/genome}"
# ──────────────────────────────────────────────────────────────────

python 12_match_ctcf_peaks.py mESC mm10 \
    ${PROJECT_ROOT}/mESC/loading_site_define/ren_mm10.motifs.txt \
    ${PROJECT_ROOT}/mESC/chip-seq/NIPBL/peak/mESC_NIPBL_chip-seq_peak.bed
