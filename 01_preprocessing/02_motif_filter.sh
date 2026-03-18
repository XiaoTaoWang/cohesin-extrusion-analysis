#!/bin/bash

# Define paths

# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(dirname "$(dirname "$SCRIPT_DIR")")}"
DATA_ROOT="${DATA_ROOT:-$(dirname "$PROJECT_ROOT")}"
GENOME_ROOT="${GENOME_ROOT:-$(dirname "$DATA_ROOT")/genome}"
# ──────────────────────────────────────────────────────────────────

MOTIF_BED="${PROJECT_ROOT}/ChIA-PET/2.CTCF_motif_filter/K562_CTCF_ENCFF221SKA.motifs_best_site_4col.bed"
SIGNAL_BW="${DATA_ROOT}/mly/ChIA-PET/bw_all/WT_CTCF_h019_153g_deep.RPGC.bw"
OUTPUT_BED="${PROJECT_ROOT}/ChIA-PET/2.CTCF_motif_filter/K562_CTCF_ENCFF221SKA.motifs_filtered_25k.bed"
SCRIPT_PATH="${PROJECT_ROOT}/ChIA-PET/2.CTCF_motif_filter/filter_ctcf_by_signal.py"

# Create the python script (content provided above)
cat << 'EOF' > "$SCRIPT_PATH"
import sys
import random
import numpy as np
import pyBigWig

def main():
    if len(sys.argv) < 4:
        print("Usage: python filter_ctcf_by_signal.py <motif_bed> <signal_bw> <output_bed> [target_count]")
        sys.exit(1)

    motif_bed = sys.argv[1]
    signal_bw = sys.argv[2]
    output_bed = sys.argv[3]
    target_count = int(sys.argv[4]) if len(sys.argv) > 4 else 25000

    print(f"Processing motifs from: {motif_bed}")
    print(f"Using signal from: {signal_bw}")
    
    # Open BigWig
    bw = pyBigWig.open(signal_bw)
    chroms = bw.chroms()
    chrom_list = list(chroms.keys())
    chrom_sizes = {k: v for k, v in chroms.items()}

    # 1. Read motifs and calculate signal
    motifs = []
    motif_len = 19 # Default, will update from data
    
    with open(motif_bed, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3: continue
            
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            # Update motif length estimate
            if len(motifs) == 0:
                motif_len = end - start

            if chrom not in chroms:
                continue
            
            # Get mean signal for the region
            # Handle potential out of bounds
            if start < 0 or end > chrom_sizes[chrom]:
                continue

            try:
                # type='mean' calculates average signal over the region
                val = bw.stats(chrom, start, end, type='mean')[0]
                if val is None: val = 0
            except RuntimeError:
                val = 0
            
            motifs.append((val, line.strip()))

    print(f"Total motifs read: {len(motifs)}")

    # 2. Sort motifs by signal descending
    motifs.sort(key=lambda x: x[0], reverse=True)

    # 3. Select top N motifs
    kept_motifs = motifs[:target_count]
    if not kept_motifs:
        print("No motifs found or kept.")
        bw.close()
        sys.exit(1)

    cutoff_signal = kept_motifs[-1][0]
    print(f"Signal cutoff for top {len(kept_motifs)} motifs: {cutoff_signal:.4f}")

    # 4. Permutation test for significance
    # Generate random regions to build a null distribution of signals
    print("Running permutation test with 100,000 random regions...")
    random_signals = []
    num_random = 100000
    
    for _ in range(num_random):
        c = random.choice(chrom_list)
        size = chrom_sizes[c]
        if size <= motif_len: continue
        
        r_start = random.randint(0, size - motif_len)
        r_end = r_start + motif_len
        
        try:
            val = bw.stats(c, r_start, r_end, type='mean')[0]
            if val is None: val = 0
            random_signals.append(val)
        except:
            pass
    
    random_signals = np.array(random_signals)
    
    # Calculate empirical p-value: fraction of random regions with signal >= cutoff
    n_above = np.sum(random_signals >= cutoff_signal)
    p_value = (n_above + 1) / (len(random_signals) + 1)
    
    print(f"Empirical P-value of the cutoff: {p_value:.6g}")
    
    if p_value > 0.05:
        print("WARNING: The cutoff signal is not significantly higher than random background (p > 0.05).")
    else:
        print("The cutoff signal is significant (p <= 0.05).")

    # 5. Write output
    with open(output_bed, 'w') as f:
        for _, line in kept_motifs:
            f.write(line + '\n')
    
    print(f"Written {len(kept_motifs)} filtered motifs to: {output_bed}")
    bw.close()

if __name__ == "__main__":
    main()
EOF

# Run the python script
# Ensure you have pyBigWig and numpy installed (pip install pyBigWig numpy)
python3 "$SCRIPT_PATH" "$MOTIF_BED" "$SIGNAL_BW" "$OUTPUT_BED" 25000
