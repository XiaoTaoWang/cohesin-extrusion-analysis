# ==============================================================================
# Script Name: batch_chiadrop_predefined_regions.R
# Function:    ChIA-Drop Analysis based on Pre-filtered BED File
#              Format: cr{ID}_S (Left), cr{ID}_E (Right), cr{ID}_M-{num} (Loading)
# Date:        2026-02-03
# ==============================================================================


# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = dirname(dirname(normalizePath(sys.frame(1)$ofile))))
DATA_ROOT    <- Sys.getenv("DATA_ROOT", unset = dirname(PROJECT_ROOT))
GENOME_ROOT  <- Sys.getenv("GENOME_ROOT", unset = file.path(dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(GenomicRanges)
  library(data.table)
  library(ggpubr)
  library(stringr)
})

# ==============================================================================
# 1. Configure file paths
# ==============================================================================

# BED with all S/E/M definitions
FILE_REGIONS <- paste0(DATA_ROOT, "/ChIA-Drop/3.anchoring_middle/GM12878-cohesin-specific-regions_20260122_filt_25-75percentile_20260122.bed")

# RDS data paths
RDS_FILES <- list(
  "Cohesin" = paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between/GSE158897_GM12878-cohesin-pooled_comp.rds"),
  "CTCF"    = paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between/GSE158897_GM12878-CTCF-pooled_comp.rds"),
  "RNAPII"  = paste0(DATA_ROOT, "/ChIA-Drop/GSE158897/GSE158897_GM12878-RNAPII-pooledv2_comp.rds")
)

# Anchor buffer (expand anchors outward; set to 0 if S/E already precise)
# Allow small buffer for boundary fragments
ANCHOR_BUFFER <- 0 
OUTPUT_PDF    <- "ChIA_Drop_Predefined_Regions_Boxplot.pdf"
INTERMEDIATE_RDS <- "ChIA_Drop_Predefined_Regions_Stats.rds"

# ==============================================================================
# 2. Parse BED (S/E/M matching)
# ==============================================================================

message(">>> 1. Parsing Region Definition File...")

raw_bed <- fread(FILE_REGIONS, header=FALSE, col.names=c("chr", "start", "end", "label"))

# Parse label: cr1_S -> loop_id=cr1, type=S
raw_bed[, c("loop_id", "tag") := tstrsplit(label, "_", fixed=TRUE, keep=1:2)]

# Split S, E, M
# S: Left Anchor
anchors_s <- raw_bed[tag == "S", .(loop_id, chr, s_start=start, s_end=end)]

# E: Right Anchor
anchors_e <- raw_bed[tag == "E", .(loop_id, e_start=start, e_end=end)]

# M: Loading Sites
loading_m <- raw_bed[grepl("^M", tag), .(loop_id, m_label=label, m_start=start, m_end=end)]

# Combine into "analysis pairs"
# Each loading site is an independent analysis row linked to loop S and E
analysis_pairs <- merge(loading_m, anchors_s, by="loop_id")
analysis_pairs <- merge(analysis_pairs, anchors_e, by="loop_id")

# Generate unique ID per row
analysis_pairs[, pair_id := 1:.N]
analysis_pairs[, m_mid := (m_start + m_end)/2]

# Compute loop search range (from S_start to E_end)
analysis_pairs[, `:=`(
  search_start = pmin(s_start, e_start),
  search_end   = pmax(s_end, e_end)
)]

message(sprintf("    Total Loading Sites (Analysis Units): %d", nrow(analysis_pairs)))
message(sprintf("    Unique Loops Covered: %d", length(unique(analysis_pairs$loop_id))))

# ==============================================================================
# 3. Processing and classification
# ==============================================================================

if (file.exists(INTERMEDIATE_RDS)) {
  message(sprintf(">>> Loading cached results from: %s", INTERMEDIATE_RDS))
  final_stats_list <- readRDS(INTERMEDIATE_RDS)
} else {
  final_stats_list <- list()

  for (sample_name in names(RDS_FILES)) {
    rds_file <- RDS_FILES[[sample_name]]
    message(sprintf("\n>>> Analyzing %s...", sample_name))
  
  if (!file.exists(rds_file)) {
    warning(paste("File missing:", rds_file))
    next
  }
  
  # 3.1 Load Data
  data_obj <- readRDS(rds_file)
  if (is(data_obj, "GRanges")) frags_dt <- as.data.table(data_obj) else frags_dt <- as.data.table(data_obj)
  setnames(frags_dt, "seqnames", "chr", skip_absent=TRUE)
  if (!"cluster_id" %in% colnames(frags_dt)) frags_dt[, cluster_id := CBMB]
  
  frags_gr <- makeGRangesFromDataFrame(frags_dt, keep.extra.columns = TRUE)
  
  # 3.2 Find Overlaps with Loop Regions
  # First find all fragments within the loop range
  # Note: pair_id used in GRanges may have overlaps; findOverlaps handles them
  target_gr <- makeGRangesFromDataFrame(
    analysis_pairs[, .(chr, start=search_start, end=search_end, pair_id)]
  )
  
  ov <- findOverlaps(frags_gr, target_gr)
  dt_ov <- data.table(frag_idx = queryHits(ov), pair_id = subjectHits(ov))
  
  # 3.3 Merge Coordinates
  # Merge loading/anchor coordinates
  dt_merged <- merge(dt_ov, analysis_pairs, by="pair_id")
  # Merge fragment coordinates
  dt_merged <- cbind(dt_merged, frags_dt[dt_merged$frag_idx, .(cluster_id, start, end)])
  dt_merged[, mid := (start + end)/2]
  
  message("    Classifying complexes...")
  
  # 3.4 State flags (relative to current pair definition)
  dt_merged[, `:=`(
    # 1. Loading Overlap (Exact M site)
    is_loading = mid >= m_start & mid <= m_end,
    
    # 2. Anchor Overlap (S or E)
    is_anchor  = (mid >= (s_start - ANCHOR_BUFFER) & mid <= (s_end + ANCHOR_BUFFER)) | 
                 (mid >= (e_start - ANCHOR_BUFFER) & mid <= (e_end + ANCHOR_BUFFER)),
    
    # 3. Strictly Left (Between S_end and M_start)
    is_left    = mid > (s_end + ANCHOR_BUFFER) & mid < m_start,
    
    # 4. Strictly Right (Between M_end and E_start)
    is_right   = mid > m_end & mid < (e_start - ANCHOR_BUFFER),
    
    # 5. Raw Sides (for simple Two-sided count requirement)
    raw_l      = mid < m_start,
    raw_r      = mid > m_end
  )]
  
  # 3.5 Aggregate to Cluster Level per Pair
  cl_summ <- dt_merged[, .(
    c_start = min(mid),
    c_end   = max(mid),
    
    hit_anch = any(is_anchor),
    hit_load = any(is_loading),
    
    cnt_strict_l = sum(is_left),
    cnt_strict_r = sum(is_right),
    
    raw_l = sum(raw_l),
    raw_r = sum(raw_r)
  ), by = .(pair_id, cluster_id)]
  
  cl_summ[, span := c_end - c_start]
  
  # 3.6 Classification Logic
  
  # A. Anchor Complex (Top Priority)
  cl_summ[, is_class_anchor := hit_anch]
  
  # B. Two-sided Candidate
  # No Anchor, No Loading, Left >=1, Right >=1 (Strictly raw counts logic first)
  cl_summ[, is_cand_ts := !hit_anch & !hit_load & raw_l >= 1 & raw_r >= 1]
  
  # C. One-sided Candidate
  # No Anchor, Hit Loading, (Strict Left > 0 OR Strict Right > 0)
  cl_summ[, is_cand_os := !hit_anch & hit_load & (cnt_strict_l > 0 | cnt_strict_r > 0)]
  
  # D. Progressive Nesting for Two-sided (Per Pair_ID)
  # Russian Doll filtering
  
  ts_data <- cl_summ[is_cand_ts == TRUE]
  
  if (nrow(ts_data) > 0) {
    # Must sort by span within each pair_id
    setorder(ts_data, pair_id, span)
    
    kept_ids <- ts_data[, {
      if (.N == 0) integer(0)
      else if (.N == 1) cluster_id
      else {
        keep <- logical(.N)
        keep[1] <- TRUE
        p_s <- c_start[1]
        p_e <- c_end[1]
        for (i in 2:.N) {
          # Nesting condition: Left is Lefter, Right is Righter
          if (c_start[i] < p_s && c_end[i] > p_e) {
            keep[i] <- TRUE
            p_s <- c_start[i]
            p_e <- c_end[i]
          }
        }
        cluster_id[keep]
      }
    }, by = pair_id]$V1
    
    # Create validation keys
    valid_keys <- paste(ts_data[cluster_id %in% kept_ids]$pair_id, 
                        ts_data[cluster_id %in% kept_ids]$cluster_id, sep="_")
    
    cl_summ[, key := paste(pair_id, cluster_id, sep="_")]
    cl_summ[, final_twosided := key %in% valid_keys]
  } else {
    cl_summ[, final_twosided := FALSE]
  }
  
  # 3.7 Summarize Counts
  pair_stats <- cl_summ[, .(
    n_anchor   = sum(is_class_anchor),
    n_twosided = sum(final_twosided),
    n_onesided = sum(is_cand_os)
  ), by = pair_id]
  
  # Fill potentially missing pairs (count = 0)
  all_pairs <- data.table(pair_id = analysis_pairs$pair_id)
  pair_stats <- merge(all_pairs, pair_stats, by="pair_id", all.x=TRUE)
  pair_stats[is.na(pair_stats)] <- 0
  
  # Melt for Plotting
  melted <- melt(pair_stats, id.vars="pair_id", measure.vars=c("n_anchor", "n_twosided", "n_onesided"),
                 variable.name="Type", value.name="Count")
  
  melted[, Type := factor(Type, levels=c("n_twosided", "n_onesided", "n_anchor"), 
                          labels=c("Two-sided", "One-sided", "Anchor"))]
  melted[, Sample := sample_name]
  
  final_stats_list[[sample_name]] <- melted
}
  saveRDS(final_stats_list, INTERMEDIATE_RDS)
}

# ==============================================================================
# 4. Plot (PDF output)
# ==============================================================================

if (length(final_stats_list) > 0) {
  plot_df <- rbindlist(final_stats_list)
  
  # Handle log10 (remove zeros for boxplot; or use log1p)
  # Original plot looked like zero values were removed
  plot_df_clean <- plot_df[Count > 0]
  plot_df_clean[, log_val := log10(Count)]
  
  # N count (Analysis Units = loading sites)
  n_stats <- plot_df[, .(N = length(unique(pair_id))), by=Sample]
  
  # Statistical tests
  comparisons <- list(c("Two-sided", "One-sided"), 
                      c("One-sided", "Anchor"), 
                      c("Two-sided", "Anchor"))
  
  # Custom colors (Cohesin=dark green, CTCF=blue, RNAPII=purple)
  sample_colors <- c("Cohesin"="#008000", "CTCF"="#0000F3", "RNAPII"="#8000b7")

  p <- ggplot(plot_df_clean, aes(x=Type, y=log_val)) +
    geom_boxplot(color="black", fill="white", outlier.shape=1, outlier.size=1, width=0.5, lwd=0.6, fatten=0, show.legend=FALSE) +
    stat_summary(fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x)), geom="errorbar", aes(color=Sample), width=0.5, lwd=0.6, show.legend=FALSE) +
    facet_wrap(~Sample, nrow=1, scales="free") +
    
    # P-values
    stat_compare_means(comparisons = comparisons, label = "p.signif", method="wilcox.test") +
    
    # N number text
    geom_text(data=n_stats, aes(x=2, y=5.2, label=paste0("n=", N)), inherit.aes=FALSE) +
    
    # Median text
    stat_summary(fun.data = function(y) {
      m <- median(y)
      data.frame(y=4.5, label=sprintf("med=%.1f", 10^m))
    }, geom="text", size=3) +
    
    scale_color_manual(values = sample_colors) +
    scale_y_continuous(limits=c(0, 5.5), breaks=0:5) + 
    labs(y = "log10(# of complexes)", x = NULL) +
    theme_classic(base_size=15) +
    theme(
      axis.text = element_text(color="black"),
      strip.background = element_blank(),
      strip.text = element_text(face="bold")
    )
  
  ggsave(OUTPUT_PDF, p, width=12, height=6)
  message(sprintf(">>> Plot saved to: %s", OUTPUT_PDF))
  
} else {
  message("No data available for plotting.")
}