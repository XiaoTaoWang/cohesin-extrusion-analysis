# ==============================================================================
# Script Name: batch_loop_stats_MEMORY_SAFE.R
# Function:    Memory-efficient batch processing with chunking
# Date:        2026-02-15
# ==============================================================================


# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = dirname(dirname(normalizePath(sys.frame(1)$ofile))))
DATA_ROOT    <- Sys.getenv("DATA_ROOT", unset = dirname(PROJECT_ROOT))
GENOME_ROOT  <- Sys.getenv("GENOME_ROOT", unset = file.path(dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(ggplot2)
  library(scales)
})

# ==============================================================================
# 1. Configuration
# ==============================================================================

DIR_MAIN  <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between")
DIR_OTHER <- paste0(DATA_ROOT, "/ChIA-Drop/GSE158897")
MOTIF_FILE <- file.path(DIR_MAIN, "CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_id_v2.bed")
LOOP_BED_FILE <- paste0(DATA_ROOT, "/ChIA-Drop/3.anchoring_middle/GM12878-cohesin-specific-regions_20260122.bed")

RDS_FILES <- list(
  "CTCF" = file.path(DIR_MAIN, "GSE158897_GM12878-CTCF-pooled_comp.rds"),
  "Cohesin" = file.path(DIR_MAIN, "GSE158897_GM12878-cohesin-pooled_comp.rds"),
  "RNAPII" = paste0(DATA_ROOT, "/ChIA-Drop/GSE158897/GSE158897_GM12878-RNAPII-pooledv2_comp.rds")
)

MIN_FRAGS_IN_ROI   <- 2
EDGE_BUFFER        <- 500
EXTEND_BP          <- 50000

# === Memory Control ===
CHUNK_SIZE <- 100  # process 100 loops per chunk (tunable)

# ==============================================================================
# 2. Load Motifs (Small, keep in memory)
# ==============================================================================
message(">>> Loading CTCF motifs...")
motifs_raw <- read.table(MOTIF_FILE, stringsAsFactors = FALSE)
colnames(motifs_raw)[1:5] <- c("seqnames","start","end","strand","motif_id")

# ==============================================================================
# 3. Parse Loop BED File
# ==============================================================================
message(">>> Parsing loop BED file...")
bed_raw <- read.table(LOOP_BED_FILE, header = FALSE, stringsAsFactors = FALSE)
colnames(bed_raw)[1:4] <- c("chr", "start", "end", "label")

loops_list <- bed_raw %>%
  mutate(
    loop_id = sub("_[SE]$", "", label),
    anchor_type = ifelse(grepl("_S$", label), "left", 
                  ifelse(grepl("_E$", label), "right", "loading"))
  ) %>%
  filter(anchor_type %in% c("left", "right")) %>%
  select(chr, start, end, loop_id, anchor_type)

loops_df <- loops_list %>%
  pivot_wider(
    id_cols = c(chr, loop_id),
    names_from = anchor_type,
    values_from = c(start, end),
    names_sep = "_"
  ) %>%
  filter(!is.na(start_left) & !is.na(start_right)) %>%
  mutate(
    region_start = start_left - EXTEND_BP,
    region_end = end_right + EXTEND_BP,
    anchor_left_start = start_left,
    anchor_left_end = end_left,
    anchor_right_start = start_right,
    anchor_right_end = end_right
  )

message(sprintf("    Found %d loops", nrow(loops_df)))

# ==============================================================================
# 4. Simple Analysis Functions (No fancy tricks, just solid code)
# ==============================================================================

process_single_loop <- function(loop, df_raw, motifs_raw) {
  # Extract region
  roi_gr <- GRanges(seqnames = loop$chr, 
                    ranges = IRanges(start = loop$region_start, end = loop$region_end))
  
  if (inherits(df_raw, "GRanges")) {
    df <- as.data.frame(subsetByOverlaps(df_raw, roi_gr))
  } else {
    df <- df_raw %>% 
      filter(seqnames == loop$chr, end >= loop$region_start, start <= loop$region_end)
  }
  
  if (nrow(df) == 0) return(NULL)
  
  # Add cluster_id
  if (!"cluster_id" %in% colnames(df)) {
    if ("CBMB" %in% colnames(df)) df$cluster_id <- df$CBMB 
    else df$cluster_id <- 1:nrow(df)
  }
  
  df$frag_mid <- (df$start + df$end) / 2
  df$display_strand <- NA
  
  # Motif overlap
  target_chr <- loop$chr
  motifs_sub <- motifs_raw
  if (grepl("chr", target_chr)) {
    if (!any(grepl("chr", motifs_sub$seqnames))) 
      motifs_sub$seqnames <- paste0("chr", motifs_sub$seqnames)
  } else {
    motifs_sub$seqnames <- gsub("chr", "", motifs_sub$seqnames)
  }
  
  motifs_roi <- motifs_sub %>% 
    filter(seqnames == target_chr, end >= loop$region_start, start <= loop$region_end)
  
  if (nrow(motifs_roi) > 0) {
    motifs_gr <- makeGRangesFromDataFrame(motifs_roi, keep.extra.columns = TRUE)
    frag_gr <- GRanges(seqnames = df$seqnames, 
                       ranges = IRanges(start = df$start, end = df$end))
    ov <- findOverlaps(frag_gr, motifs_gr)
    
    if (length(ov) > 0) {
      ov_df <- data.frame(
        f_idx = queryHits(ov),
        m_str = as.character(strand(motifs_gr[subjectHits(ov)]))
      )
      strand_summ <- ov_df %>% 
        group_by(f_idx) %>%
        summarise(has_pos = "+" %in% m_str, has_neg = "-" %in% m_str, .groups = "drop") %>%
        mutate(final_str = case_when(
          has_pos & !has_neg ~ "+",
          !has_pos & has_neg ~ "-",
          TRUE ~ "Conflict"
        ))
      df$display_strand[strand_summ$f_idx] <- strand_summ$final_str
    }
  }
  
  # Filter valid clusters
  df$in_display <- df$frag_mid >= loop$region_start & df$frag_mid <= loop$region_end
  
  valid_ids <- df %>% 
    filter(in_display) %>% 
    group_by(cluster_id) %>% 
    summarize(n = n(), .groups = "drop") %>% 
    filter(n >= MIN_FRAGS_IN_ROI) %>% 
    pull(cluster_id)
  
  df <- df %>% filter(cluster_id %in% valid_ids)
  if (nrow(df) == 0) return(NULL)
  
  # Boundary filter
  bound_ids <- df %>%
    group_by(cluster_id) %>%
    summarize(l = min(frag_mid), r = max(frag_mid), .groups = "drop") %>%
    filter(l >= (loop$region_start + EDGE_BUFFER), 
           r <= (loop$region_end - EDGE_BUFFER)) %>%
    pull(cluster_id)
  
  df <- df %>% filter(cluster_id %in% bound_ids)
  if (nrow(df) == 0) return(NULL)
  
  # Calculate stats
  stats <- df %>%
    group_by(cluster_id) %>%
    summarise(
      c_start = min(frag_mid),
      c_end = max(frag_mid),
      is_anchored_left = any(!is.na(display_strand) & display_strand == "+" & 
                             frag_mid >= loop$anchor_left_start & 
                             frag_mid <= loop$anchor_left_end),
      is_anchored_right = any(!is.na(display_strand) & display_strand == "-" & 
                              frag_mid >= loop$anchor_right_start & 
                              frag_mid <= loop$anchor_right_end),
      has_plus = any(!is.na(display_strand) & display_strand == "+"),
      has_minus = any(!is.na(display_strand) & display_strand == "-"),
      .groups = "drop"
    ) %>%
    filter(
      (is_anchored_left & c_end > loop$anchor_left_end & c_end <= loop$anchor_right_end) |
      (is_anchored_right & c_start >= loop$anchor_left_start & c_start < loop$anchor_right_start)
    )
  
  if (nrow(stats) == 0) return(NULL)
  
  n_forward <- sum(stats$is_anchored_left & 
                   stats$c_end > loop$anchor_left_end & 
                   stats$c_end <= loop$anchor_right_end)
  n_reverse <- sum(stats$is_anchored_right & 
                   stats$c_start >= loop$anchor_left_start & 
                   stats$c_start < loop$anchor_right_start)
  n_complete <- sum(stats$has_plus & stats$has_minus)
  n_ongoing <- sum((stats$has_plus & !stats$has_minus) | 
                   (!stats$has_plus & stats$has_minus))
  n_total <- n_forward + n_reverse
  
  if (n_total == 0) return(NULL)
  
  data.frame(
    loop_id = loop$loop_id,
    chr = loop$chr,
    forward_n = n_forward,
    reverse_n = n_reverse,
    complete_n = n_complete,
    ongoing_n = n_ongoing,
    total_n = n_total,
    forward_pct = n_forward / n_total * 100,
    reverse_pct = n_reverse / n_total * 100,
    complete_pct = n_complete / n_total * 100,
    ongoing_pct = n_ongoing / n_total * 100
  )
}

# ==============================================================================
# 5. Chunked Processing (Memory-Safe)
# ==============================================================================
message(sprintf("\n>>> Processing loops in chunks of %d...\n", CHUNK_SIZE))

all_results <- list()

for (sample_name in names(RDS_FILES)) {
  message(sprintf("=== Processing %s ===", sample_name))
  
  # Load data ONCE per sample
  message("    Loading RDS...")
  df_raw <- readRDS(RDS_FILES[[sample_name]])
  
  # Process in chunks
  n_loops <- nrow(loops_df)
  n_chunks <- ceiling(n_loops / CHUNK_SIZE)
  
  sample_results <- list()
  
  for (chunk_i in 1:n_chunks) {
    start_idx <- (chunk_i - 1) * CHUNK_SIZE + 1
    end_idx <- min(chunk_i * CHUNK_SIZE, n_loops)
    
    message(sprintf("  Chunk %d/%d (loops %d-%d)", chunk_i, n_chunks, start_idx, end_idx))
    
    chunk_results <- lapply(start_idx:end_idx, function(i) {
      if (i %% 50 == 0) cat(sprintf("    Loop %d/%d\n", i, n_loops))
      
      tryCatch({
        process_single_loop(loops_df[i, ], df_raw, motifs_raw)
      }, error = function(e) NULL)
    })
    
    sample_results <- c(sample_results, chunk_results[!sapply(chunk_results, is.null)])
    
    # Force garbage collection after each chunk
    gc()
  }
  
  if (length(sample_results) > 0) {
    all_results[[sample_name]] <- bind_rows(sample_results) %>%
      mutate(sample = sample_name)
  }
  
  # Clear large object
  rm(df_raw)
  gc()
}

# Combine all results
all_results_df <- bind_rows(all_results)

message(sprintf("\n>>> Analysis complete! %d valid loops", nrow(all_results_df)))

# ==============================================================================
# 6. Save Results
# ==============================================================================
output_file <- "loop_stats_summary.tsv"
write.table(all_results_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf(">>> Saved to: %s", output_file))

# ==============================================================================
# 7. Generate Boxplot
# ==============================================================================
message("\n>>> Generating boxplot...")

plot_data <- all_results_df %>%
  select(sample, loop_id, forward_pct, reverse_pct, complete_pct) %>%
  pivot_longer(cols = c(forward_pct, reverse_pct, complete_pct),
               names_to = "category",
               values_to = "percentage") %>%
  mutate(
    category = factor(category, 
                      levels = c("forward_pct", "reverse_pct", "complete_pct"),
                      labels = c("Forward", "Reverse", "Complete"))
  )

category_colors <- c(
  "Forward" = "#EF4444",
  "Reverse" = "#A78BFA",
  "Complete" = "#10B981"
)

p_box <- ggplot(plot_data, aes(x = category, y = percentage, fill = category)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 0.5) +
  facet_wrap(~ sample, nrow = 1) +
  scale_fill_manual(values = category_colors) +
  labs(x = NULL, y = "Percentage (%)", 
       title = "Distribution of Forward/Reverse/Complete Interactions Across Loops") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave("loop_stats_boxplot.pdf", p_box, width = 10, height = 4)
print(p_box)
message(">>> Saved boxplot to: loop_stats_boxplot.pdf")

# ==============================================================================
# 8. Summary Statistics
# ==============================================================================
message("\n>>> Summary Statistics:")
summary_stats <- plot_data %>%
  group_by(sample, category) %>%
  summarise(
    n_loops = n(),
    median = median(percentage, na.rm = TRUE),
    mean = mean(percentage, na.rm = TRUE),
    q25 = quantile(percentage, 0.25, na.rm = TRUE),
    q75 = quantile(percentage, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

message("\n>>> All done! ✅")