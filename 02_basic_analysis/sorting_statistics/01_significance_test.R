#!/usr/bin/env Rscript
# ==============================================================================
# Script: analyze_loading_site_distribution_v2.R
# Purpose: Classify complexes for each loading site within convergent loops
# Date: 2026-02-03
# 
# Input: Pre-processed file with loop-loading site mappings
#   Format: chr start end type
#   Type: S (left anchor), E (right anchor), M (loading site)
#   ID format: cr1_S, cr1_E, cr1_M-1, cr1_M-2, ...
# ==============================================================================


# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = dirname(dirname(normalizePath(sys.frame(1)$ofile))))
DATA_ROOT    <- Sys.getenv("DATA_ROOT", unset = dirname(PROJECT_ROOT))
GENOME_ROOT  <- Sys.getenv("GENOME_ROOT", unset = file.path(dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

select <- dplyr::select
filter <- dplyr::filter
summarise <- dplyr::summarise
mutate <- dplyr::mutate

# ==============================================================================
# 1. Configuration
# ==============================================================================

LOOP_LOADING_FILE <- paste0(DATA_ROOT, "/ChIA-Drop/3.anchoring_middle/GM12878-cohesin-specific-regions_20260122_filt_25-75percentile_20260122.bed")

FRAGMENT_FILES <- list(
  "CTCF" = paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between/GSE158897_GM12878-CTCF-pooled_comp.rds"),
  "Cohesin" = paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between/GSE158897_GM12878-cohesin-pooled_comp.rds"),
  "RNAPII" = paste0(DATA_ROOT, "/ChIA-Drop/GSE158897/GSE158897_GM12878-RNAPII-pooledv2_comp.rds")
)

OUTPUT_PREFIX <- "loading_site_distribution_v2"

# ==============================================================================
# 2. Load and Parse Loop-Loading Mappings
# ==============================================================================

cat("Loading loop-loading mappings...\n")

regions <- read.table(LOOP_LOADING_FILE, header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("chr", "start", "end", "id"))

regions <- regions %>%
  mutate(
    loop_id = sub("_.*", "", id),
    type = sub(".*_", "", id),
    region_mid = (start + end) / 2
  )

cat(sprintf("Total regions: %d\n", nrow(regions)))
cat(sprintf("  Anchors (S): %d\n", sum(regions$type == "S")))
cat(sprintf("  Anchors (E): %d\n", sum(regions$type == "E")))
cat(sprintf("  Loading (M-*): %d\n", sum(grepl("^M-", regions$type))))

anchors <- regions %>% filter(type %in% c("S", "E"))
loading_sites <- regions %>% filter(grepl("^M-", type))

cat(sprintf("\nUnique loops: %d\n", n_distinct(regions$loop_id)))
cat(sprintf("Loading sites: %d\n", nrow(loading_sites)))

loop_info <- regions %>%
  group_by(loop_id, chr) %>%
  summarise(
    anchor_left_start = min(start[type == "S"]),
    anchor_left_end = max(end[type == "S"]),
    anchor_right_start = min(start[type == "E"]),
    anchor_right_end = max(end[type == "E"]),
    n_loading = sum(grepl("^M-", type)),
    .groups = "drop"
  )

cat(sprintf("Loops with complete anchors: %d\n", nrow(loop_info)))

# ==============================================================================
# 3. Function: Classify complexes for a single loading site
# ==============================================================================

classify_loading_site <- function(loading_row, loop_row, fragments_gr) {
  chr <- loading_row$chr
  loading_start <- loading_row$start
  loading_end <- loading_row$end
  
  anchor_left_start <- loop_row$anchor_left_start
  anchor_left_end <- loop_row$anchor_left_end
  anchor_right_start <- loop_row$anchor_right_start
  anchor_right_end <- loop_row$anchor_right_end
  
  loop_start <- anchor_left_start
  loop_end <- anchor_right_end
  
  loop_gr <- GRanges(seqnames = chr, ranges = IRanges(start = loop_start, end = loop_end))
  frags_in_loop <- subsetByOverlaps(fragments_gr[seqnames(fragments_gr) == chr], loop_gr)
  
  if (length(frags_in_loop) == 0) {
    return(data.frame(anchor = 0, one_sided = 0, two_sided = 0))
  }
  
  frags_df <- as.data.frame(frags_in_loop)
  
  if (!"cluster_id" %in% colnames(frags_df)) {
    if ("CBMB" %in% colnames(frags_df)) {
      frags_df$cluster_id <- frags_df$CBMB
    } else {
      frags_df$cluster_id <- seq_len(nrow(frags_df))
    }
  }
  
  frags_df$frag_mid <- (frags_df$start + frags_df$end) / 2
  
  anchor_left_gr <- GRanges(seqnames = chr, ranges = IRanges(start = anchor_left_start, end = anchor_left_end))
  anchor_right_gr <- GRanges(seqnames = chr, ranges = IRanges(start = anchor_right_start, end = anchor_right_end))
  loading_gr <- GRanges(seqnames = chr, ranges = IRanges(start = loading_start, end = loading_end))
  
  frags_df <- frags_df %>%
    mutate(
      overlaps_anchor_left = overlapsAny(
        GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end)),
        anchor_left_gr
      ),
      overlaps_anchor_right = overlapsAny(
        GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end)),
        anchor_right_gr
      ),
      overlaps_loading = overlapsAny(
        GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end)),
        loading_gr
      ),
      is_left_of_loading = frag_mid < loading_start,
      is_right_of_loading = frag_mid > loading_end,
      is_between_loading_left = frag_mid < loading_start & frag_mid > anchor_left_end,
      is_between_loading_right = frag_mid > loading_end & frag_mid < anchor_right_start
    )
  
  complex_class <- frags_df %>%
    group_by(cluster_id) %>%
    summarise(
      has_anchor_left = any(overlaps_anchor_left),
      has_anchor_right = any(overlaps_anchor_right),
      has_loading_overlap = any(overlaps_loading),
      n_left_of_loading = sum(is_left_of_loading),
      n_right_of_loading = sum(is_right_of_loading),
      n_between_loading_left = sum(is_between_loading_left),
      n_between_loading_right = sum(is_between_loading_right),
      span_start = min(frag_mid),
      span_end = max(frag_mid),
      extends_past_left = any(frag_mid < anchor_left_start),
      extends_past_right = any(frag_mid > anchor_right_end),
      .groups = "drop"
    ) %>%
    mutate(
      category = case_when(
        has_anchor_left | has_anchor_right ~ "anchor",
        has_loading_overlap & (n_between_loading_left >= 1 | n_between_loading_right >= 1) ~ "one_sided",
        !has_loading_overlap & n_left_of_loading >= 1 & n_right_of_loading >= 1 &
          !extends_past_left & !extends_past_right ~ "two_sided_candidate",
        TRUE ~ "other"
      )
    )
  
  two_sided_candidates <- complex_class %>%
    filter(category == "two_sided_candidate") %>%
    arrange(span_end - span_start)
  
  two_sided_ids <- c()
  if (nrow(two_sided_candidates) > 0) {
    prev_left <- Inf
    prev_right <- -Inf
    
    for (i in seq_len(nrow(two_sided_candidates))) {
      current <- two_sided_candidates[i, ]
      if (i == 1) {
        two_sided_ids <- c(two_sided_ids, current$cluster_id)
        prev_left <- current$span_start
        prev_right <- current$span_end
      } else {
        if (current$span_start < prev_left && current$span_end > prev_right) {
          two_sided_ids <- c(two_sided_ids, current$cluster_id)
          prev_left <- current$span_start
          prev_right <- current$span_end
        }
      }
    }
  }
  
  complex_class <- complex_class %>%
    mutate(
      category = case_when(
        cluster_id %in% two_sided_ids ~ "two_sided",
        category == "two_sided_candidate" ~ "other",
        TRUE ~ category
      )
    )
  
  counts <- complex_class %>%
    filter(category %in% c("anchor", "one_sided", "two_sided")) %>%
    count(category) %>%
    pivot_wider(names_from = category, values_from = n, values_fill = 0)
  
  data.frame(
    anchor = if ("anchor" %in% names(counts)) counts$anchor else 0,
    one_sided = if ("one_sided" %in% names(counts)) counts$one_sided else 0,
    two_sided = if ("two_sided" %in% names(counts)) counts$two_sided else 0
  )
}

# ==============================================================================
# 4. Process each dataset
# ==============================================================================

cat(strrep("=", 78), "\n")
cat("Processing datasets\n")
cat(strrep("=", 78), "\n")

results_list <- list()

for (dataset_name in names(FRAGMENT_FILES)) {
  cat(sprintf("\n[%s] %s\n", Sys.time(), dataset_name))
  
  frag_file <- FRAGMENT_FILES[[dataset_name]]
  if (!file.exists(frag_file)) {
    cat(sprintf("  WARNING: File not found\n"))
    next
  }
  
  fragments <- readRDS(frag_file)
  
  if (!inherits(fragments, "GRanges")) {
    fragments_gr <- makeGRangesFromDataFrame(fragments, keep.extra.columns = TRUE)
  } else {
    fragments_gr <- fragments
  }
  
  cat(sprintf("  Fragments: %s\n", format(length(fragments_gr), big.mark = ",")))
  cat("  Classifying...\n")
  
  loading_results <- loading_sites %>%
    left_join(loop_info, by = c("loop_id", "chr")) %>%
    rowwise() %>%
    mutate(
      classification = list(classify_loading_site(
        data.frame(chr = chr, start = start, end = end),
        data.frame(anchor_left_start = anchor_left_start, 
                   anchor_left_end = anchor_left_end,
                   anchor_right_start = anchor_right_start, 
                   anchor_right_end = anchor_right_end),
        fragments_gr
      ))
    ) %>%
    ungroup() %>%
    unnest_wider(classification)
  
  loading_results$dataset <- dataset_name
  results_list[[dataset_name]] <- loading_results
  
  cat(sprintf("  Anchor: %d, One-sided: %d, Two-sided: %d\n",
              sum(loading_results$anchor),
              sum(loading_results$one_sided),
              sum(loading_results$two_sided)))
}

# ==============================================================================
# 5. Save results
# ==============================================================================

cat(strrep("=", 78), "\n")
cat("Saving results\n")
cat(strrep("=", 78), "\n")

all_results <- bind_rows(results_list)

output_file <- paste0(OUTPUT_PREFIX, "_detailed.tsv")
write.table(all_results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[Output] %s\n", output_file))

summary_stats <- all_results %>%
  group_by(dataset) %>%
  summarise(
    n_loading_sites = n(),
    n_unique_loops = n_distinct(loop_id),
    total_anchor = sum(anchor),
    total_one_sided = sum(one_sided),
    total_two_sided = sum(two_sided),
    median_anchor = median(anchor),
    median_one_sided = median(one_sided),
    median_two_sided = median(two_sided)
  )

output_summary <- paste0(OUTPUT_PREFIX, "_summary.tsv")
write.table(summary_stats, output_summary, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[Output] %s\n", output_summary))

cat(strrep("=", 78), "\n")
cat("SUMMARY\n")
cat(strrep("=", 78), "\n")
print(summary_stats, n = Inf)

# ==============================================================================
# 6. Visualization
# ==============================================================================

cat(strrep("=", 78), "\n")
cat("Visualization\n")
cat(strrep("=", 78), "\n")

plot_data <- all_results %>%
  select(dataset, id, loop_id, anchor, one_sided, two_sided) %>%
  pivot_longer(cols = c(anchor, one_sided, two_sided),
               names_to = "category", values_to = "count") %>%
  filter(count > 0) %>%
  mutate(log10_count = log10(count))

plot_data$category <- factor(
  plot_data$category,
  levels = c("two_sided", "one_sided", "anchor"),
  labels = c("Two-sided", "One-sided", "Anchor")
)

plot_data$dataset <- factor(
  plot_data$dataset,
  levels = c("CTCF", "Cohesin", "RNAPII")
)

stat_results <- data.frame()

for (ds in levels(plot_data$dataset)) {
  ds_data <- plot_data %>% filter(dataset == ds)
  if (nrow(ds_data) == 0) next
  
  categories <- unique(ds_data$category)
  if (length(categories) >= 2) {
    pairs <- combn(categories, 2, simplify = FALSE)
    for (pair in pairs) {
      g1 <- ds_data %>% filter(category == pair[1]) %>% pull(log10_count)
      g2 <- ds_data %>% filter(category == pair[2]) %>% pull(log10_count)
      
      if (length(g1) > 0 && length(g2) > 0) {
        wt <- wilcox.test(g1, g2)
        stat_results <- rbind(stat_results, data.frame(
          dataset = ds,
          comparison = paste(pair[1], "vs", pair[2]),
          p_value = wt$p.value,
          statistic = as.numeric(wt$statistic)
        ))
      }
    }
  }
}

stat_file <- paste0(OUTPUT_PREFIX, "_statistics.tsv")
write.table(stat_results, stat_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[Output] %s\n", stat_file))

p <- ggplot(plot_data, aes(x = category, y = log10_count)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1, width = 0.6) +
  facet_wrap(~dataset, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "log10(# of complexes)") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

pdf_file <- paste0(OUTPUT_PREFIX, "_boxplot.pdf")
ggsave(pdf_file, p, width = 10, height = 4)
cat(sprintf("[Output] %s\n", pdf_file))

cat(strrep("=", 78), "\n")
cat("Complete!\n")
cat(strrep("=", 78), "\n")