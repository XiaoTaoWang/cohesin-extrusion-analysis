# ==============================================================================
# Script Name: batch_final_parallel_v2.R
# Function:    Batch Generate ChIA-Drop Plots (Multi-threaded)
#              1. Reads loop definitions from CSVs.
#              2. Iterates through loops in parallel (6 cores).
#              3. Generates plots for each loop (CTCF & Cohesin).
# Date:        2026-02-01
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
  library(scales)
  library(GenomicRanges)
  library(patchwork)
  library(rtracklayer)
  library(parallel) 
})

slice <- dplyr::slice
select <- dplyr::select

# ==============================================================================
# 1. Load Loop Data
# ==============================================================================
csv_files <- c(
  paste0(DATA_ROOT, "/ChIA-Drop/4.correlation_NIPBL/NIPBL_vs_Cohesin_bottom_right.csv"),
  paste0(DATA_ROOT, "/ChIA-Drop/4.correlation_NIPBL/NIPBL_vs_Cohesin_top_left.csv"),
  paste0(DATA_ROOT, "/ChIA-Drop/4.correlation_NIPBL/NIPBL_vs_CTCF_bottom_right.csv"),
  paste0(DATA_ROOT, "/ChIA-Drop/4.correlation_NIPBL/NIPBL_vs_CTCF_top_left.csv")
)

# Read and combine
loop_data <- do.call(rbind, lapply(csv_files, function(f) {
  if(file.exists(f)) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    df$category_group <- tools::file_path_sans_ext(basename(f))
    return(df)
  } else NULL
}))

if (is.null(loop_data) || nrow(loop_data) == 0) {
  stop("No loop data found in CSV files.")
}

# Deduplicate
plot_queue <- loop_data %>%
  select(loop_id, chr1, start1, end1, chr2, start2, end2, category_group) %>%
  distinct()

message(sprintf(">>> Found %d unique loops to process.", nrow(plot_queue)))

# ==============================================================================
# 2. Global Configuration & Constants
# ==============================================================================

# === Visual Parameters ===
GAP_SIZE_PERCENT   <- 0.01        
MIN_FRAGS_IN_ROI   <- 2             
EDGE_BUFFER        <- 500           

FRAGMENT_SIZE      <- 1.5           
MOTIF_SHAPE_CODE   <- 15           
no_motif_size      <- 1.0           

MOTIF_COLOR_P      <- "#EF4444"           
MOTIF_COLOR_N      <- "#A78BFA"           
COLOR_LINE         <- "#acada2"
LINE_WIDTH         <- 0.1
COLOR_NIPBL        <- "#E69F00" 
LOADING_COLOR      <- "#ffd04edd"

# === Rectangle Config ===
RECT_WIDTH_BP  <- 3000    
RECT_HEIGHT    <- 2.5     

# === File Paths ===
DIR_MAIN  <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between")
DIR_OTHER <- paste0(DATA_ROOT, "/ChIA-Drop/GSE158897")
MOTIF_FILE <- file.path(DIR_MAIN, "CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_id_v2.bed")
OUT_DIR    <- file.path(DIR_MAIN, "18.batch_plots")

if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# === Tracks Info ===
TRACKS_INFO <- list(
  "CTCF_ChIA-PET"    = c(file.path(DIR_MAIN, "GM12878_CTCF_ChIA-PET_GSE158897_LHG0052H.for.BROWSER.sorted.bw"), "#0000FF"),
  "Cohesin_ChIA-PET" = c(file.path(DIR_MAIN, "GM12878_cohesin_ChIA-PET_GSE158897_LHG0051H_0104V.for.BROWSER.sorted.bw"), "#009E73"),
  "NIPBL"   = c(file.path(DIR_OTHER, "GM12878_NIPBL_chip-seq_GSE158897_CHG0030.q30.nr.sorted.bw"), COLOR_NIPBL), 
  "WAPL"    = c(file.path(DIR_OTHER, "GM12878_WAPL_chip-seq_GSE158897_CHG0032.q30.nr.sorted.bw"), "#A52A2A"),
  "RNAPII"  = c(file.path(DIR_MAIN, "GM12878_RNAPII_ChIA-PET_GSE158897_LHG0035N_0035V_0045V.for.BROWSER.sorted.bw"), "#800080"),
  "H3K27ac" = c(file.path(DIR_OTHER, "GM12878_h3k27ac_chip-seq_ENCFF340JIF.bigWig"), "#D55E00"),
  "H3K4me1" = c(file.path(DIR_OTHER, "GM12878_h3k4me1_chip-seq_ENCFF831ZHL.bigWig"), "#8B0000")
)

# ==============================================================================
# 3. Caching Setup
# ==============================================================================
if (!exists("G_CACHE_CTCF")) { G_CACHE_CTCF <- NULL }
if (!exists("G_CACHE_COHESIN")) { G_CACHE_COHESIN <- NULL }

BATCH_CONFIG <- list(
  "CTCF" = list(
    rds_file = file.path(DIR_MAIN, "GSE158897_GM12878-CTCF-pooled_comp.rds"),
    color_dot = "#0000f3",
    cache_name = "G_CACHE_CTCF"
  ),
  "Cohesin" = list(
    rds_file = file.path(DIR_MAIN, "GSE158897_GM12878-cohesin-pooled_comp.rds"),
    color_dot = "#008000",
    cache_name = "G_CACHE_COHESIN"
  )
)

# ==============================================================================
# 4. Helpers & Plotting Functions
# ==============================================================================

clean_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(color = "black", size = 10, margin = margin(t = 5)),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_blank(),           
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(color = "black", size = 9, hjust = 1, margin = margin(r = 5), face = "bold"),
    axis.title.y = element_blank(),         
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none",       
    plot.margin = margin(t = 2, b = 2, r = 20, l = 10)   
  )

round_smart <- function(x) {
  if (is.na(x) || x == 0) return(0)
  if (x <= 1) return(round(x, 2))
  if (x <= 10) return(ceiling(x * 10) / 10)     
  if (x <= 50) return(ceiling(x))                 
  return(ceiling(x / 10) * 10)                   
}

plot_motif_track <- function(bed_file, chrom, start_pos, end_pos) {
  base_plot <- ggplot() +       
    annotate("rect", xmin = LOADING_START, xmax = LOADING_END, ymin = -Inf, ymax = Inf, fill = LOADING_COLOR, alpha = 0.5) +
    coord_cartesian(xlim = c(start_pos, end_pos), ylim = c(0.8, 1.2), expand = FALSE, clip = "off") +
    scale_x_continuous(limits = c(start_pos, end_pos), expand = c(0, 0)) +
    labs(y = "CTCF Motif") + clean_theme +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), plot.margin = margin(t = 5, b = 0, r = 20, l = 10))
   
  if (!file.exists(bed_file)) return(base_plot)
  df <- tryCatch({ read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)[,1:4] }, error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(base_plot)
  colnames(df) <- c("chr", "start", "end", "strand")
   
  if (grepl("chr", chrom)) { if (!any(grepl("chr", df$chr))) df$chr <- paste0("chr", df$chr) } 
  else { df$chr <- gsub("chr", "", df$chr) }
   
  df <- df %>% filter(chr == chrom, end >= start_pos, start <= end_pos)
  df$mid <- (df$start + df$end) / 2
   
  base_plot +
    geom_segment(data = df, aes(x = ifelse(strand=="-", mid + 200, mid - 200),
                                xend = ifelse(strand=="-", mid - 200, mid + 200),
                                y = 1, yend = 1, color = strand),
                 arrow = arrow(length = unit(0.20, "cm"), type = "closed"), linewidth = 0.5) +
    scale_color_manual(values = c("+" = MOTIF_COLOR_P, "-" = MOTIF_COLOR_N))
}

plot_signal_track <- function(bw_file, track_name, fill_color, chrom, start_pos, end_pos) {
  p_void <- ggplot() + 
    annotate("rect", xmin = LOADING_START, xmax = LOADING_END, ymin = -Inf, ymax = Inf, fill = LOADING_COLOR, alpha = 0.5) +
    clean_theme + labs(y = track_name) + coord_cartesian(xlim = c(start_pos, end_pos), expand = FALSE) +
    scale_x_continuous(limits = c(start_pos, end_pos), expand = c(0, 0)) +
    theme(axis.line.x=element_blank(), axis.text.y = element_blank(), axis.ticks.x=element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 9, face = "bold", margin = margin(r = 5)),
          plot.margin = margin(t = 2, b = 2, r = 20, l = 10))
   
  if (!file.exists(bw_file)) return(p_void)
  bw_obj <- BigWigFile(bw_file)
  roi_gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start_pos, end = end_pos))
   
  binned_list <- tryCatch({ summary(bw_obj, which = roi_gr, size = 1000, type = "mean") }, error = function(e) NULL)
  if (is.null(binned_list) || length(binned_list) == 0) return(p_void)
  binned_gr <- binned_list[[1]]
  if (length(binned_gr) == 0) return(p_void)
   
  plot_df <- data.frame(pos = (start(binned_gr) + end(binned_gr)) / 2, score = binned_gr$score)
  plot_df$score[is.na(plot_df$score)] <- 0; plot_df$score[plot_df$score < 0] <- 0
   
  smooth_fit <- ksmooth(plot_df$pos, plot_df$score, kernel = "normal", bandwidth = (end_pos - start_pos) / 800)
  final_df <- data.frame(pos = smooth_fit$x, score = smooth_fit$y)
   
  raw_max <- max(final_df$score, na.rm=TRUE)
  q995 <- quantile(final_df$score, 0.995, na.rm = TRUE)
  effective_max <- if (raw_max > 4 * q995 && q995 > 1) q995 * 1.5 else raw_max * 1.1
  effective_max <- max(effective_max, 1)
  rounded_max <- round_smart(effective_max)
   
  ggplot(final_df, aes(x = pos, y = score)) +
    annotate("rect", xmin = LOADING_START, xmax = LOADING_END, ymin = 0, ymax = rounded_max, fill = LOADING_COLOR, alpha = 0.5) +
    geom_area(fill = fill_color, alpha = 0.9, outline.type = "upper") + geom_line(color = fill_color, linewidth = 0.2) +
    scale_y_continuous(breaks = c(0, rounded_max), labels = c("0", rounded_max), expand = c(0, 0)) +
    scale_x_continuous(limits = c(start_pos, end_pos), expand = c(0, 0)) +
    coord_cartesian(xlim = c(start_pos, end_pos), ylim = c(0, rounded_max), expand = FALSE, clip = "on") +
    labs(y = track_name) + clean_theme +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 9, face = "bold", margin = margin(r = 5)),
          panel.grid = element_blank(), panel.border = element_blank(), plot.margin = margin(t = 2, b = 2, r = 20, l = 10))
}

process_fragments_dynamic <- function(df_raw, left_range, right_range) {
  roi_gr <- GRanges(seqnames = CHR_BW, ranges = IRanges(start = START, end = END))
   
  if(inherits(df_raw, "GRanges")) {
    df_gr <- suppressWarnings(subsetByOverlaps(df_raw, roi_gr))
    df <- as.data.frame(df_gr)
  } else if(is.data.frame(df_raw)) {
    df <- df_raw %>% filter(seqnames == as.character(seqnames(roi_gr)), end >= start(roi_gr), start <= end(roi_gr))
  } else { return(NULL) }
   
  if(nrow(df) == 0) return(NULL)
   
  if (!"cluster_id" %in% colnames(df)) {
    if("CBMB" %in% colnames(df)) df$cluster_id <- df$CBMB else df$cluster_id <- 1:nrow(df)
  }
   
  df$frag_mid <- (df$start + df$end) / 2
  df$display_strand <- NA
   
  if (file.exists(MOTIF_FILE)) {
    motifs_raw <- read.table(MOTIF_FILE, stringsAsFactors = FALSE)
    colnames(motifs_raw)[1:5] <- c("seqnames","start","end","strand","motif_id")
    target_chr <- unique(df$seqnames)[1]
    if (grepl("chr", target_chr)) { if(!any(grepl("chr", motifs_raw$seqnames))) motifs_raw$seqnames <- paste0("chr", motifs_raw$seqnames) }
    else { motifs_raw$seqnames <- gsub("chr", "", motifs_raw$seqnames) }
     
    motifs_in_region <- motifs_raw %>% filter(seqnames == target_chr, end >= START, start <= END)
     
    if (nrow(motifs_in_region) > 0) {
      motifs_gr <- makeGRangesFromDataFrame(motifs_in_region, keep.extra.columns = TRUE)
      frag_gr <- GRanges(seqnames = df$seqnames, ranges = IRanges(start = df$start, end = df$end))
      ov <- suppressWarnings(findOverlaps(frag_gr, motifs_gr))
       
      if (length(ov) > 0) {
        ov_df <- data.frame(f_idx = queryHits(ov), m_str = as.character(strand(motifs_gr[subjectHits(ov)])))
        strand_summ <- ov_df %>% group_by(f_idx) %>%
          summarise(has_pos = "+" %in% m_str, has_neg = "-" %in% m_str) %>%
          mutate(final_str = case_when(has_pos & !has_neg ~ "+", !has_pos & has_neg ~ "-", TRUE ~ "Conflict"))
        df$display_strand[strand_summ$f_idx] <- strand_summ$final_str
      }
    }
  }

  if (!is.null(left_range) || !is.null(right_range)) {
    if (!is.null(left_range)) {
      is_left <- df$frag_mid >= left_range[1] & df$frag_mid <= left_range[2]
      df$display_strand[is_left] <- "+"
    }
    if (!is.null(right_range)) {
      is_right <- df$frag_mid >= right_range[1] & df$frag_mid <= right_range[2]
      df$display_strand[is_right] <- "-"
    }
  }
   
  df$in_display_region <- df$frag_mid >= START & df$frag_mid <= END
  df$in_loading <- df$frag_mid >= LOADING_START & df$frag_mid <= LOADING_END
  return(df)
}

classify_clusters_dynamic <- function(df, left_range, right_range) {
  if(is.null(df) || nrow(df) == 0) return(NULL)
   
  l_start <- START; l_end <- START; r_start <- END; r_end <- END
   
  if(!is.null(left_range)) { l_start <- left_range[1]; l_end <- left_range[2] }
  if(!is.null(right_range)) { r_start <- right_range[1]; r_end <- right_range[2] }

  cluster_stats <- df %>% group_by(cluster_id) %>%
    summarise(
      c_start = min(frag_mid), c_end = max(frag_mid), span = c_end - c_start,
      n_frags_left = sum(frag_mid < LOADING_START), n_frags_right = sum(frag_mid > LOADING_END), n_frags_in_loading = sum(in_loading),
      has_left_anchor = any(!is.na(display_strand) & display_strand == "+" & frag_mid >= l_start & frag_mid <= l_end),
      has_right_anchor = any(!is.na(display_strand) & display_strand == "-" & frag_mid >= r_start & frag_mid <= r_end)
    )
   
  classified <- cluster_stats %>%
    mutate(
      is_forward = has_left_anchor & c_end > l_end & c_end <= r_end,
      is_reverse = has_right_anchor & c_start >= l_start & c_start < r_start
    )
   
  forward <- classified %>% filter(is_forward) %>% mutate(category = "Forward Anchored", sort_key = c_end) %>% select(cluster_id, category, span, sort_key)
  reverse <- classified %>% filter(is_reverse) %>% mutate(category = "Reverse Anchored", sort_key = c_start) %>% select(cluster_id, category, span, sort_key)
   
  anchored_ids <- unique(c(forward$cluster_id, reverse$cluster_id))
  unanchored_stats <- classified %>% filter(!cluster_id %in% anchored_ids) # Correct variable definition
   
  # Two-sided Logic
  twosided_cand <- unanchored_stats %>% filter(n_frags_left>=1, n_frags_right>=1, n_frags_in_loading==0) %>% arrange(span)
  twosided_ids <- c()
  if(nrow(twosided_cand)>0) {
    sel <- data.frame(); 
    for(i in 1:nrow(twosided_cand)) {
      curr <- twosided_cand[i,]
      if(nrow(sel)==0) sel <- rbind(sel, curr)
      else if(curr$c_start < min(sel$c_start) && curr$c_end > max(sel$c_end)) sel <- rbind(sel, curr)
    }
    if(nrow(sel)>0) twosided_ids <- sel$cluster_id
  }
   
  # FIX: Use unanchored_stats here instead of unanchored
  unanchored <- unanchored_stats %>%
    mutate(
      category = case_when(
        cluster_id %in% twosided_ids ~ "Two-sided at loading",
        n_frags_in_loading > 0 & ((n_frags_left>0 & n_frags_right==0)|(n_frags_left==0 & n_frags_right>0)) ~ "One-sided at loading",
        TRUE ~ "Other without anchor"
      ),
      extends_right = n_frags_right > 0,
      sort_key = span
    ) %>% select(cluster_id, category, span, sort_key, extends_right)
   
  bind_rows(forward, reverse, unanchored)
}

assign_layout <- function(df, left_range, right_range) {
  if(is.null(df) || nrow(df) == 0) return(NULL)
   
  valid_ids <- df %>% filter(in_display_region) %>% group_by(cluster_id) %>% summarize(n=n()) %>% filter(n >= MIN_FRAGS_IN_ROI) %>% pull(cluster_id)
  df <- df %>% filter(cluster_id %in% valid_ids)
   
  bound_ids <- df %>% group_by(cluster_id) %>% summarize(l=min(frag_mid), r=max(frag_mid)) %>% 
    filter(l >= (START + EDGE_BUFFER) & r <= (END - EDGE_BUFFER)) %>% pull(cluster_id)
  df <- df %>% filter(cluster_id %in% bound_ids)
   
  cats <- classify_clusters_dynamic(df, left_range, right_range)
  if(is.null(cats) || nrow(cats)==0) return(NULL)
   
  fwd <- cats %>% filter(category == "Forward Anchored") %>% arrange(sort_key, span)
  rev <- cats %>% filter(category == "Reverse Anchored") %>% arrange(desc(sort_key), span)
  one <- cats %>% filter(category == "One-sided at loading") %>% arrange(desc(extends_right), span)
  two <- cats %>% filter(category == "Two-sided at loading") %>% arrange(span)
   
  total_rows <- nrow(fwd) + nrow(rev) + nrow(one) + nrow(two)
  gap <- max(2, round(total_rows * GAP_SIZE_PERCENT))
  cursor <- total_rows + (4 * gap)
   
  make_grp <- function(d, lbl, h) {
    if(nrow(d)==0) return(list(df=data.frame(), lbl=data.frame(), h=h))
    d$y_idx <- seq(h, h - nrow(d) + 1)
    list(df = d %>% select(cluster_id, y_idx, category), lbl = data.frame(label=paste0(lbl, "\nn=", nrow(d)), y=h-nrow(d)/2), h=h-nrow(d))
  }
   
  s1 <- make_grp(fwd, "Forward\nAnchored", cursor)
  s2 <- make_grp(rev, "Reverse\nAnchored", s1$h - gap)
  s3 <- make_grp(one, "One-sided\nat loading", s2$h - gap)
  s4 <- make_grp(two, "Two-sided\nat loading", s3$h - gap)
   
  layout_res <- bind_rows(s1$df, s2$df, s3$df, s4$df)
  lbl_res <- bind_rows(s1$lbl, s2$lbl, s3$lbl, s4$lbl)
  if(nrow(layout_res)==0) return(NULL)
  list(data = df %>% inner_join(layout_res, by="cluster_id", relationship="many-to-many"), labels = lbl_res)
}

plot_final_view <- function(layout_res, current_dot_color) {
  df <- layout_res$data
  labels <- layout_res$labels
   
  df_motif <- df %>% filter(!is.na(display_strand) & display_strand %in% c("+", "-"))
  df_loading_orange <- df %>% filter(in_loading, !category %in% c("Forward Anchored", "Reverse Anchored"))
  ids_special <- c(paste(df_motif$cluster_id, df_motif$frag_mid, df_motif$y_idx),
                   paste(df_loading_orange$cluster_id, df_loading_orange$frag_mid, df_loading_orange$y_idx))
  df_dot <- df %>% filter(!paste(cluster_id, frag_mid, y_idx) %in% ids_special)
   
  p <- ggplot() +
    annotate("rect", xmin = LOADING_START, xmax = LOADING_END, ymin = -Inf, ymax = Inf, fill = LOADING_COLOR, alpha = 0.5) +
    geom_line(data = df, aes(x = frag_mid, y = y_idx, group = interaction(cluster_id, category)), linewidth = LINE_WIDTH, color = COLOR_LINE) +
     
    geom_point(data = df_dot, aes(x = frag_mid, y = y_idx), size = no_motif_size, color = current_dot_color, shape = 19) +
     
    geom_point(data = df_loading_orange, aes(x = frag_mid, y = y_idx), color = COLOR_NIPBL, size = FRAGMENT_SIZE, shape = MOTIF_SHAPE_CODE) +
     
    geom_point(data = df_motif, aes(x = frag_mid, y = y_idx, color = display_strand), size = FRAGMENT_SIZE, shape = MOTIF_SHAPE_CODE) +
    scale_color_manual(values = c("+" = MOTIF_COLOR_P, "-" = MOTIF_COLOR_N)) +
     
    scale_x_continuous(limits = c(START, END), breaks = breaks_extended(n = 5), labels = comma, expand = c(0,0)) +
    labs(y = NULL, x = paste0(CHR_BW, ":", comma(START), "-", comma(END))) +       
    coord_cartesian(xlim = c(START, END), expand = FALSE, clip = "off") +
    clean_theme +
    scale_y_continuous(breaks = labels$y, labels = labels$label, expand = expansion(mult = c(0.05, 0)))
     
  return(p)
}

# ==============================================================================
# 5. Main Execution Loop (Parallel)
# ==============================================================================

# Pre-load data
message(">>> Pre-loading RDS files into memory...")
for (sample_key in names(BATCH_CONFIG)) {
  cfg <- BATCH_CONFIG[[sample_key]]
  if (!exists(cfg$cache_name) || is.null(get(cfg$cache_name))) {
    message(sprintf("    Loading %s...", cfg$rds_file))
    if(file.exists(cfg$rds_file)) {
      assign(cfg$cache_name, readRDS(cfg$rds_file), envir = .GlobalEnv)
    } else {
      warning(sprintf("    File not found: %s", cfg$rds_file))
    }
  }
}

process_single_loop <- function(i) {
  row <- plot_queue[i,]
   
  # Globals in child fork
  CHR_BW <<- row$chr1
  START  <<- row$start1 - 10000
  END    <<- row$end2 + 10000
   
  MANUAL_ANCHOR_LEFT_RANGE  <<- c(row$start1, row$end1)
  MANUAL_ANCHOR_RIGHT_RANGE <<- c(row$start2, row$end2)
   
  LOADING_START <<- 0
  LOADING_END   <<- 0
   
  loop_id_str <- paste0("Loop_", row$loop_id, "_", CHR_BW, "_", row$start1, "_", row$end2)
  message(sprintf("[Worker %d] Processing %s...", Sys.getpid(), loop_id_str))
   
  tryCatch({
    # Tracks
    p_motif <- plot_motif_track(MOTIF_FILE, CHR_BW, START, END)
     
    signal_plots <- list()
    for(name in names(TRACKS_INFO)) {
      info <- TRACKS_INFO[[name]]
      signal_plots[[name]] <- plot_signal_track(info[1], name, info[2], CHR_BW, START, END)
    }
     
    # Samples
    for (sample_name in names(BATCH_CONFIG)) {
      config <- BATCH_CONFIG[[sample_name]]
       
      if (exists(config$cache_name) && !is.null(get(config$cache_name))) {
        df_raw <- get(config$cache_name)
         
        df_tagged <- process_fragments_dynamic(df_raw, MANUAL_ANCHOR_LEFT_RANGE, MANUAL_ANCHOR_RIGHT_RANGE)
        layout_res <- assign_layout(df_tagged, MANUAL_ANCHOR_LEFT_RANGE, MANUAL_ANCHOR_RIGHT_RANGE)
         
        p_frag <- if(!is.null(layout_res)) {
          plot_final_view(layout_res, config$color_dot)
        } else {
          ggplot() + theme_void() + labs(title = paste(sample_name, "- No Clusters"))
        }
         
        h_ratios <- c(0.6, rep(0.8, length(signal_plots)), 10)
        final_plot <- wrap_plots(c(list(p_motif), signal_plots, list(p_frag)), ncol = 1) + 
          plot_layout(heights = h_ratios)
         
        cat_dir <- file.path(OUT_DIR, row$category_group)
        if(!dir.exists(cat_dir)) dir.create(cat_dir, recursive = TRUE)
         
        out_filename <- file.path(cat_dir, paste0(loop_id_str, "_", sample_name, ".png"))
        ggsave(out_filename, plot = final_plot, width = 12, height = 22, units = "in", limitsize = FALSE, dpi = 200)
         
      } else {
        message(sprintf("    Skipping %s (Data not loaded)", sample_name))
      }
    }
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", loop_id_str, e$message))
    return(FALSE)
  })
}

message(">>> Starting Parallel Batch Processing (6 Cores)...")

results <- mclapply(1:nrow(plot_queue), process_single_loop, mc.cores = 6)
