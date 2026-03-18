# ==============================================================================
# Script Name: plot_chiadrop_custom_interval_v2.R
# Function:    ChIA-Drop Visualization (Custom Interval Anchors)
#              1. Supports Custom Interval for Anchors (Left/Right Range)
#              2. Logic: Overlap with Interval -> Classified as Anchored
#              3. Color Priority: Anchor (Red/Purple) > Green (if Anchored) > Orange (only if unanchored)
# Date:        2026-01-28
# function 1. self-defined loading site and anchor site
# function 2. color the anchor site and loading site
# function 3. omit the middle motif (focus on the anchor motifs)
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
})
slice <- dplyr::slice
select <- dplyr::select

# ==============================================================================
# simple comment
# ==============================================================================

# === Region of Interest ===
# CHR_BW <- "chr3"
# START  <- 185990000
# END    <- 186590000
CHR_BW <- "chr3"
START  <- 187700000
END    <- 189000000
# CHR_BW <- "chr3"
# START  <- 187730000
# END    <- 189000000
# CHR_BW <- "chr4"
# START  <- 102615000
# END    <- 103015000
# CHR_BW <- "chr10"
# START  <- 123577500
# END    <- 123974000
# CHR_BW <- "chr2"
# START  <- 37300000
# END    <- 38000000
# chr3:119882000-120156000
# CHR_BW <- "chr3"
# START  <- 119882000
# END    <- 120156000


# === Loading Site ===
# LOADING_START  <- 102819516 
# LOADING_END    <- 102829528
# LOADING_START  <- 37663953 
# LOADING_END    <- 37673032
# LOADING_START  <- 123833000
# LOADING_END    <- 123850000
# LOADING_START <- 120090000
# LOADING_END   <- 120098000
LOADING_START <- 0
LOADING_END   <- 0

# LOADING_COLOR  <- "#ffd04edd"
LOADING_COLOR <- "#ffffff"

# === Anchor Sites (Set to NULL for Auto) ===
# MANUAL_ANCHOR_LEFT_RANGE  <- c(102616000, 102630000) 
# MANUAL_ANCHOR_RIGHT_RANGE <- c(103005308, 103013327)
# MANUAL_ANCHOR_LEFT_RANGE  <- c(37410595, 37419679) 
# MANUAL_ANCHOR_RIGHT_RANGE <- c(37873670, 37881689)

MANUAL_ANCHOR_LEFT_RANGE  <- NULL
MANUAL_ANCHOR_RIGHT_RANGE <- NULL


# === Visual Parameters ===
GAP_SIZE_PERCENT   <- 0.01        
MIN_FRAGS_IN_ROI   <- 2            
EDGE_BUFFER        <- 500          

MOTIF_COLOR_P      <- "#EF4444"          
MOTIF_COLOR_N      <- "#A78BFA"          
COLOR_DOT          <- "#0000f3"
COLOR_LINE         <- "#acada2"
LINE_WIDTH         <- 0.1

# Defined NIPBL Color for Loading Site Fragments
COLOR_NIPBL        <- "#E69F00" 

# === Rectangle Configuration ===
RECT_WIDTH_BP  <- 3000   
RECT_HEIGHT    <- 1.0    

# === File Paths ===
MOTIF_FILE   <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_id_v2.bed")
DIR_MAIN  <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between")
DIR_OTHER <- paste0(DATA_ROOT, "/ChIA-Drop/GSE158897")

# simple comment
if (!exists("G_CACHE_CTCF")) { G_CACHE_CTCF <- NULL }
if (!exists("G_CACHE_COHESIN")) { G_CACHE_COHESIN <- NULL }
if (!exists("G_CACHE_RNAPII")) { G_CACHE_RNAPII <- NULL }

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
  ),
  "RNAPII" = list(
    rds_file = paste0(DATA_ROOT, "/ChIA-Drop/GSE158897/GSE158897_GM12878-RNAPII-pooledv2_comp.rds"),
    color_dot = "#800080",
    cache_name = "G_CACHE_RNAPII"
  )
)

TRACKS_INFO <- list(
  "CTCF"    = c(file.path(DIR_MAIN, "GM12878_CTCF_ChIA-PET_GSE158897_LHG0052H.for.BROWSER.sorted.bw"), "#0000FF"),
  "Cohesin" = c(file.path(DIR_MAIN, "GM12878_cohesin_ChIA-PET_GSE158897_LHG0051H_0104V.for.BROWSER.sorted.bw"), "#009E73"),
  "NIPBL"   = c(file.path(DIR_OTHER, "GM12878_NIPBL_chip-seq_GSE158897_CHG0030.q30.nr.sorted.bw"), COLOR_NIPBL), 
  "WAPL"    = c(file.path(DIR_OTHER, "GM12878_WAPL_chip-seq_GSE158897_CHG0032.q30.nr.sorted.bw"), "#A52A2A"),
  "RNAPII"  = c(file.path(DIR_MAIN, "GM12878_RNAPII_ChIA-PET_GSE158897_LHG0035N_0035V_0045V.for.BROWSER.sorted.bw"), "#800080"),
  "H3K27ac" = c(file.path(DIR_OTHER, "GM12878_h3k27ac_chip-seq_ENCFF340JIF.bigWig"), "#D55E00"),
  "H3K4me1" = c(file.path(DIR_OTHER, "GM12878_h3k4me1_chip-seq_ENCFF831ZHL.bigWig"), "#8B0000")
)

message(sprintf(">>> Plotting Region: %s:%s-%s", CHR_BW, comma(START), comma(END)))
message(sprintf(">>> Loading Site: %s:%s-%s", CHR_BW, comma(LOADING_START), comma(LOADING_END)))

# ==============================================================================
# 2. Theme & Helpers
# ==============================================================================
COORD_XLIM <- c(START, END)
COORD_EXPAND <- FALSE

clean_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(color = "black", size = 10, margin = margin(t = 5)),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_blank(),          
    axis.line.y = element_blank(),      
    axis.ticks.y = element_blank(),
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

# ==============================================================================
# 3. Track Functions
# ==============================================================================
plot_motif_track <- function(bed_file, chrom, start_pos, end_pos) {
  base_plot <- ggplot() +      
    annotate("rect", xmin = LOADING_START, xmax = LOADING_END, 
             ymin = -Inf, ymax = Inf, fill = LOADING_COLOR, alpha = 0.5) +
    coord_cartesian(xlim = COORD_XLIM, ylim = c(0.8, 1.2), expand = COORD_EXPAND, clip = "off") +
    scale_x_continuous(limits = COORD_XLIM, expand = c(0, 0)) +
    labs(y = "CTCF Motif") + clean_theme +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), plot.margin = margin(t = 5, b = 0, r = 20, l = 10))
  
  if (!file.exists(bed_file)) return(base_plot)
  
  df <- tryCatch({ 
    df_raw <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
    df_raw[,1:4]
  }, error = function(e) NULL)
  
  if (is.null(df) || nrow(df) == 0) return(base_plot)
  colnames(df) <- c("chr", "start", "end", "strand")
  
  if (grepl("chr", chrom)) { 
    if (!any(grepl("chr", df$chr))) df$chr <- paste0("chr", df$chr) 
  } else { 
    df$chr <- gsub("chr", "", df$chr) 
  }
  
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
    annotate("rect", xmin = LOADING_START, xmax = LOADING_END, 
             ymin = -Inf, ymax = Inf, fill = LOADING_COLOR, alpha = 0.5) +
    clean_theme + labs(y = track_name) +      
    coord_cartesian(xlim = COORD_XLIM, expand = COORD_EXPAND) +
    scale_x_continuous(limits = COORD_XLIM, expand = c(0, 0)) +
    theme(axis.line.x=element_blank(), axis.text.y = element_blank(), axis.ticks.x=element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 9, face = "bold", margin = margin(r = 5)),
          plot.margin = margin(t = 2, b = 2, r = 20, l = 10))
  
  if (!file.exists(bw_file)) {
    message(sprintf("Warning: File not found for %s: %s", track_name, bw_file))
    return(p_void)
  }
  
  bw_obj <- BigWigFile(bw_file)
  roi_gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start_pos, end = end_pos))
  
  target_bins <- 1000 
  
  binned_list <- tryCatch({
    summary(bw_obj, which = roi_gr, size = target_bins, type = "mean")
  }, error = function(e) NULL)
  
  if (is.null(binned_list) || length(binned_list) == 0) return(p_void)
  
  binned_gr <- binned_list[[1]]
  if (length(binned_gr) == 0) return(p_void)
  
  plot_df <- data.frame(
    pos = (start(binned_gr) + end(binned_gr)) / 2,
    score = binned_gr$score
  )
  
  plot_df$score[is.na(plot_df$score)] <- 0
  plot_df$score[plot_df$score < 0] <- 0
  
  smooth_fit <- ksmooth(plot_df$pos, plot_df$score, kernel = "normal", bandwidth = (end_pos - start_pos) / 800)
  
  final_df <- data.frame(pos = smooth_fit$x, score = smooth_fit$y)
  
  raw_max <- max(final_df$score, na.rm=TRUE)
  q995 <- quantile(final_df$score, 0.995, na.rm = TRUE)
  effective_max <- if (raw_max > 4 * q995 && q995 > 1) q995 * 1.5 else raw_max * 1.1
  effective_max <- max(effective_max, 1)
  rounded_max <- round_smart(effective_max)
  
  ggplot(final_df, aes(x = pos, y = score)) +
    annotate("rect", xmin = LOADING_START, xmax = LOADING_END, 
             ymin = 0, ymax = rounded_max, fill = LOADING_COLOR, alpha = 0.5) +
    
    geom_area(fill = fill_color, alpha = 0.9, outline.type = "upper") + 
    geom_line(color = fill_color, linewidth = 0.2) +
    
    scale_y_continuous(breaks = c(0, rounded_max), labels = c("0", rounded_max), expand = c(0, 0)) +
    scale_x_continuous(limits = COORD_XLIM, expand = c(0, 0)) +
    coord_cartesian(xlim = COORD_XLIM, ylim = c(0, rounded_max), expand = COORD_EXPAND, clip = "on") +
    
    labs(y = track_name) + clean_theme +
    theme(
      axis.line.x = element_blank(), 
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 9, face = "bold", margin = margin(r = 5)),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 2, b = 2, r = 20, l = 10)
    )
}

# ==============================================================================
# 4. Processing Fragments (Supports DUAL MODES: Custom Interval vs Auto Motif)
# ==============================================================================
process_fragments_dynamic <- function(df_raw, left_range, right_range) {
  roi_gr <- GRanges(seqnames = CHR_BW, ranges = IRanges(start = START, end = END))
  
  if(class(df_raw)[1] == "GRanges") {
    df_gr <- suppressWarnings(subsetByOverlaps(df_raw, roi_gr))
    df <- as.data.frame(df_gr)
  } else if(is.data.frame(df_raw)) {
    df <- df_raw %>%
      filter(seqnames == as.character(seqnames(roi_gr)), end >= start(roi_gr), start <= end(roi_gr))
  } else { stop("Unknown data type!") }
  
  if(nrow(df) == 0) return(NULL)
  
  if (!"cluster_id" %in% colnames(df)) {
    if("CBMB" %in% colnames(df)) df$cluster_id <- df$CBMB else df$cluster_id <- 1:nrow(df)
  }
  
  df$frag_mid <- (df$start + df$end) / 2
  df$display_strand <- NA
  df$motif_id <- NA
  
  # --------------------------------------------------------------------------
  # MODE 1: CUSTOM INTERVAL MODE
  # --------------------------------------------------------------------------
  if (!is.null(left_range) || !is.null(right_range)) {
    message(">>> Mode: Custom Interval Overlap (Ignoring Motif File for Classification)")
    
    # Check Overlap with Left Interval -> Assign "+" (Forward)
    if (!is.null(left_range)) {
      is_left <- df$frag_mid >= left_range[1] & df$frag_mid <= left_range[2]
      df$display_strand[is_left] <- "+"
    }
    
    # Check Overlap with Right Interval -> Assign "-" (Reverse)
    if (!is.null(right_range)) {
      is_right <- df$frag_mid >= right_range[1] & df$frag_mid <= right_range[2]
      df$display_strand[is_right] <- "-"
    }
    
    # Note: If a fragment overlaps both (unlikely given distances), Right wins in this simple logic, 
    # but practically they are distinct regions.
    
  } else {
    # --------------------------------------------------------------------------
    # MODE 2: AUTO MOTIF MODE (Original Logic)
    # --------------------------------------------------------------------------
    message(">>> Mode: Auto Motif Detection (Fallback)")
    
    if (!file.exists(MOTIF_FILE)) return(NULL)
    
    motifs_raw <- read.table(MOTIF_FILE, stringsAsFactors = FALSE)
    colnames(motifs_raw)[1:5] <- c("seqnames","start","end","strand","motif_id")
    
    target_chr <- unique(df$seqnames)[1]
    if (grepl("chr", target_chr)) {      
      if(!any(grepl("chr", motifs_raw$seqnames))) motifs_raw$seqnames <- paste0("chr", motifs_raw$seqnames)    
    } else { motifs_raw$seqnames <- gsub("chr", "", motifs_raw$seqnames) }
    
    # Filter to only leftmost and rightmost motifs in region
    motifs_in_region <- motifs_raw %>%
      filter(seqnames == target_chr, end >= START, start <= END)
    
    if (nrow(motifs_in_region) > 0) {
      motifs_in_region$mid <- (motifs_in_region$start + motifs_in_region$end) / 2
      leftmost_idx <- which.min(motifs_in_region$mid)
      rightmost_idx <- which.max(motifs_in_region$mid)
      
      if (leftmost_idx != rightmost_idx) {
        motifs_in_region <- motifs_in_region[c(leftmost_idx, rightmost_idx), ]
      } else {
        motifs_in_region <- motifs_in_region[leftmost_idx, , drop = FALSE]
      }
    }
    
    motifs <- makeGRangesFromDataFrame(motifs_in_region, keep.extra.columns = TRUE)
    frag_gr <- GRanges(seqnames = df$seqnames, ranges = IRanges(start = df$start, end = df$end))
    
    ov <- suppressWarnings(findOverlaps(frag_gr, motifs))
    ov_df <- data.frame(f_idx = queryHits(ov), m_idx = subjectHits(ov), m_str = as.character(strand(motifs[subjectHits(ov)])))
    
    strand_summ <- ov_df %>%
      group_by(f_idx) %>%
      summarise(has_pos = "+" %in% m_str, has_neg = "-" %in% m_str) %>%
      mutate(final_str = case_when(has_pos & !has_neg ~ "+", !has_pos & has_neg ~ "-", TRUE ~ "Conflict"))
    
    df$display_strand[strand_summ$f_idx] <- strand_summ$final_str
  }
  
  # Flags for loading and display region
  df$in_display_region <- df$frag_mid >= START & df$frag_mid <= END
  df$in_loading <- df$frag_mid >= LOADING_START & df$frag_mid <= LOADING_END
  
  return(df)
}

# ==============================================================================
# simple comment
# ==============================================================================
classify_clusters_dynamic <- function(df, left_range, right_range) {
  if(is.null(df) || nrow(df) == 0) return(NULL)
  
  # --------------------------------------------------------------------------
  # simple comment
  # --------------------------------------------------------------------------
  # simple comment
  # simple comment
  # "+" = Hit Left Interval / Forward Motif
  # "-" = Hit Right Interval / Reverse Motif
  
  clusters_with_anchors <- df %>%
    group_by(cluster_id) %>%
    summarise(has_anchor_signal = any(!is.na(display_strand) & display_strand %in% c("+", "-"))) %>%
    filter(has_anchor_signal) %>%
    pull(cluster_id)
  
  # --------------------------------------------------------------------------
  # simple comment
  # --------------------------------------------------------------------------
  cluster_stats <- df %>%
    group_by(cluster_id) %>%
    summarise(
      c_start = min(frag_mid),
      c_end = max(frag_mid),
      span = c_end - c_start,
      
      n_frags_left = sum(frag_mid < LOADING_START),
      n_frags_right = sum(frag_mid > LOADING_END),
      n_frags_in_loading = sum(in_loading)
    )
  
  # --------------------------------------------------------------------------
  # simple comment
  # --------------------------------------------------------------------------
  anchored_analysis <- df %>%
    filter(cluster_id %in% clusters_with_anchors) %>%
    filter(!is.na(display_strand)) %>%
    group_by(cluster_id) %>%
    summarise(
      has_forward = "+" %in% display_strand, 
      has_reverse = "-" %in% display_strand  
    ) %>%
    mutate(
      is_convergent = has_forward & has_reverse
    )
  
  anchored_classified <- anchored_analysis %>%
    left_join(cluster_stats, by = "cluster_id") %>%
    mutate(
      category = case_when(
        is_convergent ~ "Convergent",
        has_forward ~ "Forward Anchored",
        has_reverse ~ "Reverse Anchored",
        TRUE ~ "Other Anchored"
      ),
      # simple comment
      # simple comment
      # simple comment
      # simple comment
      # simple comment
      sort_key = case_when(
        category == "Forward Anchored" ~ c_end,
        category == "Reverse Anchored" ~ c_start,
        TRUE ~ c_start
      )
    )
  
  # --------------------------------------------------------------------------
  # simple comment
  # --------------------------------------------------------------------------
  clusters_unanchored <- setdiff(unique(df$cluster_id), clusters_with_anchors)
  
  unanchored_stats <- cluster_stats %>%
    filter(cluster_id %in% clusters_unanchored)
  
  # simple comment
  twosided_candidates <- unanchored_stats %>%
    filter(
      n_frags_left >= 1,
      n_frags_right >= 1,
      n_frags_in_loading == 0
    ) %>%
    arrange(span)
  
  twosided_selected <- data.frame()
  if(nrow(twosided_candidates) > 0) {
    for(i in 1:nrow(twosided_candidates)) {
      current <- twosided_candidates[i, ]
      if(nrow(twosided_selected) == 0) {
        twosided_selected <- rbind(twosided_selected, current)
      } else {
        prev_max_left <- min(twosided_selected$c_start)
        prev_max_right <- max(twosided_selected$c_end)
        # simple comment
        if(current$c_start < prev_max_left && current$c_end > prev_max_right) {
          twosided_selected <- rbind(twosided_selected, current)
        }
      }
    }
  }
  
  twosided_ids <- if(nrow(twosided_selected) > 0) twosided_selected$cluster_id else c()
  
  unanchored_classified <- unanchored_stats %>%
    mutate(
      category = case_when(
        cluster_id %in% twosided_ids ~ "Two-sided at loading",
        n_frags_in_loading > 0 & ((n_frags_left > 0 & n_frags_right == 0) | (n_frags_left == 0 & n_frags_right > 0)) ~ "One-sided at loading",
        TRUE ~ "Other without anchor"
      ),
      extends_right = n_frags_right > 0,
      # simple comment
      sort_key = span 
    )
  
  # --------------------------------------------------------------------------
  # simple comment
  # --------------------------------------------------------------------------
  final_categories <- bind_rows(
    anchored_classified %>% select(cluster_id, category, span, sort_key) %>% mutate(extends_right = NA),
    unanchored_classified %>% select(cluster_id, category, span, extends_right, sort_key)
  )
  
  return(final_categories)
}

# ==============================================================================
# simple comment
# ==============================================================================
assign_layout <- function(df, left_range, right_range) {
  if(is.null(df) || nrow(df) == 0) return(NULL)
  
  # --------------------------------------------------------------------------
  # simple comment
  # --------------------------------------------------------------------------
  valid_ids <- df %>% 
    filter(in_display_region) %>% 
    group_by(cluster_id) %>% summarize(n=n()) %>% filter(n >= MIN_FRAGS_IN_ROI) %>% pull(cluster_id)
  
  df <- df %>% filter(cluster_id %in% valid_ids)
  
  # simple comment
  cluster_boundaries <- df %>%
    group_by(cluster_id) %>%
    summarize(l=min(frag_mid), r=max(frag_mid)) %>%
    filter(l >= (START + EDGE_BUFFER) & r <= (END - EDGE_BUFFER)) %>% pull(cluster_id)
  
  df <- df %>% filter(cluster_id %in% cluster_boundaries)
  
  # --------------------------------------------------------------------------
  # simple comment
  # --------------------------------------------------------------------------
  cats <- classify_clusters_dynamic(df, left_range, right_range)
  if(is.null(cats) || nrow(cats) == 0) return(NULL)
  
  # --------------------------------------------------------------------------
  # simple comment
  # --------------------------------------------------------------------------
  
  # simple comment
  # simple comment
  forward <- cats %>% 
    filter(category == "Forward Anchored") %>% 
    arrange(sort_key, span)
  
  # simple comment
  # simple comment
  reverse <- cats %>% 
    filter(category == "Reverse Anchored") %>% 
    arrange(desc(sort_key), span)
  
  # simple comment
  conv <- cats %>% 
    filter(category == "Convergent") %>% 
    arrange(span)
  
  # simple comment
  onesided <- cats %>% 
    filter(category == "One-sided at loading") %>% 
    arrange(desc(extends_right), span)
  
  # simple comment
  twosided <- cats %>% 
    filter(category == "Two-sided at loading") %>% 
    arrange(span)
  
  # --------------------------------------------------------------------------
  # simple comment
  # --------------------------------------------------------------------------
  total_clusters <- nrow(forward) + nrow(reverse) + nrow(conv) + nrow(onesided) + nrow(twosided)
  gap_size <- max(2, round(total_clusters * GAP_SIZE_PERCENT))
  n_groups <- 5
  
  # simple comment
  cursor <- total_clusters + (n_groups * gap_size)
  
  # simple comment
  add_group <- function(grp, label, h) {
    if(nrow(grp) == 0) return(list(df=data.frame(), lbl=data.frame(), h=h))
    
    # simple comment
    grp$y_idx <- seq(h, h - nrow(grp) + 1)
    
    lbl <- data.frame(label = paste0(label, "\nn=", nrow(grp)), y = h - (nrow(grp)/2))
    list(
      df = grp %>% select(cluster_id, y_idx, category), 
      lbl = lbl, 
      h = h - nrow(grp)
    )
  }
  
  res_list <- list()
  lbl_list <- list()
  
  # simple comment
  step1 <- add_group(forward, "Forward\nAnchored", cursor)
  step2 <- add_group(reverse, "Reverse\nAnchored", step1$h - gap_size)
  step3 <- add_group(conv, "Convergent", step2$h - gap_size)
  step4 <- add_group(onesided, "One-sided\nat loading", step3$h - gap_size)
  step5 <- add_group(twosided, "Two-sided\nat loading", step4$h - gap_size)
  
  layout_df <- bind_rows(step1$df, step2$df, step3$df, step4$df, step5$df)
  labels_df <- bind_rows(step1$lbl, step2$lbl, step3$lbl, step4$lbl, step5$lbl)
  
  if(nrow(layout_df) == 0) return(NULL)
  
  return(list(
    data = df %>% inner_join(layout_df, by="cluster_id"),
    labels = labels_df
  ))
}

# ==============================================================================
# 7. Plotting Final View (Revised Color Logic)
# ==============================================================================
plot_final_view <- function(layout_res, current_dot_color) {
  df <- layout_res$data
  labels <- layout_res$labels
  coord_label <- paste0(CHR_BW, ":", comma(START), "-", comma(END))
  
  SIZE_DOT    <- 1.5    
  
  # === REVISED COLOR LOGIC ===
  # Priority 1: Anchors (Left/Forward = Red, Right/Reverse = Purple)
  # Priority 2: Unanchored Loading (Orange)
  # Priority 3: Everything else (Green/Dot) -> Includes Anchored fragments passing through Loading
  
  # 1. Anchor Fragments (Rectangles)
  # Captured by 'display_strand' which is set by custom interval or motif
  df_anchor <- df %>% 
    filter(!is.na(display_strand) & display_strand %in% c("+", "-"))
  
  # 2. Loading Fragments (Orange Rectangles)
  # CRITICAL CHANGE: Only color Orange if the cluster is NOT anchored.
  # If the cluster is Forward/Reverse/Convergent, force Green (avoid Orange).
  anchored_categories <- c("Forward Anchored", "Reverse Anchored", "Convergent")
  
  df_loading_orange <- df %>%
    filter(in_loading) %>%
    filter(!category %in% anchored_categories) # Exclude anchored complexes from orange coloring
  
  # 3. Green Dots (Everything else)
  # Exclude fragments already drawn as Anchor or Orange Loading
  ids_anchor_rows <- paste(df_anchor$cluster_id, df_anchor$frag_mid)
  ids_orange_rows <- paste(df_loading_orange$cluster_id, df_loading_orange$frag_mid)
  ids_current_rows <- paste(df$cluster_id, df$frag_mid)
  
  df_dot <- df %>%
    filter(!ids_current_rows %in% c(ids_anchor_rows, ids_orange_rows))
  
  p <- ggplot() +
    # Background
    annotate("rect", xmin = LOADING_START, xmax = LOADING_END, 
             ymin = -Inf, ymax = Inf, fill = LOADING_COLOR, alpha = 0.5) +
    
    # Layer 1: Connecting Lines
    geom_line(data = df, aes(x = frag_mid, y = y_idx, group = cluster_id), 
              linewidth = LINE_WIDTH, color = COLOR_LINE) +
    
    # Layer 2: Green Dots (Includes anchored fragments that cross loading site)
    geom_point(data = df_dot,
               aes(x = frag_mid, y = y_idx), 
               size = SIZE_DOT, color = current_dot_color, shape = 19) +
    
    # Layer 3: Unanchored Loading Fragments (Orange)
    geom_point(data = df_loading_orange,
              aes(x = frag_mid, y = y_idx),
              color = COLOR_NIPBL, 
              size = SIZE_DOT, shape = 15) +
    
    # Layer 4: Anchor Fragments (Red/Purple) - Highest Priority
    geom_point(data = df_anchor,
              aes(x = frag_mid, y = y_idx, color = display_strand),
              size = 2.0, shape = 15) +
    
    # Scales
    scale_color_manual(values = c("+" = MOTIF_COLOR_P, "-" = MOTIF_COLOR_N)) +
    
    scale_x_continuous(limits = COORD_XLIM, breaks = breaks_extended(n = 5), labels = comma, expand = c(0,0)) +
    labs(y = NULL, x = coord_label) +      
    coord_cartesian(xlim = COORD_XLIM, expand = COORD_EXPAND, clip = "off") +
    clean_theme
  
  if(!is.null(labels) && nrow(labels) > 0) {
    p <- p + scale_y_continuous(
      breaks = labels$y, 
      labels = labels$label, 
      expand = expansion(mult = c(0.05, 0))
    )
  } else {
    p <- p + scale_y_continuous(breaks = NULL, expand = expansion(mult = c(0.05, 0)))
  }
  
  return(p)
}

get_interaction_stats <- function(df, left_range, right_range) {
  if (is.null(df) || nrow(df) == 0) return("Complete: 0 | Ongoing: 0")
  
  # 1. Filter valid clusters (Must match assign_layout logic)
  valid_ids <- df %>% filter(in_display_region) %>% group_by(cluster_id) %>% summarize(n=n()) %>% filter(n >= MIN_FRAGS_IN_ROI) %>% pull(cluster_id)
  df_sub <- df %>% filter(cluster_id %in% valid_ids)
  
  bound_ids <- df_sub %>% group_by(cluster_id) %>% summarize(l=min(frag_mid), r=max(frag_mid)) %>% 
    filter(l >= (START + EDGE_BUFFER) & r <= (END - EDGE_BUFFER)) %>% pull(cluster_id)
  df_sub <- df_sub %>% filter(cluster_id %in% bound_ids)
  
  if(nrow(df_sub) == 0) return("Complete: 0 | Ongoing: 0")
  
  # 2. Define Anchors (Logic copied from classify_clusters_dynamic)
  l_start <- START; l_end <- START; r_start <- END; r_end <- END
  
  # Note: In this script, anchors are defined by custom intervals or motifs in process_fragments_dynamic
  # and stored in display_strand. We can use display_strand directly.
  # However, to be consistent with "Complete" definition (has + and -), we check display_strand.
  
  # 3. Calculate Stats
  stats <- df_sub %>%
    group_by(cluster_id) %>%
    summarise(
      has_plus_total = any(!is.na(display_strand) & display_strand == "+"),
      has_minus_total = any(!is.na(display_strand) & display_strand == "-")
    ) %>%
    filter(has_plus_total | has_minus_total) %>%
    summarise(
      complete_count = sum(has_plus_total & has_minus_total),
      ongoing_count = sum((has_plus_total & !has_minus_total) | (!has_plus_total & has_minus_total))
    )
  
  return(sprintf("Complete: %d | Ongoing: %d", stats$complete_count, stats$ongoing_count))
}

# ==============================================================================
# 8. Execution
# ==============================================================================
message("1. Initializing...")

# Determine if we are in Custom Mode or Auto Mode
use_custom_intervals <- !is.null(MANUAL_ANCHOR_LEFT_RANGE) || !is.null(MANUAL_ANCHOR_RIGHT_RANGE)

if (use_custom_intervals) {
  message(">>> CUSTOM INTERVAL MODE ACTIVE")
  message(sprintf("   Left Interval: %s - %s", comma(MANUAL_ANCHOR_LEFT_RANGE[1]), comma(MANUAL_ANCHOR_LEFT_RANGE[2])))
  message(sprintf("   Right Interval: %s - %s", comma(MANUAL_ANCHOR_RIGHT_RANGE[1]), comma(MANUAL_ANCHOR_RIGHT_RANGE[2])))
} else {
  message(">>> AUTO MOTIF MODE ACTIVE (No custom intervals provided)")
}

message("2. Generating Tracks...")
p_motif <- plot_motif_track(MOTIF_FILE, CHR_BW, START, END)

signal_plots <- list()
for(name in names(TRACKS_INFO)) {
  info <- TRACKS_INFO[[name]]
  signal_plots[[name]] <- plot_signal_track(info[1], name, info[2], CHR_BW, START, END)
}

process_sample_plot <- function(sample_name, config) {
  message(sprintf("\n>>> Processing Sample: %s", sample_name))
  message(sprintf("    Cache Var: %s", config$cache_name))
  
  # 1. Cache Logic
  if (!exists(config$cache_name) || is.null(get(config$cache_name))) {
    message(sprintf("    Status: Cache miss. Loading RDS: %s", config$rds_file))
    if(!file.exists(config$rds_file)) {
      message("    ERROR: RDS file not found!")
      return(NULL)
    }
    df_temp <- readRDS(config$rds_file)
    assign(config$cache_name, df_temp, envir = .GlobalEnv)
    df_raw <- df_temp
  } else {
    message("    Status: Cache hit! Using existing data from memory.")
    df_raw <- get(config$cache_name)
  }
  
  # 2. Process
  df_tagged <- process_fragments_dynamic(df_raw, MANUAL_ANCHOR_LEFT_RANGE, MANUAL_ANCHOR_RIGHT_RANGE)
  
  # 2.1 Calculate Stats
  stats_text <- get_interaction_stats(df_tagged, MANUAL_ANCHOR_LEFT_RANGE, MANUAL_ANCHOR_RIGHT_RANGE)
  
  # 3. Layout
  layout_res <- assign_layout(df_tagged, MANUAL_ANCHOR_LEFT_RANGE, MANUAL_ANCHOR_RIGHT_RANGE)
  
  n_clusters <- 0
  if (!is.null(layout_res)) {
    n_clusters <- length(unique(layout_res$data$cluster_id))
  }

  # 4. Plot
  p_frag <- if(!is.null(layout_res)) {
    plot_final_view(layout_res, config$color_dot)
  } else {
    ggplot() + theme_void() + labs(title = paste(sample_name, "- No Clusters"))
  }
  
  p_frag <- p_frag + ggtitle(paste0(sample_name, "  [ ", stats_text, " ]")) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    
  return(list(plot = p_frag, count = n_clusters))
}

fragment_results <- list()
for (sample in names(BATCH_CONFIG)) {
  fragment_results[[sample]] <- process_sample_plot(sample, BATCH_CONFIG[[sample]])
}

fragment_plots <- lapply(fragment_results, function(x) x$plot)
fragment_counts <- sapply(fragment_results, function(x) x$count)

# 5. Assemble Final Combined Plot
message("\n>>> Assembling Combined Plot...")
plot_list <- c(list(p_motif), signal_plots, fragment_plots)

# Dynamic height calculation
track_heights <- c(0.6, rep(0.8, length(signal_plots)))
# Scale: 1 inch per 35 clusters, min 4 inches
frag_heights <- pmax(fragment_counts * (1/20), 4)
h_ratios <- c(track_heights, frag_heights)

final_plot <- wrap_plots(plot_list, ncol = 1) + plot_layout(heights = h_ratios)
print(final_plot)

total_pdf_height <- sum(h_ratios) + 1 # Add small buffer
output_file <- "ChIA_Drop_Combined_CTCF_Cohesin_RNAPII_Separate.pdf"
ggsave(output_file, plot = final_plot, width = 14, height = total_pdf_height, units = "in", limitsize = FALSE)
message(sprintf("    Saved: %s", output_file))

message("\n>>> Batch Processing Complete!")