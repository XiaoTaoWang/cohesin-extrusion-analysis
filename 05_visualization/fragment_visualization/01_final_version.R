# ==============================================================================
# Script Name: plot_chiadrop_batch_dual_cache_v2.R
# Function:    Batch Generate ChIA-Drop Plots (CTCF & Cohesin)
#              FEATURE: Dual Caching Mechanism (Keeps both datasets in memory)
#              1. Checks G_CACHE_CTCF / G_CACHE_COHESIN before loading.
#              2. Generates two separate PDFs.
# Date:        2026-01-31
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
# CHR_BW <- "chr4"
# START  <- 102615000
# END    <- 103015000
# CHR_BW <- "chr10"
# START  <- 123577500
# END    <- 123974000
# CHR_BW <- "chr2"
# START  <- 37300000
# END    <- 38000000

# === Loading Site ===
# LOADING_START  <- 102819516 
# LOADING_END    <- 102829528
# LOADING_START  <- 37663953 
# LOADING_END    <- 37673032
LOADING_START  <- 123834500
LOADING_END    <- 123849500

LOADING_COLOR  <- "#ffd04edd"
# LOADING_COLOR <- "#ffffff"

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

FRAGMENT_SIZE      <- 1.5         # Point size
MOTIF_SHAPE_CODE   <- 15          # 15 = Square
no_motif_size      <- 1.0          # Point size for fragments without motif

MOTIF_COLOR_P      <- "#EF4444"          
MOTIF_COLOR_N      <- "#A78BFA"          
COLOR_LINE         <- "#acada2"
LINE_WIDTH         <- 0.1
COLOR_NIPBL        <- "#E69F00" 

# === Rectangle Config ===
RECT_WIDTH_BP  <- 3000   
RECT_HEIGHT    <- 2.5    

# === File Paths ===
DIR_MAIN  <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between")
DIR_OTHER <- paste0(DATA_ROOT, "/ChIA-Drop/GSE158897")
MOTIF_FILE <- file.path(DIR_MAIN, "CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_id_v2.bed")

# === Tracks Info (Shared) ===
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
# 2. Batch Config & Caching Setup
# ==============================================================================
# simple comment
if (!exists("G_CACHE_CTCF")) { G_CACHE_CTCF <- NULL }
if (!exists("G_CACHE_COHESIN")) { G_CACHE_COHESIN <- NULL }

BATCH_CONFIG <- list(
  "CTCF" = list(
    rds_file = file.path(DIR_MAIN, "GSE158897_GM12878-CTCF-pooled_comp.rds"),
    color_dot = "#0000f3",
    output_name = "ChIA_Drop_CTCF_Final.pdf",
    cache_name = "G_CACHE_CTCF"  # simple comment
  ),
  "Cohesin" = list(
    rds_file = file.path(DIR_MAIN, "GSE158897_GM12878-cohesin-pooled_comp.rds"),
    color_dot = "#008000",
    output_name = "ChIA_Drop_Cohesin_Final.pdf",
    cache_name = "G_CACHE_COHESIN" # simple comment
  )
)

# ==============================================================================
# 3. Helpers
# ==============================================================================
COORD_XLIM <- c(START, END)
COORD_EXPAND <- FALSE

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

# ==============================================================================
# 4. Track Functions
# ==============================================================================
plot_motif_track <- function(bed_file, chrom, start_pos, end_pos) {
  base_plot <- ggplot() +      
    annotate("rect", xmin = LOADING_START, xmax = LOADING_END, ymin = -Inf, ymax = Inf, fill = LOADING_COLOR, alpha = 0.5) +
    coord_cartesian(xlim = COORD_XLIM, ylim = c(0.8, 1.2), expand = COORD_EXPAND, clip = "off") +
    scale_x_continuous(limits = COORD_XLIM, expand = c(0, 0)) +
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
    clean_theme + labs(y = track_name) + coord_cartesian(xlim = COORD_XLIM, expand = COORD_EXPAND) +
    scale_x_continuous(limits = COORD_XLIM, expand = c(0, 0)) +
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
    scale_x_continuous(limits = COORD_XLIM, expand = c(0, 0)) +
    coord_cartesian(xlim = COORD_XLIM, ylim = c(0, rounded_max), expand = COORD_EXPAND, clip = "on") +
    labs(y = track_name) + clean_theme +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 9, face = "bold", margin = margin(r = 5)),
          panel.grid = element_blank(), panel.border = element_blank(), plot.margin = margin(t = 2, b = 2, r = 20, l = 10))
}

# ==============================================================================
# 5. Core Logic
# ==============================================================================
process_fragments_dynamic <- function(df_raw, left_range, right_range) {
  roi_gr <- GRanges(seqnames = CHR_BW, ranges = IRanges(start = START, end = END))
  
  # Handle GRanges or DataFrame
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
  
  # 1. Motif Overlap (Always run to capture middle motifs)
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

  # 2. Custom Interval Mode (Override anchors if provided)
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
  
  if (file.exists(MOTIF_FILE)) {
    motifs_raw <- read.table(MOTIF_FILE, stringsAsFactors = FALSE)
    colnames(motifs_raw)[1:5] <- c("seqnames","start","end","strand","motif_id")
    target_chr <- unique(df$seqnames)[1]
    if (grepl("chr", target_chr)) { if(!any(grepl("chr", motifs_raw$seqnames))) motifs_raw$seqnames <- paste0("chr", motifs_raw$seqnames) }
    else { motifs_raw$seqnames <- gsub("chr", "", motifs_raw$seqnames) }
    
    motifs_plus <- motifs_raw %>% filter(seqnames == target_chr, strand == "+", end >= START, start <= END)
    motifs_minus <- motifs_raw %>% filter(seqnames == target_chr, strand == "-", end >= START, start <= END)
    
    if(!is.null(left_range)) { l_start <- left_range[1]; l_end <- left_range[2] }
    else if(nrow(motifs_plus) > 0) {
      idx <- which.min((motifs_plus$start + motifs_plus$end)/2)
      l_start <- motifs_plus$start[idx]; l_end <- motifs_plus$end[idx]
    }
    
    if(!is.null(right_range)) { r_start <- right_range[1]; r_end <- right_range[2] }
    else if(nrow(motifs_minus) > 0) {
      idx <- which.max((motifs_minus$start + motifs_minus$end)/2)
      r_start <- motifs_minus$start[idx]; r_end <- motifs_minus$end[idx]
    }
  }

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
  unanchored_stats <- classified %>% filter(!cluster_id %in% anchored_ids)
  
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
  
  if (file.exists(MOTIF_FILE)) {
    motifs_raw <- read.table(MOTIF_FILE, stringsAsFactors = FALSE)
    colnames(motifs_raw)[1:5] <- c("seqnames","start","end","strand","motif_id")
    target_chr <- unique(df$seqnames)[1]
    if (grepl("chr", target_chr)) { if(!any(grepl("chr", motifs_raw$seqnames))) motifs_raw$seqnames <- paste0("chr", motifs_raw$seqnames) }
    else { motifs_raw$seqnames <- gsub("chr", "", motifs_raw$seqnames) }
    
    motifs_plus <- motifs_raw %>% filter(seqnames == target_chr, strand == "+", end >= START, start <= END)
    motifs_minus <- motifs_raw %>% filter(seqnames == target_chr, strand == "-", end >= START, start <= END)
    
    if(!is.null(left_range)) { l_start <- left_range[1]; l_end <- left_range[2] }
    else if(nrow(motifs_plus) > 0) {
      idx <- which.min((motifs_plus$start + motifs_plus$end)/2)
      l_start <- motifs_plus$start[idx]; l_end <- motifs_plus$end[idx]
    }
    
    if(!is.null(right_range)) { r_start <- right_range[1]; r_end <- right_range[2] }
    else if(nrow(motifs_minus) > 0) {
      idx <- which.max((motifs_minus$start + motifs_minus$end)/2)
      r_start <- motifs_minus$start[idx]; r_end <- motifs_minus$end[idx]
    }
  }
  
  # 3. Calculate Stats
  stats <- df_sub %>%
    group_by(cluster_id) %>%
    summarise(
      is_anchored_left = any(!is.na(display_strand) & display_strand == "+" & frag_mid >= l_start & frag_mid <= l_end),
      is_anchored_right = any(!is.na(display_strand) & display_strand == "-" & frag_mid >= r_start & frag_mid <= r_end),
      has_plus_total = any(!is.na(display_strand) & display_strand == "+"),
      has_minus_total = any(!is.na(display_strand) & display_strand == "-")
    ) %>%
    filter(is_anchored_left | is_anchored_right) %>%
    summarise(
      complete_count = sum(has_plus_total & has_minus_total),
      ongoing_count = sum((has_plus_total & !has_minus_total) | (!has_plus_total & has_minus_total))
    )
  
  return(sprintf("Complete: %d | Ongoing: %d", stats$complete_count, stats$ongoing_count))
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
    
    scale_x_continuous(limits = COORD_XLIM, breaks = breaks_extended(n = 5), labels = comma, expand = c(0,0)) +
    labs(y = NULL, x = paste0(CHR_BW, ":", comma(START), "-", comma(END))) +      
    coord_cartesian(xlim = COORD_XLIM, expand = COORD_EXPAND, clip = "off") +
    clean_theme +
    scale_y_continuous(breaks = labels$y, labels = labels$label, expand = expansion(mult = c(0.05, 0)))
    
  return(p)
}

# ==============================================================================
# 6. Main Execution (With Dual Cache Support)
# ==============================================================================

# A. Generate Common Tracks (Run Once)
message(">>> Generating Common Tracks (Motif & BigWigs)...")
p_motif <- plot_motif_track(MOTIF_FILE, CHR_BW, START, END)
signal_plots <- list()
for(name in names(TRACKS_INFO)) {
  info <- TRACKS_INFO[[name]]
  signal_plots[[name]] <- plot_signal_track(info[1], name, info[2], CHR_BW, START, END)
}

# B. Loop Over Configurations
run_sample_analysis <- function(sample_name, config) {
  message(sprintf("\n>>> Processing Sample: %s", sample_name))
  message(sprintf("    Cache Var: %s", config$cache_name))
  
  # 1. Cache Logic (Explicit Global Variable Check)
  # simple comment
  if (!exists(config$cache_name) || is.null(get(config$cache_name))) {
    message(sprintf("    Status: Cache miss. Loading RDS: %s", config$rds_file))
    if(!file.exists(config$rds_file)) {
      message("    ERROR: RDS file not found!")
      return(NULL)
    }
    df_temp <- readRDS(config$rds_file)
    # simple comment
    assign(config$cache_name, df_temp, envir = .GlobalEnv)
    df_raw <- df_temp
    message("    Status: Loaded and cached to global environment.")
  } else {
    message("    Status: Cache hit! Using existing data from memory.")
    df_raw <- get(config$cache_name)
  }
  
  # 2. Process
  df_tagged <- process_fragments_dynamic(df_raw, MANUAL_ANCHOR_LEFT_RANGE, MANUAL_ANCHOR_RIGHT_RANGE)
  
  # 2.1 Calculate Stats (Ongoing vs Complete)
  stats_text <- get_interaction_stats(df_tagged, MANUAL_ANCHOR_LEFT_RANGE, MANUAL_ANCHOR_RIGHT_RANGE)
  
  # 3. Layout
  layout_res <- assign_layout(df_tagged, MANUAL_ANCHOR_LEFT_RANGE, MANUAL_ANCHOR_RIGHT_RANGE)
  
  # 4. Plot
  p_frag <- if(!is.null(layout_res)) {
    plot_final_view(layout_res, config$color_dot)
  } else {
    ggplot() + theme_void() + labs(title = paste(sample_name, "- No Clusters"))
  }
  
  # 5. Assemble
  h_ratios <- c(0.6, rep(0.8, length(signal_plots)), 10)
  final_plot <- wrap_plots(c(list(p_motif), signal_plots, list(p_frag)), ncol = 1) + 
    plot_layout(heights = h_ratios) +
    plot_annotation(title = paste0(sample_name, "  [ ", stats_text, " ]"),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
  
  # 6. Save
  ggsave(config$output_name, plot = final_plot, width = 12, height = 22, units = "in", limitsize = FALSE)
  message(sprintf("    Saved: %s", config$output_name))
  print(final_plot)
}

for (sample in names(BATCH_CONFIG)) {
  run_sample_analysis(sample, BATCH_CONFIG[[sample]])
}

