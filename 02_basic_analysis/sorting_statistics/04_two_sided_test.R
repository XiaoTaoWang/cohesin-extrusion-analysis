# ==============================================================================
# Script Name: plot_chiadrop_twosided_progressive.R
# Function:    ChIA-Drop Visualization - Two-sided Progressive Filtering
# Date:        2026-01-27
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
# 1. Configuration
# ==============================================================================
# chr10:123560000-123980000 (420 kb)
# chr4	101344930	101348138
# chr4	102340756	102346087
# chr4	102819516	102829528
# chr4	104487765	104495478
# chr4	107988850	107995179
# chr4	109432932	109435097
# chr4	109702244	109703823
# chr4	102616803	102624803	cr2240_S
# chr4	103005308	103013308	cr2240_E
# chr4	102819516	102829528	cr2240_M-1
# chr4	102616803	102624803	cr2241_S
# chr4	103095948	103103948	cr2241_E
# chr4	102819516	102829528	cr2241_M-1
# chr10:123560000-123980000 (420 kb)
# chr10	123390435	123398473
# chr10	123834337	123849353
# chr10	124686695	124690191
# chr10	124695557	124698107
# chr10	124717806	124724115
# chr10	124740072	124744963


CHR_BW <- "chr10"
START  <- 123565000
END    <- 123980000

# === Loading Site Configuration ===
LOADING_START  <- 123835000 
LOADING_END    <- 123849000 
LOADING_COLOR  <- "#f9ed81"

# === Anchor Sites (leftmost and rightmost motifs) ===
# These will be determined automatically from the motif file
# or can be set manually here if known
ANCHOR_LEFT  <- NULL  # Will be auto-detected
ANCHOR_RIGHT <- NULL  # Will be auto-detected

GAP_SIZE_PERCENT   <- 0.01        
MIN_FRAGS_IN_ROI   <- 2            
EDGE_BUFFER        <- 500         

MOTIF_COLOR_P      <- "#EF4444"         
MOTIF_COLOR_N      <- "#A78BFA"         
COLOR_DOT      <- "#008000"
COLOR_LINE     <- "#acada2"
LINE_WIDTH     <- 0.1

# === File Paths ===
MOTIF_FILE   <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_id_v2.bed")
RDS_FILE     <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between/GSE158897_GM12878-cohesin-pooled_comp.rds")

DIR_MAIN  <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between")
DIR_OTHER <- paste0(DATA_ROOT, "/ChIA-Drop/GSE158897")

TRACKS_INFO <- list(
  "CTCF"    = c(file.path(DIR_MAIN, "GM12878_CTCF_ChIA-PET_GSE158897_LHG0052H.for.BROWSER.sorted.bw"), "#0000FF"),
  "Cohesin" = c(file.path(DIR_MAIN, "GM12878_cohesin_ChIA-PET_GSE158897_LHG0051H_0104V.for.BROWSER.sorted.bw"), "#009E73"),
  "NIPBL"   = c(file.path(DIR_OTHER, "GM12878_NIPBL_chip-seq_GSE158897_CHG0030.q30.nr.sorted.bw"), "#E69F00"),
  "WAPL"    = c(file.path(DIR_OTHER, "GM12878_WAPL_chip-seq_GSE158897_CHG0032.q30.nr.sorted.bw"), "#A52A2A"),
  "RNAPII"  = c(file.path(DIR_MAIN, "GM12878_RNAPII_ChIA-PET_GSE158897_LHG0035N_0035V_0045V.for.BROWSER.sorted.bw"), "#800080"),
  "H3K27ac" = c(file.path(DIR_OTHER, "GM12878_h3k27ac_chip-seq_ENCFF340JIF.bigWig"), "#D55E00"),
  "H3K4me1" = c(file.path(DIR_OTHER, "GM12878_h3k4me1_chip-seq_ENCFF831ZHL.bigWig"), "#8B0000")
)

message(sprintf(">>> Plotting Region: %s:%s-%s", CHR_BW, comma(START), comma(END)))
message(sprintf(">>> Loading Site: %s:%s-%s", CHR_BW, comma(LOADING_START), comma(LOADING_END)))

if (!exists("G_FRAGMENTS_CACHE")) { G_FRAGMENTS_CACHE <- NULL }

# ==============================================================================
# 2. Theme
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
# 4. Get Anchor Positions
# ==============================================================================
get_anchor_positions <- function(motif_file, chrom, start_pos, end_pos) {
  if (!file.exists(motif_file)) {
    message("Warning: Motif file not found, cannot determine anchors")
    return(list(left = NULL, right = NULL))
  }
  
  motifs_raw <- read.table(motif_file, stringsAsFactors = FALSE)
  colnames(motifs_raw)[1:4] <- c("seqnames","start","end","strand")
  
  if (grepl("chr", chrom)) {     
    if(!any(grepl("chr", motifs_raw$seqnames))) motifs_raw$seqnames <- paste0("chr", motifs_raw$seqnames)   
  } else { 
    motifs_raw$seqnames <- gsub("chr", "", motifs_raw$seqnames) 
  }
  
  motifs_in_region <- motifs_raw %>%
    filter(seqnames == chrom, end >= start_pos, start <= end_pos)
  
  if (nrow(motifs_in_region) == 0) {
    return(list(left = NULL, right = NULL))
  }
  
  motifs_in_region$mid <- (motifs_in_region$start + motifs_in_region$end) / 2
  
  left_anchor <- min(motifs_in_region$mid)
  right_anchor <- max(motifs_in_region$mid)
  
  message(sprintf(">>> Detected Anchors: Left=%s, Right=%s", comma(left_anchor), comma(right_anchor)))
  
  return(list(left = left_anchor, right = right_anchor))
}

# ==============================================================================
# 5. Processing Fragments (with edge motifs only)
# ==============================================================================
process_fragments_and_motifs <- function(anchor_left, anchor_right) {
  roi_gr <- GRanges(seqnames = CHR_BW, ranges = IRanges(start = START, end = END))
  
  if (!is.null(G_FRAGMENTS_CACHE)) { 
    df_raw <- G_FRAGMENTS_CACHE 
  } else {     
    if(!file.exists(RDS_FILE)) stop("RDS not found")
    df_raw <- readRDS(RDS_FILE)
    G_FRAGMENTS_CACHE <<- df_raw
  }
  
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
  
  if (!file.exists(MOTIF_FILE)) return(NULL)
  
  motifs_raw <- read.table(MOTIF_FILE, stringsAsFactors = FALSE)
  colnames(motifs_raw)[1:5] <- c("seqnames","start","end","strand","motif_id")
  
  target_chr <- unique(df$seqnames)[1]
  if (grepl("chr", target_chr)) {     
    if(!any(grepl("chr", motifs_raw$seqnames))) motifs_raw$seqnames <- paste0("chr", motifs_raw$seqnames)   
  } else { motifs_raw$seqnames <- gsub("chr", "", motifs_raw$seqnames) }
  
  # Filter to only leftmost and rightmost motifs
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
  
  id_selection <- ov_df %>%     
    inner_join(strand_summ %>% dplyr::select(f_idx, final_str), by="f_idx") %>%
    filter( (final_str == "+" & m_str == "+") | (final_str == "-" & m_str == "-") ) %>%
    group_by(f_idx) %>% dplyr::slice(1) %>% dplyr::select(f_idx, valid_m_idx = m_idx)
  
  df$display_strand[strand_summ$f_idx] <- strand_summ$final_str
  df$motif_id[id_selection$f_idx] <- id_selection$valid_m_idx
  
  df$in_display_region <- df$frag_mid >= START & df$frag_mid <= END
  df$in_loading <- df$frag_mid >= LOADING_START & df$frag_mid <= LOADING_END
  
  return(df)
}

# ==============================================================================
# 6. Classification with Two-sided Progressive Filtering
# ==============================================================================
classify_all_with_twosided_progressive <- function(df, anchor_left, anchor_right) {
  if(is.null(df) || nrow(df) == 0) return(NULL)
  
  # 1. Identify clusters with motifs
  clusters_with_motifs <- df %>%
    group_by(cluster_id) %>%
    summarise(has_any_motif = any(!is.na(display_strand) & display_strand != "Conflict")) %>%
    filter(has_any_motif) %>%
    pull(cluster_id)
  
  # Motif Classification
  cluster_motif_analysis <- df %>%
    filter(cluster_id %in% clusters_with_motifs) %>%
    filter(!is.na(display_strand) & display_strand != "Conflict") %>%
    group_by(cluster_id) %>%
    summarise(
      has_forward = "+" %in% display_strand,
      has_reverse = "-" %in% display_strand,
      leftmost_plus_pos = ifelse(any(display_strand == "+"), min(frag_mid[display_strand == "+"]), NA),
      rightmost_minus_pos = ifelse(any(display_strand == "-"), max(frag_mid[display_strand == "-"]), NA)
    ) %>%
    mutate(
      is_convergent = has_forward & has_reverse & 
        !is.na(leftmost_plus_pos) & !is.na(rightmost_minus_pos) & 
        leftmost_plus_pos < rightmost_minus_pos
    )
  
  cluster_stats_motif <- df %>%
    filter(cluster_id %in% clusters_with_motifs) %>%
    group_by(cluster_id) %>%
    summarise(
      c_start = min(frag_mid),
      c_end = max(frag_mid),
      span = c_end - c_start,
      start_strand = display_strand[which.min(frag_mid)],
      end_strand = display_strand[which.max(frag_mid)],
      start_motif_id = motif_id[which.min(frag_mid)],
      end_motif_id = motif_id[which.max(frag_mid)],
      start_has_motif = !is.na(start_strand) & start_strand != "Conflict",
      end_has_motif = !is.na(end_strand) & end_strand != "Conflict"
    )
  
  motif_classified <- cluster_motif_analysis %>%
    left_join(cluster_stats_motif, by = "cluster_id") %>%
    mutate(
      category = case_when(
        is_convergent ~ "Convergent",
        start_has_motif & start_strand == "+" ~ "Forward Anchored",
        end_has_motif & end_strand == "-" ~ "Reverse Anchored",
        TRUE ~ "Other with motif"
      ),
      sort_key = case_when(
        category == "Forward Anchored" ~ ifelse(is.na(start_motif_id), c_start, start_motif_id),
        category == "Reverse Anchored" ~ ifelse(is.na(end_motif_id), c_end, end_motif_id),
        TRUE ~ c_start
      )
    )
  
  # 2. Clusters without motifs - including Loading Analysis
  clusters_without_motifs <- df %>%
    filter(!cluster_id %in% clusters_with_motifs) %>%
    pull(cluster_id) %>%
    unique()
  
  # Calculate cluster statistics for non-motif clusters
  cluster_stats_no_motif <- df %>%
    filter(cluster_id %in% clusters_without_motifs) %>%
    group_by(cluster_id) %>%
    summarise(
      c_start = min(frag_mid),
      c_end = max(frag_mid),
      span = c_end - c_start,
      
      n_frags_left = sum(frag_mid < LOADING_START),
      n_frags_right = sum(frag_mid > LOADING_END),
      n_frags_in_loading = sum(in_loading),
      
      extends_past_left_anchor = ifelse(!is.null(anchor_left), c_start < anchor_left, FALSE),
      extends_past_right_anchor = ifelse(!is.null(anchor_right), c_end > anchor_right, FALSE)
    )
  
  # Identify two-sided candidates
  twosided_candidates <- cluster_stats_no_motif %>%
    filter(
      n_frags_left >= 1,
      n_frags_right >= 1,
      n_frags_in_loading == 0,
      !extends_past_left_anchor,
      !extends_past_right_anchor
    ) %>%
    arrange(span)
  
  # Progressive filtering for two-sided
  twosided_selected <- data.frame()
  
  if(nrow(twosided_candidates) > 0) {
    message(sprintf(">>> Found %d two-sided candidates", nrow(twosided_candidates)))
    
    for(i in 1:nrow(twosided_candidates)) {
      current <- twosided_candidates[i, ]
      
      if(nrow(twosided_selected) == 0) {
        twosided_selected <- rbind(twosided_selected, current)
      } else {
        prev_max_left <- min(twosided_selected$c_start)
        prev_max_right <- max(twosided_selected$c_end)
        
        if(current$c_start < prev_max_left && current$c_end > prev_max_right) {
          twosided_selected <- rbind(twosided_selected, current)
        }
      }
    }
    message(sprintf(">>> After progressive filtering: %d two-sided retained", nrow(twosided_selected)))
  }
  
  # Mark selected two-sided clusters
  twosided_cluster_ids <- if(nrow(twosided_selected) > 0) twosided_selected$cluster_id else c()
  
  # Classify remaining non-motif clusters (including one-sided loading)
  no_motif_classified <- cluster_stats_no_motif %>%
    mutate(
      category = case_when(
        cluster_id %in% twosided_cluster_ids ~ "Two-sided at loading",
        
        # One-sided at loading
        n_frags_in_loading > 0 & 
          ((n_frags_left > 0 & n_frags_right == 0) | (n_frags_left == 0 & n_frags_right > 0)) ~ "One-sided at loading",
        
        TRUE ~ "Other without motif"
      ),
      extends_right = n_frags_right > 0,
      sort_key = span
    )
  
  # Merge all categories
  all_categories <- bind_rows(
    motif_classified %>% 
      filter(category %in% c("Convergent", "Forward Anchored", "Reverse Anchored")) %>%
      select(cluster_id, category, span, sort_key) %>%
      mutate(extends_right = NA),
    
    no_motif_classified %>% 
      filter(category %in% c("One-sided at loading", "Two-sided at loading")) %>%
      select(cluster_id, category, span, extends_right, sort_key)
  )
  
  return(all_categories)
}

# ==============================================================================
# 7. Layout Assignment with All Categories
# ==============================================================================
assign_cluster_positions_all_categories <- function(df, anchor_left, anchor_right) {
  if(is.null(df) || nrow(df) == 0) return(NULL)
  
  # Filter fragments in ROI
  valid_ids <- df %>% 
    filter(in_display_region) %>% 
    group_by(cluster_id) %>% 
    summarize(n_in_roi = n()) %>% 
    filter(n_in_roi >= MIN_FRAGS_IN_ROI) %>% 
    pull(cluster_id)
  
  df <- df %>% filter(cluster_id %in% valid_ids)
  if(nrow(df) == 0) return(NULL)
  
  # Filter boundaries
  cluster_boundaries <- df %>%
    group_by(cluster_id) %>%
    summarize(leftmost = min(frag_mid), rightmost = max(frag_mid)) %>%
    filter(leftmost >= (START + EDGE_BUFFER) & rightmost <= (END - EDGE_BUFFER)) %>%
    pull(cluster_id)
  
  df <- df %>% filter(cluster_id %in% cluster_boundaries)
  if(nrow(df) == 0) return(NULL)
  
  # Classify with progressive filtering for two-sided
  cluster_categories <- classify_all_with_twosided_progressive(df, anchor_left, anchor_right)
  if(is.null(cluster_categories) || nrow(cluster_categories) == 0) return(NULL)
  
  # Sort each category
  forward_anchored <- cluster_categories %>%
    filter(category == "Forward Anchored") %>%
    arrange(sort_key, span)
  
  reverse_anchored <- cluster_categories %>%
    filter(category == "Reverse Anchored") %>%
    arrange(desc(sort_key), span)
  
  convergent <- cluster_categories %>%
    filter(category == "Convergent") %>%
    arrange(span)
  
  one_sided_loading <- cluster_categories %>%
    filter(category == "One-sided at loading") %>%
    arrange(desc(extends_right), span)
  
  two_sided_loading <- cluster_categories %>%
    filter(category == "Two-sided at loading") %>%
    arrange(span)
  
  # Calculate layout
  n_groups <- sum(c(
    nrow(forward_anchored) > 0,
    nrow(reverse_anchored) > 0,
    nrow(convergent) > 0,
    nrow(one_sided_loading) > 0,
    nrow(two_sided_loading) > 0
  ))
  
  total_clusters <- nrow(forward_anchored) + nrow(reverse_anchored) + nrow(convergent) + 
                   nrow(one_sided_loading) + nrow(two_sided_loading)
  
  gap_size <- max(2, round(total_clusters * GAP_SIZE_PERCENT))
  
  total_height <- total_clusters + ((n_groups - 1) * gap_size)
  cursor <- total_height
  
  add_group <- function(grp_df, group_label, cur_h) {
    if(nrow(grp_df) == 0) return(list(df=data.frame(), lbl=data.frame(), h=cur_h))
    mid_point <- cur_h - (nrow(grp_df) / 2)
    lbl <- data.frame(
      label = paste0(group_label, "\nn=", nrow(grp_df)), 
      y = mid_point
    )
    grp_df$y_idx <- seq(cur_h, cur_h - nrow(grp_df) + 1)
    return(list(
      df = grp_df %>% select(cluster_id, y_idx, category), 
      lbl = lbl, 
      h = cur_h - nrow(grp_df)
    ))
  }
  
  results <- list()
  labels_list <- list()
  
  if(nrow(forward_anchored) > 0) {
    res <- add_group(forward_anchored, "Forward\nAnchored", cursor)
    results <- c(results, list(res$df))
    labels_list <- c(labels_list, list(res$lbl))
    cursor <- res$h - gap_size
  }
  
  if(nrow(reverse_anchored) > 0) {
    res <- add_group(reverse_anchored, "Reverse\nAnchored", cursor)
    results <- c(results, list(res$df))
    labels_list <- c(labels_list, list(res$lbl))
    cursor <- res$h - gap_size
  }
  
  if(nrow(convergent) > 0) {
    res <- add_group(convergent, "Convergent", cursor)
    results <- c(results, list(res$df))
    labels_list <- c(labels_list, list(res$lbl))
    cursor <- res$h - gap_size
  }
  
  if(nrow(one_sided_loading) > 0) {
    res <- add_group(one_sided_loading, "One-sided\nat loading", cursor)
    results <- c(results, list(res$df))
    labels_list <- c(labels_list, list(res$lbl))
    cursor <- res$h - gap_size
  }
  
  if(nrow(two_sided_loading) > 0) {
    res <- add_group(two_sided_loading, "Two-sided\nat loading\n(progressive)", cursor)
    results <- c(results, list(res$df))
    labels_list <- c(labels_list, list(res$lbl))
  }
  
  layout_df <- bind_rows(results)
  labels_df <- bind_rows(labels_list)
  
  if(nrow(layout_df) == 0) return(NULL)
  
  return(list(
    data = df %>% inner_join(layout_df, by="cluster_id"),
    labels = labels_df
  ))
}

# ==============================================================================
# 8. Plotting Final View
# ==============================================================================
plot_final_view <- function(layout_res, anchor_left, anchor_right) {
  df <- layout_res$data
  labels <- layout_res$labels
  coord_label <- paste0(CHR_BW, ":", comma(START), "-", comma(END))
  
  SIZE_DOT    <- 0.5   
  SIZE_SQUARE <- 1.0   
  
  p <- ggplot() +
    # Highlight Loading Site
    annotate("rect", xmin = LOADING_START, xmax = LOADING_END, 
             ymin = -Inf, ymax = Inf, fill = LOADING_COLOR, alpha = 0.5) +
    
    geom_line(data = df, aes(x = frag_mid, y = y_idx, group = cluster_id), 
              linewidth = LINE_WIDTH, color = COLOR_LINE) +
    geom_point(data = df %>% filter(is.na(display_strand) | display_strand == "Conflict"),
               aes(x = frag_mid, y = y_idx), size = SIZE_DOT, color = COLOR_DOT, shape = 19) +
    geom_point(data = df %>% filter(display_strand %in% c("+", "-")),
               aes(x = frag_mid, y = y_idx, color = display_strand), size = SIZE_SQUARE, shape = 15) +
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

# ==============================================================================
# 9. Execution
# ==============================================================================
message("1. Determining Anchor Positions...")
anchors <- get_anchor_positions(MOTIF_FILE, CHR_BW, START, END)
ANCHOR_LEFT <- anchors$left
ANCHOR_RIGHT <- anchors$right

message("2. Generating Top Motif Track...")
p_motif <- plot_motif_track(MOTIF_FILE, CHR_BW, START, END)

message("3. Generating Signal Tracks...")
signal_plots <- list()
for(name in names(TRACKS_INFO)) {
  info <- TRACKS_INFO[[name]]
  message(sprintf("   Processing %s...", name))
  signal_plots[[name]] <- plot_signal_track(info[1], name, info[2], CHR_BW, START, END)
}

message("4. Processing Fragments...")
df_tagged <- process_fragments_and_motifs(ANCHOR_LEFT, ANCHOR_RIGHT)

if(!is.null(df_tagged) && !is.null(ANCHOR_LEFT) && !is.null(ANCHOR_RIGHT)) {
  message("5. Computing Layout with All Categories (Two-sided uses Progressive Filtering)...")
  layout_res <- assign_cluster_positions_all_categories(df_tagged, ANCHOR_LEFT, ANCHOR_RIGHT)
  p_frag <- if(!is.null(layout_res)) {
    plot_final_view(layout_res, ANCHOR_LEFT, ANCHOR_RIGHT)
  } else {
    ggplot() + theme_void() + labs(title="No Two-sided Complexes Found")
  }
} else {
  p_frag <- ggplot() + theme_void() + labs(title="No Data or Anchors Not Found")
}

message("6. Assembling Final Plot...")

all_plots_list <- c(list(p_motif), signal_plots, list(p_frag))
h_ratios <- c(0.6, rep(0.8, length(signal_plots)), 10)

final_plot <- wrap_plots(all_plots_list, ncol = 1) + 
  plot_layout(heights = h_ratios)

message("Saving plot to ChIA_Drop_All_Categories_Progressive.pdf...")
ggsave("ChIA_Drop_All_Categories_Progressive.pdf", 
       plot = final_plot, 
       width = 12, 
       height = 22, 
       units = "in",
       limitsize = FALSE)

print(final_plot)
message("Done!")
message(sprintf(">>> Final output: Two-sided complexes with progressive filtering"))