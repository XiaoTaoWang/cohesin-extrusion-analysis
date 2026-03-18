# ==============================================================================
# Script Name: plot_rnapii_ep_interactions_filtered.R
# Function:    RNAPII ChIA-Drop Visualization (Filtered for P/E overlaps)
# Date:        2026-01-15
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
# Plot region settings
CHR_BW <- "chr5"
START  <- 116500000
END    <- 117000000

# Filtering and display parameters
MIN_FRAGS_IN_ROI   <- 2            
GAP_SIZE_PERCENT   <- 0.02         
EDGE_BUFFER        <- 500  

# === Fragment Overlap Settings ===
EXTEND_SIZE        <- 4000 
FRAG_SIZE_P_E      <- 1.0  
FRAG_SIZE_NONE     <- 0.5  

# === Colors ===
COLOR_PROMOTER     <- "#E31A1C" # Red
COLOR_ENHANCER     <- "#FF7F00" # Orange
COLOR_NONE         <- "#33A02C" # Green
COLOR_LINE         <- "#acada2"
LINE_WIDTH         <- 0.1

LABEL_SIGNAL       <- "RNAPII Signal"

# === File Paths ===
# [UPDATED] Gene Annotation
GENE_BED_FILE <- paste0(GENOME_ROOT, "/hg38/hg38.gene.bed")

# Annotations
ENHANCER_FILE <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between/GM12878_active_enhancers.bed")
PROMOTER_FILE <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between/GM12878_active_promoters.bed")

# Data
RDS_FILE      <- paste0(DATA_ROOT, "/ChIA-Drop/GSE158897/GSE158897_GM12878-RNAPII-pooledv2_comp.rds")
BW_FILE       <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between/GM12878_RNAPII_ChIA-PET_GSE158897_LHG0035N_0035V_0045V.for.BROWSER.sorted.bw")

message(sprintf(">>> Plotting Region: %s:%s-%s", CHR_BW, comma(START), comma(END)))

if (!exists("G_FRAGMENTS_CACHE")) { G_FRAGMENTS_CACHE <- NULL }

# ==============================================================================
# 2. Theme & Helper Functions
# ==============================================================================
clean_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(color = "black", size = 10, margin = margin(t = 5)),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_blank(),         
    axis.line.y = element_blank(),     
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(color = "black", size = 10, hjust = 1, face = "bold"),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 10, face = "bold", margin = margin(r = 10)),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none",     
    plot.margin = margin(t = 2, b = 2, r = 20, l = 20)   
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

# 3.1 Gene Track (Custom BED Parser)
plot_gene_track <- function(chrom, start_pos, end_pos) {
  message("   Generating Gene Track (Custom BED)...")
  p_void <- ggplot() + clean_theme + labs(y = "Genes") + theme(axis.line.x=element_blank(), axis.text=element_blank())

  if (!file.exists(GENE_BED_FILE)) {
    message("   [WARN] Gene file not found.")
    return(p_void)
  }

  df <- tryCatch({
    d <- read.table(GENE_BED_FILE, header=FALSE, stringsAsFactors=FALSE, sep="\t", quote="")
    # Format provided in your notes:
    # V1=chr, V2=start, V3=end, V4=strand, V5=ID, V6=Symbol, V7=Type
    data.frame(
      chr    = d[,1],
      start  = d[,2] + 1, # 0-based to 1-based
      end    = d[,3],
      strand = d[,4],     # column 4 = strand
      symbol = d[,6]      # column 6 = gene symbol
    )
  }, error = function(e) { message(e); NULL })

  if (is.null(df) || nrow(df) == 0) return(p_void)

  if (grepl("chr", chrom)) { 
    if (!any(grepl("chr", df$chr))) df$chr <- paste0("chr", df$chr) 
  } else { 
    df$chr <- gsub("chr", "", df$chr) 
  }

  df_genes <- df %>% filter(chr == chrom, end >= start_pos, start <= end_pos)
  
  if (nrow(df_genes) == 0) return(p_void)

  # Layout: avoid overlapping gene labels
  df_genes <- df_genes %>% arrange(start)
  df_genes$y <- 1
  if(nrow(df_genes) > 1) {
    for(i in 2:nrow(df_genes)) {
      if(df_genes$start[i] < df_genes$end[i-1]) {
        df_genes$y[i] <- 3 - df_genes$y[i-1] 
      } else {
        df_genes$y[i] <- 1
      }
    }
  }
  
  df_genes$mid <- (pmax(df_genes$start, start_pos) + pmin(df_genes$end, end_pos)) / 2

  ggplot() +
    # Gene body
    geom_segment(data = df_genes, aes(x = pmax(start, start_pos), xend = pmin(end, end_pos), 
                                      y = y, yend = y), color = "#377eb8", linewidth = 2) +
    # Gene label
    geom_text(data = df_genes, aes(x = mid, y = y + 0.45, label = symbol), size = 3, fontface = "italic") +
    # Direction arrow
    geom_segment(data = df_genes, 
                 aes(x = ifelse(strand == "+", pmin(end, end_pos)-200, pmax(start, start_pos)+200),
                     xend = ifelse(strand == "+", pmin(end, end_pos), pmax(start, start_pos)),
                     y = y, yend = y),
                 arrow = arrow(length = unit(0.2, "cm"), type = "open"), 
                 linewidth = 0.5, color = "black") +
    scale_y_continuous(limits = c(0.5, 3)) +
    coord_cartesian(xlim = c(start_pos, end_pos), expand = FALSE) +
    labs(y = "Genes") + clean_theme +
    theme(axis.line.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
}

# 3.2 Annotation Track
plot_bed_track <- function(bed_file, track_name, color_fill, chrom, start_pos, end_pos) {
  p_void <- ggplot() + clean_theme + labs(y = track_name) + theme(axis.line.x=element_blank(), axis.text=element_blank())
  
  if (!file.exists(bed_file)) return(p_void)
  
  df <- tryCatch({ 
    d <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
    data.frame(chr = d[,1], start = d[,2] + 1, end = d[,3])
  }, error = function(e) NULL)
  
  if (is.null(df) || nrow(df) == 0) return(p_void)
  
  if (grepl("chr", chrom)) { 
    if (!any(grepl("chr", df$chr))) df$chr <- paste0("chr", df$chr) 
  } else { 
    df$chr <- gsub("chr", "", df$chr) 
  }
  
  df <- df %>% filter(chr == chrom, end >= start_pos, start <= end_pos)
  
  ggplot(df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.8), fill = color_fill, color = NA) +
    coord_cartesian(xlim = c(start_pos, end_pos), ylim = c(0, 1), expand = FALSE) +
    labs(y = track_name) + clean_theme +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), plot.margin = margin(t = 2, b = 2, r = 20, l = 20))
}

# 3.3 Signal Track
plot_signal_track <- function(bw_file, track_name, fill_color, chrom, start_pos, end_pos) {
  p_void <- ggplot() + clean_theme + labs(y = track_name) +     
    theme(axis.line.x=element_blank(), axis.text.y = element_blank())
  
  if (!file.exists(bw_file)) return(p_void)
  
  roi_gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start_pos, end = end_pos))
  bw_gr <- tryCatch({ import(bw_file, which = roi_gr) }, error = function(e) NULL)
  if (is.null(bw_gr)) return(p_void)
  
  bw_df <- as.data.frame(bw_gr)
  if(nrow(bw_df)==0) bw_df <- data.frame(start=start_pos,end=start_pos,score=0)
  
  raw_max <- max(bw_df$score, na.rm=TRUE)
  effective_max <- if(raw_max > 0) raw_max * 1.1 else 1
  rounded_max <- round_smart(effective_max)
  
  ggplot(bw_df, aes(xmin = start, xmax = end, ymin = 0, ymax = score)) +
    geom_rect(fill = fill_color, color = fill_color, linewidth = 0.1) +     
    scale_y_continuous(breaks = c(0, rounded_max), labels = c("0", rounded_max)) +
    coord_cartesian(xlim = c(start_pos, end_pos), ylim = c(0, rounded_max), expand = FALSE) +
    labs(y = track_name) + clean_theme +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(t = 2, b = 2, r = 20, l = 20))
}

# ==============================================================================
# 4. Processing Fragments (Annotate P/E & Filter)
# ==============================================================================
process_rnapii_interactions <- function() {
  roi_gr <- GRanges(seqnames = CHR_BW, ranges = IRanges(start = START, end = END))
  
  if (!is.null(G_FRAGMENTS_CACHE)) { 
    df_raw <- G_FRAGMENTS_CACHE 
  } else {     
    if(!file.exists(RDS_FILE)) stop(paste("RDS not found:", RDS_FILE))
    df_raw <- readRDS(RDS_FILE)
    G_FRAGMENTS_CACHE <<- df_raw 
  }
  
  # Subset to ROI
  if(inherits(df_raw, "GRanges")) {
    df_gr <- suppressWarnings(subsetByOverlaps(df_raw, roi_gr))
    df <- as.data.frame(df_gr)
  } else {
    df <- df_raw %>% filter(seqnames == as.character(seqnames(roi_gr)), end >= start(roi_gr), start <= end(roi_gr))
  }
  
  if(nrow(df) == 0) return(NULL)
  
  if (!"cluster_id" %in% colnames(df)) {
    if("CBMB" %in% colnames(df)) df$cluster_id <- df$CBMB else df$cluster_id <- 1:nrow(df)
  }
  
  df$frag_mid <- (df$start + df$end) / 2
  
  # 1. Annotate Fragments
  ext_gr <- GRanges(
    seqnames = df$seqnames,
    ranges   = IRanges(start = floor(df$frag_mid - EXTEND_SIZE), 
                       end   = floor(df$frag_mid + EXTEND_SIZE))
  )
  
  read_bed_as_gr <- function(f) {
    if(!file.exists(f)) return(GRanges())
    d <- read.table(f, header=F, stringsAsFactors=F)
    GRanges(seqnames = d[,1], ranges = IRanges(start = d[,2]+1, end = d[,3]))
  }
  
  promoters_gr <- read_bed_as_gr(PROMOTER_FILE)
  enhancers_gr <- read_bed_as_gr(ENHANCER_FILE)
  
  ov_p <- findOverlaps(ext_gr, promoters_gr)
  ov_e <- findOverlaps(ext_gr, enhancers_gr)
  
  df$type <- "None"
  df$type[unique(queryHits(ov_e))] <- "Enhancer"
  df$type[unique(queryHits(ov_p))] <- "Promoter" 
  
  # 2. Filter Clusters: Must have at least one P or E
  # Identify clusters that contain at least one Promoter or Enhancer
  interesting_clusters <- df %>%
    group_by(cluster_id) %>%
    summarise(has_feature = any(type %in% c("Promoter", "Enhancer"))) %>%
    filter(has_feature) %>%
    pull(cluster_id)
  
  message(sprintf("   Clusters filtering: Total %d -> Keeping %d (must overlap P/E)", 
                  length(unique(df$cluster_id)), length(interesting_clusters)))
  
  df <- df %>% filter(cluster_id %in% interesting_clusters)
  
  if(nrow(df) == 0) return(NULL)

  # 3. Visualization Filter (Keep fragments in ROI)
  df <- df %>% filter(frag_mid >= START & frag_mid <= END)
  
  # Filter small clusters (visual)
  valid_ids <- df %>% 
    group_by(cluster_id) %>% 
    summarize(n = n()) %>% 
    filter(n >= MIN_FRAGS_IN_ROI) %>% 
    pull(cluster_id)
  
  df <- df %>% filter(cluster_id %in% valid_ids)
  
  # Edge Buffer
  cluster_boundaries <- df %>%
    group_by(cluster_id) %>%
    summarize(
      leftmost = min(frag_mid),
      rightmost = max(frag_mid)
    ) %>%
    filter(leftmost >= (START + EDGE_BUFFER) & rightmost <= (END - EDGE_BUFFER)) %>%
    pull(cluster_id)
  
  df <- df %>% filter(cluster_id %in% cluster_boundaries)
  
  if(nrow(df) == 0) return(NULL)
  return(df)
}

# ==============================================================================
# 5. Sorting & Layout
# ==============================================================================
assign_rnapii_positions <- function(df) {
  if(is.null(df) || nrow(df) == 0) return(NULL)
  
  cluster_stats <- df %>%
    group_by(cluster_id) %>%
    summarise(
      has_P = any(type == "Promoter"),
      has_E = any(type == "Enhancer"),
      min_P_pos = ifelse(any(type == "Promoter"), min(frag_mid[type == "Promoter"]), NA),
      min_E_pos = ifelse(any(type == "Enhancer"), min(frag_mid[type == "Enhancer"]), NA),
      c_start = min(frag_mid)
    ) %>%
    mutate(
      category = case_when(
        has_P & has_E ~ "E-P",
        has_P & !has_E ~ "P-P", 
        !has_P & has_E ~ "E-E", 
        TRUE ~ "Other" # Should be empty due to filter, but kept for safety
      )
    )
  
  # Sort Priority: P position > E position > Start
  cluster_stats <- cluster_stats %>%
    mutate(sort_key = coalesce(min_P_pos, min_E_pos, c_start))
  
  df_ep <- cluster_stats %>% filter(category == "E-P") %>% arrange(sort_key)
  df_pp <- cluster_stats %>% filter(category == "P-P") %>% arrange(sort_key)
  df_ee <- cluster_stats %>% filter(category == "E-E") %>% arrange(sort_key)
  df_ot <- cluster_stats %>% filter(category == "Other") %>% arrange(sort_key)
  
  message(sprintf("   Classification: E-P: %d, P-P: %d, E-E: %d, Other: %d",
                  nrow(df_ep), nrow(df_pp), nrow(df_ee), nrow(df_ot)))
  
  gap_size <- max(2, round(nrow(cluster_stats) * GAP_SIZE_PERCENT))
  cursor <- nrow(cluster_stats) + (3 * gap_size) 
  
  layout_df <- data.frame()
  labels_df <- data.frame()
  
  process_group <- function(sub_df, label_text) {
    if(nrow(sub_df) == 0) return(NULL)
    y_idx <- seq(cursor, cursor - nrow(sub_df) + 1)
    mid_y <- mean(y_idx)
    res <- sub_df %>% select(cluster_id) %>% mutate(y_idx = y_idx)
    return(list(res = res, label = data.frame(label=label_text, y=mid_y), count = nrow(sub_df)))
  }
  
  if(nrow(df_ep) > 0) {
    res <- process_group(df_ep, "E-P")
    layout_df <- rbind(layout_df, res$res)
    labels_df <- rbind(labels_df, res$label)
    cursor <- cursor - res$count - gap_size
  }
  
  if(nrow(df_pp) > 0) {
    res <- process_group(df_pp, "P-P")
    layout_df <- rbind(layout_df, res$res)
    labels_df <- rbind(labels_df, res$label)
    cursor <- cursor - res$count - gap_size
  }
  
  if(nrow(df_ee) > 0) {
    res <- process_group(df_ee, "E-E")
    layout_df <- rbind(layout_df, res$res)
    labels_df <- rbind(labels_df, res$label)
    cursor <- cursor - res$count - gap_size
  }
  
  if(nrow(df_ot) > 0) {
    res <- process_group(df_ot, "Other")
    layout_df <- rbind(layout_df, res$res)
  }
  
  if(nrow(layout_df) == 0) return(NULL)
  
  df_final <- df %>% inner_join(layout_df, by="cluster_id")
  return(list(data = df_final, labels = labels_df))
}

# ==============================================================================
# 6. Plotting
# ==============================================================================
plot_interactions_view <- function(layout_res) {
  df <- layout_res$data
  labels <- layout_res$labels
  
  coord_label <- paste0(CHR_BW, ":", comma(START), "-", comma(END))
  
  p <- ggplot() +
    geom_line(data = df, aes(x = frag_mid, y = y_idx, group = cluster_id),
              linewidth = LINE_WIDTH, color = COLOR_LINE) +
    
    geom_point(data = df %>% filter(type == "None"),
               aes(x = frag_mid, y = y_idx),
               size = FRAG_SIZE_NONE, color = COLOR_NONE, shape = 19) + 
    
    geom_point(data = df %>% filter(type == "Promoter"),
               aes(x = frag_mid, y = y_idx),
               size = FRAG_SIZE_P_E, color = COLOR_PROMOTER, shape = 15) + 
    
    geom_point(data = df %>% filter(type == "Enhancer"),
               aes(x = frag_mid, y = y_idx),
               size = FRAG_SIZE_P_E, color = COLOR_ENHANCER, shape = 15) + 
    
    scale_x_continuous(limits = c(START, END), breaks = breaks_extended(n = 5), labels = comma, expand = c(0,0)) +
    labs(y = NULL, x = coord_label) +     
    coord_cartesian(clip = "off") +     
    clean_theme
  
  if(!is.null(labels) && nrow(labels) > 0) {
    p <- p + scale_y_continuous(breaks = labels$y, labels = labels$label, expand = expansion(mult = c(0.05, 0)))
  } else {
    p <- p + scale_y_continuous(breaks = NULL, expand = expansion(mult = c(0.05, 0)))
  }
  
  return(p)
}

# ==============================================================================
# 7. Execution
# ==============================================================================
message("1. Generating Gene Track...")
p_gene <- plot_gene_track(CHR_BW, START, END)

message("2. Generating Annotation Tracks...")
p_promoter <- plot_bed_track(PROMOTER_FILE, "Promoters", COLOR_PROMOTER, CHR_BW, START, END)
p_enhancer <- plot_bed_track(ENHANCER_FILE, "Enhancers", COLOR_ENHANCER, CHR_BW, START, END)

message("3. Generating Signal Track...")
p_signal <- plot_signal_track(BW_FILE, LABEL_SIGNAL, "#1f77b4", CHR_BW, START, END)

message("4. Processing RNAPII Interactions...")
df_tagged <- process_rnapii_interactions()

if(!is.null(df_tagged)) {
  message("5. Computing Layout...")
  layout_res <- assign_rnapii_positions(df_tagged)
  
  if(!is.null(layout_res)) {
    p_inter <- plot_interactions_view(layout_res)
  } else {
    p_inter <- ggplot() + theme_void() + labs(title="No Interactions Found")
  }
} else {
  p_inter <- ggplot() + theme_void() + labs(title="No Valid Data")
}

message("6. Assembling Final Plot...")
final_plot <- p_gene / p_promoter / p_enhancer / p_signal / p_inter + 
  plot_layout(heights = c(0.5, 0.2, 0.2, 0.6, 5))

ggsave("RNAPII_ChIA_Drop_EP_View_Filtered.pdf", 
       plot = final_plot,
       width = 11,
       height = 15,
       units = "in") 

print(final_plot)