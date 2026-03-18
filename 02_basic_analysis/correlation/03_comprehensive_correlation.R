#!/usr/bin/env Rscript

# Set working directory

# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = dirname(dirname(normalizePath(sys.frame(1)$ofile))))
DATA_ROOT    <- Sys.getenv("DATA_ROOT", unset = dirname(PROJECT_ROOT))
GENOME_ROOT  <- Sys.getenv("GENOME_ROOT", unset = file.path(dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

setwd(paste0(DATA_ROOT, "/ChIA-Drop/4.correlation_NIPBL"))
cat(sprintf("Working directory: %s\n", getwd()))

# Load required libraries
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
  library(parallel)
  library(data.table)
  library(gridExtra)
  library(MASS)
  library(viridis)
})

# ==============================================================================
# 1. Configuration
# ==============================================================================

# Number of cores for parallel processing

n_cores <- 6

# Directories
DIR_MAIN  <- paste0(DATA_ROOT, "/ChIA-Drop/2.CTCF_motif_fragment_between")
DIR_OTHER <- paste0(DATA_ROOT, "/ChIA-Drop/GSE158897")
DIR_CHIP  <- paste0(DATA_ROOT, "/public_databases/minji/minji_2024_chip-seq")

# Output Directories
OUTPUT_DIR <- "correlation_plots"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

CACHE_DIR <- "correlation_plots/cache"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# Define all tracks
TRACKS_INFO <- list(
  "CTCF_ChIA-PET"     = file.path(DIR_MAIN, "GM12878_CTCF_ChIA-PET_GSE158897_LHG0052H.for.BROWSER.sorted.bw"),
  "CTCF_ChIP-seq"     = file.path(DIR_CHIP, "GM12878-CTCF-ChIP-seq_GSM5983423_CHG0031.q30.nr.sorted.bw"),
  "CTCF_ChIA-Drop"    = file.path(DIR_MAIN, "ChIA-Drop_GSE158897_GM12878-CTCF-pooled.bw"),
  "Cohesin_ChIA-PET"  = file.path(DIR_MAIN, "GM12878_cohesin_ChIA-PET_GSE158897_LHG0051H_0104V.for.BROWSER.sorted.bw"),
  "Cohesin_ChIP-seq"  = file.path(DIR_CHIP, "GM12878-RAD21-ChIP-seq_GSM5983424_CHG0034.q30.nr.sorted.bw"),
  "Cohesin_ChIA-Drop" = file.path(DIR_MAIN, "ChIA-Drop_GSE158897_GM12878-cohesin-pooled.bw"),
  "NIPBL"             = file.path(DIR_OTHER, "GM12878_NIPBL_chip-seq_GSE158897_CHG0030.q30.nr.sorted.bw"),
  "WAPL"              = file.path(DIR_OTHER, "GM12878_WAPL_chip-seq_GSE158897_CHG0032.q30.nr.sorted.bw"),
  "RNAPII"            = file.path(DIR_MAIN, "GM12878_RNAPII_ChIA-PET_GSE158897_LHG0035N_0035V_0045V.for.BROWSER.sorted.bw"),
  "H3K27ac"           = file.path(DIR_OTHER, "GM12878_h3k27ac_chip-seq_ENCFF340JIF.bigWig"),
  "H3K4me1"           = file.path(DIR_OTHER, "GM12878_h3k4me1_chip-seq_ENCFF831ZHL.bigWig")
)

# Reference tracks to use as X-axis
X_AXIS_TRACKS <- c("NIPBL", "WAPL", "RNAPII", "H3K27ac", "H3K4me1")

# Input bedpe file
bedpe_file <- paste0(DATA_ROOT, "/ChIA-Drop/3.anchoring_middle/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot.sorted.convergent.loops")

# Metrics to calculate
METRICS <- c("mean", "max", "sum")

# ==============================================================================
# 2. Data Loading and Processing
# ==============================================================================

cat("\nReading BEDPE file...\n")
bedpe <- fread(bedpe_file, header = FALSE)
colnames(bedpe)[1:6] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
cat(sprintf("Total interactions: %d\n", nrow(bedpe)))

# Calculate Loop Metadata (Distance & PET Count)
bedpe$dist <- abs((bedpe$start1 + bedpe$end1)/2 - (bedpe$start2 + bedpe$end2)/2)
# Assuming V7 is PET count/Score based on standard BEDPE and filename "PETcnt"
if(ncol(bedpe) >= 7) bedpe$score <- bedpe$V7 else bedpe$score <- NA

cat(sprintf("Calculated loop distances and extracted scores (if available).\n"))

# Create GRanges for anchor1 and anchor2
anchor1 <- GRanges(seqnames = bedpe$chr1, ranges = IRanges(start = bedpe$start1, end = bedpe$end1))
anchor2 <- GRanges(seqnames = bedpe$chr2, ranges = IRanges(start = bedpe$start2, end = bedpe$end2))

# Combine both anchors
all_regions <- c(anchor1, anchor2)
cat(sprintf("Total regions (both anchors): %d\n", length(all_regions)))

# Function to calculate signal
calculate_signal <- function(bw_file, regions, method = "mean", track_name = "", use_parallel = TRUE) {
  safe_track_name <- gsub("[^[:alnum:]]", "_", track_name)
  cache_file <- file.path(CACHE_DIR, sprintf("%s_%s.rds", safe_track_name, method))
  
  if (file.exists(cache_file)) {
    cat(sprintf("    [Cache Hit] %s (%s)\n", track_name, method))
    return(readRDS(cache_file))
  }
  
  cat(sprintf("    [Computing] %s (%s)...\n", track_name, method))
  
  if (!file.exists(bw_file)) {
    warning(sprintf("File not found: %s", bw_file))
    signal <- rep(NA, length(regions))
    saveRDS(signal, cache_file)
    return(signal)
  }
  
  tryCatch({
    bw <- import(bw_file, format = "BigWig")
    overlaps <- findOverlaps(regions, bw)
    
    if (use_parallel && n_cores > 1) {
      chunk_size <- ceiling(length(regions) / (n_cores * 4))
      chunks <- split(1:length(regions), ceiling(seq_along(1:length(regions)) / chunk_size))
      
      chunk_results <- mclapply(chunks, function(chunk_indices) {
        sapply(chunk_indices, function(i) {
          idx <- queryHits(overlaps) == i
          if (sum(idx) == 0) return(0)
          
          region_width <- width(regions[i])
          overlapping_ranges <- bw[subjectHits(overlaps)[idx]]
          scores <- score(overlapping_ranges)
          widths <- width(overlapping_ranges)
          
          if (method == "mean") return(sum(scores * widths) / region_width)
          else if (method == "max") return(max(scores))
          else if (method == "sum") return(sum(scores * widths))
        })

      }, mc.cores = n_cores, mc.preschedule = FALSE)
      

      # Check if all chunks returned results
      chunk_lengths <- sapply(chunk_results, length)
      expected_lengths <- sapply(chunks, length)
      
      if (!all(chunk_lengths == expected_lengths)) {
        warning(sprintf("Parallel processing failed for %s, falling back to serial processing", track_name))
        use_parallel <- FALSE
      } else {
        signal <- unlist(chunk_results)
      }
    } 
    
    if (!use_parallel || n_cores <= 1) {
      # Serial processing
      signal <- sapply(1:length(regions), function(i) {
        idx <- queryHits(overlaps) == i
        if (sum(idx) == 0) return(0)
        region_width <- width(regions[i])
        overlapping_ranges <- bw[subjectHits(overlaps)[idx]]
        scores <- score(overlapping_ranges)
        widths <- width(overlapping_ranges)
        
        if (method == "mean") return(sum(scores * widths) / region_width)
        else if (method == "max") return(max(scores))
        else if (method == "sum") return(sum(scores * widths))
      })
    }
    
    if (length(signal) != length(regions)) {
      stop(sprintf("Signal length mismatch for %s: got %d, expected %d", track_name, length(signal), length(regions)))
    }
    
    saveRDS(signal, cache_file)
    return(signal)
  }, error = function(e) {
    warning(sprintf("Error processing %s: %s", bw_file, e$message))
    return(rep(NA, length(regions)))
  })
}

# ==============================================================================
# 3. Plotting Functions
# ==============================================================================

get_density <- function(x, y, n = 100) {
  valid <- is.finite(x) & is.finite(y)
  x_clean <- x[valid]
  y_clean <- y[valid]
  
  if(length(x_clean) < 10) return(rep(1, length(x)))
  
  dens <- MASS::kde2d(x = x_clean, y = y_clean, n = n)
  ix <- findInterval(x_clean, dens$x)
  iy <- findInterval(y_clean, dens$y)
  ii <- cbind(ix, iy)
  
  res <- rep(NA, length(x))
  res[valid] <- dens$z[ii]
  return(res)
}

create_scatter_plot <- function(df, x_col, y_col, metric, title_prefix) {
  if (!all(c(x_col, y_col) %in% names(df))) return(NULL)
  
  sub_df <- na.omit(df[, c(x_col, y_col)])
  names(sub_df) <- c("X", "Y")
  
  # Outlier Removal (Quantile 0.01% - 99.99%)
  qx <- quantile(sub_df$X, probs = c(0.0001, 0.9999), na.rm = TRUE)
  qy <- quantile(sub_df$Y, probs = c(0.0001, 0.9999), na.rm = TRUE)
  
  keep_idx <- sub_df$X >= qx[1] & sub_df$X <= qx[2] &
              sub_df$Y >= qy[1] & sub_df$Y <= qy[2]
  
  sub_df <- sub_df[keep_idx, ]
  
  n_points <- nrow(sub_df)
  if(n_points < 10) return(NULL)
  
  sub_df$density <- get_density(sub_df$X, sub_df$Y)
  
  r_p <- round(cor(sub_df$X, sub_df$Y, method = "pearson"), 2)
  r_s <- round(cor(sub_df$X, sub_df$Y, method = "spearman"), 2)
  
  p <- ggplot(sub_df, aes(x = X, y = Y)) +
    geom_point(aes(color = density), size = 0.5, alpha = 0.6) + 
    scale_color_gradientn(colours = c("#00008B", "#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000")) +
    labs(x = paste0(x_col, " (", metric, ")"), 
         y = paste0(y_col, " (", metric, ")"),
         title = paste0(title_prefix, "\nn=", n_points, " | P: ", r_p, " | S: ", r_s)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(color = "black", size = 8)
    )
  return(p)
}

create_validation_plot <- function(df, x_col, y_col, color_col, metric) {
  if (!all(c(x_col, y_col, color_col) %in% names(df))) return(NULL)
  
  sub_df <- na.omit(df[, c(x_col, y_col, color_col)])
  names(sub_df) <- c("X", "Y", "Color")
  
  # Outlier Removal (Same as main plot)
  qx <- quantile(sub_df$X, probs = c(0.0001, 0.9999), na.rm = TRUE)
  qy <- quantile(sub_df$Y, probs = c(0.0001, 0.9999), na.rm = TRUE)
  keep_idx <- sub_df$X >= qx[1] & sub_df$X <= qx[2] &
              sub_df$Y >= qy[1] & sub_df$Y <= qy[2]
  sub_df <- sub_df[keep_idx, ]
  
  if(nrow(sub_df) < 10) return(NULL)
  
  p <- ggplot(sub_df, aes(x = X, y = Y)) +
    geom_point(aes(color = Color), size = 0.5, alpha = 0.6) + 
    scale_color_viridis(option = "magma", name = color_col) +
    labs(x = paste0(x_col, " (", metric, ")"), 
         y = paste0(y_col, " (", metric, ")"),
         title = paste0(x_col, " vs ", y_col, "\nColored by ", color_col)) +
    theme_classic() +
    theme(legend.position = "right")
    
  return(p)
}

# ==============================================================================
# 4. Main Loop
# ==============================================================================

for (method in METRICS) {
  cat(sprintf("\n=== Processing Metric: %s ===\n", method))
  
  # 1. Calculate/Load Signals
  signal_data <- list()
  for (track_name in names(TRACKS_INFO)) {
    signal_data[[track_name]] <- calculate_signal(
      TRACKS_INFO[[track_name]], 
      all_regions, 
      method = method,
      track_name = track_name,
      use_parallel = TRUE
    )
  }
  
  df <- as.data.frame(signal_data)
  
  # Add Loop Metadata to df (Repeat for anchor1 and anchor2)
  # Note: all_regions = c(anchor1, anchor2), so we concatenate the metadata similarly
  # df$Log10_Distance <- log10(c(bedpe$dist, bedpe$dist) + 1)
  # df$PET_Count <- c(bedpe$score, bedpe$score)
  
  # 2. Generate Plots for each X-axis track
  for (x_track in X_AXIS_TRACKS) {
    if (!x_track %in% names(df)) {
      warning(sprintf("Reference track %s not found in data, skipping...", x_track))
      next
    }
    
    cat(sprintf("  Generating plots with X-axis: %s\n", x_track))
    
    plot_list <- list()
    y_tracks <- setdiff(names(df), x_track)
    
    for (y_track in y_tracks) {
      # Skip if all NA
      if (all(is.na(df[[y_track]]))) next
      
      p <- create_scatter_plot(
        df, x_track, y_track, method, 
        paste(x_track, "vs", y_track)
      )
      
      if (!is.null(p)) {
        plot_list[[y_track]] <- p
      }
    }
    
    if (length(plot_list) > 0) {
      # Arrange and Save
      ncol_set <- if(length(plot_list) > 6) 4 else 3
      
      output_pdf <- file.path(OUTPUT_DIR, sprintf("Correlation_%s_vs_Others_%s.pdf", x_track, method))
      
      pdf_width <- 4 * ncol_set
      pdf_height <- 3.5 * ceiling(length(plot_list) / ncol_set)
      
      ggsave(
        filename = output_pdf, 
        plot = arrangeGrob(grobs = plot_list, ncol = ncol_set, 
                           top = paste0("Correlation: ", x_track, " vs Others [", method, "]")), 
        width = pdf_width, height = pdf_height, limitsize = FALSE
      )
      cat(sprintf("    Saved: %s\n", output_pdf))
    }
    
    # Special Validation: NIPBL vs CTCF colored by H3K27ac/RNAPII to explain the "split"
    if (x_track == "NIPBL" && "CTCF_ChIP-seq" %in% names(df)) {
      val_plots <- list()
      
      # Define variables to use for coloring
      color_vars <- c("H3K27ac", "RNAPII", "Log10_Distance", "PET_Count")
      
      for (c_var in color_vars) {
        if (c_var %in% names(df)) 
          val_plots[[c_var]] <- create_validation_plot(df, "NIPBL", "CTCF_ChIP-seq", c_var, method)
      }
        
      if (length(val_plots) > 0) {
        val_pdf <- file.path(OUTPUT_DIR, sprintf("Validation_NIPBL_CTCF_Split_%s.pdf", method))
        ggsave(filename = val_pdf, plot = arrangeGrob(grobs = val_plots, ncol = 2),
               width = 12, height = 5)
        cat(sprintf("    Saved Validation Plot: %s\n", val_pdf))
      }
    }
  }
}
