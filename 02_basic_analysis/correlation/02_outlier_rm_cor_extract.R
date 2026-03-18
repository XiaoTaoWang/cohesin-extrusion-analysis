# ==============================================================================
# 1. Load required packages
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(MASS)    # 2D density
  library(viridis)
})

# ==============================================================================
# 2. Configure paths and parameters
# ==============================================================================
# [Important] Set this to your cache directory
CACHE_DIR <- "./correlation_plots/cache" 

METRICS <- c("mean", "max", "sum")

TRACK_NAMES <- c(
  "NIPBL",
  "CTCF_ChIP_seq",
  "Cohesin_ChIP_seq",
  "CTCF_ChIA_PET",
  "Cohesin_ChIA_PET",
  "CTCF_ChIA_Drop",
  "Cohesin_ChIA_Drop"
)

# ==============================================================================
# simple comment
# ==============================================================================

# simple comment
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

# simple comment
create_scatter_plot <- function(df, x_col, y_col, metric, title_prefix) {
  if (!all(c(x_col, y_col) %in% names(df))) return(NULL)
  
  sub_df <- na.omit(df[, c(x_col, y_col)])
  names(sub_df) <- c("X", "Y")
  
  # --- NEW: Outlier Removal (Quantile 1% - 99%) ---
  # Remove extreme values to prevent axis compression
  qx <- quantile(sub_df$X, probs = c(0.0001, 0.9999), na.rm = TRUE)
  qy <- quantile(sub_df$Y, probs = c(0.0001, 0.9999), na.rm = TRUE)
  
  # Keep data within the [1%, 99%] range for both axes
  keep_idx <- sub_df$X >= qx[1] & sub_df$X <= qx[2] &
              sub_df$Y >= qy[1] & sub_df$Y <= qy[2]
  
  sub_df <- sub_df[keep_idx, ]
  
  # Check if enough data remains
  n_points <- nrow(sub_df)
  if(n_points < 10) return(NULL)
  # ------------------------------------------------
  
  # Recalculate density on filtered data
  sub_df$density <- get_density(sub_df$X, sub_df$Y)
  
  # Calculate correlations
  r_p <- round(cor(sub_df$X, sub_df$Y, method = "pearson"), 2)
  r_s <- round(cor(sub_df$X, sub_df$Y, method = "spearman"), 2)
  
  # Plotting
  p <- ggplot(sub_df, aes(x = X, y = Y)) +
    geom_point(aes(color = density), size = 0.5, alpha = 0.6) + 
    # Custom color palette: DarkBlue -> Yellow -> Red
    scale_color_gradientn(colours = c("#00008B", "#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000")) +
    labs(x = paste0(x_col, " (", metric, ")"), 
         y = paste0(y_col, " (", metric, ")"),
         title = paste0(title_prefix, "\nn=", n_points, " | (Filtered) P: ", r_p, " | S: ", r_s)) +
    theme_classic() + # Ensures white background
    theme(
      legend.position = "none",
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(color = "black", size = 8)
    )
  return(p)
}

analyze_one_metric <- function(metric) {
  message(sprintf("\n>>> Processing Metric: %s (Outliers Removed) <<<", metric))
  
  # --- Step 1: Read Data ---
  data_list <- list()
  for (track in TRACK_NAMES) {
    file_path <- file.path(CACHE_DIR, paste0(track, "_", metric, ".rds"))
    if (file.exists(file_path)) {
      data_list[[track]] <- readRDS(file_path)
    }
  }
  
  if (length(data_list) == 0) {
    warning(paste("No cache files found for", metric))
    return(NULL)
  }
  
  df <- as.data.frame(data_list)
  
  # --- Step 2: Generate Plots ---
  plot_list <- list()
  
  # 1. NIPBL vs Others
  if ("NIPBL" %in% names(df)) {
    others <- setdiff(names(df), "NIPBL")
    for (target in others) {
      plot_list[[paste0("NIPBL_", target)]] <- create_scatter_plot(
        df, "NIPBL", target, metric, paste("NIPBL vs", target)
      )
    }
  }
  
  # 2. Pairs
  plot_list[["ChIP_Seq_Pair"]] <- create_scatter_plot(
    df, "CTCF_ChIP_seq", "Cohesin_ChIP_seq", metric, "CTCF vs Cohesin (ChIP-seq)"
  )
  plot_list[["ChIA_PET_Pair"]] <- create_scatter_plot(
    df, "CTCF_ChIA_PET", "Cohesin_ChIA_PET", metric, "CTCF vs Cohesin (ChIA-PET)"
  )
  plot_list[["ChIA_Drop_Pair"]] <- create_scatter_plot(
    df, "CTCF_ChIA_Drop", "Cohesin_ChIA_Drop", metric, "CTCF vs Cohesin (ChIA-Drop)"
  )
  
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  return(plot_list)
}

# ==============================================================================
# 4. Loop over metrics and output
# ==============================================================================

for (met in METRICS) {
  plots <- analyze_one_metric(met)
  
  if (!is.null(plots) && length(plots) > 0) {
    # Determine layout columns
    ncol_set <- if(length(plots) > 6) 4 else 3
    
    # Print to screen
    grid.arrange(grobs = plots, ncol = ncol_set, 
                 top = paste0("Signal Correlation [", met, "] - Top/Bottom 1% Filtered"))
    
    # Save to file with NEW filename
    ggsave(
      filename = paste0("Correlation_", met, "_OutlierRm.pdf"), 
      plot = arrangeGrob(grobs = plots, ncol = ncol_set), 
      width = 15, height = 10
    )
  }
}

message("\nDone.")