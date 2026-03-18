# ==============================================================================
# 1. Load required packages
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(MASS)    # for 2D density
  library(viridis)
})

# ==============================================================================
# 2. Configure paths and parameters
# ==============================================================================
# [Important] Set this to your cache directory
CACHE_DIR <- "./correlation_plots/cache" 

# All metrics to analyze
METRICS <- c("mean", "max", "sum")

# Track names (ensure they match file names)
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
# 3. Core functions
# ==============================================================================

# 3.1 Helper: compute density colors
get_density <- function(x, y, n = 100) {
  valid <- is.finite(x) & is.finite(y)
  x_clean <- x[valid]
  y_clean <- y[valid]
  
  # If too few points, return all 1s
  if(length(x_clean) < 10) return(rep(1, length(x)))
  
  dens <- MASS::kde2d(x = x_clean, y = y_clean, n = n)
  ix <- findInterval(x_clean, dens$x)
  iy <- findInterval(y_clean, dens$y)
  ii <- cbind(ix, iy)
  
  # Fill back to original length
  res <- rep(NA, length(x))
  res[valid] <- dens$z[ii]
  return(res)
}

# 3.2 Plot unit function
create_scatter_plot <- function(df, x_col, y_col, metric, title_prefix) {
  # Verify columns exist
  if (!all(c(x_col, y_col) %in% names(df))) return(NULL)
  
  sub_df <- na.omit(df[, c(x_col, y_col)])
  names(sub_df) <- c("X", "Y")
  
  # Compute density
  sub_df$density <- get_density(sub_df$X, sub_df$Y)
  
  # Compute correlations
  r_p <- round(cor(sub_df$X, sub_df$Y, method = "pearson"), 2)
  r_s <- round(cor(sub_df$X, sub_df$Y, method = "spearman"), 2)
  
  # Plot
  p <- ggplot(sub_df, aes(x = X, y = Y)) +
    geom_point(aes(color = density), size = 0.5, alpha = 0.6) + 
    # Reproduce specified DeepBlue -> Yellow -> Red palette
    scale_color_gradientn(colours = c("#00008B", "#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000")) +
    labs(x = paste0(x_col, " (", metric, ")"), 
         y = paste0(y_col, " (", metric, ")"),
         title = paste0(title_prefix, "\nP: ", r_p, " | S: ", r_s)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(color = "black", size = 8)
    )
  return(p)
}

# 3.3 Main analysis function (per metric)
analyze_one_metric <- function(metric) {
  message(sprintf("\n>>> Processing metric: %s <<<", metric))

  # --- Step 1: read data ---
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

  # --- Step 2: build plot list ---
  plot_list <- list()
  
  # 1. NIPBL vs all other tracks
  if ("NIPBL" %in% names(df)) {
    others <- setdiff(names(df), "NIPBL")
    for (target in others) {
      plot_list[[paste0("NIPBL_", target)]] <- create_scatter_plot(
        df, "NIPBL", target, metric, paste("NIPBL vs", target)
      )
    }
  }
  
  # 2. Special pair: CTCF ChIP vs Cohesin ChIP
  plot_list[["ChIP_Seq_Pair"]] <- create_scatter_plot(
    df, "CTCF_ChIP_seq", "Cohesin_ChIP_seq", metric, "CTCF vs Cohesin (ChIP-seq)"
  )
  
  # 3. Special pair: CTCF ChIA-PET vs Cohesin ChIA-PET (new request)
  plot_list[["ChIA_PET_Pair"]] <- create_scatter_plot(
    df, "CTCF_ChIA_PET", "Cohesin_ChIA_PET", metric, "CTCF vs Cohesin (ChIA-PET)"
  )
  
  # 4. Special pair: CTCF ChIA-Drop vs Cohesin ChIA-Drop (new request)
  plot_list[["ChIA_Drop_Pair"]] <- create_scatter_plot(
    df, "CTCF_ChIA_Drop", "Cohesin_ChIA_Drop", metric, "CTCF vs Cohesin (ChIA-Drop)"
  )
  
  # Drop empty plots (if some files missing)
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  return(plot_list)
}

# ==============================================================================
# 4. Loop over metrics and output
# ==============================================================================

for (met in METRICS) {
  # Generate all plots for current metric
  plots <- analyze_one_metric(met)
  
  if (!is.null(plots) && length(plots) > 0) {
    # Auto layout: if >6 plots use 4 cols, else 3
    ncol_set <- if(length(plots) > 6) 4 else 3
    
    grid.arrange(grobs = plots, ncol = ncol_set, 
                 top = paste0("Signal Correlation Analysis [Metric: ", met, "]"))
    
    # Hint: uncomment to save PDF
    ggsave(paste0("Correlation_", met, ".pdf"), arrangeGrob(grobs = plots, ncol = ncol_set), width = 15, height = 10)
  }
}

