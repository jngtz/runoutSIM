library(lhs)
library(future.apply)

# --- Setup parallel plan ---
plan(multisession, workers = parallel::detectCores() - 2)

# --- Generate random samples ---
n_samples <- 50


# Generate Latin Hypercube Samples for 5 parameters: 3 RW + 2 PCM
lhs_combined <- randomLHS(n_samples, 5)

combined_samples <- data.frame(
  # RW parameters
  rw_slp = qunif(lhs_combined[,1], 20, 40),
  rw_ex  = qunif(lhs_combined[,2], 1.3, 3),
  rw_per = qunif(lhs_combined[,3], 1.5, 2),
  
  # PCM parameters
  pcm_mu = qunif(lhs_combined[,4], 0.05, 0.4),
  pcm_md = qunif(lhs_combined[,5], 20, 120)
)
# --- Evaluation function ---
evaluate_combined <- function(rw, pcm) {
  n_slides <- nrow(runout_polygons)
  
  # Pre-allocate results matrix
  results <- matrix(NA, nrow = 3, ncol = n_slides,
                    dimnames = list(c("roc", "length_error", "length_relerr"), NULL))
  
  for (i in 1:n_slides) {
    out <- tryCatch({
      pcmPerformance(
        dem            = dem,
        slide_plys     = runout_polygons,
        slide_src      = source_points,
        slide_id       = i,
        rw_slp         = rw$rw_slp,
        rw_ex          = rw$rw_ex,
        rw_per         = rw$rw_per,
        pcm_mu         = pcm$pcm_mu,
        pcm_md         = pcm$pcm_md,
        gpp_iter       = 1000,
        buffer_source  = NULL,
        plot_eval      = FALSE,
        return_features= FALSE
      )
    }, error = function(e) NULL)
    
    if (is.null(out)) {
      results[, i] <- c(roc = 0, length_error = Inf, length_relerr = Inf)
    } else {
      results[, i] <- c(roc = out$roc, 
                        length_error = out$length.error, 
                        length_relerr = out$length.relerr)
    }
  }
  
  # Compute mean and median over all slides
  mean_roc <- mean(results["roc", ], na.rm = TRUE)
  mean_len <- mean(results["length_error", ], na.rm = TRUE)
  mean_relerr <- mean(results["length_relerr", ], na.rm = TRUE)
  
  median_roc <- median(results["roc", ], na.rm = TRUE)
  median_len <- median(results["length_error", ], na.rm = TRUE)
  median_relerr <- median(results["length_relerr", ], na.rm = TRUE)
  
  return(list(
    mean_roc = mean_roc,
    mean_length_error = mean_len,
    mean_relerr = mean_relerr,
    median_roc = median_roc,
    median_length_error = median_len,
    median_relerr = median_relerr
  ))
}


# --- Evaluate all combinations ---
# --- Evaluate all combined samples in parallel ---
combined_results <- future_lapply(1:nrow(combined_samples), function(k) {
  row <- as.list(combined_samples[k, ])
  
  res <- evaluate_combined(row, row)  # rw and pcm are from the same row
  
  c(sample_idx = k,
    rw_slp = row$rw_slp,
    rw_ex  = row$rw_ex,
    rw_per = row$rw_per,
    pcm_mu = row$pcm_mu,
    pcm_md = row$pcm_md,
    mean_roc = res$mean_roc,
    mean_length_error = res$mean_length_error,
    mean_relerr = res$mean_relerr,
    median_roc = res$median_roc,
    median_length_error = res$median_length_error,
    median_relerr = res$median_relerr)
})

# --- Convert results to data frame ---
combined_df <- do.call(rbind, lapply(combined_results, function(x) as.data.frame(as.list(x))))

# --- Optional: compute a weighted score ---
weight_len <- 2
combined_df$score <- combined_df$mean_roc - weight_len * combined_df$mean_length_error / max(combined_df$mean_length_error, 1e-6)

# --- Find best combination ---
best_row <- combined_df[which.max(combined_df$score), ]
best_params <- combined_df[which.max(combined_df$score), c("rw_slp","rw_ex","rw_per","pcm_mu","pcm_md")]


# --- Example: weighted score (length error twice as important as ROC) ---
weight_len <- 2
combined_df$score <- combined_df$mean_roc - weight_len * combined_df$mean_length_error / max(combined_df$mean_length_error,1e-6)

# --- Find best combination ---
best_row <- combined_df[which.max(combined_df$score), ]
best_rw <- rw_samples[best_row$rw_idx, ]
best_pcm <- pcm_samples[best_row$pcm_idx, ]

print(list(best_rw = best_rw, best_pcm = best_pcm, score = best_row$score))


future:::ClusterRegistry("stop")
##########################################################

# This works
library(lhs)
library(future.apply)

# --- Setup parallel plan ---
plan(multisession, workers = parallel::detectCores() - 2)

# --- Generate random samples for RW parameters ---
n_samples <- 50
lhs_design <- randomLHS(n_samples, 3)

rw_samples <- data.frame(
  rw_slp = qunif(lhs_design[,1], 20, 40),   # slope threshold (deg)
  rw_ex  = qunif(lhs_design[,2], 1.3, 3),     # exponent
  rw_per = qunif(lhs_design[,3], 1.5, 2)    # persistence
)

# --- Evaluation function (maximize ROC) ---
evaluate_rw <- function(x) {
  rocs <- sapply(seq_len(nrow(runout_polygons)), function(i) {
    out <- tryCatch({
      pcmPerformance(
        dem            = dem,
        slide_plys     = runout_polygons,
        slide_src      = source_points,
        slide_id       = i,
        rw_slp         = x$rw_slp,
        rw_ex          = x$rw_ex,
        rw_per         = x$rw_per,
        pcm_mu         = 0.0001,   # fixed PCM params
        pcm_md         = 40,
        gpp_iter       = 1000,
        buffer_source  = NULL,
        plot_eval      = FALSE,
        return_features= FALSE
      )
    }, error = function(e) NULL)
    
    if (is.null(out)) return(0) else return(out$roc)
  })
  
  mean(rocs, na.rm = TRUE)
}

# --- Run evaluations in parallel across random parameter sets ---
rw_scores <- future_sapply(1:n_samples, function(i) {
  evaluate_rw(as.list(rw_samples[i, ]))
})

# --- Combine parameters with scores ---
rw_results <- cbind(rw_samples, score = rw_scores)

# --- Pick best set ---
best_rw <- rw_results[which.max(rw_results$score), ]

print(best_rw)

### PCM #######

# --- Setup parallel plan ---
plan(multisession, workers = parallel::detectCores() - 2)

# --- Generate random samples for PCM parameters ---
n_samples <- 50
lhs_design <- randomLHS(n_samples, 2)

pcm_samples <- data.frame(
  pcm_mu = qunif(lhs_design[,1], 0.05, 0.4),   # friction coefficient
  pcm_md = qunif(lhs_design[,2], 20, 120)      # mass-to-drag ratio
)

# --- Evaluation function (minimize length error) ---
evaluate_pcm <- function(x) {
  errs <- sapply(seq_len(nrow(runout_polygons)), function(i) {
    out <- tryCatch({
      pcmPerformance(
        dem            = dem,
        slide_plys     = runout_polygons,
        slide_src      = source_points,
        slide_id       = i,
        rw_slp         = 40,   # fixed RW params (can be best_rw)
        rw_ex          = 2,
        rw_per         = 2,
        pcm_mu         = x$pcm_mu,
        pcm_md         = x$pcm_md,
        gpp_iter       = 1000,
        buffer_source  = 20,
        plot_eval      = FALSE,
        return_features= FALSE
      )
    }, error = function(e) NULL)
    
    if (is.null(out)) return(Inf) else return(out$length.error)
  })
  
  mean(errs, na.rm = TRUE)
}

# --- Run evaluations in parallel across random parameter sets ---
pcm_scores <- future_sapply(1:n_samples, function(i) {
  evaluate_pcm(as.list(pcm_samples[i, ]))
})

# --- Combine parameters with scores ---
pcm_results <- cbind(pcm_samples, score = pcm_scores)

# --- Pick best set (lowest length error) ---
best_pcm <- pcm_results[which.min(pcm_results$score), ]

print(best_pcm)
