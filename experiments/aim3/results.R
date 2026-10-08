aim3_add_relative_se <- function(results) {
  if ("aumvc_relative_mc_se" %in% names(results)) return(results)
  relative <- rep(NA_real_, nrow(results))
  valid <- is.finite(results$aumvc_normalized) &
    results$aumvc_normalized > 0 &
    is.finite(results$aumvc_normalized_mc_se)
  relative[valid] <- results$aumvc_normalized_mc_se[valid] /
    results$aumvc_normalized[valid]
  results$aumvc_relative_mc_se <- relative
  results
}

aim3_add_ambient_deltas <- function(results) {
  ambient <- results[results$representation == "ambient", , drop = FALSE]
  if (nrow(ambient) != 1L) {
    stop("Each Aim 3 case must contain one ambient result.", call. = FALSE)
  }
  results$roc_delta_ambient <- results$roc_auc - ambient$roc_auc
  results$pr_delta_ambient <- results$pr_auc - ambient$pr_auc
  aim3_add_relative_se(results)
}

aim3_select_columns <- function(results, columns) {
  missing <- setdiff(columns, names(results))
  if (length(missing) > 0L) {
    stop("Missing Aim 3 result columns: ", paste(missing, collapse = ", "))
  }
  results[, columns, drop = FALSE]
}

aim3_synthetic_all_columns <- c(
  "run", "truth", "intrinsic_dim", "original_dim", "effective_dim",
  "snr", "sample_size", "embedding_dim", "representation", "aumvc",
  "aumvc_normalized", "aumvc_mc_se", "aumvc_relative_mc_se",
  "aumvc_subset_se", "minimum_hit_count", "zero_occupancy",
  "low_occupancy", "box_log_volume", "roc_auc", "pr_auc",
  "roc_distributional", "pr_distributional",
  "roc_structural", "pr_structural", "roc_delta_ambient",
  "pr_delta_ambient", "ocsvm_kkt_gap",
  "embedding_oos_stress", "runtime_seconds"
)

aim3_real_all_columns <- c(
  "dataset", "run", "sample_size", "original_dim", "effective_dim",
  "embedding_dim", "representation", "aumvc", "aumvc_normalized",
  "aumvc_mc_se", "aumvc_relative_mc_se", "aumvc_subset_se",
  "minimum_hit_count", "zero_occupancy", "low_occupancy",
  "box_log_volume", "roc_auc", "pr_auc", "roc_delta_ambient",
  "pr_delta_ambient", "ocsvm_kkt_gap",
  "embedding_oos_stress", "runtime_seconds"
)

aim3_synthetic_key_columns <- c(
  "run", "truth", "intrinsic_dim", "original_dim", "snr",
  "sample_size", "embedding_dim", "representation",
  "aumvc_normalized", "aumvc_relative_mc_se", "minimum_hit_count",
  "zero_occupancy", "low_occupancy",
  "roc_auc", "pr_auc", "roc_distributional", "pr_distributional",
  "roc_structural", "pr_structural", "roc_delta_ambient",
  "pr_delta_ambient", "embedding_oos_stress", "runtime_seconds"
)

aim3_real_key_columns <- c(
  "dataset", "run", "sample_size", "original_dim", "embedding_dim",
  "representation", "aumvc_normalized", "aumvc_relative_mc_se",
  "minimum_hit_count", "zero_occupancy", "low_occupancy",
  "roc_auc", "pr_auc", "roc_delta_ambient", "pr_delta_ambient",
  "embedding_oos_stress", "runtime_seconds"
)

aim3_save_result_set <- function(results, directory, prefix, synthetic = FALSE) {
  results <- aim3_add_relative_se(results)
  if (synthetic) {
    all <- aim3_select_columns(results, aim3_synthetic_all_columns)
    key <- aim3_select_columns(results, aim3_synthetic_key_columns)
  } else {
    all <- aim3_select_columns(results, aim3_real_all_columns)
    key <- aim3_select_columns(results, aim3_real_key_columns)
  }
  all_path <- file.path(directory, paste0(prefix, "_all_results.csv"))
  key_path <- file.path(directory, paste0(prefix, "_key_results.csv"))
  aim3_save_csv(all, all_path)
  aim3_save_csv(key, key_path)
  c(all = all_path, key = key_path)
}
