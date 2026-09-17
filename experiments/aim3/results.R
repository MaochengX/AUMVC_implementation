aim3_open_run_batch <- function(root, settings, runner, n_runs) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  state <- file.path(dirname(root), "state")
  dir.create(state, recursive = TRUE, showWarnings = FALSE)
  n_runs <- as.integer(n_runs)
  if (length(n_runs) != 1L || is.na(n_runs) || n_runs < 1L) {
    stop("n_runs must be a positive integer.")
  }

  entries <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  numbers <- suppressWarnings(as.integer(sub("^exp", "", basename(entries))))
  starts <- entries[!is.na(numbers) & grepl("^exp[0-9]+$", basename(entries))]
  starts <- starts[order(as.integer(sub("^exp", "", basename(starts))))]

  for (first in starts) {
    manifest <- file.path(state, paste0(basename(first), ".rds"))
    done <- file.path(state, paste0(basename(first), "_", runner, ".done"))
    if (!file.exists(manifest) || file.exists(done)) next
    batch <- readRDS(manifest)
    if (!identical(batch$settings, settings) ||
        !identical(batch$n_runs, n_runs) ||
        !identical(batch$layout_version, 2L)) next
    paths <- file.path(root, paste0("exp", batch$numbers))
    if (all(dir.exists(paths))) return(paths)
    stop("An existing experiment batch is missing run directories.")
  }

  first <- if (length(numbers) && any(!is.na(numbers))) {
    max(numbers, na.rm = TRUE) + 1L
  } else 1L
  indices <- seq.int(first, length.out = n_runs)
  paths <- file.path(root, paste0("exp", indices))
  for (path in paths) dir.create(path)
  saveRDS(list(settings = settings, n_runs = n_runs, numbers = indices,
               layout_version = 2L),
          file.path(state, paste0("exp", first, ".rds")))
  paths
}

aim3_finish_run_batch <- function(paths, runner) {
  state <- file.path(dirname(dirname(paths[1L])), "state")
  done <- file.path(state, paste0(basename(paths[1L]), "_", runner, ".done"))
  writeLines("complete", done)
}

aim3_run_state_dir <- function(run_dir) {
  path <- file.path(dirname(dirname(run_dir)), "state", basename(run_dir))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

aim3_save_csv <- function(result, path) {
  temporary <- paste0(path, ".tmp")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(result, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Could not save result CSV: ", path)
}

aim3_full_columns <- c(
  "dataset", "run", "truth", "intrinsic_dim", "ambient_dim", "snr",
  "sample_size", "embedding_dim", "representation", "aumvc",
  "aumvc_normalized", "aumvc_mc_se", "aumvc_normalized_mc_se",
  "zero_occupancy", "low_occupancy", "box_log_volume", "roc_auc",
  "pr_auc", "roc_distributional", "pr_distributional",
  "roc_structural", "pr_structural", "runtime_seconds",
  "roc_delta_ambient", "pr_delta_ambient"
)

aim3_full_results <- function(results, columns = aim3_full_columns) {
  for (name in setdiff(columns, names(results))) {
    results[[name]] <- rep(NA, nrow(results))
  }
  results[, columns, drop = FALSE]
}

# Label performance tests detection quality; occupancy, MC error and runtime
# test whether the criterion remains usable as ambient dimension increases.
aim3_key_results <- function(results) {
  relative_se <- rep(NA_real_, nrow(results))
  valid <- is.finite(results$aumvc_normalized) &
    results$aumvc_normalized > 0 &
    is.finite(results$aumvc_normalized_mc_se)
  relative_se[valid] <- results$aumvc_normalized_mc_se[valid] /
    results$aumvc_normalized[valid]
  results$aumvc_relative_mc_se <- relative_se
  columns <- c(
    "dataset", "run", "sample_size", "ambient_dim", "embedding_dim",
    "representation", "roc_auc", "pr_auc",
    "roc_delta_ambient", "pr_delta_ambient", "aumvc_normalized",
    "aumvc_normalized_mc_se", "aumvc_relative_mc_se",
    "zero_occupancy", "low_occupancy", "runtime_seconds"
  )
  results[, columns, drop = FALSE]
}

aim3_save_synthetic <- function(run_dir, results) {
  results$dataset <- "synthetic"
  aim3_save_csv(aim3_full_results(results),
                file.path(run_dir, "synthetic.csv"))
}

aim3_synthetic_complete <- function(run_dir, expected_rows) {
  path <- file.path(run_dir, "synthetic.csv")
  if (!file.exists(path)) return(FALSE)
  results <- utils::read.csv(path, na.strings = c("", "NA"))
  nrow(results) == expected_rows &&
    "dataset" %in% names(results) &&
    all(results$dataset == "synthetic")
}

aim3_save_dataset <- function(run_dir, dataset, results) {
  if (dataset == "synthetic") stop("Use aim3_save_synthetic for synthetic data.")
  full_path <- file.path(run_dir, "all_results.csv")
  key_path <- file.path(run_dir, "key_results.csv")
  results$dataset <- dataset
  previous <- NULL
  if (file.exists(full_path)) {
    previous <- utils::read.csv(full_path, na.strings = c("", "NA"))
    if (any(previous$dataset == "synthetic")) {
      stop("This run uses the older mixed synthetic/real output layout.")
    }
    previous <- previous[previous$dataset != dataset, , drop = FALSE]
  }
  columns <- c(aim3_full_columns,
               setdiff(union(names(previous), names(results)), aim3_full_columns))
  current <- aim3_full_results(results, columns)
  if (!is.null(previous)) {
    current <- rbind(aim3_full_results(previous, columns), current)
  }
  dataset_order <- c("synthetic", "ECG", "fashion", "shuttle")
  current <- current[order(match(current$dataset, dataset_order)), , drop = FALSE]
  rownames(current) <- NULL
  aim3_save_csv(current, full_path)
  aim3_save_csv(aim3_key_results(current), key_path)
}

aim3_dataset_complete <- function(run_dir, dataset, expected_rows) {
  full_path <- file.path(run_dir, "all_results.csv")
  if (!file.exists(full_path)) return(FALSE)
  results <- utils::read.csv(full_path, na.strings = c("", "NA"))
  if (!"dataset" %in% names(results)) {
    stop("This run uses the older per-dataset output layout.")
  }
  if (any(results$dataset == "synthetic")) {
    stop("This run uses the older mixed synthetic/real output layout.")
  }
  if (sum(results$dataset == dataset) != expected_rows) return(FALSE)
  # Repair a missing or stale key file after an interrupted two-file write.
  aim3_save_csv(aim3_key_results(aim3_full_results(results)),
                file.path(run_dir, "key_results.csv"))
  TRUE
}

aim3_add_ambient_deltas <- function(results) {
  ambient <- results[results$representation == "ambient", ]
  results$roc_delta_ambient <- results$roc_auc - ambient$roc_auc
  results$pr_delta_ambient <- results$pr_auc - ambient$pr_auc
  results
}
