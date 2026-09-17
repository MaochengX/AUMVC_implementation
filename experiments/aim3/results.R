aim3_open_run_batch <- function(root, settings, runner, n_runs) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  n_runs <- as.integer(n_runs)
  if (length(n_runs) != 1L || is.na(n_runs) || n_runs < 1L) {
    stop("n_runs must be a positive integer.")
  }

  entries <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  numbers <- suppressWarnings(as.integer(sub("^exp", "", basename(entries))))
  starts <- entries[!is.na(numbers) & grepl("^exp[0-9]+$", basename(entries))]
  starts <- starts[order(as.integer(sub("^exp", "", basename(starts))))]

  for (first in starts) {
    manifest <- file.path(first, "batch.rds")
    if (!file.exists(manifest)) next
    batch <- readRDS(manifest)
    if (!identical(batch$settings, settings) ||
        !identical(batch$n_runs, n_runs) ||
        file.exists(file.path(first, paste0(runner, ".done")))) next

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
  saveRDS(list(settings = settings, n_runs = n_runs, numbers = indices),
          file.path(paths[1L], "batch.rds"))
  for (path in paths) saveRDS(settings, file.path(path, "settings.rds"))
  paths
}

aim3_finish_run_batch <- function(paths, runner) {
  writeLines("complete", file.path(paths[1L], paste0(runner, ".done")))
}

aim3_save_result_csv <- function(result, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  utils::write.csv(result, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop("Could not save result CSV: ", path)
}

# Each run gets a full CSV and a compact view of the same measurements.
aim3_key_results <- function(results) {
  columns <- c(
    "run", "truth", "intrinsic_dim", "ambient_dim", "snr",
    "sample_size", "embedding_dim", "representation",
    "roc_auc", "pr_auc", "aumvc_normalized",
    "aumvc_normalized_mc_se", "zero_occupancy", "low_occupancy",
    "runtime_seconds"
  )
  results[, intersect(columns, names(results)), drop = FALSE]
}

aim3_save_results <- function(dataset_dir, results) {
  aim3_save_result_csv(results, file.path(dataset_dir, "all_results.csv"))
  aim3_save_result_csv(aim3_key_results(results),
                  file.path(dataset_dir, "key_results.csv"))
}

aim3_results_complete <- function(dataset_dir) {
  full_path <- file.path(dataset_dir, "all_results.csv")
  key_path <- file.path(dataset_dir, "key_results.csv")
  if (!file.exists(full_path)) return(FALSE)
  if (!file.exists(key_path)) {
    aim3_save_result_csv(aim3_key_results(utils::read.csv(full_path)), key_path)
  }
  TRUE
}

aim3_add_ambient_deltas <- function(results) {
  ambient <- results[results$representation == "ambient", ]
  results$roc_delta_ambient <- results$roc_auc - ambient$roc_auc
  results$pr_delta_ambient <- results$pr_auc - ambient$pr_auc
  results
}
