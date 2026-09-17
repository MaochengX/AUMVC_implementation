source("experiments/aim3/settings.R")
source("experiments/aim3/synthetic.R")
source("experiments/aim3/representation_comparison.R")
source("experiments/aim3/results.R")

settings <- AIM3_SETTINGS
run_dirs <- aim3_open_run_batch(
  "experiments/aim3/result", settings, "synthetic", settings$n_runs
)

grid <- expand.grid(
  truth = settings$synthetic$truth_functions,
  intrinsic_dim = settings$synthetic$intrinsic_dims,
  ambient_dim = settings$synthetic$ambient_dims,
  snr = settings$synthetic$snr_levels,
  stringsAsFactors = FALSE
)
total <- settings$n_runs * nrow(grid)
done <- 0L
start <- proc.time()[[3L]]

for (run in seq_len(settings$n_runs)) {
  dataset_dir <- file.path(run_dirs[run], "synthetic")
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- file.path(dataset_dir, "checkpoints")
  dir.create(checkpoint_dir, showWarnings = FALSE)
  if (aim3_results_complete(dataset_dir)) {
    done <- done + nrow(grid)
    next
  }
  cases <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    checkpoint <- file.path(checkpoint_dir, sprintf("case%d.rds", i))
    if (file.exists(checkpoint)) {
      cases[[i]] <- readRDS(checkpoint)
    } else {
      row <- grid[i, ]
      # The same seed is reused across d so the smaller coordinate set is paired.
      case_seed <- settings$seed + run * 100000L +
        match(row$truth, settings$synthetic$truth_functions) * 10000L +
        row$intrinsic_dim * 100L + as.integer(row$snr)
      case <- aim3_generate_synthetic(
        row$intrinsic_dim, row$snr, row$truth, case_seed,
        settings$synthetic, row$ambient_dim
      )
      comparison <- aim3_compare_representations(
        case, case_seed + 50000L, settings
      )
      cases[[i]] <- data.frame(row, run = run,
                               aim3_add_ambient_deltas(comparison))
      temporary <- paste0(checkpoint, ".tmp")
      saveRDS(cases[[i]], temporary)
      if (!file.rename(temporary, checkpoint)) stop("Could not save checkpoint.")
    }
    done <- done + 1L
    elapsed <- proc.time()[[3L]] - start
    eta <- if (done > 0L) elapsed * (total - done) / done else NA_real_
    cat(sprintf("\rSynthetic %d/%d | elapsed %.0fs | ETA %.0fs",
                done, total, elapsed, eta))
    flush.console()
  }
  aim3_save_results(dataset_dir, do.call(rbind, cases))
}
cat("\n")
aim3_finish_run_batch(run_dirs, "synthetic")
cat("Saved: ", run_dirs[1L], " through ", tail(run_dirs, 1L),
    "/synthetic/{all_results,key_results}.csv\n", sep = "")
