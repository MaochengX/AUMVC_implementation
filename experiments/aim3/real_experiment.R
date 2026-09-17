source("experiments/aim3/settings.R")
source("experiments/aim3/representation_comparison.R")
source("experiments/aim3/results.R")
source("experiments/aim3/real_data.R")

settings <- AIM3_SETTINGS
real_settings <- settings$real
loaders <- list(
  ecg200 = aim3_load_ecg200,
  fashion_mnist = aim3_load_fashion_mnist,
  shuttle = aim3_load_shuttle
)
samplers <- list(
  ecg200 = aim3_sample_ecg200,
  fashion_mnist = aim3_sample_fashion_mnist,
  shuttle = aim3_sample_shuttle
)

if (any(!real_settings$datasets %in% names(loaders))) {
  stop("Unknown real dataset in settings.R.")
}
run_dirs <- aim3_open_run_batch(
  "experiments/aim3/result", settings, "real", settings$n_runs
)

for (dataset in real_settings$datasets) {
  dataset_settings <- real_settings[[dataset]]
  data <- loaders[[dataset]](dataset_settings$data_dir)
  output_name <- switch(dataset, ecg200 = "ECG",
                        fashion_mnist = "fashion", shuttle = "shuttle")
  current_md5 <- aim3_real_data_md5(data$files)
  start <- proc.time()[[3L]]
  for (run in seq_len(settings$n_runs)) {
    run_dir <- run_dirs[run]
    state_dir <- file.path(aim3_run_state_dir(run_dir), "real")
    dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
    md5_path <- file.path(state_dir, paste0(output_name, "_md5.rds"))
    if (file.exists(md5_path)) {
      if (!identical(readRDS(md5_path), current_md5)) {
        stop("Dataset files changed while resuming ", dataset)
      }
    } else {
      saveRDS(current_md5, md5_path)
    }
    if (aim3_dataset_complete(run_dir, output_name, 3L)) next
    seed <- settings$seed + run * 100000L
    sampled <- samplers[[dataset]](data, seed, dataset_settings)
    analysis <- settings
    analysis$split_counts <- aim3_real_split_counts(nrow(sampled$x))
    case <- list(x = sampled$x, labels = sampled$labels,
                 intrinsic_dim = as.integer(real_settings$embedding_dim))
    comparison <- aim3_compare_representations(case, seed + 50000L, analysis)
    result <- data.frame(
      run = run, sample_size = nrow(sampled$x), ambient_dim = ncol(sampled$x),
      embedding_dim = case$intrinsic_dim, aim3_add_ambient_deltas(comparison)
    )
    aim3_save_dataset(run_dir, output_name, result)
    elapsed <- proc.time()[[3L]] - start
    cat(sprintf("\r%s %d/%d | elapsed %.0fs | ETA %.0fs",
                dataset, run, settings$n_runs, elapsed,
                elapsed * (settings$n_runs - run) / run))
    flush.console()
  }
  cat("\n")
  cat("Completed ", dataset, "\n", sep = "")
}
aim3_finish_run_batch(run_dirs, "real")
cat("Saved: ", run_dirs[1L], " through ", tail(run_dirs, 1L),
    " (all_results.csv and key_results.csv for real data in each run)\n",
    sep = "")
