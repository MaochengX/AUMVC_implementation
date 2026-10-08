source("experiments/aim3/settings.R")
source("experiments/aim3/synthetic.R")
source("experiments/aim3/representation_comparison.R")
source("experiments/aim3/results.R")

settings <- AIM3_SETTINGS
grid <- expand.grid(
  truth = settings$synthetic$truth_functions,
  intrinsic_dim = settings$synthetic$intrinsic_dims,
  ambient_dim = settings$synthetic$ambient_dims,
  snr = settings$synthetic$snr_levels,
  stringsAsFactors = FALSE
)
total <- settings$n_runs * nrow(grid)
results <- vector("list", total)
done <- 0L
start <- proc.time()[[3L]]

for (run in seq_len(settings$n_runs)) {
  for (index in seq_len(nrow(grid))) {
    row <- grid[index, ]
    case_seed <- settings$seed + run * 100000L +
      match(row$truth, settings$synthetic$truth_functions) * 10000L +
      row$intrinsic_dim * 100L + as.integer(row$snr)
    case <- aim3_generate_synthetic(
      row$intrinsic_dim,
      row$snr,
      row$truth,
      case_seed,
      settings$synthetic,
      row$ambient_dim
    )
    comparison <- aim3_compare_representations(
      case,
      case_seed + 50000L,
      settings
    )
    effective_dim <- attr(comparison, "effective_dim")
    done <- done + 1L
    results[[done]] <- data.frame(
      run = run,
      truth = row$truth,
      intrinsic_dim = row$intrinsic_dim,
      original_dim = row$ambient_dim,
      effective_dim = effective_dim,
      snr = row$snr,
      sample_size = nrow(case$x),
      embedding_dim = row$intrinsic_dim,
      aim3_add_ambient_deltas(comparison)
    )
    elapsed <- proc.time()[[3L]] - start
    eta <- elapsed * (total - done) / done
    cat(sprintf(
      "\rSynthetic %d/%d | elapsed %.0fs | ETA %.0fs",
      done,
      total,
      elapsed,
      eta
    ))
    flush.console()
  }
}

cat("\n")
directory <- aim3_next_experiment("experiments/aim3/result/synthetic")
paths <- aim3_save_result_set(
  do.call(rbind, results),
  directory,
  "synthetic",
  TRUE
)
cat("Saved: ", paths["all"], "\n", sep = "")
cat("Saved: ", paths["key"], "\n", sep = "")
