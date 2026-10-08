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
display_names <- c(
  ecg200 = "ECG200",
  fashion_mnist = "Fashion-MNIST",
  shuttle = "Shuttle"
)

if (
  length(real_settings$datasets) == 0L ||
  any(!real_settings$datasets %in% names(loaders))
) {
  stop("Select valid real datasets in Aim 3 settings.", call. = FALSE)
}

total <- settings$n_runs * length(real_settings$datasets)
results <- vector("list", total)
done <- 0L
start <- proc.time()[[3L]]
directory <- aim3_next_experiment("experiments/aim3/result/real")

for (dataset in real_settings$datasets) {
  dataset_settings <- real_settings[[dataset]]
  data <- loaders[[dataset]](dataset_settings$data_dir)
  for (run in seq_len(settings$n_runs)) {
    seed <- settings$seed + run * 100000L +
      match(dataset, real_settings$datasets) * 10000L
    sampled <- samplers[[dataset]](data, seed, dataset_settings)
    analysis <- settings
    analysis$split_counts <- aim3_real_split_counts(nrow(sampled$x))
    case <- list(
      x = sampled$x,
      labels = sampled$labels,
      intrinsic_dim = as.integer(real_settings$embedding_dim)
    )
    comparison <- aim3_compare_representations(
      case,
      seed + 50000L,
      analysis
    )
    effective_dim <- attr(comparison, "effective_dim")
    done <- done + 1L
    results[[done]] <- data.frame(
      dataset = display_names[[dataset]],
      run = run,
      sample_size = nrow(sampled$x),
      original_dim = ncol(sampled$x),
      effective_dim = effective_dim,
      embedding_dim = case$intrinsic_dim,
      aim3_add_ambient_deltas(comparison)
    )
    paths <- aim3_save_result_set(
      do.call(rbind, results[seq_len(done)]),
      directory,
      "real",
      FALSE
    )
    elapsed <- proc.time()[[3L]] - start
    eta <- elapsed * (total - done) / done
    cat(sprintf(
      "\rReal data %d/%d | elapsed %.0fs | ETA %.0fs",
      done,
      total,
      elapsed,
      eta
    ))
    flush.console()
  }
}

cat("\n")
cat("Saved: ", paths["all"], "\n", sep = "")
cat("Saved: ", paths["key"], "\n", sep = "")
