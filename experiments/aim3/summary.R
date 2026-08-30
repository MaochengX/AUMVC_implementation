results_root <- "experiments/aim3/results"
latest_file <- file.path(results_root, "latest_experiment.txt")

if (!file.exists(latest_file)) {
  stop(
    "No completed Aim 3 experiment was found.",
    call. = FALSE
  )
}

experiment_id <- trimws(
  readLines(latest_file, warn = FALSE)[1L]
)

experiment_dir <- file.path(
  results_root,
  experiment_id
)

path <- file.path(
  experiment_dir,
  "synthetic.rds"
)

if (!file.exists(path)) {
  stop(
    "The latest Aim 3 experiment has no synthetic.rds file.",
    call. = FALSE
  )
}

saved <- readRDS(path)
results <- saved$results
runs <- sort(unique(results$run))

for (run in runs) {
  run_results <- results[
    results$run == run,
    ,
    drop = FALSE
  ]

  output <- run_results[, c(
    "truth",
    "intrinsic_dim",
    "snr",
    "representation",
    "aumvc",
    "aumvc_normalized",
    "aumvc_normalized_mc_se",
    "zero_occupancy",
    "runtime_seconds"
  )]

  names(output) <- c(
    "truth",
    "intrinsic_dim",
    "snr",
    "representation",
    "AUMVC",
    "AUMVC_normalized",
    "MC_SE_normalized",
    "zero_occupancy",
    "runtime_seconds"
  )

  output_file <- file.path(
    experiment_dir,
    sprintf("run_%02d.csv", run)
  )

  write.csv(
    output,
    output_file,
    row.names = FALSE
  )

  cat("Saved: ", output_file, "\n", sep = "")
}
