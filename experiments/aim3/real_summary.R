source("experiments/settings.R")

settings <- AIM3_SETTINGS
results_root <- "experiments/aim3/results/real"
latest_file <- file.path(results_root, "latest_experiment.txt")

if (!file.exists(latest_file)) {
  stop(
    "No completed Aim 3 real-data experiment was found.",
    call. = FALSE
  )
}

experiment_id <- trimws(
  readLines(latest_file, warn = FALSE)[1L]
)

datasets <- settings$real$datasets

for (dataset in datasets) {
  experiment_dir <- file.path(
    results_root,
    dataset,
    experiment_id
  )

  path <- file.path(
    experiment_dir,
    "real.rds"
  )

  if (!file.exists(path)) {
    stop(
      "Missing real.rds for ",
      dataset,
      ".",
      call. = FALSE
    )
  }

  results <- readRDS(path)$results
  runs <- sort(unique(results$run))

  for (run in runs) {
    output <- results[
      results$run == run,
      c(
        "representation",
        "aumvc",
        "aumvc_normalized",
        "aumvc_normalized_mc_se",
        "zero_occupancy",
        "runtime_seconds"
      ),
      drop = FALSE
    ]

    names(output) <- c(
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
}
