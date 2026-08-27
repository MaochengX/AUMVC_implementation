source("experiments/settings.R")

settings <- AIM3_SETTINGS
path <- "experiments/aim3/results/synthetic.rds"

if (!file.exists(path)) {
  stop(
    "Run experiments/aim3/synthetic_experiment.R first.",
    call. = FALSE
  )
}

saved <- readRDS(path)

if (!identical(saved$settings, settings)) {
  stop("Saved Aim 3 results use different settings.", call. = FALSE)
}

results <- saved$results
summary_table <- results[, c(
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

names(summary_table) <- c(
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

write.csv(
  summary_table,
  "experiments/aim3/results/summary.csv",
  row.names = FALSE
)

print(summary_table, row.names = FALSE, digits = 5)
