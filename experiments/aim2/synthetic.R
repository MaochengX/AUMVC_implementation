source("experiments/aim2/evaluation.R")
settings <- AIM2_SETTINGS
set.seed(settings$seed)

x_normal <- matrix(rnorm(settings$synthetic$n_normal * 2), ncol = 2)
angle <- runif(settings$synthetic$n_anomaly, 0, 2 * pi)
radius <- runif(
  settings$synthetic$n_anomaly,
  settings$synthetic$radius_min,
  settings$synthetic$radius_max
)
x_anomaly <- cbind(radius * cos(angle), radius * sin(angle))
x <- rbind(x_normal, x_anomaly)
labels <- c(
  rep(0L, settings$synthetic$n_normal),
  rep(1L, settings$synthetic$n_anomaly)
)

output <- run_dataset(
  x,
  labels,
  "Synthetic",
  settings
)
directory <- next_experiment("experiments/aim2/result/synthetic")
path <- file.path(directory, "synthetic_results.csv")
save_csv(report_rows(output), path)
cat("Saved: ", path, "\n", sep = "")
