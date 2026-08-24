source("experiments/aim2/validation_agianst_labels.R")

set.seed(AIM2_SETTINGS$seed)

cfg <- AIM2_SETTINGS$synthetic

x_normal <- matrix(
  rnorm(cfg$n_normal * 2),
  ncol = 2
)

angle <- runif(
  cfg$n_anomaly,
  0,
  2 * pi
)

radius <- runif(
  cfg$n_anomaly,
  cfg$radius_min,
  cfg$radius_max
)

x_anomaly <- cbind(
  radius * cos(angle),
  radius * sin(angle)
)

x <- rbind(
  x_normal,
  x_anomaly
)

labels <- c(
  rep(0L, cfg$n_normal),
  rep(1L, cfg$n_anomaly)
)

run_aim2_dataset(
  x = x,
  labels = labels,
  dataset = "Synthetic",
  counts = AIM2_SETTINGS$split_counts$synthetic
)
