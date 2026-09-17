AIM1_SETTINGS <- list(
  seed = 1234L,
  n_reference = 20000L,
  n_mc_repetitions = 5L,
  aumvc_alpha_grid = seq(0.9, 0.999, by = 0.0001),
  one_cluster = list(
    n_train = 200L,
    n_eval = 200L,
    n_anomaly = 2L,
    anomaly_mean = 5,
    anomaly_sd = 0.3
  ),
  two_clusters = list(
    n_per_cluster = 100L,
    cluster_mean = 2,
    cluster_sd = 0.7
  )
)

experiment_seed <- function(settings, offset = 0L) {
  settings$seed + as.integer(offset)
}
