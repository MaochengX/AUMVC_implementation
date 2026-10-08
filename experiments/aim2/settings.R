AIM2_SETTINGS <- list(
  seed = 1234L,
  n_runs = 50L,
  n_reference = 100000L,
  n_mc_repetitions = 1L,
  reference_chunk_size = 5000L,
  aumvc_alpha_grid = seq(0.9, 0.999, by = 0.0001),
  real_datasets = c("adult", "http", "pima", "smtp", "wilt"),
  max_anomaly_fraction = 0.10,
  minimum_hit_count = 10L,
  maximum_relative_mc_se = 0.10,
  confidence_level = 0.95,
  split_fractions = c(
    detector_train = 0.40,
    reference = 0.20,
    aumvc = 0.20,
    label_eval = 0.20
  ),
  split_maximums = c(
    detector_train = 10000L,
    reference = 20000L,
    aumvc = 20000L,
    label_eval = 50000L
  ),
  detectors = list(
    ocsvm = list(nu = 0.5),
    lof = list(k = 20L),
    iforest = list(ntrees = 100L, sample_size = 256L)
  ),
  synthetic = list(
    n_normal = 1800L,
    n_anomaly = 200L,
    radius_min = 1.8,
    radius_max = 3.2
  )
)

experiment_run_seed <- function(settings, run, offset = 0L) {
  settings$seed + (as.integer(run) - 1L) * 1000L + as.integer(offset)
}
