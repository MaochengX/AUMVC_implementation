AIM2_SETTINGS <- list(
  seed = 1234L,
  n_runs = 2L,
  n_reference = 20000L,
  n_mc_repetitions = 5L,
  aumvc_alpha_grid = seq(0.9, 0.999, by = 0.0001),
  real_datasets = c("adult", "http", "pima", "smtp", "wilt"),
  concordance_tolerance = 1e-8,
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
  ),
  split_counts = list(
    synthetic = c(500L, 400L, 500L, 600L),
    adult = c(1500L, 2000L, 2000L, 20000L),
    http = c(1500L, 2000L, 2000L, 50000L),
    pima = c(200L, 150L, 250L, 168L),
    smtp = c(1500L, 2000L, 2000L, 89656L),
    wilt = c(1200L, 1200L, 1200L, 1239L)
  )
)

experiment_run_seed <- function(settings, run, offset = 0L) {
  settings$seed + (as.integer(run) - 1L) * 1000L + as.integer(offset)
}
