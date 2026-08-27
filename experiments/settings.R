AIM1_SETTINGS <- list(
  seed = 1111L,
  n_reference = 20000L,
  n_mc_repetitions = 5L,
  aumvc_alpha_grid = seq(0.9, 0.999, by = 0.0001),
  auemc_tau_grid = c(0, 10^seq(-6, 6, length.out = 800)),
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

AIM2_SETTINGS <- list(
  seed = 1111L,
  n_runs = 2L,
  n_reference = 20000L,
  n_mc_repetitions = 5L,
  aumvc_alpha_grid = seq(0.9, 0.999, by = 0.0001),
  auemc_tau_grid = c(0, 10^seq(-6, 6, length.out = 800)),
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

AIM3_SETTINGS <- list(
  seed = 1111L,
  n_runs = 1L,
  synthetic = list(
    n = 1000L,
    ambient_dim = 200L,
    intrinsic_dims = c(3L, 5L, 10L),
    snr_levels = c(1, 5, 10),
    truth_functions = c(
      "polynomial_interaction",
      "oscillatory_local"
    ),
    n_distributional = 50L,
    n_structural = 50L,
    distributional_tail = c(0.995, 0.9995),
    structural_shift = 3
  ),
  split_counts = c(
    embedding = 250L,
    detector_train = 250L,
    reference = 250L,
    evaluation = 250L
  ),
  detector = list(
    name = "OCSVM",
    nu = 0.5
  ),
  goix_subsampling = list(
    n_subsets = 50L,
    subset_dim = 5L
  ),
  mds = list(
    interpolation_neighbors = 20L,
    interpolation_ridge = 1e-3
  ),
  reference = list(
    points_per_volume = 1000,
    min_points = 20000L,
    max_points = 100000L,
    n_mc_repetitions = 5L
  ),
  aumvc_alpha_grid = seq(0.9, 0.999, by = 0.0001)
)

AIM4_SETTINGS <- list(
  seed = 1111L
)

aim_seed <- function(settings, offset = 0L) {
  settings$seed + as.integer(offset)
}

aim_run_seed <- function(settings, run) {
  settings$seed + (as.integer(run) - 1L) * 1000L
}
