AIM3_SETTINGS <- list(
  seed = 1234L,
  n_runs = 2L,
  synthetic = list(
    n = 1000L,
    ambient_dims = c(80L, 200L, 400L),
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
    embedding = 200L,
    detector_train = 200L,
    reference = 200L,
    evaluation = 200L,
    label_eval = 200L
  ),
  detector = list(nu = 0.5),
  embedding = list(method = "mds"),
  goix_subsampling = list(
    n_subsets = 50L,
    subset_dim = 5L
  ),
  real = list(
    datasets = c("ecg200", "fashion_mnist", "shuttle"),
    embedding_dim = 5L,
    ecg200 = list(
      data_dir = "dataset/ECG200"
    ),
    fashion_mnist = list(
      data_dir = "dataset/fashion",
      sample_size = 1000L,
      normal_class = 0L,
      anomaly_fraction = 0.10
    ),
    shuttle = list(
      data_dir = "dataset/shuttle",
      sample_size = 1000L,
      excluded_class = 4L
    )
  ),
  n_reference = 100000L,
  n_mc_repetitions = 5L,
  aumvc_alpha_grid = seq(0.9, 0.999, by = 0.0001)
)
