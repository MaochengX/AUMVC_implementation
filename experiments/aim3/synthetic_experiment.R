source("experiments/aim3/settings.R")
source("experiments/aim3/utils.R")
source("experiments/aim3/representation_comparison.R")
source("experiments/aim3/analysis.R")

unit_rows <- function(x) {
  norms <- sqrt(rowSums(x^2))
  norms[norms == 0] <- 1
  x / norms
}

next_coordinate <- function(q) c(seq.int(2L, q), 1L)

polynomial_basis <- function(z) {
  next_index <- next_coordinate(ncol(z))
  cbind(z, z^2, z^3, z * z[, next_index, drop = FALSE])
}

oscillatory_basis <- function(z) {
  next_index <- next_coordinate(ncol(z))
  cbind(
    sin(z), cos(2 * z), tanh(1.5 * z), exp(-(z^2) / 2),
    sin(z * z[, next_index, drop = FALSE])
  )
}

basis <- function(z, truth) {
  switch(
    truth,
    polynomial_interaction = polynomial_basis(z),
    oscillatory_local = oscillatory_basis(z),
    stop("Unknown truth function.", call. = FALSE)
  )
}

fit_basis_scaler <- function(phi) {
  scale <- apply(phi, 2L, sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(center = colMeans(phi), scale = scale)
}

apply_basis_scaler <- function(phi, scaler) {
  phi <- sweep(phi, 2L, scaler$center, "-")
  sweep(phi, 2L, scaler$scale, "/")
}

distributional_latent <- function(n, q, tail, seed) {
  set.seed(seed)
  direction <- unit_rows(matrix(rnorm(n * q), nrow = n))
  radius <- sqrt(qchisq(runif(n, tail[1L], tail[2L]), df = q))
  direction * radius
}

structural_basis <- function(z, w, truth, severity) {
  if (!identical(dim(z), dim(w)) || !is.finite(severity) || severity < 0 || severity > 1) {
    stop("Invalid structural settings.", call. = FALSE)
  }
  phi <- basis(z, truth)
  alternative <- basis(w, truth)
  q <- ncol(z)
  columns <- seq.int(ncol(phi) - q + 1L, ncol(phi))
  phi[, columns] <- (1 - severity) * phi[, columns, drop = FALSE] +
    severity * alternative[, columns, drop = FALSE]
  phi
}

generate_synthetic <- function(intrinsic_dim, ambient_dim, snr, truth, seed, settings) {
  q <- as.integer(intrinsic_dim)
  d <- as.integer(ambient_dim)
  n_distributional <- as.integer(settings$n_distributional)
  n_structural <- as.integer(settings$n_structural)
  n_normal <- as.integer(settings$n) - n_distributional - n_structural
  if (n_normal < 1L || d <= 5L * q) stop("Invalid synthetic dimensions.", call. = FALSE)
  if (!q %in% settings$intrinsic_dims || !d %in% settings$ambient_dims ||
      !snr %in% settings$snr_levels || !truth %in% settings$truth_functions) {
    stop("Unsupported synthetic setting.", call. = FALSE)
  }
  set.seed(seed)
  z_normal <- matrix(rnorm(n_normal * q), nrow = n_normal)
  z_distributional <- distributional_latent(
    n_distributional, q, settings$distributional_tail, seed + 1L
  )
  set.seed(seed + 2L)
  z_structural <- matrix(rnorm(n_structural * q), nrow = n_structural)
  set.seed(seed + 3L)
  w_structural <- matrix(rnorm(n_structural * q), nrow = n_structural)
  normal_basis <- basis(z_normal, truth)
  scaler <- fit_basis_scaler(normal_basis)
  normal_basis <- apply_basis_scaler(normal_basis, scaler)
  distributional_basis <- apply_basis_scaler(basis(z_distributional, truth), scaler)
  structural <- apply_basis_scaler(
    structural_basis(z_structural, w_structural, truth, settings$structural_severity),
    scaler
  )
  set.seed(seed + 4L)
  coefficients <- matrix(rnorm(ncol(normal_basis) * d), nrow = ncol(normal_basis)) /
    sqrt(ncol(normal_basis))
  signal_normal <- normal_basis %*% coefficients
  signal_distributional <- distributional_basis %*% coefficients
  signal_structural <- structural %*% coefficients
  noise_sd <- sqrt(mean(apply(signal_normal, 2L, var)) / snr)
  set.seed(seed + 5L)
  noise <- matrix(rnorm(settings$n * d, sd = noise_sd), nrow = settings$n)
  types <- c(
    rep("normal", n_normal),
    rep("distributional", n_distributional),
    rep("structural", n_structural)
  )
  list(
    x = rbind(signal_normal, signal_distributional, signal_structural) + noise,
    labels = as.integer(types != "normal"),
    outlier_type = types,
    intrinsic_dim = q
  )
}

synthetic_grid <- function(settings) {
  expand.grid(
    truth = settings$truth_functions,
    intrinsic_dim = settings$intrinsic_dims,
    original_dim = settings$ambient_dims,
    snr = settings$snr_levels,
    stringsAsFactors = FALSE
  )
}

synthetic_seed <- function(row, run, settings, base_seed) {
  base_seed + run * 100000L +
    match(row$truth, settings$truth_functions) * 10000L +
    as.integer(row$intrinsic_dim) * 100L +
    match(row$snr, settings$snr_levels) * 10L
}

synthetic_result_columns <- c(
  "run", "truth", "intrinsic_dim", "snr", "sample_size", "original_dim",
  "embedding_dim", "detector_dim", "representation", "gamma_scale",
  "aumvc_normalized", "aumvc_relative_mc_se", "min_hit_count",
  "zero_occupancy", "low_occupancy", "box_log_volume", "roc_auc", "pr_auc",
  "roc_distributional", "pr_distributional", "roc_structural", "pr_structural",
  "embedding_oos_stress", "runtime_sec", "aumvc_reliable"
)

synthetic_bandwidth_columns <- c(
  "run", "truth", "intrinsic_dim", "snr", "sample_size", "original_dim",
  "embedding_dim", "detector_dim", "representation", "compared_pairs_count",
  "unreliable_pairs_rate", "selected_gamma_scale", "agreement_pairs_rank_rate",
  "selected_gamma_scale_reliable", "reliable_agreement_pairs_rank_rate", "tie_rate"
)

evaluate_synthetic_task <- function(task) {
  row <- task$row
  run <- task$run
  settings <- task$settings
  seed <- synthetic_seed(row, run, settings$synthetic, settings$seed)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  case <- generate_synthetic(
    row$intrinsic_dim, row$original_dim, row$snr, row$truth,
    seed, settings$synthetic
  )
  evaluated <- evaluate_case(case, seed + 50000L, settings)
  metadata <- data.frame(
    run = run,
    truth = row$truth,
    intrinsic_dim = row$intrinsic_dim,
    snr = row$snr,
    sample_size = nrow(case$x),
    original_dim = row$original_dim,
    stringsAsFactors = FALSE
  )
  list(
    results = attach_metadata(metadata, evaluated$results),
    bandwidth = attach_metadata(metadata, evaluated$bandwidth)
  )
}

run_synthetic_experiment <- function(settings = AIM3_SETTINGS) {
  n_runs <- validate_runs(settings$n_runs)
  grid <- synthetic_grid(settings$synthetic)
  total <- n_runs * nrow(grid)
  results <- vector("list", total)
  bandwidth <- vector("list", total)
  dirs <- result_dirs("synthetic")
  done <- 0L
  start <- proc.time()[[3L]]
  cat("\naim3 synthetic\n")
  parallel <- isTRUE(settings$parallel$enabled)
  workers <- if (parallel) as.integer(settings$parallel$workers) else 1L
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("parallel$workers must be a positive integer.", call. = FALSE)
  }
  if (parallel && workers > 1L) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Install 'future' and 'future.apply' for parallel execution.", call. = FALSE)
    }
    previous_plan <- future::plan()
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(previous_plan), add = TRUE)
  }
  for (run in seq_len(n_runs)) {
    jobs <- lapply(seq_len(nrow(grid)), function(i) {
      list(row = grid[i, , drop = FALSE], run = run, settings = settings)
    })
    cases <- if (parallel && workers > 1L) {
      future.apply::future_lapply(jobs, evaluate_synthetic_task, future.seed = TRUE)
    } else {
      lapply(jobs, evaluate_synthetic_task)
    }
    index <- done + seq_along(cases)
    results[index] <- lapply(cases, `[[`, "results")
    bandwidth[index] <- lapply(cases, `[[`, "bandwidth")
    done <- done + length(cases)
    elapsed <- proc.time()[[3L]] - start
    cat(sprintf("\r%d/%d | elapsed %.0fs | eta %.0fs", done, total, elapsed, elapsed * (total - done) / done))
    flush.console()
    current_results <- do.call(rbind, results[seq_len(done)])
    current_bandwidth <- do.call(rbind, bandwidth[seq_len(done)])
    current_results <- current_results[, synthetic_result_columns, drop = FALSE]
    current_bandwidth <- current_bandwidth[, synthetic_bandwidth_columns, drop = FALSE]
    write_report(current_results, file.path(dirs$experiment, "syn_results.csv"), settings$output_digits)
    write_report(current_bandwidth, file.path(dirs$experiment, "syn_bandwidth_exp.csv"), settings$output_digits)
  }
  results <- do.call(rbind, results)
  bandwidth <- do.call(rbind, bandwidth)
  results <- results[, synthetic_result_columns, drop = FALSE]
  bandwidth <- bandwidth[, synthetic_bandwidth_columns, drop = FALSE]
  write_report(results, file.path(dirs$experiment, "syn_results.csv"), settings$output_digits)
  write_report(bandwidth, file.path(dirs$experiment, "syn_bandwidth_exp.csv"), settings$output_digits)
  analyze_results(results, bandwidth, "synthetic", dirs$analysis, settings$output_digits)
  cat("\ncomplete: ", dirname(dirs$experiment), "\n", sep = "")
  invisible(list(results = results, bandwidth = bandwidth))
}

if (sys.nframe() == 0L) run_synthetic_experiment()
