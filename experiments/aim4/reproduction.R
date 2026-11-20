source("detectors/ocsvm.R")
source("aumvc/input_validation.R")
source("aumvc/level_set.R")

REPRODUCTION_SETTINGS <- list(
  seed = 1234L,
  nu = 0.4,
  alpha = 0.95,
  mass_grid = seq(0.91, 0.99, length.out = 10L),
  n_mc = 10000L,
  train_fraction = 0.8,
  sigma_boston = seq(0.01, 4, length.out = 30L),
  sigma_gaussian_mixture = seq(0.01, 3, length.out = 30L),
  boston_models = 25L,
  boston_calibration_tolerance = 0.02,
  solver_tolerance = 1e-3,
  gaussian_mixture_models = 10L,
  gaussian_mixture_n = 1000L,
  truth_n = 1000000L,
  true_mass_n = 100000L,
  output = "experiments/aim4/reproduction"
)

reproduction_normality <- function(model, x) {
  -score_ocsvm(model, x)
}

reproduction_cutoff <- function(scores, mass) {
  scores <- sort(as.numeric(scores), decreasing = TRUE)
  k <- max(1L, min(length(scores), ceiling(mass * length(scores))))
  if (k == length(scores) || scores[k] == scores[k + 1L]) return(scores[k])
  (scores[k] + scores[k + 1L]) / 2
}

reproduction_splits <- function(n, count, fraction, seed) {
  set.seed(seed)
  replicate(count, {
    train <- sample.int(n, floor(n * fraction))
    list(train = train, calibration = setdiff(seq_len(n), train))
  }, simplify = FALSE)
}

reproduction_reference <- function(x, n_mc, seed) {
  box <- fit_reference_box(x)
  set.seed(seed)
  u <- matrix(runif(n_mc * ncol(x)), ncol = ncol(x))
  draws <- sweep(sweep(u, 2L, box$width, "*"), 2L, box$lower, "+")
  list(points = draws, box = box)
}

reproduction_fit <- function(x, sigma, splits, reference_points, masses, nu,
                             keep_models = FALSE, evaluate_on_data = FALSE,
                             solver_tolerance = REPRODUCTION_SETTINGS$solver_tolerance) {
  n_models <- length(splits)
  n_masses <- length(masses)
  thresholds <- matrix(NA_real_, nrow = n_models, ncol = n_masses)
  achieved <- matrix(NA_real_, nrow = n_models, ncol = n_masses)
  reference_mean <- numeric(nrow(reference_points))
  reference_centered <- matrix(0, nrow(reference_points), n_masses)
  data_centered <- if (evaluate_on_data) {
    matrix(0, nrow(x), n_masses)
  } else NULL
  models <- if (keep_models) vector("list", n_models) else NULL
  gamma <- 1 / (2 * sigma^2)
  for (b in seq_along(splits)) {
    split <- splits[[b]]
    model <- fit_ocsvm(x[split$train, , drop = FALSE], nu = nu, gamma = gamma,
                       tolerance = solver_tolerance)
    scores <- reproduction_normality(model, x[split$calibration, , drop = FALSE])
    thresholds[b, ] <- vapply(masses, function(a) reproduction_cutoff(scores, a), numeric(1))
    achieved[b, ] <- vapply(thresholds[b, ], function(t) mean(scores >= t), numeric(1))
    ref_score <- reproduction_normality(model, reference_points)
    reference_mean <- reference_mean + ref_score / n_models
    if (evaluate_on_data) data_score <- reproduction_normality(model, x)
    for (j in seq_along(masses)) {
      reference_centered[, j] <- reference_centered[, j] +
        (ref_score - thresholds[b, j]) / n_models
      if (evaluate_on_data) {
        data_centered[, j] <- data_centered[, j] +
          (data_score - thresholds[b, j]) / n_models
      }
    }
    if (keep_models) models[[b]] <- model
  }
  cutoff <- colMeans(thresholds)
  occupancy <- colMeans(reference_centered >= 0)
  target_error <- abs(sweep(achieved, 2L, masses, "-"))
  scale <- max(1, abs(cutoff))
  rounded_reference <- signif(reference_mean, digits = 10L)
  reference_flat_rate <- max(tabulate(match(rounded_reference,
    unique(rounded_reference)))) / length(rounded_reference)
  list(
    sigma = sigma,
    gamma = gamma,
    masses = masses,
    cutoff = cutoff,
    occupancy = occupancy,
    reference_scores = reference_mean,
    reference_inside_95 = reference_centered[, match(0.95, masses)] >= 0,
    models = models,
    calibration_mass = colMeans(achieved),
    calibration_max_error = max(target_error),
    calibration_mean_error = mean(target_error),
    calibration_mass_90 = mean(achieved[, match(0.90, masses)]),
    calibration_mass_95 = mean(achieved[, match(0.95, masses)]),
    calibration_failure_fraction = mean(target_error > 0.02),
    threshold_collapsed = diff(range(cutoff)) <= 1e-10 * scale,
    reference_flat_rate = reference_flat_rate,
    reference_collapsed = reference_flat_rate >= 0.95,
    full_box = any(occupancy >= 1 - 1 / nrow(reference_points)),
    empirical_mass = if (evaluate_on_data) colMeans(data_centered >= 0) else NULL,
    nested = all(diff(cutoff) <= 1e-10)
  )
}

reproduction_predict <- function(fitted, x) {
  scores <- numeric(nrow(x))
  for (model in fitted$models) {
    scores <- scores + reproduction_normality(model, x) / length(fitted$models)
  }
  scores
}

reproduction_tune <- function(x, sigmas, splits, reference, masses, nu,
                              solver_tolerance) {
  all_masses <- sort(unique(c(masses, 0.90, 0.95)))
  index <- match(masses, all_masses)
  rows <- lapply(seq_along(sigmas), function(i) {
    fitted <- reproduction_fit(x, sigmas[i], splits, reference$points,
                               all_masses, nu, solver_tolerance = solver_tolerance)
    valid <- is.finite(fitted$calibration_max_error) &&
      fitted$calibration_max_error <= 0.02 &&
      !fitted$threshold_collapsed && !fitted$reference_collapsed && !fitted$full_box
    data.frame(
      sigma = sigmas[i], gamma = fitted$gamma,
      mv_area = trapezoid_area(masses,
                               fitted$occupancy[index] * reference$box$volume),
      volume_95 = fitted$occupancy[match(0.95, all_masses)] * reference$box$volume,
      calibration_max_error = fitted$calibration_max_error,
      reference_collapsed = fitted$reference_collapsed,
      threshold_collapsed = fitted$threshold_collapsed,
      full_box = fitted$full_box,
      valid_calibration = valid
    )
  })
  table <- do.call(rbind, rows)
  if (any(!is.finite(table$mv_area))) stop("Nonfinite MV area.")
  chosen <- which.min(table$mv_area)
  table$selected <- seq_len(nrow(table)) == chosen
  selected <- reproduction_fit(x, sigmas[chosen], splits, reference$points,
                               all_masses, nu, keep_models = TRUE,
                               solver_tolerance = solver_tolerance)
  list(table = table, model = selected)
}

reproduction_gaussian_mixture <- function(n) {
  group <- sample.int(2L, n, replace = TRUE)
  center <- ifelse(group == 1L, 2.5, 7.5)
  cbind(rnorm(n, center), rnorm(n, center))
}

reproduction_gaussian_mixture_density <- function(x) {
  0.5 * dnorm(x[, 1L], 2.5) * dnorm(x[, 2L], 2.5) +
    0.5 * dnorm(x[, 1L], 7.5) * dnorm(x[, 2L], 7.5)
}

reproduction_boston <- function(path = "dataset/boston/Boston_housing.data") {
  if (!file.exists(path)) stop("Boston Housing data not found: ", path)
  values <- scan(path, what = numeric(), quiet = TRUE)
  if (length(values) != 506L * 14L || any(!is.finite(values))) {
    stop("Boston Housing file must have 506 rows and 14 numeric columns.")
  }
  raw <- matrix(values, ncol = 14L, byrow = TRUE)
  x <- raw[, c(6L, 13L), drop = FALSE]
  center <- colMeans(x)
  scale <- sqrt(colMeans(sweep(x, 2L, center, "-")^2))
  if (any(scale <= 0)) stop("Zero-variance Boston feature.")
  sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
}

reproduction_true_mass <- function(fitted, x, chunk_size = 5000L) {
  n <- nrow(x)
  inside <- logical(n)
  cutoff <- fitted$cutoff[match(0.95, fitted$masses)]
  for (start in seq.int(1L, n, by = chunk_size)) {
    index <- start:min(start + chunk_size - 1L, n)
    inside[index] <- reproduction_predict(fitted, x[index, , drop = FALSE]) >= cutoff
  }
  mean(inside)
}

run_boston_reproduction <- function(settings = REPRODUCTION_SETTINGS) {
  x <- reproduction_boston()
  splits <- reproduction_splits(nrow(x), settings$boston_models,
                                settings$train_fraction, settings$seed + 1L)
  reference <- reproduction_reference(x, settings$n_mc, settings$seed + 2L)
  fitted <- reproduction_tune(x, settings$sigma_boston, splits, reference,
                              settings$mass_grid, settings$nu, settings$solver_tolerance)
  folder <- file.path(settings$output, "boston")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(fitted$table, file.path(folder, "bandwidth.csv"), row.names = FALSE)
  cat(sprintf("Boston: %d sigmas, %d models; selected sigma %.5f\n",
              length(settings$sigma_boston), length(splits),
              fitted$table$sigma[fitted$table$selected]))
  invisible(fitted$table)
}

run_gaussian_mixture_reproduction <- function(settings = REPRODUCTION_SETTINGS) {
  set.seed(settings$seed + 10L)
  x <- reproduction_gaussian_mixture(settings$gaussian_mixture_n)
  splits <- reproduction_splits(nrow(x), settings$gaussian_mixture_models,
                                settings$train_fraction, settings$seed + 11L)
  reference <- reproduction_reference(x, settings$n_mc, settings$seed + 12L)
  fitted <- reproduction_tune(x, settings$sigma_gaussian_mixture, splits,
                              reference, settings$mass_grid, settings$nu,
                              settings$solver_tolerance)
  folder <- file.path(settings$output, "gaussian_mixture")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(fitted$table, file.path(folder, "bandwidth.csv"), row.names = FALSE)

  set.seed(settings$seed + 13L)
  cutoff_sample <- reproduction_gaussian_mixture(settings$truth_n)
  density_cutoff <- as.numeric(stats::quantile(
    reproduction_gaussian_mixture_density(cutoff_sample),
    probs = 1 - settings$alpha, names = FALSE
  ))
  truth <- reproduction_gaussian_mixture_density(reference$points) >= density_cutoff
  estimated <- fitted$model$reference_inside_95
  mismatch <- xor(truth, estimated)
  mismatch_rate <- mean(mismatch)
  true_volume <- reference$box$volume * mean(truth)
  if (true_volume <= 0) stop("True set has no Monte Carlo hits.")

  set.seed(settings$seed + 14L)
  mass_sample <- reproduction_gaussian_mixture(settings$true_mass_n)
  estimated_mass <- reproduction_true_mass(fitted$model, mass_sample)
  box <- reference$box
  inside_box <- rowSums(sweep(mass_sample, 2L, box$lower, ">=") &
                        sweep(mass_sample, 2L, box$upper, "<=")) == ncol(mass_sample)
  selected_sigma <- fitted$table$sigma[fitted$table$selected]
  comparison <- data.frame(
    sigma = selected_sigma, target_mass = settings$alpha,
    true_mass = estimated_mass,
    true_mass_se = sqrt(estimated_mass * (1 - estimated_mass) / nrow(mass_sample)),
    symmetric_difference_volume_in_box = reference$box$volume * mismatch_rate,
    symmetric_difference_mc_se = reference$box$volume *
      sqrt(mismatch_rate * (1 - mismatch_rate) / length(mismatch)),
    relative_symmetric_difference = mismatch_rate / mean(truth),
    reference_box_mass = mean(inside_box)
  )
  utils::write.csv(comparison, file.path(folder, "set_comparison.csv"), row.names = FALSE)
  cat(sprintf("Gaussian mixture: selected sigma %.5f; true mass %.4f\n",
              selected_sigma, estimated_mass))
  invisible(comparison)
}

run_reproduction <- function(settings = REPRODUCTION_SETTINGS) {
  studies <- trimws(strsplit(Sys.getenv("REPRODUCTION_STUDIES", "boston,gaussian_mixture"),
                             ",", fixed = TRUE)[[1L]])
  if (any(!studies %in% c("boston", "gaussian_mixture"))) stop("Unknown study.")
  if ("boston" %in% studies) run_boston_reproduction(settings)
  if ("gaussian_mixture" %in% studies) run_gaussian_mixture_reproduction(settings)
}

if (sys.nframe() == 0L) run_reproduction()
