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

reproduction_tune <- function(x, sigmas, splits, reference, masses, nu, name,
                              truth = NULL, check_calibration = FALSE,
                              calibration_tolerance = 0.02,
                              diagnostic_path = NULL,
                              solver_tolerance = REPRODUCTION_SETTINGS$solver_tolerance) {
  all_masses <- sort(unique(c(masses, 0.90, 0.95)))
  target <- match(0.95, all_masses)
  grid_indices <- match(masses, all_masses)
  results <- lapply(seq_along(sigmas), function(i) {
    cat(sprintf("%s: bandwidth %d/%d (sigma %.5f)\n", name, i, length(sigmas), sigmas[i]))
    fitted <- reproduction_fit(x, sigmas[i], splits, reference$points, all_masses, nu,
                               evaluate_on_data = check_calibration,
                               solver_tolerance = solver_tolerance)
    area <- fitted$occupancy[grid_indices] * reference$box$volume
    difference <- if (is.null(truth)) NA_real_ else {
      mean(xor(fitted$reference_inside_95, truth)) / mean(truth)
    }
    valid <- is.finite(fitted$calibration_max_error) &&
      fitted$calibration_max_error <= calibration_tolerance &&
      !fitted$threshold_collapsed && !fitted$reference_collapsed &&
      !fitted$full_box
    data.frame(
      sigma = sigmas[i], gamma = fitted$gamma,
      mv_area = trapezoid_area(masses, area),
      mv_area_normalized = trapezoid_area(masses, fitted$occupancy[grid_indices]),
      volume_95 = fitted$occupancy[target] * reference$box$volume,
      calibrated_symmetric_difference = difference,
      nested = as.numeric(fitted$nested),
      calibration_max_error = fitted$calibration_max_error,
      calibration_mean_error = fitted$calibration_mean_error,
      calibration_failure_fraction = fitted$calibration_failure_fraction,
      calibration_mass_90 = fitted$calibration_mass_90,
      calibration_mass_95 = fitted$calibration_mass_95,
      threshold_collapsed = fitted$threshold_collapsed,
      reference_flat_rate = fitted$reference_flat_rate,
      reference_collapsed = fitted$reference_collapsed,
      full_box = fitted$full_box,
      valid_calibration = valid,
      empirical_mass_90 = if (check_calibration) fitted$empirical_mass[match(0.90, all_masses)] else NA_real_,
      empirical_mass_95 = if (check_calibration) fitted$empirical_mass[match(0.95, all_masses)] else NA_real_
    )
  })
  results <- do.call(rbind, results)
  rownames(results) <- NULL
  if (any(!is.finite(results$mv_area))) stop("Nonfinite MV area in reproduction.")
  ranking <- order(results$mv_area, results$sigma)
  chosen <- ranking[1L]
  results$rank <- match(seq_len(nrow(results)), ranking)
  results$selected <- seq_len(nrow(results)) == chosen
  if (!is.null(diagnostic_path)) {
    dir.create(dirname(diagnostic_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(results, diagnostic_path, row.names = FALSE)
  }
  if (check_calibration && !results$valid_calibration[chosen]) {
    warning(sprintf("%s: paper-selected sigma %.4f fails diagnostic checks; retained for faithful reproduction.",
                    name, sigmas[chosen]), call. = FALSE)
  }
  selected <- reproduction_fit(x, sigmas[chosen], splits, reference$points,
                               all_masses, nu, keep_models = TRUE,
                               evaluate_on_data = check_calibration,
                               solver_tolerance = solver_tolerance)
  list(table = results, chosen = chosen, model = selected)
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
  if (!file.exists(path)) {
    stop(sprintf("Boston Housing data not found at %s", path), call. = FALSE)
  }
  values <- scan(path, what = numeric(), quiet = TRUE)
  if (length(values) != 506L * 14L || anyNA(values) ||
      !all(is.finite(values))) {
    stop("Boston Housing file must contain 506 observations with 14 numeric columns.", call. = FALSE)
  }
  raw <- matrix(values, ncol = 14L, byrow = TRUE)
  x <- raw[, c(6L, 13L), drop = FALSE]
  means <- colMeans(x)
  scales <- sqrt(colMeans(sweep(x, 2L, means, "-")^2))
  if (any(scales <= 0)) stop("Boston features have zero variance.")
  x <- sweep(sweep(x, 2L, means, "-"), 2L, scales, "/")
  colnames(x) <- c("standardized RM", "standardized LSTAT")
  list(x = x, source = path)
}

reproduction_plot_sigma <- function(table, title, target = NULL,
                                    selected = NULL) {
  plot(table$sigma, table$mv_area, type = "l", lwd = 1.4,
       col = "#2878B5", xlab = expression(sigma), ylab = "Area under MV curve",
       main = title)
  if ("valid_calibration" %in% names(table)) {
    points(table$sigma, table$mv_area,
           pch = ifelse(table$valid_calibration, 16L, 1L), cex = 0.7,
           col = ifelse(table$valid_calibration, "#2878B5", "#A44A3F"))
  } else {
    points(table$sigma, table$mv_area, pch = 16L, cex = 0.6, col = "#2878B5")
  }
  j <- if (is.null(selected)) which.min(table$mv_area) else selected
  abline(v = table$sigma[j], lty = 2, col = "#D55E00")
  legend("topright", sprintf("selected sigma = %.4f", table$sigma[j]),
         col = "#D55E00", lty = 2, bty = "n")
  if (!is.null(target)) abline(v = target, col = "gray45", lty = 3)
}

reproduction_plot_sets <- function(x, fitted, masses, title, truth_cutoff = NULL,
                                   standard_model = NULL, standard = FALSE) {
  xr <- range(x[, 1L]); yr <- range(x[, 2L])
  padding_x <- diff(xr) * 0.05; padding_y <- diff(yr) * 0.05
  xx <- seq(xr[1L] - padding_x, xr[2L] + padding_x, length.out = 140L)
  yy <- seq(yr[1L] - padding_y, yr[2L] + padding_y, length.out = 140L)
  grid <- as.matrix(expand.grid(xx, yy))
  score <- if (standard) reproduction_normality(standard_model, grid) else reproduction_predict(fitted, grid)
  z <- matrix(score, nrow = length(xx))
  plot(x, pch = 16, cex = 0.25, col = "gray45", asp = 1,
       xlim = range(xx), ylim = range(yy), xlab = "Coordinate 1", ylab = "Coordinate 2", main = title)
  for (i in seq_along(masses)) {
    cutoff <- if (standard) 0 else fitted$cutoff[which.min(abs(fitted$masses - masses[i]))]
    if (min(z) < cutoff && max(z) > cutoff) {
      contour(xx, yy, z, levels = cutoff, add = TRUE, drawlabels = FALSE,
              col = c("#2878B5", "#D55E00")[i], lwd = 2)
    }
  }
  if (!is.null(truth_cutoff)) {
    density <- matrix(reproduction_gaussian_mixture_density(grid), nrow = length(xx))
    contour(xx, yy, density, levels = truth_cutoff, add = TRUE,
            drawlabels = FALSE, col = "darkgreen", lwd = 2, lty = 2)
  }
  legend("topright", legend = c(paste0("estimated mass ", masses),
           if (!is.null(truth_cutoff)) "true mass 0.95"),
         lty = c(rep(1, length(masses)), if (!is.null(truth_cutoff)) 2),
         col = c(c("#2878B5", "#D55E00")[seq_along(masses)],
                 if (!is.null(truth_cutoff)) "darkgreen"), bty = "n", cex = 0.8)
}

reproduction_top3 <- function(table, study, n, models, standard_oracle = NULL) {
  top <- table[order(table$rank)[seq_len(min(3L, nrow(table)))], , drop = FALSE]
  data.frame(
    study = study, n = n, models = models, rank = top$rank,
    sigma = top$sigma, gamma = top$gamma, gamma_scale = 2 * top$gamma,
    mv_area = top$mv_area, valid_calibration = top$valid_calibration,
    calibration_max_error = top$calibration_max_error,
    empirical_mass_90 = top$empirical_mass_90,
    empirical_mass_95 = top$empirical_mass_95,
    nested = as.logical(top$nested),
    calibrated_symmetric_difference = top$calibrated_symmetric_difference,
    standard_oracle_sigma = if (is.null(standard_oracle)) NA_real_ else standard_oracle$sigma,
    standard_oracle_symmetric_difference = if (is.null(standard_oracle)) NA_real_ else standard_oracle$error
  )
}

run_boston_reproduction <- function(settings = REPRODUCTION_SETTINGS) {
  data <- reproduction_boston()
  x <- data$x
  splits <- reproduction_splits(nrow(x), settings$boston_models,
                                settings$train_fraction, settings$seed + 1L)
  ref <- reproduction_reference(x, settings$n_mc, settings$seed + 2L)
  dir <- file.path(settings$output, "boston")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  fit <- reproduction_tune(x, settings$sigma_boston, splits, ref,
                           settings$mass_grid, settings$nu, "Boston",
                           check_calibration = TRUE,
                           calibration_tolerance = settings$boston_calibration_tolerance,
                           diagnostic_path = file.path(dir, "bandwidth.csv"),
                           solver_tolerance = settings$solver_tolerance)
  selected <- fit$model
  indices <- match(c(0.90, 0.95), selected$masses)
  empirical <- selected$empirical_mass[indices]
  if (any(abs(empirical - c(0.90, 0.95)) > 0.05)) {
    warning("Boston aggregated sets are not near their target masses; inspect the diagnostics.")
  }
  grDevices::pdf(file.path(dir, "figures.pdf"), width = 7.5, height = 6)
  tryCatch({
    reproduction_plot_sigma(fit$table, "Boston: calibrated aggregate MV area", 0.42,
                            selected = fit$chosen)
    reproduction_plot_sets(x, selected, c(0.90, 0.95), "Boston: nested MV sets")
  }, finally = grDevices::dev.off())
  reproduction_top3(fit$table, "boston", nrow(x), length(splits))
}

run_gaussian_mixture_reproduction <- function(settings = REPRODUCTION_SETTINGS) {
  set.seed(settings$seed + 10L)
  x <- reproduction_gaussian_mixture(settings$gaussian_mixture_n)
  splits <- reproduction_splits(nrow(x), settings$gaussian_mixture_models,
                                settings$train_fraction, settings$seed + 11L)
  ref <- reproduction_reference(x, settings$n_mc, settings$seed + 12L)
  set.seed(settings$seed + 13L)
  truth_cutoff <- as.numeric(stats::quantile(reproduction_gaussian_mixture_density(reproduction_gaussian_mixture(settings$truth_n)),
                                               probs = 1 - settings$alpha, names = FALSE))
  truth <- reproduction_gaussian_mixture_density(ref$points) >= truth_cutoff
  truth_volume <- mean(truth)
  if (truth_volume <= 0) stop("True MV set has zero reference occupancy.")
  fit <- reproduction_tune(x, settings$sigma_gaussian_mixture, splits, ref,
                           settings$mass_grid, settings$nu, "Gaussian mixture", truth,
                           solver_tolerance = settings$solver_tolerance)
  errors <- lapply(seq_along(settings$sigma_gaussian_mixture), function(i) {
    sigma <- settings$sigma_gaussian_mixture[i]
    cat(sprintf("Standard OCSVM: bandwidth %d/%d\n", i, length(settings$sigma_gaussian_mixture)))
    ordinary <- fit_ocsvm(x, nu = 1 - settings$alpha, gamma = 1 / (2 * sigma^2),
                          tolerance = settings$solver_tolerance)
    error <- mean(xor(reproduction_normality(ordinary, ref$points) >= 0, truth)) / truth_volume
    data.frame(sigma = sigma, standard_symmetric_difference = error)
  })
  errors <- do.call(rbind, errors)
  table <- merge(fit$table, errors, by = "sigma", sort = FALSE)
  table <- table[match(settings$sigma_gaussian_mixture, table$sigma), , drop = FALSE]
  standard_best <- which.min(table$standard_symmetric_difference)
  standard <- fit_ocsvm(x, nu = 1 - settings$alpha,
                        gamma = 1 / (2 * table$sigma[standard_best]^2),
                        tolerance = settings$solver_tolerance)
  dir <- file.path(settings$output, "gaussian_mixture")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(table, file.path(dir, "bandwidth.csv"), row.names = FALSE)
  grDevices::pdf(file.path(dir, "figures.pdf"), width = 7.5, height = 6)
  tryCatch({
    plot(table$sigma, table$standard_symmetric_difference, type = "b", pch = 16,
         col = "gray35", cex = 0.55, xlab = expression(sigma),
         ylab = "Relative symmetric-difference volume",
         ylim = range(c(table$standard_symmetric_difference,
                        table$calibrated_symmetric_difference), finite = TRUE),
         main = "Gaussian mixture: true-set estimation error")
    lines(table$sigma, table$calibrated_symmetric_difference, type = "b", pch = 16,
          col = "#2878B5", cex = 0.55)
    legend("topright", c("Standard OCSVM", "Calibrated aggregate"),
           col = c("gray35", "#2878B5"), lty = 1, bty = "n")
    reproduction_plot_sets(x, fit$model, settings$alpha,
                           "Gaussian mixture: calibrated aggregate", truth_cutoff)
    reproduction_plot_sets(x, NULL, settings$alpha,
                           "Gaussian mixture: oracle-selected standard OCSVM",
                           truth_cutoff, standard, standard = TRUE)
  }, finally = grDevices::dev.off())
  reproduction_top3(table, "gaussian_mixture", nrow(x), length(splits),
                    standard_oracle = list(
                      sigma = table$sigma[standard_best],
                      error = table$standard_symmetric_difference[standard_best]))
}

run_reproduction <- function(studies = c("boston", "gaussian_mixture"),
                             settings = REPRODUCTION_SETTINGS) {
  if (!length(studies) || !all(studies %in% c("boston", "gaussian_mixture"))) {
    stop("Studies must be boston and/or gaussian_mixture.")
  }
  summaries <- list()
  if ("boston" %in% studies) summaries$boston <- run_boston_reproduction(settings)
  if ("gaussian_mixture" %in% studies) summaries$gaussian_mixture <- run_gaussian_mixture_reproduction(settings)
  all_columns <- unique(unlist(lapply(summaries, names)))
  summaries <- lapply(summaries, function(x) {
    x[setdiff(all_columns, names(x))] <- NA
    x[all_columns]
  })
  summary <- do.call(rbind, summaries)
  summary_path <- file.path(settings$output, "summary.csv")
  if (length(studies) == 1L && file.exists(summary_path)) {
    old <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
    if (all(c("study", "rank", "sigma", "mv_area") %in% names(old))) {
      old <- old[!old$study %in% studies, , drop = FALSE]
      summary <- rbind(summary, old[, names(summary), drop = FALSE])
    } else {
      warning("Old summary.csv used the previous one-row format; run both studies to recreate the top-three summary.",
              call. = FALSE)
    }
  }
  rownames(summary) <- NULL
  utils::write.csv(summary, summary_path, row.names = FALSE)
  print(summary)
  invisible(summary)
}

if (sys.nframe() == 0L) {
  studies <- strsplit(Sys.getenv("REPRODUCTION_STUDIES", "boston,gaussian_mixture"), ",", fixed = TRUE)[[1L]]
  run_reproduction(trimws(studies))
}
