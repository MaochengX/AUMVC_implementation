auemc_from_scores <- function(
    evaluation_scores,
    reference_scores,
    box_volume,
    score_direction = "normality",
    tau_grid = c(0, 10^seq(-6, 6, length.out = 800))
) {
  tau_grid <- validate_tau_grid(tau_grid)

  evaluation_scores <- orient_scores(
    validate_scores(evaluation_scores),
    score_direction
  )

  reference_scores <- orient_scores(
    validate_scores(reference_scores),
    score_direction
  )

  sets <- level_set_table(
    evaluation_scores,
    reference_scores,
    box_volume
  )

  mass <- c(0, sets$mass)
  occupancy <- c(0, sets$occupancy)

  em <- vapply(
    tau_grid,
    function(tau) max(mass - tau * occupancy),
    numeric(1)
  )

  crossing <- which(em <= 0.9)[1L]

  if (is.na(crossing)) {
    stop("tau_grid does not reach EM = 0.9.", call. = FALSE)
  }

  if (crossing == 1L) {
    tau_cutoff <- 0
    area_tau <- 0
  } else {
    i <- crossing - 1L

    tau_cutoff <- tau_grid[i] +
      (0.9 - em[i]) *
      (tau_grid[crossing] - tau_grid[i]) /
      (em[crossing] - em[i])

    tau_area <- c(tau_grid[seq_len(i)], tau_cutoff)
    em_area <- c(em[seq_len(i)], 0.9)

    area_tau <- trapezoid_area(tau_area, em_area)
  }

  list(
    em_curve = data.frame(
      tau = tau_grid,
      t = tau_grid / box_volume,
      em = em
    ),
    tau_cutoff = tau_cutoff,
    t_cutoff = tau_cutoff / box_volume,
    auemc = area_tau / box_volume,
    auemc_normalized = area_tau
  )
}

auemc <- function(
    x_eval,
    reference,
    score_fun,
    score_direction = "normality",
    tau_grid = c(0, 10^seq(-6, 6, length.out = 800))
) {
  x_eval <- validate_matrix(x_eval, "x_eval")

  if (!is.function(score_fun)) stop("score_fun must be a function", call. = FALSE)

  evaluation_scores <- validate_scores(
    score_fun(x_eval),
    nrow(x_eval),
    "evaluation_scores"
  )

  reference_scores <- score_reference_repetitions(reference, score_fun)

  repetitions <- lapply(
    reference_scores,
    function(scores) {
      auemc_from_scores(
        evaluation_scores,
        scores,
        reference$box$volume,
        score_direction,
        tau_grid
      )
    }
  )

  values <- vapply(repetitions, function(x) x$auemc, numeric(1))
  normalized <- vapply(
    repetitions,
    function(x) x$auemc_normalized,
    numeric(1)
  )

  em_matrix <- do.call(
    cbind,
    lapply(repetitions, function(x) x$em_curve$em)
  )

  curve <- repetitions[[1L]]$em_curve
  curve$em <- rowMeans(em_matrix)

  if (length(values) > 1L) {
    curve$em_mc_se <- apply(em_matrix, 1L, sd) /
      sqrt(length(values))

    mc_sd <- sd(values)
    mc_se <- mc_sd / sqrt(length(values))
  } else {
    curve$em_mc_se <- NA_real_
    mc_sd <- NA_real_
    mc_se <- NA_real_
  }

  list(
    em_curve = curve,
    auemc = mean(values),
    auemc_normalized = mean(normalized),
    auemc_mc_sd = mc_sd,
    auemc_mc_se = mc_se,
    t_cutoff = mean(
      vapply(repetitions, function(x) x$t_cutoff, numeric(1))
    )
  )
}
