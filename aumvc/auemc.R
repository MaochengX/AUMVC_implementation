trapezoid_area_em <- function(x, y) {
  sum(
    diff(x) *
      (head(y, -1L) + tail(y, -1L)) / 2
  )
}

auemc_from_scores <- function(
    evaluation_scores,
    reference_scores,
    box_volume,
    score_direction = "normality",
    penalty_grid = seq(0, 100, by = 0.01),
    penalty_scale = "normalized",
    em_cutoff = 0.9
) {
  evaluation_scores <- orient_scores(
    evaluation_scores,
    score_direction
  )

  reference_scores <- orient_scores(
    reference_scores,
    score_direction
  )

  level_sets <- level_set_table(
    evaluation_scores,
    reference_scores,
    box_volume
  )

  if (penalty_scale == "normalized") {
    tau <- penalty_grid
    t <- tau / box_volume
  } else {
    t <- penalty_grid
    tau <- t * box_volume
  }

  em <- vapply(
    tau,
    function(current_tau) {
      max(
        0,
        level_sets$mass -
          current_tau * level_sets$occupancy
      )
    },
    numeric(1L)
  )

  em_curve <- data.frame(
    tau = tau,
    t = t,
    em = em
  )

  crossing <- which(em_curve$em <= em_cutoff)

  cutoff <- if (length(crossing) > 0L) {
    crossing[1L]
  } else {
    nrow(em_curve)
  }

  used <- seq_len(cutoff)

  list(
    em_curve = em_curve,
    auemc = trapezoid_area_em(
      em_curve$t[used],
      em_curve$em[used]
    ),
    auemc_normalized = trapezoid_area_em(
      em_curve$tau[used],
      em_curve$em[used]
    )
  )
}

auemc <- function(
    x_eval,
    reference,
    score_fun,
    score_direction = "normality",
    penalty_grid = seq(0, 100, by = 0.01),
    penalty_scale = "normalized",
    em_cutoff = 0.9
) {
  evaluation_scores <- score_fun(x_eval)

  repetitions <- lapply(
    seq_len(reference$n_mc_repetitions),
    function(r) {
      reference_points <- sample_reference_points(
        reference$box,
        reference$n_reference,
        reference$seeds[r]
      )

      reference_scores <- score_fun(reference_points)

      auemc_from_scores(
        evaluation_scores,
        reference_scores,
        reference$box$volume,
        score_direction,
        penalty_grid,
        penalty_scale,
        em_cutoff
      )
    }
  )

  auemc_values <- vapply(
    repetitions,
    function(x) x$auemc,
    numeric(1L)
  )

  normalized_values <- vapply(
    repetitions,
    function(x) x$auemc_normalized,
    numeric(1L)
  )

  em_matrix <- do.call(
    cbind,
    lapply(
      repetitions,
      function(x) x$em_curve$em
    )
  )

  em_curve <- repetitions[[1L]]$em_curve
  em_curve$em <- rowMeans(em_matrix)

  if (reference$n_mc_repetitions > 1L) {
    em_curve$em_mc_sd <- apply(
      em_matrix,
      1L,
      sd
    )

    em_curve$em_mc_se <- em_curve$em_mc_sd /
      sqrt(reference$n_mc_repetitions)

    auemc_mc_sd <- sd(auemc_values)
    auemc_mc_se <- auemc_mc_sd /
      sqrt(reference$n_mc_repetitions)
  } else {
    em_curve$em_mc_sd <- NA_real_
    em_curve$em_mc_se <- NA_real_
    auemc_mc_sd <- NA_real_
    auemc_mc_se <- NA_real_
  }

  list(
    em_curve = em_curve,
    auemc = mean(auemc_values),
    auemc_normalized = mean(normalized_values),
    auemc_mc_sd = auemc_mc_sd,
    auemc_mc_se = auemc_mc_se,
    replicates = data.frame(
      repetition = seq_along(auemc_values),
      seed = reference$seeds,
      auemc = auemc_values,
      auemc_normalized = normalized_values
    )
  )
}
