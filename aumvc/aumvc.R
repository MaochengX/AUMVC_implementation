trapezoid_area <- function(x, y) {
  sum(
    diff(x) *
      (head(y, -1L) + tail(y, -1L)) / 2
  )
}

aumvc_from_scores <- function(
    evaluation_scores,
    reference_scores,
    box_volume,
    score_direction = "normality",
    alpha_grid = seq(0.9, 0.999, length.out = 100L)
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

  mv_curve <- do.call(
    rbind,
    lapply(
      alpha_grid,
      function(alpha) {
        candidates <- which(level_sets$mass >= alpha)

        chosen <- candidates[
          which.min(level_sets$volume[candidates])
        ]

        data.frame(
          alpha = alpha,
          threshold = level_sets$threshold[chosen],
          empirical_mass = level_sets$mass[chosen],
          volume = level_sets$volume[chosen],
          volume_se = level_sets$volume_se[chosen]
        )
      }
    )
  )

  mv_curve$volume_normalized <- mv_curve$volume / box_volume

  list(
    mv_curve = mv_curve,
    level_sets = level_sets,
    aumvc = trapezoid_area(
      mv_curve$alpha,
      mv_curve$volume
    ),
    aumvc_normalized = trapezoid_area(
      mv_curve$alpha,
      mv_curve$volume_normalized
    )
  )
}

aumvc <- function(
    x_eval,
    reference,
    score_fun,
    score_direction = "normality",
    alpha_grid = seq(0.9, 0.999, length.out = 100L)
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

      aumvc_from_scores(
        evaluation_scores,
        reference_scores,
        reference$box$volume,
        score_direction,
        alpha_grid
      )
    }
  )

  aumvc_values <- vapply(
    repetitions,
    function(x) x$aumvc,
    numeric(1L)
  )

  normalized_values <- vapply(
    repetitions,
    function(x) x$aumvc_normalized,
    numeric(1L)
  )

  volume_matrix <- do.call(
    cbind,
    lapply(
      repetitions,
      function(x) x$mv_curve$volume
    )
  )

  normalized_matrix <- do.call(
    cbind,
    lapply(
      repetitions,
      function(x) x$mv_curve$volume_normalized
    )
  )

  mv_curve <- repetitions[[1L]]$mv_curve
  mv_curve$volume <- rowMeans(volume_matrix)
  mv_curve$volume_normalized <- rowMeans(normalized_matrix)

  if (reference$n_mc_repetitions > 1L) {
    mv_curve$volume_mc_sd <- apply(
      volume_matrix,
      1L,
      sd
    )

    mv_curve$volume_mc_se <- mv_curve$volume_mc_sd /
      sqrt(reference$n_mc_repetitions)

    aumvc_mc_sd <- sd(aumvc_values)
    aumvc_mc_se <- aumvc_mc_sd /
      sqrt(reference$n_mc_repetitions)
  } else {
    mv_curve$volume_mc_sd <- NA_real_
    mv_curve$volume_mc_se <- NA_real_
    aumvc_mc_sd <- NA_real_
    aumvc_mc_se <- NA_real_
  }

  list(
    mv_curve = mv_curve,
    aumvc = mean(aumvc_values),
    aumvc_normalized = mean(normalized_values),
    aumvc_mc_sd = aumvc_mc_sd,
    aumvc_mc_se = aumvc_mc_se,
    replicates = data.frame(
      repetition = seq_along(aumvc_values),
      seed = reference$seeds,
      aumvc = aumvc_values,
      aumvc_normalized = normalized_values
    )
  )
}
