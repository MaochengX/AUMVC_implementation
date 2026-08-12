aumvc_from_scores <- function(
    evaluation_scores,
    reference_scores,
    box_volume,
    score_direction = "normality",
    alpha_grid = seq(0.9, 0.999, by = 0.0001)
) {
  alpha_grid <- validate_alpha_grid(alpha_grid)

  evaluation_scores <- orient_scores(
    validate_scores(evaluation_scores),
    score_direction
  )

  reference_scores <- orient_scores(
    validate_scores(reference_scores),
    score_direction
  )

  n <- length(evaluation_scores)
  ordered_scores <- sort(evaluation_scores, decreasing = TRUE)

  k <- pmin(n, ceiling(alpha_grid * n))
  threshold <- ordered_scores[k]

  occupancy <- vapply(
    threshold,
    function(u) mean(reference_scores >= u),
    numeric(1)
  )

  curve <- data.frame(
    alpha = alpha_grid,
    threshold = threshold,
    empirical_mass = vapply(
      threshold,
      function(u) mean(evaluation_scores >= u),
      numeric(1)
    ),
    volume = box_volume * occupancy,
    volume_normalized = occupancy
  )

  list(
    mv_curve = curve,
    aumvc = trapezoid_area(curve$alpha, curve$volume),
    aumvc_normalized = trapezoid_area(
      curve$alpha,
      curve$volume_normalized
    )
  )
}

aumvc <- function(
    x_eval,
    reference,
    score_fun,
    score_direction = "normality",
    alpha_grid = seq(0.9, 0.999, by = 0.0001)
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
      aumvc_from_scores(
        evaluation_scores,
        scores,
        reference$box$volume,
        score_direction,
        alpha_grid
      )
    }
  )

  values <- vapply(repetitions, function(x) x$aumvc, numeric(1))
  normalized <- vapply(
    repetitions,
    function(x) x$aumvc_normalized,
    numeric(1)
  )

  volume_matrix <- do.call(
    cbind,
    lapply(repetitions, function(x) x$mv_curve$volume)
  )

  curve <- repetitions[[1L]]$mv_curve
  curve$volume <- rowMeans(volume_matrix)
  curve$volume_normalized <- curve$volume / reference$box$volume

  if (length(values) > 1L) {
    curve$volume_mc_se <- apply(volume_matrix, 1L, sd) /
      sqrt(length(values))

    mc_sd <- sd(values)
    mc_se <- mc_sd / sqrt(length(values))
  } else {
    curve$volume_mc_se <- NA_real_
    mc_sd <- NA_real_
    mc_se <- NA_real_
  }

  list(
    mv_curve = curve,
    aumvc = mean(values),
    aumvc_normalized = mean(normalized),
    aumvc_mc_sd = mc_sd,
    aumvc_mc_se = mc_se
  )
}
