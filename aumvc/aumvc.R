aumvc_from_scores <- function(
    evaluation_scores,
    reference_scores,
    box_volume,
    score_direction = "normality",
    alpha_grid = seq(0.9, 0.999, by = 0.0001),
    box_log_volume = NULL
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
  threshold <- ordered_scores[pmin(n, ceiling(alpha_grid * n))]

  count_below <- function(sorted_values, value) {
    left <- 1L
    right <- length(sorted_values) + 1L
    while (left < right) {
      middle <- floor((left + right) / 2)
      if (middle <= length(sorted_values) && sorted_values[middle] < value) {
        left <- middle + 1L
      } else {
        right <- middle
      }
    }
    as.integer(left - 1L)
  }
  sorted_reference <- sort(reference_scores)

  hit_count <- length(sorted_reference) - vapply(
    threshold,
    function(value) count_below(sorted_reference, value),
    integer(1)
  )
  occupancy <- hit_count / length(reference_scores)

  volume <- scale_occupancy_volume(
    occupancy,
    box_volume,
    box_log_volume
  )

  if (is.null(box_log_volume)) box_log_volume <- log(box_volume)
  log_volume <- rep(-Inf, length(occupancy))
  positive <- occupancy > 0
  log_volume[positive] <- box_log_volume + log(occupancy[positive])

  aumvc_log <- log_trapezoid_area(alpha_grid, log_volume)

  widths <- diff(alpha_grid)
  weights <- c(
    widths[1L] / 2,
    (head(widths, -1L) + tail(widths, -1L)) / 2,
    tail(widths, 1L) / 2
  )

  threshold_order <- order(threshold)
  ordered_threshold <- threshold[threshold_order]
  cumulative_weight <- cumsum(weights[threshold_order])
  contribution_index <- findInterval(reference_scores, ordered_threshold)
  mc_contribution <- numeric(length(reference_scores))
  positive_index <- contribution_index > 0L
  mc_contribution[positive_index] <- cumulative_weight[
    contribution_index[positive_index]
  ]

  curve <- data.frame(
    alpha = alpha_grid,
    threshold = threshold,
    empirical_mass = vapply(
      threshold,
      function(value) mean(evaluation_scores >= value),
      numeric(1)
    ),
    hit_count = hit_count,
    volume = volume,
    volume_normalized = occupancy
  )

  list(
    mv_curve = curve,
    aumvc = exp_if_representable(aumvc_log),
    aumvc_log = aumvc_log,
    aumvc_normalized = mean(mc_contribution),
    mc_contribution = mc_contribution
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
  combined_scores <- unlist(reference_scores, use.names = FALSE)

  result <- aumvc_from_scores(
    evaluation_scores,
    combined_scores,
    reference$box$volume,
    score_direction,
    alpha_grid,
    reference$box$log_volume
  )
  curve <- result$mv_curve
  occupancy_se <- sqrt(
    curve$volume_normalized * (1 - curve$volume_normalized) /
      length(combined_scores)
  )
  curve$volume_normalized_mc_se <- occupancy_se
  curve$volume_mc_se <- vapply(occupancy_se, function(value) {
    if (value == 0) return(0)
    exp_if_representable(reference$box$log_volume + log(value))
  }, numeric(1))

  normalized_mc_sd <- sd(result$mc_contribution)
  normalized_mc_se <- normalized_mc_sd / sqrt(length(result$mc_contribution))
  mc_sd <- if (normalized_mc_sd == 0) {
    0
  } else {
    exp_if_representable(reference$box$log_volume + log(normalized_mc_sd))
  }
  mc_se <- if (normalized_mc_se == 0) {
    0
  } else {
    exp_if_representable(reference$box$log_volume + log(normalized_mc_se))
  }
  relative_mc_se <- if (result$aumvc_normalized > 0) {
    normalized_mc_se / result$aumvc_normalized
  } else {
    NA_real_
  }
  hit_count_matrix <- matrix(curve$hit_count, ncol = 1L)
  list(
    mv_curve = curve,
    hit_count_matrix = hit_count_matrix,
    mc_contribution = result$mc_contribution,
    aumvc = result$aumvc,
    aumvc_log = result$aumvc_log,
    aumvc_normalized = result$aumvc_normalized,
    aumvc_mc_sd = mc_sd,
    aumvc_mc_se = mc_se,
    aumvc_normalized_mc_sd = normalized_mc_sd,
    aumvc_normalized_mc_se = normalized_mc_se,
    aumvc_relative_mc_se = relative_mc_se,
    minimum_hit_count = min(curve$hit_count),
    zero_occupancy = mean(curve$hit_count == 0),
    low_occupancy = mean(curve$hit_count < 10)
  )
}
