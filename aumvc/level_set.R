orient_scores <- function(scores, direction) {
  if (direction == "normality") return(as.numeric(scores))
  if (direction == "anomaly") return(-as.numeric(scores))
  stop("score_direction must be 'normality' or 'anomaly'", call. = FALSE)
}

trapezoid_area <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

fit_reference_box <- function(x_reference) {
  x_reference <- validate_matrix(x_reference, "x_reference")
  lower <- apply(x_reference, 2L, min)
  upper <- apply(x_reference, 2L, max)
  width <- upper - lower

  if (any(width <= 0)) {
    stop("The reference data have a constant coordinate.", call. = FALSE)
  }

  volume <- prod(width)
  log_volume <- sum(log(width))
  if (!is.finite(volume) && log_volume <= log(.Machine$double.xmax)) {
    volume <- exp(log_volume)
  }

  list(
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    width = as.numeric(width),
    volume = volume,
    log_volume = log_volume,
    dimension = ncol(x_reference)
  )
}

sample_reference_points <- function(box, n_reference, seed) {
  set.seed(seed)
  x <- matrix(runif(n_reference * box$dimension), nrow = n_reference)
  x <- sweep(x, 2L, box$width, "*")
  sweep(x, 2L, box$lower, "+")
}

make_reference <- function(
    x_reference,
    n_reference = 20000L,
    n_mc_repetitions = 5L,
    seed = 2030L
) {
  n_reference <- as.integer(n_reference)
  n_mc_repetitions <- as.integer(n_mc_repetitions)
  if (n_reference < 1L || n_mc_repetitions < 1L) {
    stop("Reference sample sizes must be positive.", call. = FALSE)
  }

  list(
    box = fit_reference_box(x_reference),
    n_reference = n_reference,
    n_mc_repetitions = n_mc_repetitions,
    seeds = seed + seq_len(n_mc_repetitions) - 1L
  )
}

score_reference_repetitions <- function(reference, score_fun) {
  lapply(seq_len(reference$n_mc_repetitions), function(r) {
    points <- sample_reference_points(
      reference$box,
      reference$n_reference,
      reference$seeds[r]
    )
    validate_scores(score_fun(points), reference$n_reference, "reference_scores")
  })
}

scale_occupancy_volume <- function(occupancy, box_volume, box_log_volume = NULL) {
  occupancy <- as.numeric(occupancy)
  if (is.finite(box_volume)) return(box_volume * occupancy)
  if (is.null(box_log_volume)) box_log_volume <- log(box_volume)

  volume <- numeric(length(occupancy))
  positive <- occupancy > 0
  if (!any(positive)) return(volume)

  log_volume <- box_log_volume + log(occupancy[positive])
  finite <- log_volume <= log(.Machine$double.xmax)
  indices <- which(positive)
  volume[indices[finite]] <- exp(log_volume[finite])
  volume[indices[!finite]] <- Inf
  volume
}

level_set_table <- function(
    evaluation_scores,
    reference_scores,
    box_volume,
    box_log_volume = NULL
) {
  thresholds <- sort(unique(c(evaluation_scores, reference_scores)), decreasing = TRUE)
  eval_bin <- tabulate(match(evaluation_scores, thresholds), nbins = length(thresholds))
  ref_bin <- tabulate(match(reference_scores, thresholds), nbins = length(thresholds))
  mass <- cumsum(eval_bin) / length(evaluation_scores)
  occupancy <- cumsum(ref_bin) / length(reference_scores)

  data.frame(
    threshold = thresholds,
    mass = mass,
    occupancy = occupancy,
    volume = scale_occupancy_volume(occupancy, box_volume, box_log_volume)
  )
}
