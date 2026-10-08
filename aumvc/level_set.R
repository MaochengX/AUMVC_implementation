orient_scores <- function(scores, direction) {
  if (direction == "normality") return(as.numeric(scores))
  if (direction == "anomaly") return(-as.numeric(scores))
  stop("score_direction must be 'normality' or 'anomaly'", call. = FALSE)
}

trapezoid_area <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

log_sum_exp <- function(values) {
  values <- as.numeric(values)
  if (any(values == Inf, na.rm = TRUE)) return(Inf)
  finite <- is.finite(values)
  if (!any(finite)) return(-Inf)
  maximum <- max(values[finite])
  maximum + log(sum(exp(values[finite] - maximum)))
}

log_mean_exp <- function(values) {
  log_sum_exp(values) - log(length(values))
}

log_trapezoid_area <- function(x, log_y) {
  widths <- diff(x) / 2
  terms <- vapply(seq_along(widths), function(index) {
    if (widths[index] <= 0) return(-Inf)
    log(widths[index]) + log_sum_exp(log_y[c(index, index + 1L)])
  }, numeric(1))
  log_sum_exp(terms)
}

exp_if_representable <- function(log_value) {
  if (is.na(log_value)) return(NA_real_)
  if (is.infinite(log_value) && log_value < 0) return(0)
  if (!is.finite(log_value) || log_value > log(.Machine$double.xmax)) return(Inf)
  exp(log_value)
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
  x <- matrix(
    runif(n_reference * box$dimension),
    nrow = n_reference,
    byrow = TRUE
  )
  x <- sweep(x, 2L, box$width, "*")
  sweep(x, 2L, box$lower, "+")
}

make_reference <- function(
    x_reference,
    n_reference = 100000L,
    n_mc_repetitions = 1L,
    seed = 1234L,
    chunk_size = 5000L
) {
  n_reference <- as.integer(n_reference)
  n_mc_repetitions <- as.integer(n_mc_repetitions)
  chunk_size <- as.integer(chunk_size)
  if (
    n_reference < 1L ||
    n_mc_repetitions < 1L ||
    chunk_size < 1L
  ) {
    stop("Reference sample sizes must be positive.", call. = FALSE)
  }

  list(
    box = fit_reference_box(x_reference),
    n_reference = n_reference,
    n_mc_repetitions = n_mc_repetitions,
    seeds = seed + seq_len(n_mc_repetitions) - 1L,
    chunk_size = min(chunk_size, n_reference)
  )
}

score_reference_repetitions <- function(reference, score_fun) {
  lapply(seq_len(reference$n_mc_repetitions), function(r) {
    set.seed(reference$seeds[r])
    scores <- numeric(reference$n_reference)
    starts <- seq(1L, reference$n_reference, by = reference$chunk_size)
    for (start in starts) {
      end <- min(start + reference$chunk_size - 1L, reference$n_reference)
      size <- end - start + 1L
      points <- matrix(
        runif(size * reference$box$dimension),
        nrow = size,
        byrow = TRUE
      )
      points <- sweep(points, 2L, reference$box$width, "*")
      points <- sweep(points, 2L, reference$box$lower, "+")
      scores[start:end] <- validate_scores(
        score_fun(points),
        size,
        "reference_scores"
      )
    }
    scores
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
