orient_scores <- function(scores, direction) {
  if (direction == "normality") return(as.numeric(scores))
  if (direction == "anomaly") return(-as.numeric(scores))

  stop(
    "score_direction must be 'normality' or 'anomaly'",
    call. = FALSE
  )
}

trapezoid_area <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

fit_reference_box <- function(x_train) {
  x_train <- validate_matrix(x_train, "x_train")

  lower <- apply(x_train, 2L, min)
  upper <- apply(x_train, 2L, max)
  width <- upper - lower

  if (any(width <= 0)) {
    stop("The reference split has a constant coordinate", call. = FALSE)
  }

  list(
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    width = as.numeric(width),
    volume = prod(width),
    dimension = ncol(x_train)
  )
}

sample_reference_points <- function(box, n_reference, seed) {
  set.seed(seed)

  x <- matrix(
    runif(n_reference * box$dimension),
    nrow = n_reference
  )

  x <- sweep(x, 2L, box$width, "*")
  sweep(x, 2L, box$lower, "+")
}

make_reference <- function(
    x_train,
    n_reference = 20000L,
    n_mc_repetitions = 5L,
    seed = 2030L
) {
  x_train <- validate_matrix(x_train, "x_train")

  list(
    box = fit_reference_box(x_train),
    n_reference = as.integer(n_reference),
    n_mc_repetitions = as.integer(n_mc_repetitions),
    seeds = seed + seq_len(n_mc_repetitions) - 1L,
    dimension = ncol(x_train)
  )
}

score_reference_repetitions <- function(reference, score_fun) {
  lapply(
    seq_len(reference$n_mc_repetitions),
    function(r) {
      points <- sample_reference_points(
        reference$box,
        reference$n_reference,
        reference$seeds[r]
      )

      validate_scores(
        score_fun(points),
        reference$n_reference,
        "reference_scores"
      )
    }
  )
}

level_set_table <- function(
    evaluation_scores,
    reference_scores,
    box_volume
) {
  thresholds <- sort(
    unique(c(evaluation_scores, reference_scores)),
    decreasing = TRUE
  )

  eval_bin <- tabulate(
    match(evaluation_scores, thresholds),
    nbins = length(thresholds)
  )

  ref_bin <- tabulate(
    match(reference_scores, thresholds),
    nbins = length(thresholds)
  )

  mass <- cumsum(eval_bin) / length(evaluation_scores)
  occupancy <- cumsum(ref_bin) / length(reference_scores)

  data.frame(
    threshold = thresholds,
    mass = mass,
    occupancy = occupancy,
    volume = box_volume * occupancy
  )
}
