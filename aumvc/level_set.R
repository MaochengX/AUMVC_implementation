orient_scores <- function(scores, direction) {
  if (direction == "anomaly") {
    -as.numeric(scores)
  } else {
    as.numeric(scores)
  }
}

fit_reference_box <- function(x_train) {
  lower <- apply(x_train, 2L, min)
  upper <- apply(x_train, 2L, max)
  width <- upper - lower

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

  columns <- lapply(
    seq_len(box$dimension),
    function(j) {
      runif(
        n_reference,
        min = box$lower[j],
        max = box$upper[j]
      )
    }
  )

  matrix(
    unlist(columns, use.names = FALSE),
    nrow = n_reference,
    ncol = box$dimension
  )
}

make_reference <- function(
    x_train,
    n_reference = 100000L,
    n_mc_repetitions = 5L,
    seed = 2026L
) {
  box <- fit_reference_box(x_train)

  list(
    box = box,
    n_reference = as.integer(n_reference),
    n_mc_repetitions = as.integer(n_mc_repetitions),
    seeds = as.integer(seed) + seq_len(n_mc_repetitions) - 1L,
    dimension = ncol(x_train)
  )
}

level_set_table <- function(
    evaluation_scores,
    reference_scores,
    box_volume
) {
  thresholds <- sort(
    unique(evaluation_scores),
    decreasing = TRUE
  )

  mass <- vapply(
    thresholds,
    function(u) mean(evaluation_scores >= u),
    numeric(1L)
  )

  occupancy <- vapply(
    thresholds,
    function(u) mean(reference_scores >= u),
    numeric(1L)
  )

  volume <- box_volume * occupancy

  volume_se <- box_volume * sqrt(
    occupancy * (1 - occupancy) / length(reference_scores)
  )

  data.frame(
    threshold = thresholds,
    mass = mass,
    occupancy = occupancy,
    volume = volume,
    volume_se = volume_se
  )
}
