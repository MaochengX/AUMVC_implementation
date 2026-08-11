validate_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  } else if (is.numeric(x) && is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  } else {
    x <- as.matrix(x)
  }

  if (!is.numeric(x)) {
    stop(name, "must be numeric", call. = FALSE)
  }

  storage.mode(x) <- "double"

  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(name, "cannot be empty", call. = FALSE)
  }

  if (anyNA(x) || !all(is.finite(x))) {
    stop(name, "cannot contain NA or infinite values", call. = FALSE)
  }

  x
}

validate_scores <- function(scores, n_expected = NULL, name = "scores") {
  scores <- as.numeric(scores)

  if (length(scores) < 1L) {
    stop(name, "cannot be empty", call. = FALSE)
  }

  if (anyNA(scores) || !all(is.finite(scores))) {
    stop(name, "cannot contain NA or infinite values", call. = FALSE)
  }

  if (!is.null(n_expected) && length(scores) != n_expected) {
    stop(
      name,
      "must contain",
      n_expected,
      "values, one per observation",
      call. = FALSE
    )
  }

  scores
}

validate_reference_inputs <- function(
    x_train,
    n_reference,
    n_mc_repetitions,
    seed
) {
  x_train <- validate_matrix(x_train, "x_train")

  widths <- apply(x_train, 2L, max) - apply(x_train, 2L, min)

  if (any(widths <= 0)) {
    stop("x_train contains a constant coordinate", call. = FALSE)
  }

  if (
    length(n_reference) != 1L ||
    !is.finite(n_reference) ||
    n_reference < 1 ||
    n_reference != as.integer(n_reference)
  ) {
    stop("n_reference must be a positive integer", call. = FALSE)
  }

  if (
    length(n_mc_repetitions) != 1L ||
    !is.finite(n_mc_repetitions) ||
    n_mc_repetitions < 1 ||
    n_mc_repetitions != as.integer(n_mc_repetitions)
  ) {
    stop("n_mc_repetitions must be a positive integer", call. = FALSE)
  }

  if (length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be one finite number", call. = FALSE)
  }

  x_train
}

validate_aumvc_inputs <- function(
    x_eval,
    reference,
    score_fun,
    score_direction,
    alpha_grid
) {
  x_eval <- validate_matrix(x_eval, "x_eval")

  if (!is.list(reference)) {
    stop("reference must be created by make_reference()", call. = FALSE)
  }

  if (ncol(x_eval) != reference$dimension) {
    stop("x_eval and reference must have the same dimension", call. = FALSE)
  }

  if (!is.function(score_fun)) {
    stop("score_fun must be a function", call. = FALSE)
  }

  if (!(score_direction %in% c("normality", "anomaly"))) {
    stop("score_direction must be 'normality' or 'anomaly'", call. = FALSE)
  }

  alpha_grid <- as.numeric(alpha_grid)

  if (
    anyNA(alpha_grid) ||
    !all(is.finite(alpha_grid)) ||
    any(alpha_grid <= 0 | alpha_grid > 1) ||
    any(diff(alpha_grid) <= 0)
  ) {
    stop("alpha_grid must be strictly increasing in (0, 1]", call. = FALSE)
  }

  x_eval
}

validate_auemc_inputs <- function(
    x_eval,
    reference,
    score_fun,
    score_direction,
    penalty_grid
) {
  x_eval <- validate_matrix(x_eval, "x_eval")

  if (!is.list(reference)) {
    stop("reference must be created by make_reference()", call. = FALSE)
  }

  if (ncol(x_eval) != reference$dimension) {
    stop("x_eval and reference must have the same dimension", call. = FALSE)
  }

  if (!is.function(score_fun)) {
    stop("score_fun must be a function", call. = FALSE)
  }

  if (!(score_direction %in% c("normality", "anomaly"))) {
    stop("score_direction must be 'normality' or 'anomaly'", call. = FALSE)
  }

  penalty_grid <- as.numeric(penalty_grid)

  if (
    anyNA(penalty_grid) ||
    !all(is.finite(penalty_grid)) ||
    penalty_grid[1L] != 0 ||
    any(penalty_grid < 0) ||
    any(diff(penalty_grid) <= 0)
  ) {
    stop(
      "penalty_grid must start at 0 and be strictly increasing",
      call. = FALSE
    )
  }

  x_eval
}
