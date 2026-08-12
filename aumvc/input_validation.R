validate_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  } else if (is.numeric(x) && is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  } else {
    x <- as.matrix(x)
  }

  if (!is.numeric(x)) stop(name, " must be numeric", call. = FALSE)
  storage.mode(x) <- "double"

  if (
    nrow(x) < 1L ||
    ncol(x) < 1L ||
    anyNA(x) ||
    !all(is.finite(x))
  ) {
    stop(name, " must be a finite numeric matrix", call. = FALSE)
  }

  x
}

validate_scores <- function(scores, n_expected = NULL, name = "scores") {
  scores <- as.numeric(scores)

  if (length(scores) < 1L || anyNA(scores) || !all(is.finite(scores))) {
    stop(name, " must contain finite numeric values", call. = FALSE)
  }

  if (!is.null(n_expected) && length(scores) != n_expected) {
    stop(name, " has the wrong length", call. = FALSE)
  }

  scores
}

validate_alpha_grid <- function(alpha_grid) {
  alpha_grid <- as.numeric(alpha_grid)

  if (
    length(alpha_grid) < 2L ||
    any(alpha_grid <= 0 | alpha_grid >= 1) ||
    any(diff(alpha_grid) <= 0)
  ) {
    stop("Invalid alpha_grid", call. = FALSE)
  }

  alpha_grid
}

validate_tau_grid <- function(tau_grid) {
  tau_grid <- as.numeric(tau_grid)

  if (
    length(tau_grid) < 2L ||
    tau_grid[1L] != 0 ||
    any(tau_grid < 0) ||
    any(diff(tau_grid) <= 0)
  ) {
    stop("Invalid tau_grid", call. = FALSE)
  }

  tau_grid
}
