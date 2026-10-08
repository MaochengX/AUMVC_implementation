ocsvm_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (
    nrow(x) < 1L ||
    ncol(x) < 1L ||
    anyNA(x) ||
    !all(is.finite(x))
  ) {
    stop("Data must be a finite numeric matrix", call. = FALSE)
  }

  x
}

gaussian_kernel <- function(x, y, gamma, chunk_size = 500L) {
  chunk_size <- as.integer(chunk_size)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be a positive integer", call. = FALSE)
  }
  result <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  y_norm <- rowSums(y^2)
  for (start in seq(1L, nrow(x), by = chunk_size)) {
    end <- min(start + chunk_size - 1L, nrow(x))
    block <- x[start:end, , drop = FALSE]
    d2 <- outer(rowSums(block^2), y_norm, "+") -
      2 * tcrossprod(block, y)
    result[start:end, ] <- exp(-gamma * pmax(d2, 0))
  }
  result
}

solve_ocsvm_dual <- function(K, nu, tolerance = 1e-6, max_iter = 100000L) {
  n <- nrow(K)
  cap <- 1 / (nu * n)
  alpha <- rep(1 / n, n)
  gradient <- as.numeric(K %*% alpha)
  gap <- Inf
  converged <- FALSE
  iteration <- 0L

  for (iteration in seq_len(max_iter)) {
    increase <- which(alpha < cap - 1e-12)
    decrease <- which(alpha > 1e-12)

    if (length(increase) == 0L || length(decrease) == 0L) {
      gap <- 0
      converged <- TRUE
      break
    }

    i <- increase[which.min(gradient[increase])]
    j <- decrease[which.max(gradient[decrease])]

    gap <- gradient[j] - gradient[i]

    if (gap <= tolerance) {
      converged <- TRUE
      break
    }

    curvature <- K[i, i] + K[j, j] - 2 * K[i, j]
    max_step <- min(cap - alpha[i], alpha[j])

    step <- if (curvature > 1e-14) {
      min(max_step, gap / curvature)
    } else {
      max_step
    }

    if (step <= 1e-15) break

    alpha[i] <- alpha[i] + step
    alpha[j] <- alpha[j] - step

    gradient <- gradient +
      step * (K[, i] - K[, j])
  }

  list(
    alpha = alpha,
    gradient = gradient,
    cap = cap,
    gap = gap,
    converged = converged,
    iterations = iteration
  )
}

fit_ocsvm <- function(
    x_train,
    nu = 0.5,
    gamma = 1 / ncol(x_train),
    tolerance = 1e-6,
    max_iter = 100000L
) {
  x_train <- ocsvm_matrix(x_train)

  if (
    length(nu) != 1L || !is.finite(nu) || nu <= 0 || nu > 1 ||
    length(gamma) != 1L || !is.finite(gamma) || gamma <= 0 ||
    length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
    length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1L
  ) {
    stop("Invalid OCSVM parameters", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)

  K <- gaussian_kernel(x_train, x_train, gamma)
  solution <- solve_ocsvm_dual(K, nu, tolerance, max_iter)

  alpha <- solution$alpha
  cap <- solution$cap
  f_train <- solution$gradient

  eps <- 1e-8
  free <- which(alpha > eps & alpha < cap - eps)

  if (length(free) > 0L) {
    rho <- mean(f_train[free])
  } else {
    upper <- which(alpha >= cap - eps)
    zero <- which(alpha <= eps)

    rho <- if (length(upper) > 0L && length(zero) > 0L) {
      (max(f_train[upper]) + min(f_train[zero])) / 2
    } else {
      median(f_train[alpha > eps])
    }
  }

  support <- which(alpha > eps)
  sum_error <- abs(sum(alpha) - 1)
  bound_violation <- max(c(0, -alpha, alpha - cap))

  if (length(support) == 0L || !is.finite(rho)) {
    stop("OCSVM fitting produced an invalid model", call. = FALSE)
  }

  list(
    support_vectors = x_train[support, , drop = FALSE],
    support_alpha = alpha[support],
    rho = rho,
    gamma = gamma,
    dimension = ncol(x_train),
    kkt_gap = solution$gap,
    sum_alpha_error = sum_error,
    bound_violation = bound_violation,
    iterations = solution$iterations,
    converged = solution$converged &&
      sum_error <= max(1e-8, 10 * tolerance) &&
      bound_violation <= max(1e-10, tolerance)
  )
}

score_ocsvm <- function(model, newdata, chunk_size = 2000L) {
  newdata <- ocsvm_matrix(newdata)
  chunk_size <- as.integer(chunk_size)

  if (ncol(newdata) != model$dimension) {
    stop("newdata has the wrong dimension", call. = FALSE)
  }
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be a positive integer", call. = FALSE)
  }

  scores <- numeric(nrow(newdata))

  for (start in seq(1L, nrow(newdata), by = chunk_size)) {
    end <- min(start + chunk_size - 1L, nrow(newdata))

    K <- gaussian_kernel(
      newdata[start:end, , drop = FALSE],
      model$support_vectors,
      model$gamma
    )

    scores[start:end] <- model$rho -
      as.numeric(K %*% model$support_alpha)
  }

  scores
}
