ocsvm_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (
    nrow(x) < 1 ||
    ncol(x) < 1 ||
    anyNA(x) ||
    !all(is.finite(x))
  ) {
    stop("Data must be a finite numeric matrix", call. = FALSE)
  }

  x
}

gaussian_kernel <- function(x, y, gamma) {
  d2 <- outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y)
  exp(-gamma * pmax(d2, 0))
}

solve_ocsvm_dual <- function(K, nu, tolerance = 1e-7, max_iter = 100000L) {
  n <- nrow(K)
  cap <- 1 / (nu * n)
  alpha <- rep(1 / n, n)
  gradient <- as.numeric(K %*% alpha)
  gap <- Inf
  converged <- FALSE

  if (cap <= 1 / n + 1e-14) {
    return(
      list(
        alpha = alpha,
        gradient = gradient,
        cap = cap,
        gap = 0,
        converged = TRUE
      )
    )
  }

  for (iteration in seq_len(max_iter)) {
    increase <- which(alpha < cap - 1e-12)
    decrease <- which(alpha > 1e-12)

    i <- increase[which.min(gradient[increase])]
    j <- decrease[which.max(gradient[decrease])]
    gap <- gradient[j] - gradient[i]

    if (gap <= tolerance) {
      converged <- TRUE
      break
    }

    curvature <- K[i, i] + K[j, j] - 2 * K[i, j]
    max_step <- min(cap - alpha[i], alpha[j])

    if (curvature > 1e-14) {
      step <- min(max_step, gap / curvature)
    } else {
      step <- max_step
    }

    if (step <= 1e-15) {
      break
    }

    alpha[i] <- alpha[i] + step
    alpha[j] <- alpha[j] - step
    gradient <- gradient + step * (K[, i] - K[, j])
  }

  list(
    alpha = alpha,
    gradient = gradient,
    cap = cap,
    gap = gap,
    converged = converged
  )
}

fit_ocsvm <- function(
    x_train,
    nu = 0.05,
    gamma = 1 / ncol(x_train),
    tolerance = 1e-7,
    max_iter = 100000L
) {
  x_train <- ocsvm_matrix(x_train)

  if (
    length(nu) != 1 ||
    !is.finite(nu) ||
    nu <= 0 ||
    nu > 1
  ) {
    stop("nu must lie in (0, 1]", call. = FALSE)
  }

  if (
    length(gamma) != 1 ||
    !is.finite(gamma) ||
    gamma <= 0
  ) {
    stop("gamma must be positive", call. = FALSE)
  }

  K <- gaussian_kernel(x_train, x_train, gamma)
  solution <- solve_ocsvm_dual(K, nu, tolerance, max_iter)

  alpha <- solution$alpha
  cap <- solution$cap
  f_train <- solution$gradient
  eps <- 1e-8

  free <- which(alpha > eps & alpha < cap - eps)

  if (length(free) > 0) {
    rho <- mean(f_train[free])
  } else {
    upper <- which(alpha >= cap - eps)
    zero <- which(alpha <= eps)

    if (length(upper) > 0 && length(zero) > 0) {
      rho <- (max(f_train[upper]) + min(f_train[zero])) / 2
    } else {
      rho <- median(f_train[alpha > eps])
    }
  }

  support <- which(alpha > eps)

  list(
    support_vectors = x_train[support, , drop = FALSE],
    support_alpha = alpha[support],
    rho = rho,
    nu = nu,
    gamma = gamma,
    dimension = ncol(x_train),
    kkt_gap = solution$gap,
    converged = solution$converged
  )
}

score_ocsvm <- function(model, newdata, chunk_size = 2000L) {
  newdata <- ocsvm_matrix(newdata)

  if (ncol(newdata) != model$dimension) {
    stop("newdata has the wrong dimension", call. = FALSE)
  }

  scores <- numeric(nrow(newdata))

  for (start in seq(1L, nrow(newdata), by = chunk_size)) {
    end <- min(start + chunk_size - 1L, nrow(newdata))

    K <- gaussian_kernel(
      newdata[start:end, , drop = FALSE],
      model$support_vectors,
      model$gamma
    )

    scores[start:end] <- model$rho - as.numeric(K %*% model$support_alpha)
  }

  scores
}
