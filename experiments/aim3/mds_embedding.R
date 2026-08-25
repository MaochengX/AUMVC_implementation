aim3_require_manifun <- function() {
  if (!requireNamespace("manifun", quietly = TRUE)) {
    stop(
      paste(
        "Package 'manifun' is required.",
        "Run experiments/aim3/install_manifun.R first."
      ),
      call. = FALSE
    )
  }
}

aim3_cross_distance <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)

  squared <- outer(rowSums(x^2), rowSums(y^2), "+") -
    2 * tcrossprod(x, y)

  sqrt(pmax(squared, 0))
}

aim3_fit_mds <- function(landmarks, ndim, settings) {
  aim3_require_manifun()

  landmarks <- as.matrix(landmarks)
  ndim <- as.integer(ndim)

  if (ndim < 1L || ndim >= nrow(landmarks)) {
    stop("Invalid MDS dimension.", call. = FALSE)
  }

  embed_fun <- manifun::embed
  helper_env <- new.env(parent = environment(embed_fun))
  helper_env$cmdscale <- stats::cmdscale
  environment(embed_fun) <- helper_env

  embedding <- embed_fun(
    stats::dist(landmarks),
    method = "mds",
    k = ndim
  )

  points <- as.matrix(
    manifun::extract_points(embedding, ndim = ndim)
  )

  if (nrow(points) != nrow(landmarks) || ncol(points) != ndim) {
    stop("manifun returned an unexpected MDS size.", call. = FALSE)
  }

  list(
    ndim = ndim,
    landmarks = landmarks,
    points = points,
    interpolation_neighbors = min(
      as.integer(settings$interpolation_neighbors),
      nrow(landmarks)
    ),
    interpolation_ridge = settings$interpolation_ridge
  )
}

aim3_landmark_weights <- function(x, landmarks, ridge) {
  local <- sweep(landmarks, 2L, x, "-")
  covariance <- tcrossprod(local)

  scale <- sum(diag(covariance))
  regularizer <- ridge * if (scale > 0) scale else 1
  covariance <- covariance + diag(regularizer, nrow(covariance))

  weights <- tryCatch(
    solve(covariance, rep(1, nrow(covariance))),
    error = function(e) rep(1, nrow(covariance))
  )

  total <- sum(weights)

  if (!is.finite(total) || abs(total) < 1e-12) {
    return(rep(1 / length(weights), length(weights)))
  }

  weights / total
}

aim3_project_mds <- function(model, newdata) {
  newdata <- as.matrix(newdata)

  if (ncol(newdata) != ncol(model$landmarks)) {
    stop("newdata has the wrong ambient dimension.", call. = FALSE)
  }

  distances <- aim3_cross_distance(newdata, model$landmarks)
  output <- matrix(0, nrow(newdata), model$ndim)

  for (i in seq_len(nrow(newdata))) {
    nearest <- order(distances[i, ])[
      seq_len(model$interpolation_neighbors)
    ]

    exact <- nearest[distances[i, nearest] <= 1e-12]

    if (length(exact) > 0L) {
      output[i, ] <- model$points[exact[1L], ]
      next
    }

    weights <- aim3_landmark_weights(
      newdata[i, ],
      model$landmarks[nearest, , drop = FALSE],
      model$interpolation_ridge
    )

    output[i, ] <- colSums(
      model$points[nearest, , drop = FALSE] * weights
    )
  }

  output
}
