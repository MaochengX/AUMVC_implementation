aim3_require_manifun <- function() {
  if (!requireNamespace("manifun", quietly = TRUE)) {
    stop(
      "Package 'manifun' is required. Run experiments/install_manifun.R first.",
      call. = FALSE
    )
  }
}

aim3_fit_mds <- function(landmarks, ndim) {
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

  distance <- stats::dist(landmarks)
  distance_squared <- as.matrix(distance)^2

  embedding <- embed_fun(
    distance,
    method = "mds",
    k = ndim
  )

  points <- as.matrix(
    manifun::extract_points(embedding, ndim = ndim)
  )

  if (!identical(dim(points), c(nrow(landmarks), ndim))) {
    stop("manifun returned an unexpected MDS size.", call. = FALSE)
  }

  gram <- crossprod(points)

  if (qr(gram)$rank < ndim) {
    stop("MDS coordinates are rank deficient.", call. = FALSE)
  }

  list(
    landmarks = landmarks,
    points = points,
    landmark_norm = rowMeans(distance_squared) - mean(distance_squared) / 2,
    gram_inverse = solve(gram)
  )
}

aim3_squared_cross_distance <- function(x, y) {
  squared <- outer(rowSums(x^2), rowSums(y^2), "+") -
    2 * tcrossprod(x, y)
  pmax(squared, 0)
}

aim3_project_mds <- function(model, newdata) {
  newdata <- as.matrix(newdata)

  if (ncol(newdata) != ncol(model$landmarks)) {
    stop("newdata has the wrong ambient dimension.", call. = FALSE)
  }

  distance_squared <- aim3_squared_cross_distance(
    newdata,
    model$landmarks
  )

  distance_centered <- sweep(
    distance_squared,
    1L,
    rowMeans(distance_squared),
    "-"
  )

  inner_products <- -0.5 * sweep(
    distance_centered,
    2L,
    model$landmark_norm - mean(model$landmark_norm),
    "-"
  )

  inner_products %*% model$points %*% model$gram_inverse
}
