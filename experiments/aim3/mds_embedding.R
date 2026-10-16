fit_mds <- function(landmarks, ndim) {
  require_manifun()

  landmarks <- as.matrix(landmarks)
  storage.mode(landmarks) <- "double"
  ndim <- as.integer(ndim)

  if (
    nrow(landmarks) < 2L ||
    ncol(landmarks) < 1L ||
    any(!is.finite(landmarks)) ||
    length(ndim) != 1L ||
    is.na(ndim) ||
    ndim < 1L ||
    ndim >= nrow(landmarks)
  ) {
    stop("Invalid MDS settings or landmark data.", call. = FALSE)
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

  if (
    !identical(dim(points), c(nrow(landmarks), ndim)) ||
    any(!is.finite(points))
  ) {
    stop("manifun returned invalid MDS coordinates.", call. = FALSE)
  }

  gram <- crossprod(points)
  if (qr(gram)$rank < ndim) {
    stop("MDS coordinates are rank deficient.", call. = FALSE)
  }

  model <- list(
    landmarks = landmarks,
    points = points,
    landmark_norm = rowMeans(distance_squared) - mean(distance_squared) / 2,
    gram_inverse = solve(gram)
  )

  reconstructed <- project_mds(model, landmarks)
  error <- max(abs(reconstructed - points))
  scale <- max(1, max(abs(points)))
  if (!is.finite(error) || error > 1e-6 * scale) {
    stop("MDS out-of-sample projection failed its landmark check.", call. = FALSE)
  }

  model
}

squared_cross_distance <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (
    ncol(x) != ncol(y) ||
    any(!is.finite(x)) ||
    any(!is.finite(y))
  ) {
    stop("Invalid cross-distance inputs.", call. = FALSE)
  }

  squared <- outer(rowSums(x^2), rowSums(y^2), "+") -
    2 * tcrossprod(x, y)
  pmax(squared, 0)
}

project_mds <- function(model, newdata) {
  newdata <- as.matrix(newdata)

  if (
    ncol(newdata) != ncol(model$landmarks) ||
    any(!is.finite(newdata))
  ) {
    stop("Invalid MDS projection data.", call. = FALSE)
  }

  distance_squared <- squared_cross_distance(
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
