fit_isomap <- function(landmarks, ndim, k = 10L) {
  require_manifun()
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for Isomap.", call. = FALSE)
  }

  landmarks <- as.matrix(landmarks)
  storage.mode(landmarks) <- "double"
  ndim <- as.integer(ndim)
  k <- as.integer(k)
  n <- nrow(landmarks)

  if (
    n < 3L ||
    ncol(landmarks) < 1L ||
    any(!is.finite(landmarks)) ||
    length(ndim) != 1L ||
    is.na(ndim) ||
    ndim < 1L ||
    ndim >= n ||
    length(k) != 1L ||
    is.na(k) ||
    k < 1L ||
    k >= n
  ) {
    stop("Invalid Isomap settings or landmark data.", call. = FALSE)
  }

  distance <- stats::dist(landmarks)
  selected_k <- k

  repeat {
    geodesic <- tryCatch(
      vegan::isomapdist(
        distance,
        k = selected_k,
        fragmentedOK = FALSE
      ),
      error = function(error) NULL
    )

    if (!is.null(geodesic)) {
      embedding <- tryCatch(
        suppressWarnings(
          manifun::embed(
            distance,
            method = "isomap",
            ndim = ndim,
            k = selected_k,
            fragmentedOK = FALSE
          )
        ),
        error = function(error) NULL
      )

      points <- if (is.null(embedding)) {
        NULL
      } else {
        tryCatch(
          as.matrix(
            manifun::extract_points(
              embedding,
              ndim = ndim
            )
          ),
          error = function(error) NULL
        )
      }

      if (
        !is.null(points) &&
        identical(dim(points), c(n, ndim)) &&
        all(is.finite(points)) &&
        qr(crossprod(points))$rank == ndim
      ) {
        break
      }
    }

    selected_k <- selected_k + 1L
    if (selected_k >= n) {
      stop(
        "Isomap could not produce a connected, full-rank embedding.",
        call. = FALSE
      )
    }
  }

  geodesic <- as.matrix(geodesic)
  distance_squared <- geodesic^2
  point_gram <- crossprod(points)

  model <- list(
    landmarks = landmarks,
    points = points,
    geodesic = geodesic,
    k = selected_k,
    landmark_norm = rowMeans(distance_squared) -
      mean(distance_squared) / 2,
    gram_inverse = solve(point_gram)
  )

  reconstructed <- project_isomap(model, landmarks)
  error <- max(abs(reconstructed - points))
  scale <- max(1, max(abs(points)))

  if (!is.finite(error) || error > 1e-6 * scale) {
    stop(
      "Isomap out-of-sample projection failed its landmark check.",
      call. = FALSE
    )
  }

  model
}

isomap_oos_distances <- function(model, newdata) {
  newdata <- as.matrix(newdata)

  if (
    ncol(newdata) != ncol(model$landmarks) ||
    any(!is.finite(newdata))
  ) {
    stop("Invalid Isomap projection data.", call. = FALSE)
  }

  direct <- sqrt(
    squared_cross_distance(
      newdata,
      model$landmarks
    )
  )

  geodesic <- matrix(
    0,
    nrow = nrow(newdata),
    ncol = nrow(model$landmarks)
  )

  for (index in seq_len(nrow(newdata))) {
    neighbors <- order(direct[index, ])[seq_len(model$k)]
    candidates <- sweep(
      model$geodesic[neighbors, , drop = FALSE],
      1L,
      direct[index, neighbors],
      "+"
    )
    geodesic[index, ] <- apply(candidates, 2L, min)
  }

  geodesic
}

project_isomap <- function(model, newdata) {
  distance_squared <- isomap_oos_distances(model, newdata)^2

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
