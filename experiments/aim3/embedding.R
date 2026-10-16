require_manifun <- function() {
  if (!requireNamespace("manifun", quietly = TRUE)) {
    stop(
      "Package 'manifun' is required. Run experiments/install_manifun.R first.",
      call. = FALSE
    )
  }
}

source("experiments/aim3/mds_embedding.R")
source("experiments/aim3/isomap_embedding.R")

fit_embedding <- function(
    landmarks,
    ndim,
    method,
    embedding_settings
) {
  if (method == "mds") {
    model <- fit_mds(landmarks, ndim)
    project <- function(newdata) project_mds(model, newdata)
    target_distance <- function(newdata) {
      sqrt(squared_cross_distance(newdata, model$landmarks))
    }
  } else if (method == "isomap") {
    model <- fit_isomap(
      landmarks,
      ndim,
      embedding_settings$isomap_k
    )
    project <- function(newdata) project_isomap(model, newdata)
    target_distance <- function(newdata) {
      isomap_oos_distances(model, newdata)
    }
  } else {
    stop("Unknown embedding method: ", method, call. = FALSE)
  }

  points <- as.matrix(model$points)
  if (
    !identical(dim(points), c(nrow(landmarks), as.integer(ndim))) ||
    any(!is.finite(points))
  ) {
    stop(
      "Embedding returned invalid landmark coordinates.",
      call. = FALSE
    )
  }

  list(
    method = method,
    points = points,
    project = project,
    target_distance = target_distance
  )
}

project_embedding <- function(embedding, newdata) {
  points <- as.matrix(embedding$project(newdata))
  if (
    nrow(points) != nrow(newdata) ||
    ncol(points) != ncol(embedding$points) ||
    any(!is.finite(points))
  ) {
    stop("Embedding projection returned invalid coordinates.", call. = FALSE)
  }
  points
}

embedding_oos_stress <- function(embedding, newdata) {
  projected <- project_embedding(embedding, newdata)
  target_distance <- as.matrix(embedding$target_distance(newdata))
  embedded_distance <- sqrt(
    squared_cross_distance(projected, embedding$points)
  )

  if (
    !identical(dim(target_distance), dim(embedded_distance)) ||
    any(!is.finite(target_distance)) ||
    any(!is.finite(embedded_distance))
  ) {
    stop("Invalid embedding target distances.", call. = FALSE)
  }

  embedded_sum <- sum(embedded_distance^2)
  target_sum <- sum(target_distance^2)
  if (
    !is.finite(embedded_sum) || embedded_sum <= 0 ||
    !is.finite(target_sum) || target_sum <= 0
  ) {
    return(NA_real_)
  }

  scale <- sum(target_distance * embedded_distance) / embedded_sum
  residual <- target_distance - scale * embedded_distance
  sqrt(sum(residual^2) / target_sum)
}
