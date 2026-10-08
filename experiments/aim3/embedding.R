source("experiments/aim3/mds_embedding.R")
source("experiments/aim3/isomap_embedding.R")

aim3_fit_embedding <- function(
    landmarks,
    ndim,
    method,
    embedding_settings
) {
  if (method == "mds") {
    model <- aim3_fit_mds(landmarks, ndim)
    project <- function(newdata) aim3_project_mds(model, newdata)
    target_distance <- function(newdata) {
      sqrt(aim3_squared_cross_distance(newdata, model$landmarks))
    }
  } else if (method == "isomap") {
    model <- aim3_fit_isomap(
      landmarks,
      ndim,
      embedding_settings$isomap_k
    )
    project <- function(newdata) aim3_project_isomap(model, newdata)
    target_distance <- function(newdata) {
      aim3_isomap_oos_distances(model, newdata)
    }
  } else {
    stop("Unknown embedding method: ", method, call. = FALSE)
  }

  points <- as.matrix(model$points)
  if (!identical(dim(points), c(nrow(landmarks), as.integer(ndim)))) {
    stop("Embedding returned an unexpected number of landmark coordinates.",
         call. = FALSE)
  }
  list(
    method = method,
    points = points,
    project = project,
    target_distance = target_distance
  )
}

aim3_project_embedding <- function(embedding, newdata) {
  points <- as.matrix(embedding$project(newdata))
  if (ncol(points) != ncol(embedding$points) || any(!is.finite(points))) {
    stop("Embedding projection returned invalid coordinates.", call. = FALSE)
  }
  points
}

aim3_embedding_oos_stress <- function(embedding, landmarks, newdata) {
  projected <- aim3_project_embedding(embedding, newdata)
  target_distance <- embedding$target_distance(newdata)
  embedded_distance <- sqrt(
    aim3_squared_cross_distance(projected, embedding$points)
  )
  denominator <- sum(embedded_distance^2)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  scale <- sum(target_distance * embedded_distance) / denominator
  residual <- target_distance - scale * embedded_distance
  sqrt(sum(residual^2) / sum(target_distance^2))
}
