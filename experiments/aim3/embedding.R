# An embedding must fit once and project new rows without refitting.
source("experiments/aim3/mds_embedding.R")

aim3_fit_embedding <- function(landmarks, ndim, method = "mds") {
  if (method == "mds") {
    model <- aim3_fit_mds(landmarks, ndim)
    project <- function(newdata) aim3_project_mds(model, newdata)
  } else if (method == "isomap") {
    path <- "experiments/aim3/isomap_embedding.R"
    if (!file.exists(path)) {
      stop("Implement ", path, " with aim3_fit_isomap() and ",
           "aim3_project_isomap() before selecting isomap.", call. = FALSE)
    }
    source(path, local = TRUE)
    model <- aim3_fit_isomap(landmarks, ndim)
    project <- function(newdata) aim3_project_isomap(model, newdata)
  } else {
    stop("Unknown embedding method: ", method, call. = FALSE)
  }

  points <- as.matrix(model$points)
  if (!identical(dim(points), c(nrow(landmarks), as.integer(ndim)))) {
    stop("Embedding returned an unexpected number of landmark coordinates.",
         call. = FALSE)
  }
  list(method = method, points = points, project = project)
}

aim3_project_embedding <- function(embedding, newdata) {
  embedding$project(newdata)
}
