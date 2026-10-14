lof_matrix <- function(x) {
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

fit_lof <- function(x_train, k = 20L) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required", call. = FALSE)
  }

  x_train <- lof_matrix(x_train)
  k <- as.integer(k)

  if (length(k) != 1L || is.na(k) || k < 1L || k >= nrow(x_train)) {
    stop("k must be between 1 and nrow(x_train) - 1", call. = FALSE)
  }

  neighbors <- dbscan::kNN(
    x_train,
    k = k,
    sort = TRUE
  )
  neighbor_index <- as.matrix(neighbors$id)
  neighbor_distance <- as.matrix(neighbors$dist)
  k_distance <- neighbor_distance[, k]
  neighbor_k_distance <- matrix(
    k_distance[neighbor_index],
    nrow = nrow(x_train),
    ncol = k
  )
  reachability <- pmax(
    neighbor_distance,
    neighbor_k_distance
  )
  mean_reachability <- pmax(
    rowMeans(reachability),
    .Machine$double.eps
  )
  lrd <- 1 / mean_reachability

  list(
    x_train = x_train,
    k = k,
    k_distance = k_distance,
    lrd = lrd,
    dimension = ncol(x_train)
  )
}

score_lof <- function(model, newdata) {
  newdata <- lof_matrix(newdata)

  if (ncol(newdata) != model$dimension) {
    stop("newdata has the wrong dimension", call. = FALSE)
  }

  neighbors <- dbscan::kNN(
    model$x_train,
    k = model$k,
    query = newdata,
    sort = TRUE
  )
  neighbor_index <- as.matrix(neighbors$id)
  neighbor_distance <- as.matrix(neighbors$dist)
  neighbor_k_distance <- matrix(
    model$k_distance[neighbor_index],
    nrow = nrow(newdata),
    ncol = model$k
  )
  reachability <- pmax(
    neighbor_distance,
    neighbor_k_distance
  )
  mean_reachability <- pmax(
    rowMeans(reachability),
    .Machine$double.eps
  )
  query_lrd <- 1 / mean_reachability
  neighbor_lrd <- matrix(
    model$lrd[neighbor_index],
    nrow = nrow(newdata),
    ncol = model$k
  )

  rowMeans(neighbor_lrd) / query_lrd
}
