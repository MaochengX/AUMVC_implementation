fit_lof <- function(
    x_train,
    k = 20L
) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required.", call. = FALSE)
  }

  x_train <- as.matrix(x_train)
  storage.mode(x_train) <- "double"
  k <- as.integer(k)

  if (k < 1L || k >= nrow(x_train)) {
    stop("k must be between 1 and nrow(x_train) - 1.", call. = FALSE)
  }

  neighbors <- FNN::get.knn(
    x_train,
    k = k
  )

  neighbor_index <- neighbors$nn.index
  neighbor_distance <- neighbors$nn.dist

  k_distance <- neighbor_distance[, k]

  reachability <- matrix(
    0,
    nrow = nrow(x_train),
    ncol = k
  )

  for (j in seq_len(k)) {
    neighbor_j <- neighbor_index[, j]

    reachability[, j] <- pmax(
      neighbor_distance[, j],
      k_distance[neighbor_j]
    )
  }

  lrd <- 1 / rowMeans(reachability)

  list(
    x_train = x_train,
    k = k,
    k_distance = k_distance,
    lrd = lrd,
    dimension = ncol(x_train)
  )
}

score_lof <- function(model, newdata) {
  newdata <- as.matrix(newdata)
  storage.mode(newdata) <- "double"

  neighbors <- FNN::get.knnx(
    data = model$x_train,
    query = newdata,
    k = model$k
  )

  neighbor_index <- neighbors$nn.index
  neighbor_distance <- neighbors$nn.dist

  reachability <- matrix(
    0,
    nrow = nrow(newdata),
    ncol = model$k
  )

  neighbor_lrd <- matrix(
    0,
    nrow = nrow(newdata),
    ncol = model$k
  )

  for (j in seq_len(model$k)) {
    neighbor_j <- neighbor_index[, j]

    reachability[, j] <- pmax(
      neighbor_distance[, j],
      model$k_distance[neighbor_j]
    )

    neighbor_lrd[, j] <- model$lrd[neighbor_j]
  }

  query_lrd <- 1 / rowMeans(reachability)

  rowMeans(neighbor_lrd) / query_lrd
}
