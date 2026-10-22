lof_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (nrow(x) < 1L || ncol(x) < 1L || anyNA(x) || !all(is.finite(x))) {
    stop("Data must be a finite numeric matrix", call. = FALSE)
  }

  x
}

fit_lof <- function(x_train, k = 20L) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required", call. = FALSE)
  }
  if (utils::packageVersion("dbscan") < "1.2.0") {
    stop("Package 'dbscan' >= 1.2.0 is required", call. = FALSE)
  }

  x_train <- lof_matrix(x_train)
  if (length(k) != 1L || !is.finite(k) || k != as.integer(k)) {
    stop("k must be an integer", call. = FALSE)
  }
  k <- as.integer(k)
  if (k < 1L || k >= nrow(x_train)) {
    stop("k must be between 1 and nrow(x_train) - 1", call. = FALSE)
  }

  training_scores <- as.numeric(dbscan::lof(x_train, minPts = k + 1L))
  neighbors <- dbscan::kNN(x_train, k = k, sort = TRUE)
  neighbor_index <- as.matrix(neighbors$id)
  neighbor_distance <- as.matrix(neighbors$dist)
  k_distance <- neighbor_distance[, k]

  reachability <- pmax(
    neighbor_distance,
    matrix(k_distance[neighbor_index], nrow = nrow(x_train), ncol = k)
  )
  positive_distances <- k_distance[k_distance > 0 & is.finite(k_distance)]
  typical_distance <- if (length(positive_distances)) median(positive_distances) else 1
  distance_floor <- max(1e-12, 1e-10 * typical_distance)
  training_lrd <- 1 / pmax(rowMeans(reachability), distance_floor)

  list(
    x_train = x_train,
    k = k,
    minPts = k + 1L,
    dimension = ncol(x_train),
    training_scores = training_scores,
    k_distance = k_distance,
    training_lrd = training_lrd,
    distance_floor = distance_floor
  )
}

score_lof <- function(model, newdata) {
  newdata <- lof_matrix(newdata)
  if (ncol(newdata) != model$dimension) {
    stop("newdata has the wrong dimension", call. = FALSE)
  }

  if (identical(unname(newdata), unname(model$x_train))) {
    return(model$training_scores)
  }

  neighbors <- dbscan::kNN(model$x_train, k = model$k,
                           query = newdata, sort = TRUE)
  neighbor_index <- as.matrix(neighbors$id)
  neighbor_distance <- as.matrix(neighbors$dist)

  reachability <- pmax(
    neighbor_distance,
    matrix(model$k_distance[neighbor_index],
           nrow = nrow(newdata), ncol = model$k)
  )
  query_lrd <- 1 / pmax(rowMeans(reachability), model$distance_floor)
  neighbor_lrd <- matrix(
    model$training_lrd[neighbor_index],
    nrow = nrow(newdata), ncol = model$k
  )

  scores <- rowMeans(neighbor_lrd) / query_lrd
  dense_duplicates <- rowSums(neighbor_distance == 0) == model$k
  scores[dense_duplicates] <- 1
  if (!all(is.finite(scores))) {
    stop("LOF produced non-finite out-of-sample scores", call. = FALSE)
  }
  as.numeric(scores)
}
