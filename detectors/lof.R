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

euclidean_distances <- function(x, y) {
  d2 <- outer(rowSums(x^2), rowSums(y^2), "+") -
    2 * tcrossprod(x, y)

  sqrt(pmax(d2, 0))
}

k_neighbors <- function(distances, k) {
  kth <- sort(distances, partial = k)[k]
  which(distances <= kth)
}

fit_lof <- function(x_train, k = 20L, chunk_size = 250L) {
  x_train <- lof_matrix(x_train)

  k <- as.integer(k)
  chunk_size <- as.integer(chunk_size)
  n <- nrow(x_train)

  if (
    length(k) != 1L || is.na(k) || k < 1L || k >= n ||
    length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L
  ) {
    stop("k must be between 1 and nrow(x_train) - 1", call. = FALSE)
  }

  neighbors <- vector("list", n)
  k_distance <- numeric(n)
  for (start in seq(1L, n, by = chunk_size)) {
    end <- min(start + chunk_size - 1L, n)
    rows <- start:end
    distances <- euclidean_distances(
      x_train[rows, , drop = FALSE],
      x_train
    )
    distances[cbind(seq_along(rows), rows)] <- Inf
    for (position in seq_along(rows)) {
      index <- k_neighbors(distances[position, ], k)
      neighbors[[rows[position]]] <- index
      k_distance[rows[position]] <- max(distances[position, index])
    }
  }

  lrd <- vapply(
    seq_len(n),
    function(i) {
      index <- neighbors[[i]]
      difference <- sweep(
        x_train[index, , drop = FALSE],
        2L,
        x_train[i, ],
        "-"
      )
      distance <- sqrt(rowSums(difference^2))

      reachability <- pmax(
        k_distance[index],
        distance
      )

      1 / (mean(reachability) + 1e-10)
    },
    numeric(1)
  )

  list(
    x_train = x_train,
    k = k,
    k_distance = k_distance,
    lrd = lrd,
    dimension = ncol(x_train)
  )
}

score_lof <- function(model, newdata, chunk_size = 500L) {
  newdata <- lof_matrix(newdata)
  chunk_size <- as.integer(chunk_size)

  if (ncol(newdata) != model$dimension) {
    stop("newdata has the wrong dimension", call. = FALSE)
  }
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be a positive integer", call. = FALSE)
  }

  scores <- numeric(nrow(newdata))

  for (start in seq(1L, nrow(newdata), by = chunk_size)) {
    end <- min(start + chunk_size - 1L, nrow(newdata))

    distances <- euclidean_distances(
      newdata[start:end, , drop = FALSE],
      model$x_train
    )

    scores[start:end] <- vapply(
      seq_len(nrow(distances)),
      function(i) {
        index <- k_neighbors(distances[i, ], model$k)

        reachability <- pmax(
          model$k_distance[index],
          distances[i, index]
        )

        query_lrd <- 1 / (mean(reachability) + 1e-10)
        mean(model$lrd[index]) / query_lrd
      },
      numeric(1)
    )
  }

  scores
}
