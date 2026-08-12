lof_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (
    nrow(x) < 1L ||
    ncol(x) < 1L ||
    anyNA(x) ||
    !all(is.finite(x))
  ) {
    stop(
      "Data must be a finite numeric matrix", call. = FALSE)
  }

  x
}

euclidean_distances <- function(x, y) {
  d2 <- outer(
    rowSums(x^2),
    rowSums(y^2),
    "+"
  ) - 2 * tcrossprod(x, y)

  sqrt(
    pmax(d2, 0)
  )
}

nearest_indices <- function(distances, k) {
  order(
    distances,
    decreasing = FALSE
  )[seq_len(k)]
}

fit_lof <- function(
    x_train,
    k = 20L
) {
  x_train <- lof_matrix(
    x_train
  )

  k <- as.integer(k)
  n <- nrow(x_train)

  if (
    k < 1L ||
    k >= n
  ) {
    stop(
      "k must be between 1 and nrow(x_train) - 1", call. = FALSE)
  }

  distances <- euclidean_distances(
    x_train,
    x_train
  )

  diag(distances) <- Inf

  neighbors <- lapply(
    seq_len(n),
    function(i) {
      nearest_indices(
        distances[i, ],
        k
      )
    }
  )

  k_distance <- vapply(
    seq_len(n),
    function(i) {
      max(
        distances[
          i,
          neighbors[[i]]
        ]
      )
    },
    numeric(1)
  )

  lrd <- vapply(
    seq_len(n),
    function(i) {
      index <- neighbors[[i]]

      reachability <- pmax(
        k_distance[index],
        distances[i, index]
      )

      1 / (
        mean(reachability) +
          1e-10
      )
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

score_lof <- function(
    model,
    newdata,
    chunk_size = 500L
) {
  newdata <- lof_matrix(
    newdata
  )

  if (
    ncol(newdata) != model$dimension
  ) {
    stop(
      "newdata has the wrong dimension", call. = FALSE)
  }

  scores <- numeric(
    nrow(newdata)
  )

  for (
    start in seq(
      1L,
      nrow(newdata),
      by = chunk_size
    )
  ) {
    end <- min(
      start + chunk_size - 1L,
      nrow(newdata)
    )

    query <- newdata[
      start:end,
      ,
      drop = FALSE
    ]

    distances <- euclidean_distances(
      query,
      model$x_train
    )

    scores[start:end] <- vapply(
      seq_len(nrow(query)),
      function(i) {
        index <- nearest_indices(
          distances[i, ],
          model$k
        )

        reachability <- pmax(
          model$k_distance[index],
          distances[i, index]
        )

        query_lrd <- 1 / (
          mean(reachability) +
            1e-10
        )

        mean(
          model$lrd[index]
        ) / query_lrd
      },
      numeric(1)
    )
  }

  scores
}
