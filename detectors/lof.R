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

kth_distance <- function(x, k) {
  sort(
    x,
    partial = k
  )[k]
}

k_distinct_distance <- function(
    distances,
    k,
    zero_tolerance = 1e-12
) {
  positive <- distances[
    distances > zero_tolerance
  ]

  if (length(positive) < k) {
    stop(
      "Not enough distinct spatial coordinates for this value of k", call. = FALSE)
  }

  kth_distance(
    positive,
    k
  )
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

  unique_x <- unique(
    x_train
  )

  if (nrow(unique_x) <= k) {
    stop(
      "LOF requires at least k + 1 distinct spatial coordinates", call. = FALSE)
  }

  distances <- euclidean_distances(
    x_train,
    x_train
  )

  diag(distances) <- Inf

  distinct_distances <- euclidean_distances(
    x_train,
    unique_x
  )

  k_distance <- vapply(
    seq_len(n),
    function(i) {
      k_distinct_distance(
        distinct_distances[i, ],
        k
      )
    },
    numeric(1)
  )

  neighborhoods <- lapply(
    seq_len(n),
    function(i) {
      tolerance <- 1e-12 * max(
        1,
        k_distance[i]
      )

      which(
        distances[i, ] <=
          k_distance[i] + tolerance
      )
    }
  )

  lrd <- vapply(
    seq_len(n),
    function(i) {
      neighbors <- neighborhoods[[i]]

      reachability <- pmax(
        k_distance[neighbors],
        distances[i, neighbors]
      )

      1 / mean(
        reachability
      )
    },
    numeric(1)
  )

  list(
    x_train = x_train,
    unique_x = unique_x,
    k = k,
    k_distance = k_distance,
    lrd = lrd,
    dimension = ncol(x_train)
  )
}

score_lof <- function(
    model,
    newdata,
    chunk_size = 1000L
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

    distinct_distances <- euclidean_distances(
      query,
      model$unique_x
    )

    scores[start:end] <- vapply(
      seq_len(nrow(query)),
      function(i) {
        d <- distances[i, ]

        query_k <- k_distinct_distance(
          distinct_distances[i, ],
          model$k
        )

        tolerance <- 1e-12 * max(
          1,
          query_k
        )

        neighbors <- which(
          d <= query_k + tolerance
        )

        reachability <- pmax(
          model$k_distance[neighbors],
          d[neighbors]
        )

        query_lrd <- 1 / mean(
          reachability
        )

        mean(
          model$lrd[neighbors]
        ) / query_lrd
      },
      numeric(1)
    )
  }

  scores
}
