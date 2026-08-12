iforest_matrix <- function(x) {
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

iforest_c <- function(n) {
  if (n <= 1L) return(0)

  2 * sum(1 / seq_len(n - 1L)) -
    2 * (n - 1) / n
}

build_itree <- function(x, height, height_limit) {
  n <- nrow(x)

  if (height >= height_limit || n <= 1L) {
    return(list(external = TRUE, size = n))
  }

  minimum <- apply(x, 2L, min)
  maximum <- apply(x, 2L, max)
  varying <- which(maximum > minimum)

  if (length(varying) == 0L) {
    return(list(external = TRUE, size = n))
  }

  attribute <- sample(varying, 1L)

  split_value <- runif(
    1L,
    minimum[attribute],
    maximum[attribute]
  )

  left <- x[, attribute] < split_value

  list(
    external = FALSE,
    split_attribute = attribute,
    split_value = split_value,
    left = build_itree(
      x[left, , drop = FALSE],
      height + 1L,
      height_limit
    ),
    right = build_itree(
      x[!left, , drop = FALSE],
      height + 1L,
      height_limit
    )
  )
}

tree_path_length <- function(tree, x, height = 0L) {
  if (nrow(x) == 0L) return(numeric(0))

  if (tree$external) {
    return(rep(height + iforest_c(tree$size), nrow(x)))
  }

  left <- x[, tree$split_attribute] < tree$split_value
  result <- numeric(nrow(x))

  if (any(left)) {
    result[left] <- tree_path_length(
      tree$left,
      x[left, , drop = FALSE],
      height + 1L
    )
  }

  if (any(!left)) {
    result[!left] <- tree_path_length(
      tree$right,
      x[!left, , drop = FALSE],
      height + 1L
    )
  }

  result
}

fit_isolation_forest <- function(
    x_train,
    ntrees = 100L,
    sample_size = min(256L, nrow(x_train)),
    seed = 2030L
) {
  x_train <- iforest_matrix(x_train)

  sample_size <- min(
    as.integer(sample_size),
    nrow(x_train)
  )

  height_limit <- ceiling(log2(sample_size))
  set.seed(seed)

  trees <- lapply(
    seq_len(ntrees),
    function(i) {
      rows <- sample.int(
        nrow(x_train),
        sample_size,
        replace = FALSE
      )

      build_itree(
        x_train[rows, , drop = FALSE],
        0L,
        height_limit
      )
    }
  )

  list(
    trees = trees,
    ntrees = as.integer(ntrees),
    sample_size = sample_size,
    dimension = ncol(x_train)
  )
}

score_isolation_forest <- function(model, newdata) {
  newdata <- iforest_matrix(newdata)

  if (ncol(newdata) != model$dimension) {
    stop("newdata has the wrong dimension", call. = FALSE)
  }

  path_sum <- numeric(nrow(newdata))

  for (tree in model$trees) {
    path_sum <- path_sum + tree_path_length(tree, newdata)
  }

  mean_path <- path_sum / model$ntrees

  2^(
    -mean_path /
      iforest_c(model$sample_size)
  )
}
