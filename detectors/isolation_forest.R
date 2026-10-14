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

fit_isolation_forest <- function(
    x_train,
    ntrees = 100L,
    sample_size = min(256L, nrow(x_train)),
    seed = 1234L
) {
  if (!requireNamespace("isotree", quietly = TRUE)) {
    stop("Package 'isotree' is required", call. = FALSE)
  }

  x_train <- iforest_matrix(x_train)
  ntrees <- as.integer(ntrees)
  sample_size <- as.integer(sample_size)
  seed <- as.integer(seed)

  if (length(ntrees) != 1L || is.na(ntrees) || ntrees < 1L) {
    stop("ntrees must be a positive integer", call. = FALSE)
  }

  sample_size <- min(sample_size, nrow(x_train))
  if (
    length(sample_size) != 1L ||
    is.na(sample_size) ||
    sample_size < 2L
  ) {
    stop("sample_size must be at least 2", call. = FALSE)
  }

  model <- isotree::isolation.forest(
    data = x_train,
    sample_size = sample_size,
    ntrees = ntrees,
    ndim = 1L,
    max_depth = ceiling(log2(sample_size)),
    sample_with_replacement = FALSE,
    missing_action = "fail",
    scoring_metric = "depth",
    seed = seed,
    nthreads = 1L
  )

  list(
    model = model,
    dimension = ncol(x_train),
    ntrees = ntrees,
    sample_size = sample_size,
    seed = seed
  )
}

score_isolation_forest <- function(model, newdata) {
  newdata <- iforest_matrix(newdata)

  if (ncol(newdata) != model$dimension) {
    stop("newdata has the wrong dimension", call. = FALSE)
  }

  as.numeric(
    predict(
      model$model,
      newdata,
      type = "score",
      nthreads = 1L
    )
  )
}
