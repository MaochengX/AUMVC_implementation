fit_isolation_forest <- function(
    x_train,
    ntrees = 100L,
    sample_size = min(256L, nrow(x_train)),
    seed = 2026L
) {
  if (!requireNamespace("isotree", quietly = TRUE)) {
    stop("Package 'isotree' is required.", call. = FALSE)
  }

  x_train <- as.matrix(x_train)
  storage.mode(x_train) <- "double"

  model <- isotree::isolation.forest(
    data = x_train,
    sample_size = as.integer(sample_size),
    ntrees = as.integer(ntrees),
    ndim = 1L,
    missing_action = "fail",
    seed = as.integer(seed),
    nthreads = 1L
  )

  list(
    model = model,
    dimension = ncol(x_train),
    ntrees = as.integer(ntrees),
    sample_size = as.integer(sample_size),
    seed = as.integer(seed)
  )
}

score_isolation_forest <- function(
    model,
    newdata
) {
  newdata <- as.matrix(newdata)
  storage.mode(newdata) <- "double"

  as.numeric(
    predict(
      model$model,
      newdata,
      type = "score",
      nthreads = 1L
    )
  )
}
