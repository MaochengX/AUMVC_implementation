ocsvm_matrix <- function(x) {
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

fit_ocsvm <- function(
    x_train,
    nu = 0.5,
    gamma = 1 / ncol(x_train),
    tolerance = 1e-6
) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("Package 'e1071' is required", call. = FALSE)
  }

  x_train <- ocsvm_matrix(x_train)

  if (
    length(nu) != 1L || !is.finite(nu) || nu <= 0 || nu > 1 ||
    length(gamma) != 1L || !is.finite(gamma) || gamma <= 0 ||
    length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0
  ) {
    stop("Invalid OCSVM parameters", call. = FALSE)
  }

  model <- e1071::svm(
    x = x_train,
    y = NULL,
    type = "one-classification",
    kernel = "radial",
    gamma = gamma,
    nu = nu,
    scale = FALSE,
    tolerance = tolerance,
    fitted = FALSE
  )

  list(
    model = model,
    dimension = ncol(x_train),
    nu = nu,
    gamma = gamma,
    tolerance = tolerance
  )
}

score_ocsvm <- function(model, newdata) {
  newdata <- ocsvm_matrix(newdata)

  if (ncol(newdata) != model$dimension) {
    stop("newdata has the wrong dimension", call. = FALSE)
  }

  prediction <- predict(
    model$model,
    newdata,
    decision.values = TRUE
  )
  decision <- attr(prediction, "decision.values")

  if (is.null(decision)) {
    stop("OCSVM decision values are unavailable", call. = FALSE)
  }

  -as.numeric(decision)
}
