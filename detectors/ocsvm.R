fit_ocsvm <- function(
    x_train,
    nu = 0.5,
    gamma = 1 / ncol(x_train)
) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("Package 'e1071' is required.", call. = FALSE)
  }

  x_train <- as.matrix(x_train)
  storage.mode(x_train) <- "double"

  model <- e1071::svm(
    x = x_train,
    y = NULL,
    type = "one-classification",
    kernel = "radial",
    nu = nu,
    gamma = gamma,
    scale = TRUE
  )

  list(
    model = model,
    dimension = ncol(x_train),
    nu = nu,
    gamma = gamma
  )
}

score_ocsvm <- function(model, newdata) {
  newdata <- as.matrix(newdata)
  storage.mode(newdata) <- "double"

  prediction <- predict(
    model$model,
    newdata,
    decision.values = TRUE
  )

  decision <- attr(
    prediction,
    "decision.values"
  )

  -as.numeric(decision)
}
