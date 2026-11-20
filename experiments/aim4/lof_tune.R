source("experiments/aim4/tuning_common.R")
source("detectors/lof.R")

LOF_K_GRID <- c(2L, 3L, 5L, 7L, 10L, 12L, 15L, 17L, 20L, 25L, 30L, 40L)
LOF_DEFAULT_K <- 20L

score_lof_grid_shared <- function(x_train, x_query, models) {
  k_max <- max(vapply(models, `[[`, integer(1), "k"))
  neighbors <- dbscan::kNN(x_train, k = k_max, query = x_query, sort = TRUE)
  index <- as.matrix(neighbors$id)
  distance <- as.matrix(neighbors$dist)
  scores <- vapply(models, function(model) {
    k <- model$k
    ids <- index[, seq_len(k), drop = FALSE]
    d <- distance[, seq_len(k), drop = FALSE]
    reachability <- pmax(d, matrix(model$k_distance[ids], nrow = nrow(d), ncol = k))
    query_lrd <- 1 / pmax(rowMeans(reachability), model$distance_floor)
    neighbor_lrd <- matrix(model$training_lrd[ids], nrow = nrow(d), ncol = k)
    result <- rowMeans(neighbor_lrd) / query_lrd
    result[rowSums(d == 0) == k] <- 1
    if (any(!is.finite(result))) stop("Non-finite LOF candidate scores.", call. = FALSE)
    result
  }, numeric(nrow(x_query)))
  matrix(scores, nrow = nrow(x_query), ncol = length(models))
}

score_lof_reference_grid <- function(reference, x_train, models) {
  n <- reference$n_reference
  scores <- matrix(NA_real_, n * reference$n_mc_repetitions, length(models))
  for (r in seq_len(reference$n_mc_repetitions)) {
    set.seed(reference$seeds[r])
    offset <- (r - 1L) * n
    for (start in seq(1L, n, by = reference$chunk_size)) {
      end <- min(start + reference$chunk_size - 1L, n)
      count <- end - start + 1L
      points <- matrix(runif(count * reference$box$dimension), nrow = count, byrow = TRUE)
      points <- sweep(points, 2L, reference$box$width, "*")
      points <- sweep(points, 2L, reference$box$lower, "+")
      scores[offset + seq.int(start, end), ] <- score_lof_grid_shared(x_train, points, models)
    }
  }
  scores
}

score_lof_grid <- function(prepared, ks, settings) {
  train <- prepared$detector_train
  if (max(ks) >= nrow(train)) stop("The largest LOF k must be below training size.", call. = FALSE)
  models <- lapply(ks, function(k) fit_lof(train, k = k))
  evaluation <- score_lof_grid_shared(train, prepared$evaluation, models)
  validation <- score_lof_grid_shared(train, prepared$validation, models)
  test <- score_lof_grid_shared(train, prepared$test, models)

  check_rows <- seq_len(min(20L, nrow(prepared$evaluation)))
  for (i in seq_along(models)) {
    expected <- score_lof(models[[i]], prepared$evaluation[check_rows, , drop = FALSE])
    if (!isTRUE(all.equal(as.numeric(evaluation[check_rows, i]), expected,
                          tolerance = 1e-8, check.attributes = FALSE))) {
      stop("Shared neighbor scoring does not match existing LOF scoring.", call. = FALSE)
    }
  }

  mc <- score_lof_reference_grid(prepared$reference_mc, train, models)
  results <- score_aim4b_candidates(
    prepared, ks, evaluation, validation, test, mc, settings
  )
  results$minPts <- ks + 1L
  results
}

run_lof_tune <- function(settings = AIM4B_SETTINGS) {
  run_aim4b("lof", settings, score_lof_grid, LOF_K_GRID, LOF_DEFAULT_K, min(LOF_K_GRID))
}

if (sys.nframe() == 0L) run_lof_tune()
