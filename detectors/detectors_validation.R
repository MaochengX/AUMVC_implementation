args <- commandArgs(trailingOnly = FALSE)
script <- grep("^--file=", args, value = TRUE)
if (!length(script)) stop("Run with Rscript detectors/detectors_validation.R")
root <- dirname(dirname(normalizePath(sub("^--file=", "", script[1L]))))

for (package in c("dbscan", "e1071", "isotree")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install required package: ", package, call. = FALSE)
  }
}

for (file in c(
  "aumvc/input_validation.R", "aumvc/level_set.R", "aumvc/aumvc.R",
  "detectors/ocsvm.R", "detectors/lof.R", "detectors/isolation_forest.R"
)) source(file.path(root, file), local = FALSE)

same_scores <- function(x, y, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance)))
}

set.seed(1234L)
x_train <- rbind(
  matrix(rnorm(160, 0, 0.55), ncol = 2L),
  matrix(rnorm(40, 3, 0.45), ncol = 2L)
)
x_eval <- rbind(
  matrix(rnorm(90, 0, 0.55), ncol = 2L),
  matrix(rnorm(30, 3, 0.45), ncol = 2L)
)
x_reference <- matrix(rnorm(100), ncol = 2L)

ocsvm <- fit_ocsvm(x_train, nu = 0.5, gamma = 0.5)
package_ocsvm <- e1071::svm(
  x = x_train, y = NULL, type = "one-classification", kernel = "radial",
  gamma = 0.5, nu = 0.5, scale = FALSE, tolerance = 1e-6, fitted = FALSE
)
package_decisions <- attr(predict(package_ocsvm, x_eval, decision.values = TRUE),
                          "decision.values")
same_scores(score_ocsvm(ocsvm, x_eval), -package_decisions)
cat("OCSVM: identical to e1071 decision values (anomaly-score orientation).\n")

for (k in c(4L, 10L)) {
  lof <- fit_lof(x_train, k = k)
  package_scores <- dbscan::lof(x_train, minPts = k + 1L)
  same_scores(lof$training_scores, package_scores, tolerance = 1e-12)
  same_scores(score_lof(lof, x_train), package_scores, tolerance = 1e-12)
  stopifnot(all(is.finite(score_lof(lof, x_eval))))

  distance <- as.matrix(dist(x_train))
  diag(distance) <- Inf
  indices <- t(apply(distance, 1L, order))[, seq_len(k), drop = FALSE]
  kth <- distance[cbind(seq_len(nrow(x_train)), indices[, k])]
  densities <- vapply(seq_len(nrow(x_train)), function(i) {
    id <- indices[i, ]
    1 / max(mean(pmax(distance[i, id], kth[id])), lof$distance_floor)
  }, numeric(1))
  manual <- vapply(seq_len(5L), function(j) {
    ds <- sqrt(rowSums(sweep(x_train, 2L, x_eval[j, ], "-")^2))
    id <- order(ds)[seq_len(k)]
    density <- 1 / max(mean(pmax(ds[id], kth[id])), lof$distance_floor)
    mean(densities[id]) / density
  }, numeric(1))
  same_scores(score_lof(lof, x_eval[1:5, , drop = FALSE]), manual,
              tolerance = 1e-8)
}

x_duplicated <- rbind(
  matrix(rep(c(-8, -8), 8L), ncol = 2L, byrow = TRUE),
  matrix(rnorm(100, 2, 0.5), ncol = 2L)
)
lof_duplicates <- fit_lof(x_duplicated, k = 4L)
same_scores(lof_duplicates$training_scores,
            dbscan::lof(x_duplicated, minPts = 5L), tolerance = 1e-12)
stopifnot(all(lof_duplicates$training_scores[1:8] == 1))
stopifnot(all(is.finite(score_lof(lof_duplicates,
                                rbind(c(-8, -8), c(2, 2), c(10, 10))))))
cat("LOF: matches dbscan for training data; out-of-sample scores match the LOF formula.\n")

iforest <- fit_isolation_forest(x_train, ntrees = 50L,
                                sample_size = 64L, seed = 1234L)
package_iforest <- isotree::isolation.forest(
  data = x_train, sample_size = 64L, ntrees = 50L, ndim = 1L,
  max_depth = ceiling(log2(64L)), sample_with_replacement = FALSE,
  missing_action = "fail", scoring_metric = "depth", seed = 1234L,
  nthreads = 1L
)
same_scores(score_isolation_forest(iforest, x_eval),
            predict(package_iforest, x_eval, type = "score", nthreads = 1L))
cat("Isolation Forest: identical to isotree scores with matching settings.\n")

alpha <- seq(0.9, 0.999, by = 0.0001)
evaluation_scores <- seq(1, 0, length.out = 1000L)
reference_scores <- seq(1, 0, length.out = 100000L)
normality <- aumvc_from_scores(evaluation_scores, reference_scores, 1,
                                "normality", alpha)
anomaly <- aumvc_from_scores(-evaluation_scores, -reference_scores, 1,
                             "anomaly", alpha)
stopifnot(length(alpha) == 991L)
same_scores(normality$aumvc, anomaly$aumvc, tolerance = 1e-12)
same_scores(normality$aumvc_normalized, anomaly$aumvc_normalized,
            tolerance = 1e-12)

reference <- make_reference(x_reference, n_reference = 1500L,
                            n_mc_repetitions = 1L, seed = 1234L,
                            chunk_size = 250L)
lof_for_aumvc <- fit_lof(x_train, k = 10L)
score_funs <- list(
  OCSVM = function(x) score_ocsvm(ocsvm, x),
  LOF = function(x) score_lof(lof_for_aumvc, x),
  Isolation_Forest = function(x) score_isolation_forest(iforest, x)
)
for (name in names(score_funs)) {
  score_fun <- score_funs[[name]]
  scores <- score_fun(x_eval)
  stopifnot(length(scores) == nrow(x_eval), all(is.finite(scores)))
  mv <- aumvc(x_eval, reference, score_fun, score_direction = "anomaly",
              alpha_grid = seq(0.8, 0.95, length.out = 31L))
  stopifnot(is.finite(mv$aumvc_normalized), mv$aumvc_normalized >= 0,
            mv$aumvc_normalized <= 0.15 + 1e-8,
            all(is.finite(mv$mv_curve$volume_normalized)))
  flipped <- aumvc(x_eval, reference, function(x) -score_fun(x),
                   score_direction = "normality",
                   alpha_grid = seq(0.8, 0.95, length.out = 31L))
  same_scores(mv$aumvc_normalized, flipped$aumvc_normalized,
              tolerance = 1e-12)
  cat(sprintf("%-17s AUMVC = %.5f, zero occupancy = %.1f%%\n",
              name, mv$aumvc_normalized, 100 * mv$zero_occupancy))
}
cat("PASS: all three detectors agree with their package backends and support AUMVC.\n")
