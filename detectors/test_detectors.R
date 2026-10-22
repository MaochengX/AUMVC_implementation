args <- commandArgs(trailingOnly = FALSE)
script <- grep("^--file=", args, value = TRUE)
if (!length(script)) stop("Run with Rscript detectors/test_detectors.R")
root <- dirname(dirname(normalizePath(sub("^--file=", "", script[1L]))))

for (name in c("dbscan", "e1071", "isotree")) {
  if (!requireNamespace(name, quietly = TRUE)) {
    stop("Install required package: ", name, call. = FALSE)
  }
}

for (file in c(
  "aumvc/input_validation.R", "aumvc/level_set.R", "aumvc/aumvc.R",
  "detectors/ocsvm.R", "detectors/lof.R", "detectors/isolation_forest.R"
)) source(file.path(root, file), local = FALSE)

set.seed(1234)
x_train <- rbind(
  matrix(rnorm(160, 0, 0.55), ncol = 2),
  matrix(rnorm(40, 3, 0.45), ncol = 2)
)
x_eval <- rbind(
  matrix(rnorm(90, 0, 0.55), ncol = 2),
  matrix(rnorm(30, 3, 0.45), ncol = 2)
)
x_reference <- matrix(rnorm(100), ncol = 2)

for (k in c(4L, 10L)) {
  lof <- fit_lof(x_train, k = k)
  package_scores <- dbscan::lof(x_train, minPts = k + 1L)
  stopifnot(isTRUE(all.equal(lof$training_scores, as.numeric(package_scores),
                             tolerance = 1e-12)))
  stopifnot(isTRUE(all.equal(score_lof(lof, x_train), as.numeric(package_scores),
                             tolerance = 1e-12)))
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
  stopifnot(isTRUE(all.equal(score_lof(lof, x_eval[1:5, ]), manual,
                             tolerance = 1e-8)))
}

x_duplicated <- rbind(
  matrix(rep(c(-8, -8), 8L), ncol = 2L, byrow = TRUE),
  matrix(rnorm(100, 2, 0.5), ncol = 2L)
)
lof_duplicates <- fit_lof(x_duplicated, k = 4L)
reference_duplicates <- dbscan::lof(x_duplicated, minPts = 5L)
stopifnot(isTRUE(all.equal(lof_duplicates$training_scores,
                           as.numeric(reference_duplicates), tolerance = 1e-12)))
stopifnot(all(lof_duplicates$training_scores[1:8] == 1))
duplicate_query <- rbind(c(-8, -8), c(2, 2), c(10, 10))
stopifnot(all(is.finite(score_lof(lof_duplicates, duplicate_query))))
cat("LOF: dbscan training scores agree, including duplicated rows.\n")

models <- list(
  OCSVM = fit_ocsvm(x_train, nu = 0.5, gamma = 0.5),
  LOF = fit_lof(x_train, k = 10L),
  Isolation_Forest = fit_isolation_forest(
    x_train, ntrees = 50L, sample_size = 64L, seed = 1234L
  )
)
score_funs <- list(
  OCSVM = function(x) score_ocsvm(models$OCSVM, x),
  LOF = function(x) score_lof(models$LOF, x),
  Isolation_Forest = function(x) score_isolation_forest(models$Isolation_Forest, x)
)

reference <- make_reference(x_reference, n_reference = 1500L,
                            n_mc_repetitions = 1L, seed = 1234L,
                            chunk_size = 250L)
alpha <- seq(0.8, 0.95, length.out = 31L)
for (name in names(score_funs)) {
  scores <- score_funs[[name]](x_eval)
  stopifnot(length(scores) == nrow(x_eval), all(is.finite(scores)))
  mv <- aumvc(x_eval, reference, score_funs[[name]],
              score_direction = "anomaly", alpha_grid = alpha)
  stopifnot(
    length(mv$mv_curve$alpha) == length(alpha),
    is.finite(mv$aumvc_normalized),
    mv$aumvc_normalized >= 0,
    mv$aumvc_normalized <= diff(range(alpha)) + 1e-8,
    all(is.finite(mv$mv_curve$volume_normalized))
  )
  mv_normality <- aumvc(x_eval, reference,
                        function(x) -score_funs[[name]](x),
                        score_direction = "normality", alpha_grid = alpha)
  stopifnot(isTRUE(all.equal(mv$aumvc_normalized,
                             mv_normality$aumvc_normalized,
                             tolerance = 1e-12)))
  cat(sprintf("%-17s AUMVC normalized = %.5f, MC zero levels = %.1f%%\n",
              name, mv$aumvc_normalized, 100 * mv$zero_occupancy))
}
cat("PASS: all three package-backed detectors score new points for AUMVC.\n")
