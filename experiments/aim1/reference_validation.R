source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")
source("detectors/lof.R")
source("detectors/isolation_forest.R")

alpha <- seq(0.9, 0.999, by = 0.0001)
evaluation_scores <- seq(1, 0, length.out = 1000L)
reference_scores <- seq(1, 0, length.out = 100000L)
normality <- aumvc_from_scores(
  evaluation_scores,
  reference_scores,
  1,
  "normality",
  alpha
)
anomaly <- aumvc_from_scores(
  -evaluation_scores,
  -reference_scores,
  1,
  "anomaly",
  alpha
)
stopifnot(
  length(alpha) == 991L,
  isTRUE(all.equal(normality$aumvc, anomaly$aumvc, tolerance = 1e-12)),
  isTRUE(all.equal(
    normality$aumvc_normalized,
    anomaly$aumvc_normalized,
    tolerance = 1e-12
  ))
)

set.seed(1234L)
x <- rbind(
  matrix(rnorm(160L), ncol = 2L),
  matrix(rnorm(40L, 4), ncol = 2L)
)
ocsvm <- fit_ocsvm(x, nu = 0.5, gamma = 1 / ncol(x))
lof <- fit_lof(x, k = 20L)
iforest <- fit_isolation_forest(x, ntrees = 100L, sample_size = 100L, seed = 1234L)
scores <- list(
  OCSVM = score_ocsvm(ocsvm, x),
  LOF = score_lof(lof, x),
  Isolation_Forest = score_isolation_forest(iforest, x)
)
stopifnot(
  ocsvm$converged,
  all(vapply(scores, length, integer(1)) == nrow(x)),
  all(vapply(scores, function(value) all(is.finite(value)), logical(1)))
)

cat("AUMVC and detector reference validation passed.\n")
