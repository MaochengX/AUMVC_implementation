source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")
source("detectors/lof.R")
source("detectors/isolation_forest.R")

evaluation <- c(3, 2, 1, 0)
reference <- c(3, 2, 1, 0)
mv <- aumvc_from_scores(evaluation, reference, 4,
                        alpha_grid = c(0.5, 0.75))
stopifnot(abs(mv$aumvc_normalized - 0.15625) < 1e-12,
          abs(mv$aumvc - 0.625) < 1e-12)

flipped <- aumvc_from_scores(-evaluation, -reference, 4,
                             score_direction = "anomaly",
                             alpha_grid = c(0.5, 0.75))
stopifnot(isTRUE(all.equal(flipped$aumvc, mv$aumvc)))

set.seed(1234)
train <- matrix(rnorm(120), ncol = 2L)
new_points <- rbind(c(0, 0), c(7, 7))
ocsvm <- fit_ocsvm(train, nu = 0.2, gamma = 0.5)
lof <- fit_lof(train, k = 10L)
forest <- fit_isolation_forest(train, ntrees = 50L,
                               sample_size = 50L, seed = 1234L)
stopifnot(ocsvm$converged,
          score_ocsvm(ocsvm, new_points)[2L] >
            score_ocsvm(ocsvm, new_points)[1L],
          score_lof(lof, new_points)[2L] >
            score_lof(lof, new_points)[1L],
          score_isolation_forest(forest, new_points)[2L] >
            score_isolation_forest(forest, new_points)[1L])

# Optional external comparison, with no Python dependency.
if (requireNamespace("e1071", quietly = TRUE)) {
  standard <- e1071::svm(train, type = "one-classification",
                         kernel = "radial", gamma = 0.5,
                         nu = 0.2, scale = FALSE)
  predictions <- predict(standard, new_points, decision.values = TRUE)
  external <- as.numeric(attr(predictions, "decision.values"))
  stopifnot(identical(order(-score_ocsvm(ocsvm, new_points)),
                      order(external)))
}
cat("Core checks passed\n")
