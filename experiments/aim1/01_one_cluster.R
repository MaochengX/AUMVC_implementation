source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")

set.seed(2026)

n_train <- 200
n_eval <- 200
n_anomaly <- 2

x_train <- matrix(
  rnorm(n_train * 2),
  ncol = 2
)

x_eval <- rbind(
  matrix(
    rnorm((n_eval - n_anomaly) * 2),
    ncol = 2
  ),
  matrix(
    rnorm(
      n_anomaly * 2,
      mean = 6,
      sd = 0.3
    ),
    ncol = 2
  )
)

x_train <- validate_reference_inputs(
  x_train,
  n_reference = 20000L,
  n_mc_repetitions = 5L,
  seed = 2026L
)

center <- colMeans(x_train)
covariance <- cov(x_train)

score_fun <- function(newdata) {
  -mahalanobis(
    newdata,
    center = center,
    cov = covariance
  )
}

reference <- make_reference(
  x_train,
  n_reference = 20000L,
  n_mc_repetitions = 5L,
  seed = 2026L
)

x_eval <- validate_aumvc_inputs(
  x_eval,
  reference,
  score_fun,
  "normality",
  seq(0.9, 0.999, length.out = 100L)
)

validate_scores(
  score_fun(x_eval),
  n_expected = nrow(x_eval),
  name = "evaluation scores"
)

result <- aumvc(
  x_eval,
  reference,
  score_fun,
  score_direction = "normality"
)

print(
  data.frame(
    aumvc = result$aumvc,
    aumvc_normalized = result$aumvc_normalized,
    aumvc_mc_se = result$aumvc_mc_se
  )
)
