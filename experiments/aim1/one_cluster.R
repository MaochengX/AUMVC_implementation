source("experiments/aim1/settings.R")
source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")

settings <- AIM1_SETTINGS
set.seed(experiment_seed(settings, 1L))

x_train <- matrix(rnorm(settings$one_cluster$n_train * 2), ncol = 2)
x_normal <- matrix(
  rnorm((settings$one_cluster$n_eval - settings$one_cluster$n_anomaly) * 2),
  ncol = 2
)
x_anomaly <- matrix(
  rnorm(
    settings$one_cluster$n_anomaly * 2,
    mean = settings$one_cluster$anomaly_mean,
    sd = settings$one_cluster$anomaly_sd
  ),
  ncol = 2
)
x_eval <- rbind(x_normal, x_anomaly)

center <- colMeans(x_train)
covariance <- cov(x_train)
score_fun <- function(x) -mahalanobis(x, center = center, cov = covariance)

reference <- make_reference(
  x_train,
  n_reference = settings$n_reference,
  n_mc_repetitions = settings$n_mc_repetitions,
  seed = experiment_seed(settings, 101L)
)

mv <- aumvc(
  x_eval,
  reference,
  score_fun,
  score_direction = "normality",
  alpha_grid = settings$aumvc_alpha_grid
)

result <- data.frame(
  AUMVC = mv$aumvc,
  MC_SE = mv$aumvc_mc_se
)
print(result)
