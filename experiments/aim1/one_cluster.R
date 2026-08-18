source("experiments/settings.R")

source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("aumvc/auemc.R")

set.seed(experiment_seed(1L))

cfg <- AIM1_ONE_CLUSTER

x_train <- matrix(
  rnorm(cfg$n_train * 2),
  ncol = 2
)

x_normal <- matrix(
  rnorm((cfg$n_eval - cfg$n_anomaly) * 2),
  ncol = 2
)

x_anomaly <- matrix(
  rnorm(
    cfg$n_anomaly * 2,
    mean = cfg$anomaly_mean,
    sd = cfg$anomaly_sd
  ),
  ncol = 2
)

x_eval <- rbind(
  x_normal,
  x_anomaly
)

center <- colMeans(x_train)
covariance <- cov(x_train)

score_fun <- function(x) {
  -mahalanobis(
    x,
    center = center,
    cov = covariance
  )
}

reference <- make_reference(
  x_train,
  n_reference = N_REFERENCE,
  n_mc_repetitions = N_MC_REPETITIONS,
  seed = experiment_seed(101L)
)

mv <- aumvc(
  x_eval = x_eval,
  reference = reference,
  score_fun = score_fun,
  score_direction = "normality",
  alpha_grid = AUMVC_ALPHA_GRID
)

em <- auemc(
  x_eval = x_eval,
  reference = reference,
  score_fun = score_fun,
  score_direction = "normality",
  tau_grid = AUEMC_TAU_GRID
)

print(
  data.frame(
    metric = c("AUMVC", "AUEMC"),
    value = c(mv$aumvc, em$auemc),
    mc_se = c(mv$aumvc_mc_se, em$auemc_mc_se)
  )
)
