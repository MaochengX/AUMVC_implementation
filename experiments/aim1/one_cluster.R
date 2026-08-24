source("experiments/settings.R")

source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("aumvc/auemc.R")

set.seed(aim_seed(AIM1_SETTINGS, 1L))

cfg <- AIM1_SETTINGS$one_cluster

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
  n_reference = AIM1_SETTINGS$n_reference,
  n_mc_repetitions = AIM1_SETTINGS$n_mc_repetitions,
  seed = aim_seed(AIM1_SETTINGS, 101L)
)

mv <- aumvc(
  x_eval = x_eval,
  reference = reference,
  score_fun = score_fun,
  score_direction = "normality",
  alpha_grid = AIM1_SETTINGS$aumvc_alpha_grid
)

em <- auemc(
  x_eval = x_eval,
  reference = reference,
  score_fun = score_fun,
  score_direction = "normality",
  tau_grid = AIM1_SETTINGS$auemc_tau_grid
)

print(
  data.frame(
    metric = c("AUMVC", "AUEMC"),
    value = c(mv$aumvc, em$auemc),
    mc_se = c(mv$aumvc_mc_se, em$auemc_mc_se)
  )
)
