source("experiments/settings.R")

source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("aumvc/auemc.R")

set.seed(experiment_seed(2L))

cfg <- AIM1_TWO_CLUSTERS

mean_left <- -cfg$cluster_mean
mean_right <- cfg$cluster_mean

x_train <- rbind(
  matrix(
    rnorm(
      cfg$n_per_cluster * 2,
      mean = mean_left,
      sd = cfg$cluster_sd
    ),
    ncol = 2
  ),
  matrix(
    rnorm(
      cfg$n_per_cluster * 2,
      mean = mean_right,
      sd = cfg$cluster_sd
    ),
    ncol = 2
  )
)

x_eval <- rbind(
  matrix(
    rnorm(
      (cfg$n_per_cluster - 1L) * 2,
      mean = mean_left,
      sd = cfg$cluster_sd
    ),
    ncol = 2
  ),
  matrix(
    rnorm(
      (cfg$n_per_cluster - 1L) * 2,
      mean = mean_right,
      sd = cfg$cluster_sd
    ),
    ncol = 2
  ),
  matrix(
    c(
      0, 7,
      7, 0
    ),
    ncol = 2,
    byrow = TRUE
  )
)

centers <- rbind(
  c(mean_left, mean_left),
  c(mean_right, mean_right)
)

squared_distance <- function(x, center) {
  rowSums(
    sweep(x, 2L, center, "-")^2
  )
}

normality_score <- function(x) {
  pmax(
    -squared_distance(x, centers[1L, ]),
    -squared_distance(x, centers[2L, ])
  )
}

anomaly_score <- function(x) {
  -normality_score(x)
}

reference <- make_reference(
  x_train,
  n_reference = N_REFERENCE,
  n_mc_repetitions = N_MC_REPETITIONS,
  seed = experiment_seed(102L)
)

mv_normality <- aumvc(
  x_eval = x_eval,
  reference = reference,
  score_fun = normality_score,
  score_direction = "normality",
  alpha_grid = AUMVC_ALPHA_GRID
)

mv_anomaly <- aumvc(
  x_eval = x_eval,
  reference = reference,
  score_fun = anomaly_score,
  score_direction = "anomaly",
  alpha_grid = AUMVC_ALPHA_GRID
)

em_normality <- auemc(
  x_eval = x_eval,
  reference = reference,
  score_fun = normality_score,
  score_direction = "normality",
  tau_grid = AUEMC_TAU_GRID
)

em_anomaly <- auemc(
  x_eval = x_eval,
  reference = reference,
  score_fun = anomaly_score,
  score_direction = "anomaly",
  tau_grid = AUEMC_TAU_GRID
)

print(
  data.frame(
    metric = c("AUMVC", "AUEMC"),
    normality_input = c(
      mv_normality$aumvc,
      em_normality$auemc
    ),
    anomaly_input = c(
      mv_anomaly$aumvc,
      em_anomaly$auemc
    )
  )
)

stopifnot(
  isTRUE(
    all.equal(
      mv_normality$aumvc,
      mv_anomaly$aumvc,
      tolerance = 1e-12
    )
  )
)

stopifnot(
  isTRUE(
    all.equal(
      em_normality$auemc,
      em_anomaly$auemc,
      tolerance = 1e-12
    )
  )
)
