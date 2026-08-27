source("experiments/settings.R")
source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("aumvc/auemc.R")

settings <- AIM1_SETTINGS
set.seed(experiment_seed(settings, 2L))

left <- -settings$two_clusters$cluster_mean
right <- settings$two_clusters$cluster_mean
n_cluster <- settings$two_clusters$n_per_cluster
cluster_sd <- settings$two_clusters$cluster_sd

x_train <- rbind(
  matrix(rnorm(n_cluster * 2, left, cluster_sd), ncol = 2),
  matrix(rnorm(n_cluster * 2, right, cluster_sd), ncol = 2)
)

x_eval <- rbind(
  matrix(rnorm((n_cluster - 1L) * 2, left, cluster_sd), ncol = 2),
  matrix(rnorm((n_cluster - 1L) * 2, right, cluster_sd), ncol = 2),
  matrix(c(0, 7, 7, 0), ncol = 2, byrow = TRUE)
)

centers <- rbind(c(left, left), c(right, right))
squared_distance <- function(x, center) rowSums(sweep(x, 2L, center, "-")^2)
normality_score <- function(x) pmax(
  -squared_distance(x, centers[1L, ]),
  -squared_distance(x, centers[2L, ])
)
anomaly_score <- function(x) -normality_score(x)

reference <- make_reference(
  x_train,
  n_reference = settings$n_reference,
  n_mc_repetitions = settings$n_mc_repetitions,
  seed = experiment_seed(settings, 102L)
)

mv_normality <- aumvc(
  x_eval,
  reference,
  normality_score,
  score_direction = "normality",
  alpha_grid = settings$aumvc_alpha_grid
)
mv_anomaly <- aumvc(
  x_eval,
  reference,
  anomaly_score,
  score_direction = "anomaly",
  alpha_grid = settings$aumvc_alpha_grid
)
em_normality <- auemc(
  x_eval,
  reference,
  normality_score,
  score_direction = "normality",
  tau_grid = settings$auemc_tau_grid
)
em_anomaly <- auemc(
  x_eval,
  reference,
  anomaly_score,
  score_direction = "anomaly",
  tau_grid = settings$auemc_tau_grid
)

stopifnot(
  isTRUE(all.equal(mv_normality$aumvc, mv_anomaly$aumvc, tolerance = 1e-12)),
  isTRUE(all.equal(em_normality$auemc, em_anomaly$auemc, tolerance = 1e-12))
)

print(data.frame(
  metric = c("AUMVC", "AUEMC"),
  value = c(mv_normality$aumvc, em_normality$auemc)
))
cat("Score-direction test passed\n")
