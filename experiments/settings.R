SEED <- 1234L
N_RUNS <- 10L

N_REFERENCE <- 20000L
N_MC_REPETITIONS <- 5L

AUMVC_ALPHA_GRID <- seq(0.9, 0.999, by = 0.0001)
AUEMC_TAU_GRID <- c(0, 10^seq(-6, 6, length.out = 800))

CONCORDANCE_TOLERANCE <- 1e-8

OCSVM_NU <- 0.5
LOF_K <- 20L

IFOREST_NTREES <- 100L
IFOREST_SAMPLE_SIZE <- 256L

AIM1_ONE_CLUSTER <- list(
  n_train = 200L,
  n_eval = 200L,
  n_anomaly = 2L,
  anomaly_mean = 5,
  anomaly_sd = 0.3
)

AIM1_TWO_CLUSTERS <- list(
  n_per_cluster = 100L,
  cluster_mean = 2,
  cluster_sd = 0.7
)

AIM2_SYNTHETIC <- list(
  n_normal = 1800L,
  n_anomaly = 200L,
  radius_min = 1.8,
  radius_max = 3.2
)

AIM2_SPLIT_COUNTS <- list(
  synthetic = c(500L, 400L, 500L, 600L),
  adult = c(1500L, 2000L, 2000L, 20000L),
  http = c(1500L, 2000L, 2000L, 50000L),
  pima = c(200L, 150L, 250L, 168L),
  smtp = c(1500L, 2000L, 2000L, 89656L),
  wilt = c(1200L, 1200L, 1200L, 1239L)
)

run_seed <- function(run) {
  SEED + (as.integer(run) - 1L) * 1000L
}

experiment_seed <- function(offset = 0L) {
  SEED + as.integer(offset)
}
