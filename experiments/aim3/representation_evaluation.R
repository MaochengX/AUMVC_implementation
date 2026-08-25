source("experiments/settings.R")

source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")

source("experiments/aim3/mds_embedding.R")
source("experiments/aim3/reference_sampling.R")
source("experiments/aim3/goix_subsampling.R")

aim3_make_splits <- function(n, seed, settings) {
  counts <- as.integer(settings$split_counts)

  if (sum(counts) != n || any(counts < 1L)) {
    stop("Aim 3 split counts must partition the sample.", call. = FALSE)
  }

  set.seed(seed)
  rows <- sample.int(n)
  ends <- cumsum(counts)
  starts <- c(1L, head(ends, -1L) + 1L)

  split <- Map(
    function(first, last) rows[first:last],
    starts,
    ends
  )

  names(split) <- names(settings$split_counts)
  split
}

aim3_fit_standardizer <- function(x) {
  center <- colMeans(x)
  scale <- apply(x, 2L, sd)
  scale[!is.finite(scale) | scale == 0] <- 1

  list(center = center, scale = scale)
}

aim3_apply_standardizer <- function(x, standardizer) {
  x <- sweep(x, 2L, standardizer$center, "-")
  sweep(x, 2L, standardizer$scale, "/")
}

aim3_evaluate_representation <- function(
    representation,
    x_train,
    x_reference,
    x_eval,
    seed,
    settings
) {
  start <- proc.time()[[3L]]

  model <- fit_ocsvm(
    x_train,
    nu = settings$detector$nu,
    gamma = 1 / ncol(x_train)
  )

  reference <- aim3_make_reference(
    x_reference,
    seed = seed + 100L,
    settings = settings$reference
  )

  score_fun <- function(x) score_ocsvm(model, x)

  mv <- aumvc(
    x_eval = x_eval,
    reference = reference,
    score_fun = score_fun,
    score_direction = "anomaly",
    alpha_grid = settings$aumvc_alpha_grid
  )

  data.frame(
    representation = representation,
    aumvc = mv$aumvc,
    aumvc_normalized = mv$aumvc_normalized,
    aumvc_mc_se = mv$aumvc_mc_se,
    aumvc_normalized_mc_se = mv$aumvc_normalized_mc_se,
    mean_occupancy = mean(mv$mv_curve$volume_normalized),
    zero_occupancy_fraction = mean(
      mv$mv_curve$volume_normalized == 0
    ),
    log_box_volume = reference$box$log_volume,
    n_reference = reference$n_reference,
    reference_capped = reference$density_capped,
    log_reference_density = reference$log_density,
    runtime_seconds = proc.time()[[3L]] - start
  )
}

aim3_evaluate_representations <- function(case, seed, settings) {
  split <- aim3_make_splits(
    nrow(case$x),
    seed,
    settings
  )

  standardizer <- aim3_fit_standardizer(
    case$x[split$embedding, , drop = FALSE]
  )

  x_embedding <- aim3_apply_standardizer(
    case$x[split$embedding, , drop = FALSE],
    standardizer
  )

  x_train <- aim3_apply_standardizer(
    case$x[split$detector_train, , drop = FALSE],
    standardizer
  )

  x_reference <- aim3_apply_standardizer(
    case$x[split$reference, , drop = FALSE],
    standardizer
  )

  x_eval <- aim3_apply_standardizer(
    case$x[split$evaluation, , drop = FALSE],
    standardizer
  )

  ambient <- aim3_evaluate_representation(
    "ambient",
    x_train,
    x_reference,
    x_eval,
    seed + 1000L,
    settings
  )

  goix <- aim3_evaluate_goix_subsampling(
    x_train,
    x_reference,
    x_eval,
    seed + 2000L,
    settings
  )

  mds <- aim3_fit_mds(
    x_embedding,
    ndim = case$intrinsic_dim,
    settings = settings$mds
  )

  mds_train <- aim3_project_mds(mds, x_train)
  mds_reference <- aim3_project_mds(mds, x_reference)
  mds_eval <- aim3_project_mds(mds, x_eval)

  mds_result <- aim3_evaluate_representation(
    "mds",
    mds_train,
    mds_reference,
    mds_eval,
    seed + 3000L,
    settings
  )

  rbind(
    ambient,
    goix,
    mds_result
  )
}
