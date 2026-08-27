source("experiments/utils.R")
source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")
source("experiments/aim3/mds_embedding.R")

aim3_make_reference <- function(x_reference, seed, settings) {
  make_reference(
    x_reference,
    n_reference = settings$n_reference,
    n_mc_repetitions = settings$n_mc_repetitions,
    seed = seed
  )
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

  if (!model$converged) {
    warning(representation, " OCSVM solver did not converge.", call. = FALSE)
  }

  reference <- aim3_make_reference(
    x_reference,
    seed + 100L,
    settings
  )

  mv <- aumvc(
    x_eval,
    reference,
    function(x) score_ocsvm(model, x),
    score_direction = "anomaly",
    alpha_grid = settings$aumvc_alpha_grid
  )

  data.frame(
    representation = representation,
    aumvc = mv$aumvc,
    aumvc_normalized = mv$aumvc_normalized,
    aumvc_mc_se = mv$aumvc_mc_se,
    aumvc_normalized_mc_se = mv$aumvc_normalized_mc_se,
    zero_occupancy = mean(mv$mv_curve$volume_normalized == 0),
    runtime_seconds = proc.time()[[3L]] - start
  )
}

source("experiments/aim3/goix_subsampling.R")

aim3_compare_representations <- function(case, seed, settings) {
  split <- make_splits(
    nrow(case$x),
    settings$split_counts,
    names(settings$split_counts),
    seed,
    require_full = TRUE
  )

  standardizer <- fit_standardizer(
    case$x[split$embedding, , drop = FALSE]
  )

  x_embedding <- apply_standardizer(
    case$x[split$embedding, , drop = FALSE],
    standardizer
  )
  x_train <- apply_standardizer(
    case$x[split$detector_train, , drop = FALSE],
    standardizer
  )
  x_reference <- apply_standardizer(
    case$x[split$reference, , drop = FALSE],
    standardizer
  )
  x_eval <- apply_standardizer(
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
    case$intrinsic_dim
  )
  mds_standardizer <- fit_standardizer(mds$points)

  mds_result <- aim3_evaluate_representation(
    "mds",
    apply_standardizer(aim3_project_mds(mds, x_train), mds_standardizer),
    apply_standardizer(aim3_project_mds(mds, x_reference), mds_standardizer),
    apply_standardizer(aim3_project_mds(mds, x_eval), mds_standardizer),
    seed + 3000L,
    settings
  )

  rbind(ambient, goix, mds_result)
}
