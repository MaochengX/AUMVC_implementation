source("experiments/utils.R")
source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")
source("experiments/aim3/mds_embedding.R")

aim3_variable_reference_columns <- function(x_reference) {
  lower <- apply(x_reference, 2L, min)
  upper <- apply(x_reference, 2L, max)

  is.finite(lower) &
    is.finite(upper) &
    upper > lower
}

aim3_filter_reference_constants <- function(
    x_train,
    x_reference,
    x_eval,
    representation
) {
  keep <- aim3_variable_reference_columns(x_reference)

  if (!any(keep)) {
    stop(
      representation,
      " has no variable reference coordinates.",
      call. = FALSE
    )
  }

  removed <- sum(!keep)

  if (removed > 0L) {
    message(
      representation,
      ": removed ",
      removed,
      " constant reference coordinates; using ",
      sum(keep),
      "."
    )
  }

  list(
    train = x_train[, keep, drop = FALSE],
    reference = x_reference[, keep, drop = FALSE],
    eval = x_eval[, keep, drop = FALSE]
  )
}


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

  filtered <- aim3_filter_reference_constants(
    x_train,
    x_reference,
    x_eval,
    representation
  )

  x_train <- filtered$train
  x_reference <- filtered$reference
  x_eval <- filtered$eval

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

  x <- as.matrix(case$x)

  keep <- aim3_variable_reference_columns(
    x[split$reference, , drop = FALSE]
  )

  if (!any(keep)) {
    stop(
      "The reference split has no variable coordinates.",
      call. = FALSE
    )
  }

  removed <- sum(!keep)

  if (removed > 0L) {
    message(
      "Input: removed ",
      removed,
      " constant reference coordinates; using ",
      sum(keep),
      "."
    )
  }

  x <- x[, keep, drop = FALSE]

  minimum_dim <- max(
    as.integer(settings$goix_subsampling$subset_dim),
    as.integer(case$intrinsic_dim)
  )

  if (ncol(x) < minimum_dim) {
    stop(
      "Too few variable coordinates remain for Aim 3.",
      call. = FALSE
    )
  }

  standardizer <- fit_standardizer(
    x[split$embedding, , drop = FALSE]
  )

  x_embedding <- apply_standardizer(
    x[split$embedding, , drop = FALSE],
    standardizer
  )
  x_train <- apply_standardizer(
    x[split$detector_train, , drop = FALSE],
    standardizer
  )
  x_reference <- apply_standardizer(
    x[split$reference, , drop = FALSE],
    standardizer
  )
  x_eval <- apply_standardizer(
    x[split$evaluation, , drop = FALSE],
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
