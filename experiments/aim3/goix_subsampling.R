aim3_goix_subsets <- function(ambient_dim, seed, settings) {
  subset_dim <- as.integer(settings$subset_dim)
  n_subsets <- as.integer(settings$n_subsets)

  if (subset_dim > ambient_dim) {
    stop("Goix subset dimension exceeds ambient dimension.", call. = FALSE)
  }

  set.seed(seed)
  replicate(
    n_subsets,
    sample.int(ambient_dim, subset_dim, replace = FALSE),
    simplify = FALSE
  )
}

aim3_evaluate_goix_subsampling <- function(
    x_train,
    x_reference,
    x_eval,
    seed,
    settings
) {
  subsets <- aim3_goix_subsets(
    ncol(x_train),
    seed,
    settings$goix_subsampling
  )

  results <- vector("list", length(subsets))
  start <- proc.time()[[3L]]

  for (i in seq_along(subsets)) {
    columns <- subsets[[i]]
    train_subset <- x_train[, columns, drop = FALSE]
    reference_subset <- x_reference[, columns, drop = FALSE]
    eval_subset <- x_eval[, columns, drop = FALSE]

    model <- fit_ocsvm(
      train_subset,
      nu = settings$detector$nu,
      gamma = 1 / ncol(train_subset)
    )

    reference <- aim3_make_reference(
      reference_subset,
      seed + 1000L + i,
      settings
    )

    mv <- aumvc(
      eval_subset,
      reference,
      function(x) score_ocsvm(model, x),
      score_direction = "anomaly",
      alpha_grid = settings$aumvc_alpha_grid
    )

    results[[i]] <- data.frame(
      aumvc = mv$aumvc,
      aumvc_normalized = mv$aumvc_normalized,
      aumvc_mc_se = mv$aumvc_mc_se,
      aumvc_normalized_mc_se = mv$aumvc_normalized_mc_se,
      zero_occupancy = mean(mv$mv_curve$volume_normalized == 0)
    )
  }

  results <- do.call(rbind, results)

  raw_mc_se <- if (all(is.finite(results$aumvc_mc_se))) {
    sqrt(sum(results$aumvc_mc_se^2)) / nrow(results)
  } else {
    NA_real_
  }

  data.frame(
    representation = "goix_subsampling",
    aumvc = if (all(is.finite(results$aumvc))) mean(results$aumvc) else Inf,
    aumvc_normalized = mean(results$aumvc_normalized),
    aumvc_mc_se = raw_mc_se,
    aumvc_normalized_mc_se = sqrt(
      sum(results$aumvc_normalized_mc_se^2)
    ) / nrow(results),
    zero_occupancy = mean(results$zero_occupancy),
    runtime_seconds = proc.time()[[3L]] - start
  )
}
