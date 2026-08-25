aim3_make_goix_subsets <- function(
    ambient_dim,
    seed,
    settings
) {
  subset_dim <- as.integer(settings$subset_dim)
  n_subsets <- as.integer(settings$n_subsets)

  if (subset_dim > ambient_dim) {
    stop("Goix subset dimension exceeds ambient dimension.", call. = FALSE)
  }

  set.seed(seed)

  replicate(
    n_subsets,
    sample.int(
      ambient_dim,
      subset_dim,
      replace = FALSE
    ),
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
  subsets <- aim3_make_goix_subsets(
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
      seed = seed + 1000L + i,
      settings = settings$reference
    )

    score_fun <- function(x) score_ocsvm(model, x)

    mv <- aumvc(
      x_eval = eval_subset,
      reference = reference,
      score_fun = score_fun,
      score_direction = "anomaly",
      alpha_grid = settings$aumvc_alpha_grid
    )

    results[[i]] <- data.frame(
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
      log_reference_density = reference$log_density
    )
  }

  elapsed <- proc.time()[[3L]] - start
  results <- do.call(rbind, results)

  raw_aumvc <- if (all(is.finite(results$aumvc))) {
    mean(results$aumvc)
  } else {
    Inf
  }

  raw_mc_se <- if (all(is.finite(results$aumvc_mc_se))) {
    sqrt(sum(results$aumvc_mc_se^2)) / nrow(results)
  } else {
    NA_real_
  }

  normalized_mc_se <- sqrt(
    sum(results$aumvc_normalized_mc_se^2)
  ) / nrow(results)

  data.frame(
    representation = "goix_subsampling",
    aumvc = raw_aumvc,
    aumvc_normalized = mean(results$aumvc_normalized),
    aumvc_mc_se = raw_mc_se,
    aumvc_normalized_mc_se = normalized_mc_se,
    mean_occupancy = mean(results$mean_occupancy),
    zero_occupancy_fraction = mean(results$zero_occupancy_fraction),
    log_box_volume = mean(results$log_box_volume),
    n_reference = mean(results$n_reference),
    reference_capped = any(results$reference_capped),
    log_reference_density = mean(results$log_reference_density),
    runtime_seconds = elapsed
  )
}
