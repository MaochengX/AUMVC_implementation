aim3_goix_subsets <- function(ambient_dim, seed, settings) {
  subset_dim <- as.integer(settings$subset_dim)
  n_subsets <- as.integer(settings$n_subsets)
  if (
    length(subset_dim) != 1L || is.na(subset_dim) || subset_dim < 1L ||
    length(n_subsets) != 1L || is.na(n_subsets) || n_subsets < 1L ||
    subset_dim > ambient_dim
  ) {
    stop("Invalid Goix feature-subsampling settings.", call. = FALSE)
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
    x_label,
    labels_label,
    types_label,
    seed,
    settings
) {
  subsets <- aim3_goix_subsets(
    ncol(x_train),
    seed,
    settings$goix_subsampling
  )
  results <- vector("list", length(subsets))
  ensemble_scores <- matrix(0, nrow = nrow(x_label), ncol = length(subsets))
  contribution_count <- settings$n_reference * settings$n_mc_repetitions
  normalized_contribution <- numeric(contribution_count)
  raw_contribution <- numeric(contribution_count)
  start <- proc.time()[[3L]]

  for (index in seq_along(subsets)) {
    columns <- subsets[[index]]
    train_subset <- x_train[, columns, drop = FALSE]
    reference_subset <- x_reference[, columns, drop = FALSE]
    eval_subset <- x_eval[, columns, drop = FALSE]
    label_subset <- x_label[, columns, drop = FALSE]
    model <- fit_ocsvm(
      train_subset,
      nu = settings$detector$nu,
      gamma = 1 / ncol(train_subset)
    )
    if (!model$converged) {
      stop("OCSVM solver did not converge in Goix subsampling.", call. = FALSE)
    }

    reference <- aim3_make_reference(
      reference_subset,
      seed + 1000L,
      settings
    )
    score_fun <- function(x) score_ocsvm(model, x)
    mv <- aumvc(
      eval_subset,
      reference,
      score_fun,
      score_direction = "anomaly",
      alpha_grid = settings$aumvc_alpha_grid
    )
    normalized_contribution <- normalized_contribution +
      mv$mc_contribution / length(subsets)
    raw_contribution <- raw_contribution +
      mv$mc_contribution * reference$box$volume / length(subsets)
    eval_scores <- score_fun(eval_subset)
    label_scores <- score_fun(label_subset)
    ensemble_scores[, index] <- aim3_percentile_scores(
      eval_scores,
      label_scores
    )
    results[[index]] <- data.frame(
      aumvc = mv$aumvc,
      aumvc_log = mv$aumvc_log,
      aumvc_normalized = mv$aumvc_normalized,
      aumvc_mc_se = mv$aumvc_mc_se,
      aumvc_normalized_mc_se = mv$aumvc_normalized_mc_se,
      aumvc_relative_mc_se = mv$aumvc_relative_mc_se,
      minimum_hit_count = mv$minimum_hit_count,
      zero_occupancy = mv$zero_occupancy,
      low_occupancy = mv$low_occupancy,
      box_log_volume = reference$box$log_volume,
      ocsvm_kkt_gap = model$kkt_gap
    )
  }

  results <- do.call(rbind, results)
  label_scores <- rowMeans(ensemble_scores)
  type_metrics <- aim3_type_metrics(labels_label, label_scores, types_label)
  aumvc_log <- log_mean_exp(results$aumvc_log)
  raw_mc_se <- if (all(is.finite(raw_contribution))) {
    sd(raw_contribution) / sqrt(length(raw_contribution))
  } else {
    NA_real_
  }
  normalized_mc_se <- sd(normalized_contribution) /
    sqrt(length(normalized_contribution))

  data.frame(
    representation = "goix_subsampling",
    aumvc = exp_if_representable(aumvc_log),
    aumvc_normalized = mean(results$aumvc_normalized),
    aumvc_mc_se = raw_mc_se,
    aumvc_normalized_mc_se = normalized_mc_se,
    aumvc_relative_mc_se = if (mean(results$aumvc_normalized) > 0) {
      normalized_mc_se / mean(results$aumvc_normalized)
    } else {
      NA_real_
    },
    aumvc_subset_se = if (nrow(results) > 1L) {
      sd(results$aumvc_normalized) / sqrt(nrow(results))
    } else {
      NA_real_
    },
    minimum_hit_count = min(results$minimum_hit_count),
    zero_occupancy = mean(results$zero_occupancy),
    low_occupancy = mean(results$low_occupancy),
    reliable = min(results$minimum_hit_count) >= settings$minimum_hit_count &&
      all(is.finite(results$aumvc_relative_mc_se)) &&
      max(results$aumvc_relative_mc_se) <= settings$maximum_relative_mc_se,
    box_log_volume = mean(results$box_log_volume),
    roc_auc = roc_auc_score(labels_label, label_scores),
    pr_auc = pr_auc_score(labels_label, label_scores),
    roc_distributional = unname(type_metrics["roc_distributional"]),
    pr_distributional = unname(type_metrics["pr_distributional"]),
    roc_structural = unname(type_metrics["roc_structural"]),
    pr_structural = unname(type_metrics["pr_structural"]),
    ocsvm_converged = TRUE,
    ocsvm_kkt_gap = max(results$ocsvm_kkt_gap),
    embedding_oos_stress = NA_real_,
    runtime_seconds = proc.time()[[3L]] - start
  )
}
