source("experiments/aim2/settings.R")
source("experiments/aim2/utils.R")
source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")
source("detectors/lof.R")
source("detectors/isolation_forest.R")

comparison_order <- function(difference, tolerance, smaller_is_better = FALSE) {
  if (!is.finite(difference) || abs(difference) <= tolerance) return(NA_integer_)
  direction <- sign(difference)
  if (smaller_is_better) direction <- -direction
  as.integer(direction)
}

pairwise_concordance <- function(
    detectors,
    aumvc_values,
    target_values,
    tolerance = 0,
    contributions = NULL,
    confidence_level = 0.95
) {
  pairs <- combn(seq_along(detectors), 2L)
  do.call(rbind, lapply(seq_len(ncol(pairs)), function(index) {
    first <- pairs[1L, index]
    second <- pairs[2L, index]
    aumvc_order <- comparison_order(
      aumvc_values[first] - aumvc_values[second],
      tolerance,
      TRUE
    )
    target_order <- comparison_order(
      target_values[first] - target_values[second],
      tolerance
    )
    uncertainty_order <- aumvc_order
    if (!is.null(contributions)) {
      paired <- contributions[[first]] - contributions[[second]]
      difference <- mean(paired)
      standard_error <- sd(paired) / sqrt(length(paired))
      critical <- qnorm(0.5 + confidence_level / 2)
      interval <- difference + c(-1, 1) * critical * standard_error
      uncertainty_order <- if (
        any(!is.finite(interval)) || interval[1L] <= 0 && interval[2L] >= 0
      ) {
        NA_integer_
      } else {
        comparison_order(difference, 0, TRUE)
      }
    }
    data.frame(
      detector_1 = detectors[first],
      detector_2 = detectors[second],
      aumvc_order = aumvc_order,
      aumvc_uncertainty_order = uncertainty_order,
      target_order = target_order,
      agree = if (is.na(aumvc_order) || is.na(target_order)) {
        NA
      } else {
        aumvc_order == target_order
      },
      uncertainty_agree = if (
        is.na(uncertainty_order) || is.na(target_order)
      ) {
        NA
      } else {
        uncertainty_order == target_order
      }
    )
  }))
}

consensus_concordance <- function(roc_pairs, pr_pairs) {
  if (
    !identical(roc_pairs$detector_1, pr_pairs$detector_1) ||
    !identical(roc_pairs$detector_2, pr_pairs$detector_2)
  ) {
    stop("ROC and PR pair tables do not match.", call. = FALSE)
  }

  consensus <- !is.na(roc_pairs$target_order) &
    !is.na(pr_pairs$target_order) &
    roc_pairs$target_order == pr_pairs$target_order

  compare <- function(aumvc_order) {
    result <- rep(NA, length(aumvc_order))
    keep <- consensus & !is.na(aumvc_order)
    result[keep] <- aumvc_order[keep] == roc_pairs$target_order[keep]
    result
  }

  data.frame(
    detector_1 = roc_pairs$detector_1,
    detector_2 = roc_pairs$detector_2,
    agree = compare(roc_pairs$aumvc_order),
    uncertainty_agree = compare(roc_pairs$aumvc_uncertainty_order)
  )
}
fit_models <- function(x_train, seed, settings) {
  models <- list(
    OCSVM = fit_ocsvm(
      x_train,
      nu = settings$detectors$ocsvm$nu,
      gamma = 1 / ncol(x_train)
    ),
    LOF = fit_lof(x_train, k = settings$detectors$lof$k),
    Isolation_Forest = fit_isolation_forest(
      x_train,
      ntrees = settings$detectors$iforest$ntrees,
      sample_size = min(
        settings$detectors$iforest$sample_size,
        nrow(x_train)
      ),
      seed = seed + 200L
    )
  )
  models
}

make_score_functions <- function(models) {
  list(
    OCSVM = function(x) score_ocsvm(models$OCSVM, x),
    LOF = function(x) score_lof(models$LOF, x),
    Isolation_Forest = function(x) score_isolation_forest(models$Isolation_Forest, x)
  )
}

run_once <- function(x, labels, counts, seed, settings) {
  split <- make_splits(
    labels,
    counts,
    c("detector_train", "reference", "aumvc", "label_eval"),
    seed
  )
  x_train <- x[split$detector_train, , drop = FALSE]
  x_reference <- x[split$reference, , drop = FALSE]
  x_aumvc <- x[split$aumvc, , drop = FALSE]
  x_label <- x[split$label_eval, , drop = FALSE]
  labels_label <- labels[split$label_eval]

  if (
    nrow(x_train) <= settings$detectors$lof$k ||
    nrow(x_reference) < 2L ||
    nrow(x_aumvc) < 2L ||
    length(unique(labels_label)) != 2L
  ) {
    stop("The Aim 2 split is not usable.", call. = FALSE)
  }

  train_sd <- apply(x_train, 2L, sd)
  keep <- is.finite(train_sd) & train_sd > 0
  if (!any(keep)) stop("No variable training features in Aim 2.", call. = FALSE)
  x_train <- x_train[, keep, drop = FALSE]
  x_reference <- x_reference[, keep, drop = FALSE]
  x_aumvc <- x_aumvc[, keep, drop = FALSE]
  x_label <- x_label[, keep, drop = FALSE]

  standardizer <- fit_standardizer(x_train)
  x_train <- apply_standardizer(x_train, standardizer)
  x_reference <- apply_standardizer(x_reference, standardizer)
  x_aumvc <- apply_standardizer(x_aumvc, standardizer)
  x_label <- apply_standardizer(x_label, standardizer)
  reference <- make_reference(
    x_reference,
    n_reference = settings$n_reference,
    n_mc_repetitions = settings$n_mc_repetitions,
    seed = seed + 100L,
    chunk_size = settings$reference_chunk_size
  )

  score_funs <- make_score_functions(
    fit_models(x_train, seed, settings)
  )
  contributions <- vector("list", length(score_funs))
  names(contributions) <- names(score_funs)
  results <- do.call(rbind, lapply(names(score_funs), function(detector) {
    score_fun <- score_funs[[detector]]
    mv <- aumvc(
      x_aumvc,
      reference,
      score_fun,
      score_direction = "anomaly",
      alpha_grid = settings$aumvc_alpha_grid
    )
    contributions[[detector]] <<- mv$mc_contribution
    label_scores <- score_fun(x_label)
    reliable <- mv$minimum_hit_count >= settings$minimum_hit_count &&
      is.finite(mv$aumvc_relative_mc_se) &&
      mv$aumvc_relative_mc_se <= settings$maximum_relative_mc_se
    data.frame(
      detector = detector,
      aumvc = mv$aumvc,
      aumvc_log = mv$aumvc_log,
      aumvc_normalized = mv$aumvc_normalized,
      aumvc_normalized_mc_se = mv$aumvc_normalized_mc_se,
      aumvc_relative_mc_se = mv$aumvc_relative_mc_se,
      minimum_hit_count = mv$minimum_hit_count,
      zero_occupancy = mv$zero_occupancy,
      low_occupancy = mv$low_occupancy,
      reliable = reliable,
      roc_auc = roc_auc_score(labels_label, label_scores),
      pr_auc = pr_auc_score(labels_label, label_scores)
    )
  }))

  roc_pairs <- pairwise_concordance(
    results$detector,
    results$aumvc_normalized,
    results$roc_auc,
    contributions = contributions,
    confidence_level = settings$confidence_level
  )
  pr_pairs <- pairwise_concordance(
    results$detector,
    results$aumvc_normalized,
    results$pr_auc,
    contributions = contributions,
    confidence_level = settings$confidence_level
  )
  list(
    results = results,
    consensus_pairs = consensus_concordance(roc_pairs, pr_pairs)
  )
}

summarize_detectors <- function(run_results) {
  combined <- do.call(rbind, lapply(seq_along(run_results), function(run) {
    data.frame(run = run, run_results[[run]]$results)
  }))
  do.call(rbind, lapply(unique(combined$detector), function(detector) {
    values <- combined[combined$detector == detector, , drop = FALSE]
    data.frame(
      detector = detector,
      aumvc_mean = mean(values$aumvc),
      aumvc_sd = sd(values$aumvc),
      aumvc_normalized_mean = mean(values$aumvc_normalized),
      aumvc_normalized_sd = sd(values$aumvc_normalized),
      aumvc_normalized_mc_se_mean = mean(values$aumvc_normalized_mc_se),
      aumvc_relative_mc_se_mean = mean(values$aumvc_relative_mc_se),
      minimum_hit_count = min(values$minimum_hit_count),
      reliable_runs = sum(values$reliable),
      roc_mean = mean(values$roc_auc),
      roc_sd = sd(values$roc_auc),
      pr_mean = mean(values$pr_auc),
      pr_sd = sd(values$pr_auc)
    )
  }))
}

summarize_concordance <- function(run_results) {
  pairs <- do.call(
    rbind,
    lapply(run_results, function(result) result$consensus_pairs)
  )
  total_pairs <- nrow(pairs)

  summarize <- function(column, comparison) {
    compared_pairs <- sum(!is.na(pairs[[column]]))
    matching_pairs <- sum(pairs[[column]], na.rm = TRUE)
    data.frame(
      comparison = comparison,
      metric = "ROC/PR consensus",
      matching_pairs = matching_pairs,
      compared_pairs = compared_pairs,
      total_pairs = total_pairs,
      percentage = if (compared_pairs) {
        100 * matching_pairs / compared_pairs
      } else {
        NA_real_
      }
    )
  }

  rbind(
    summarize("agree", "strict"),
    summarize("uncertainty_agree", "uncertainty_aware")
  )
}

total_concordance <- function(outputs) {
  combined <- do.call(rbind, lapply(outputs, function(output) output$concordance))
  groups <- split(seq_len(nrow(combined)), combined$comparison)

  do.call(rbind, lapply(groups, function(indices) {
    part <- combined[indices, , drop = FALSE]
    matching_pairs <- sum(part$matching_pairs)
    compared_pairs <- sum(part$compared_pairs)
    total_pairs <- sum(part$total_pairs)
    data.frame(
      comparison = part$comparison[1L],
      metric = "ROC/PR consensus",
      matching_pairs = matching_pairs,
      compared_pairs = compared_pairs,
      total_pairs = total_pairs,
      percentage = if (compared_pairs) {
        100 * matching_pairs / compared_pairs
      } else {
        NA_real_
      }
    )
  }))
}
run_dataset <- function(x, labels, dataset, settings) {
  x <- validate_matrix(x, "x")
  labels <- as.integer(labels)
  if (
    length(labels) != nrow(x) || anyNA(labels) ||
    !all(labels %in% c(0L, 1L))
  ) {
    stop("labels must be binary 0/1.", call. = FALSE)
  }

  run_results <- lapply(seq_len(settings$n_runs), function(run) {
    seed <- experiment_run_seed(settings, run)
    pool <- limit_contamination(
      x,
      labels,
      settings$max_anomaly_fraction,
      seed + 10L
    )
    counts <- split_counts(nrow(pool$x), settings)
    run_once(
      pool$x,
      pool$labels,
      counts,
      seed,
      settings
    )
  })
  output <- list(
    dataset = dataset,
    n_runs = settings$n_runs,
    detector_summary = summarize_detectors(run_results),
    concordance = summarize_concordance(run_results)
  )

  display <- data.frame(
    detector = output$detector_summary$detector,
    AUMVC = mapply(
      format_mean_sd,
      output$detector_summary$aumvc_mean,
      output$detector_summary$aumvc_sd
    ),
    ROC_AUC = mapply(
      format_mean_sd,
      output$detector_summary$roc_mean,
      output$detector_summary$roc_sd
    ),
    PR_AUC = mapply(
      format_mean_sd,
      output$detector_summary$pr_mean,
      output$detector_summary$pr_sd
    )
  )
  cat(dataset, " - ", settings$n_runs, " runs\n\n", sep = "")
  print(display, row.names = FALSE)
  cat("\nROC/PR consensus comparison\n")
  print(output$concordance, row.names = FALSE)
  invisible(output)
}

comparison_rows <- function(dataset, comparison) {
  data.frame(
    dataset = rep(dataset, nrow(comparison)),
    detector = "",
    AUMVC = "",
    AUMVC_normalized = "",
    AUMVC_normalized_MC_SE = "",
    AUMVC_relative_MC_SE = "",
    minimum_hit_count = NA_real_,
    reliable_runs = NA_real_,
    ROC_AUC = "",
    PR_AUC = "",
    comparison = comparison$comparison,
    metric = comparison$metric,
    matching_pairs = comparison$matching_pairs,
    compared_pairs = comparison$compared_pairs,
    total_pairs = comparison$total_pairs,
    percentage = comparison$percentage
  )
}

report_rows <- function(output) {
  summary <- output$detector_summary
  dataset <- paste0(output$dataset, " - ", output$n_runs, " runs")
  detector_rows <- data.frame(
    dataset = rep(dataset, nrow(summary)),
    detector = summary$detector,
    AUMVC = mapply(format_mean_sd, summary$aumvc_mean, summary$aumvc_sd),
    AUMVC_normalized = mapply(
      format_mean_sd,
      summary$aumvc_normalized_mean,
      summary$aumvc_normalized_sd
    ),
    AUMVC_normalized_MC_SE = vapply(
      summary$aumvc_normalized_mc_se_mean,
      format_number,
      character(1)
    ),
    AUMVC_relative_MC_SE = vapply(
      summary$aumvc_relative_mc_se_mean,
      format_number,
      character(1)
    ),
    minimum_hit_count = summary$minimum_hit_count,
    reliable_runs = summary$reliable_runs,
    ROC_AUC = mapply(format_mean_sd, summary$roc_mean, summary$roc_sd),
    PR_AUC = mapply(format_mean_sd, summary$pr_mean, summary$pr_sd),
    comparison = "",
    metric = "",
    matching_pairs = NA_real_,
    compared_pairs = NA_real_,
    total_pairs = NA_real_,
    percentage = NA_real_
  )
  rbind(
    detector_rows,
    comparison_rows(dataset, output$concordance)
  )
}
