source("experiments/aim2/settings.R")
source("experiments/aim2/utils.R")
source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")
source("detectors/lof.R")
source("detectors/isolation_forest.R")

aim2_order <- function(difference, tolerance, smaller_is_better = FALSE) {
  if (!is.finite(difference) || abs(difference) <= tolerance) return(NA_integer_)
  direction <- sign(difference)
  if (smaller_is_better) direction <- -direction
  as.integer(direction)
}

aim2_pairwise_concordance <- function(
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
    aumvc_order <- aim2_order(
      aumvc_values[first] - aumvc_values[second],
      tolerance,
      TRUE
    )
    target_order <- aim2_order(
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
        aim2_order(difference, 0, TRUE)
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

aim2_consensus_concordance <- function(roc_pairs, pr_pairs) {
  if (
    !identical(roc_pairs$detector_1, pr_pairs$detector_1) ||
    !identical(roc_pairs$detector_2, pr_pairs$detector_2)
  ) {
    stop("ROC and PR pair tables do not match.", call. = FALSE)
  }

  agree <- mapply(function(aumvc_order, roc_order, pr_order) {
    if (
      is.na(aumvc_order) || is.na(roc_order) || is.na(pr_order) ||
      roc_order != pr_order
    ) {
      return(NA)
    }
    aumvc_order == roc_order
  }, roc_pairs$aumvc_order, roc_pairs$target_order, pr_pairs$target_order)
  uncertainty_agree <- mapply(function(aumvc_order, roc_order, pr_order) {
    if (
      is.na(aumvc_order) || is.na(roc_order) || is.na(pr_order) ||
      roc_order != pr_order
    ) {
      return(NA)
    }
    aumvc_order == roc_order
  }, roc_pairs$aumvc_uncertainty_order,
  roc_pairs$target_order,
  pr_pairs$target_order)

  data.frame(
    detector_1 = roc_pairs$detector_1,
    detector_2 = roc_pairs$detector_2,
    agree = agree,
    uncertainty_agree = uncertainty_agree
  )
}

aim2_fit_models <- function(x_train, seed, settings) {
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
  if (!models$OCSVM$converged) {
    stop("OCSVM solver did not converge in Aim 2.", call. = FALSE)
  }
  models
}

aim2_score_functions <- function(models) {
  list(
    OCSVM = function(x) score_ocsvm(models$OCSVM, x),
    LOF = function(x) score_lof(models$LOF, x),
    Isolation_Forest = function(x) score_isolation_forest(models$Isolation_Forest, x)
  )
}

aim2_run_once <- function(x, labels, counts, seed, settings) {
  split <- aim2_make_splits(
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

  score_functions <- aim2_score_functions(
    aim2_fit_models(x_train, seed, settings)
  )
  contributions <- vector("list", length(score_functions))
  names(contributions) <- names(score_functions)
  results <- do.call(rbind, lapply(names(score_functions), function(detector) {
    score_fun <- score_functions[[detector]]
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

  roc_pairs <- aim2_pairwise_concordance(
    results$detector,
    results$aumvc_normalized,
    results$roc_auc,
    contributions = contributions,
    confidence_level = settings$confidence_level
  )
  pr_pairs <- aim2_pairwise_concordance(
    results$detector,
    results$aumvc_normalized,
    results$pr_auc,
    contributions = contributions,
    confidence_level = settings$confidence_level
  )
  list(
    results = results,
    roc_pairs = roc_pairs,
    pr_pairs = pr_pairs,
    consensus_pairs = aim2_consensus_concordance(roc_pairs, pr_pairs)
  )
}

aim2_summarize_detectors <- function(run_results) {
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

aim2_summarize_concordance <- function(run_results) {
  metrics <- c(
    "AUMVC vs ROC-AUC",
    "AUMVC vs PR-AUC",
    "AUMVC vs ROC/PR consensus"
  )
  pair_names <- c("roc_pairs", "pr_pairs", "consensus_pairs")
  do.call(rbind, lapply(seq_along(metrics), function(index) {
    pairs <- do.call(
      rbind,
      lapply(run_results, function(result) result[[pair_names[index]]])
    )
    strict_matches <- sum(pairs$agree, na.rm = TRUE)
    strict_compared <- sum(!is.na(pairs$agree))
    uncertainty_matches <- sum(pairs$uncertainty_agree, na.rm = TRUE)
    uncertainty_compared <- sum(!is.na(pairs$uncertainty_agree))
    rbind(
      data.frame(
        comparison = "strict",
        metric = metrics[index],
        matches = strict_matches,
        compared = strict_compared,
        percentage = if (strict_compared > 0L) {
          100 * strict_matches / strict_compared
        } else {
          NA_real_
        }
      ),
      data.frame(
        comparison = "uncertainty_aware",
        metric = metrics[index],
        matches = uncertainty_matches,
        compared = uncertainty_compared,
        percentage = if (uncertainty_compared > 0L) {
          100 * uncertainty_matches / uncertainty_compared
        } else {
          NA_real_
        }
      )
    )
  }))
}

aim2_dataset_concordance <- function(run_results) {
  combined <- do.call(rbind, lapply(run_results, function(result) result$results))
  detectors <- unique(combined$detector)
  means <- do.call(rbind, lapply(detectors, function(detector) {
    values <- combined[combined$detector == detector, , drop = FALSE]
    data.frame(
      detector = detector,
      aumvc = mean(values$aumvc_normalized),
      roc = mean(values$roc_auc),
      pr = mean(values$pr_auc)
    )
  }))
  roc_pairs <- aim2_pairwise_concordance(
    means$detector,
    means$aumvc,
    means$roc
  )
  pr_pairs <- aim2_pairwise_concordance(
    means$detector,
    means$aumvc,
    means$pr
  )
  consensus <- aim2_consensus_concordance(roc_pairs, pr_pairs)
  tables <- list(roc_pairs, pr_pairs, consensus)
  metrics <- c(
    "AUMVC vs ROC-AUC",
    "AUMVC vs PR-AUC",
    "AUMVC vs ROC/PR consensus"
  )
  do.call(rbind, lapply(seq_along(tables), function(index) {
    agree <- tables[[index]]$agree
    matches <- sum(agree, na.rm = TRUE)
    compared <- sum(!is.na(agree))
    data.frame(
      comparison = "dataset_mean",
      metric = metrics[index],
      matches = matches,
      compared = compared,
      percentage = if (compared > 0L) 100 * matches / compared else NA_real_
    )
  }))
}

aim2_total_concordance <- function(outputs) {
  combined <- do.call(rbind, lapply(outputs, function(output) {
    rbind(output$concordance, output$dataset_concordance)
  }))
  groups <- unique(combined[, c("comparison", "metric"), drop = FALSE])
  do.call(rbind, lapply(seq_len(nrow(groups)), function(index) {
    keep <- combined$comparison == groups$comparison[index] &
      combined$metric == groups$metric[index]
    matches <- sum(combined$matches[keep])
    compared <- sum(combined$compared[keep])
    data.frame(
      comparison = groups$comparison[index],
      metric = groups$metric[index],
      matches = matches,
      compared = compared,
      percentage = if (compared > 0L) 100 * matches / compared else NA_real_
    )
  }))
}

aim2_run_dataset <- function(x, labels, dataset, settings) {
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
    pool <- aim2_limit_contamination(
      x,
      labels,
      settings$max_anomaly_fraction,
      seed + 10L
    )
    counts <- aim2_split_counts(nrow(pool$x), settings)
    aim2_run_once(
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
    detector_summary = aim2_summarize_detectors(run_results),
    concordance = aim2_summarize_concordance(run_results),
    dataset_concordance = aim2_dataset_concordance(run_results)
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
  cat("\nComparisons across all runs\n")
  print(output$concordance, row.names = FALSE)
  cat("\nComparison of detector means\n")
  print(output$dataset_concordance, row.names = FALSE)
  invisible(output)
}

aim2_comparison_rows <- function(dataset, comparison) {
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
    matches = comparison$matches,
    compared = comparison$compared,
    percentage = comparison$percentage
  )
}

aim2_report_rows <- function(output) {
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
      aim2_format_number,
      character(1)
    ),
    AUMVC_relative_MC_SE = vapply(
      summary$aumvc_relative_mc_se_mean,
      aim2_format_number,
      character(1)
    ),
    minimum_hit_count = summary$minimum_hit_count,
    reliable_runs = summary$reliable_runs,
    ROC_AUC = mapply(format_mean_sd, summary$roc_mean, summary$roc_sd),
    PR_AUC = mapply(format_mean_sd, summary$pr_mean, summary$pr_sd),
    comparison = "",
    metric = "",
    matches = NA_real_,
    compared = NA_real_,
    percentage = NA_real_
  )
  rbind(
    detector_rows,
    aim2_comparison_rows(dataset, output$concordance),
    aim2_comparison_rows(dataset, output$dataset_concordance)
  )
}
