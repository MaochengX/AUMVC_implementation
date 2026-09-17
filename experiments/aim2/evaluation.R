source("experiments/aim2/settings.R")
source("experiments/utils.R")
source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")
source("detectors/lof.R")
source("detectors/isolation_forest.R")

aim2_compare_pair <- function(aumvc_difference, target_difference, tolerance) {
  if (abs(target_difference) <= tolerance) return(NA)
  if (abs(aumvc_difference) <= tolerance) return(FALSE)
  sign(aumvc_difference) == -sign(target_difference)
}

aim2_pairwise_concordance <- function(detectors, aumvc_values, target_values, tolerance) {
  pairs <- combn(seq_along(detectors), 2L)

  do.call(rbind, lapply(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1L, k]
    j <- pairs[2L, k]
    data.frame(
      detector_1 = detectors[i],
      detector_2 = detectors[j],
      agree = aim2_compare_pair(
        aumvc_values[i] - aumvc_values[j],
        target_values[i] - target_values[j],
        tolerance
      )
    )
  }))
}

aim2_either_concordance <- function(roc_pairs, pr_pairs) {
  if (
    !identical(roc_pairs$detector_1, pr_pairs$detector_1) ||
    !identical(roc_pairs$detector_2, pr_pairs$detector_2)
  ) {
    stop("ROC and PR pair tables do not match.", call. = FALSE)
  }

  agree <- mapply(function(roc, pr) {
    values <- c(roc, pr)
    values <- values[!is.na(values)]
    if (length(values) == 0L) NA else any(values)
  }, roc_pairs$agree, pr_pairs$agree)

  data.frame(
    detector_1 = roc_pairs$detector_1,
    detector_2 = roc_pairs$detector_2,
    agree = agree
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
      sample_size = min(settings$detectors$iforest$sample_size, nrow(x_train)),
      seed = seed + 200L
    )
  )

  if (!models$OCSVM$converged) warning("OCSVM solver did not converge.", call. = FALSE)
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
  split <- make_splits(
    nrow(x),
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

  lower <- apply(x_reference, 2L, min)
  upper <- apply(x_reference, 2L, max)
  keep <- upper > lower
  if (!any(keep)) stop("No variable reference features in Aim 2.")
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
    seed = seed + 100L
  )

  score_functions <- aim2_score_functions(aim2_fit_models(x_train, seed, settings))
  results <- do.call(rbind, lapply(names(score_functions), function(detector) {
    score_fun <- score_functions[[detector]]
    label_scores <- score_fun(x_label)
    mv <- aumvc(
      x_aumvc,
      reference,
      score_fun,
      score_direction = "anomaly",
      alpha_grid = settings$aumvc_alpha_grid
    )
    data.frame(
      detector = detector,
      aumvc = mv$aumvc,
      aumvc_normalized = mv$aumvc_normalized,
      aumvc_mc_se = mv$aumvc_mc_se,
      roc_auc = roc_auc_score(labels_label, label_scores),
      pr_auc = pr_auc_score(labels_label, label_scores)
    )
  }))

  roc_pairs <- aim2_pairwise_concordance(
    results$detector,
    results$aumvc,
    results$roc_auc,
    settings$concordance_tolerance
  )
  pr_pairs <- aim2_pairwise_concordance(
    results$detector,
    results$aumvc,
    results$pr_auc,
    settings$concordance_tolerance
  )
  list(
    results = results,
    roc_pairs = roc_pairs,
    pr_pairs = pr_pairs,
    either_pairs = aim2_either_concordance(roc_pairs, pr_pairs)
  )
}

aim2_summarize_detectors <- function(run_results) {
  combined <- do.call(rbind, lapply(seq_along(run_results), function(run) {
    data.frame(run = run, run_results[[run]]$results)
  }))

  summary <- do.call(rbind, lapply(unique(combined$detector), function(detector) {
    x <- combined[combined$detector == detector, , drop = FALSE]
    data.frame(
      detector = detector,
      aumvc_mean = mean(x$aumvc),
      aumvc_sd = sd(x$aumvc),
      aumvc_normalized_mean = mean(x$aumvc_normalized),
      aumvc_normalized_sd = sd(x$aumvc_normalized),
      roc_mean = mean(x$roc_auc),
      roc_sd = sd(x$roc_auc),
      pr_mean = mean(x$pr_auc),
      pr_sd = sd(x$pr_auc)
    )
  }))

  list(combined = combined, summary = summary)
}

aim2_summarize_concordance <- function(run_results) {
  metrics <- c("AUMVC vs ROC-AUC", "AUMVC vs PR-AUC", "AUMVC vs either")
  pair_names <- c("roc_pairs", "pr_pairs", "either_pairs")

  do.call(rbind, lapply(seq_along(metrics), function(i) {
    pairs <- do.call(rbind, lapply(run_results, function(x) x[[pair_names[i]]]))
    matches <- sum(pairs$agree, na.rm = TRUE)
    compared <- sum(!is.na(pairs$agree))
    data.frame(
      metric = metrics[i],
      matches = matches,
      compared = compared,
      percentage = if (compared > 0L) 100 * matches / compared else NA_real_
    )
  }))
}

aim2_total_concordance <- function(outputs) {
  metrics <- outputs[[1L]]$concordance$metric
  do.call(rbind, lapply(seq_along(metrics), function(i) {
    matches <- sum(vapply(outputs, function(x) x$concordance$matches[i], numeric(1)))
    compared <- sum(vapply(outputs, function(x) x$concordance$compared[i], numeric(1)))
    data.frame(
      metric = metrics[i], matches = matches, compared = compared,
      percentage = if (compared > 0L) 100 * matches / compared else NA_real_
    )
  }))
}

aim2_run_dataset <- function(x, labels, dataset, counts, settings) {
  x <- validate_matrix(x, "x")
  labels <- as.integer(labels)

  if (length(labels) != nrow(x) || anyNA(labels) || !all(labels %in% c(0L, 1L))) {
    stop("labels must be binary 0/1.", call. = FALSE)
  }

  run_results <- lapply(seq_len(settings$n_runs), function(run) {
    aim2_run_once(x, labels, counts,
                  experiment_run_seed(settings, run), settings)
  })

  detector_runs <- aim2_summarize_detectors(run_results)
  concordance <- aim2_summarize_concordance(run_results)
  output <- list(
    protocol = "unsupervised_disjoint_mv",
    dataset = dataset,
    n_runs = settings$n_runs,
    base_seed = settings$seed,
    run_results = run_results,
    detector_results = detector_runs$combined,
    detector_summary = detector_runs$summary,
    concordance = concordance
  )

  display <- data.frame(
    detector = output$detector_summary$detector,
    AUMVC = format_mean_sd(output$detector_summary$aumvc_mean, output$detector_summary$aumvc_sd),
    ROC_AUC = format_mean_sd(output$detector_summary$roc_mean, output$detector_summary$roc_sd),
    PR_AUC = format_mean_sd(output$detector_summary$pr_mean, output$detector_summary$pr_sd)
  )

  cat(dataset, " - ", settings$n_runs, " runs\n\n", sep = "")
  print(display, row.names = FALSE)
  cat("\nComparisons across all runs\n")
  print(concordance, row.names = FALSE)
  invisible(output)
}

aim2_format_result <- function(mean, sd) {
  if (is.na(sd)) return(formatC(mean, digits = 4L, format = "f"))
  format_mean_sd(mean, sd)
}

aim2_comparison_rows <- function(dataset, comparison) {
  data.frame(
    dataset = rep(dataset, nrow(comparison)),
    detector = "",
    AUMVC = "",
    ROC_AUC = "",
    PR_AUC = "",
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
    AUMVC = mapply(aim2_format_result, summary$aumvc_mean, summary$aumvc_sd),
    ROC_AUC = mapply(aim2_format_result, summary$roc_mean, summary$roc_sd),
    PR_AUC = mapply(aim2_format_result, summary$pr_mean, summary$pr_sd),
    metric = "",
    matches = NA_integer_,
    compared = NA_integer_,
    percentage = NA_real_
  )
  rbind(detector_rows, aim2_comparison_rows(dataset, output$concordance))
}

aim2_save_report <- function(rows, name) {
  directory <- "experiments/aim2/result"
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, paste0(name, ".csv"))
  temporary <- paste0(path, ".tmp")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(rows, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Could not save Aim 2 CSV.")
  cat("Saved: ", path, "\n", sep = "")
}
