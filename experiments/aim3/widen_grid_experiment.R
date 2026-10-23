source("experiments/aim3/synthetic_experiment.R")
source("experiments/aim3/real_experiment.R")

original_write_report <- write_report
original_prepare_candidate_scores <- prepare_candidate_scores
original_evaluate_candidate <- evaluate_candidate
original_evaluate_goix <- evaluate_goix

score_diagnostic_columns <- c(
  "label_score_sd", "label_score_range", "label_score_distinct",
  "label_near_constant", "eval_near_constant",
  "reference_near_constant", "goix_constant_subsets"
)
synthetic_result_columns <- c(synthetic_result_columns, score_diagnostic_columns)
real_result_columns <- c(real_result_columns, score_diagnostic_columns)

score_is_constant <- function(scores) {
  if (!length(scores) || any(!is.finite(scores))) return(NA)
  diff(range(scores)) <= 1e-10 * max(1, max(abs(scores)))
}

score_diagnostics <- function(label, evaluation, reference = NULL) {
  c(
    label_score_sd = stats::sd(label),
    label_score_range = diff(range(label)),
    label_score_distinct = length(unique(signif(label, 12L))),
    label_near_constant = as.numeric(score_is_constant(label)),
    eval_near_constant = as.numeric(score_is_constant(evaluation)),
    reference_near_constant = if (is.null(reference)) NA_real_ else
      as.numeric(score_is_constant(reference)),
    goix_constant_subsets = NA_real_
  )
}

prepare_candidate_scores <- function(x_train, x_eval, x_label, reference, scales, settings) {
  scored <- original_prepare_candidate_scores(
    x_train, x_eval, x_label, reference, scales, settings
  )
  extremes <- range(scales)
  for (i in seq_along(scales)) {
    if (!scales[i] %in% extremes) next
    candidate <- scored$candidates[[i]]
    scored$candidates[[i]]$score_diagnostics <- score_diagnostics(
      candidate$label_scores, candidate$evaluation_scores,
      candidate$reference_scores
    )
  }
  scored
}

evaluate_candidate <- function(
    representation, x_train, x_eval, x_label, labels_label, types_label,
    reference, gamma_scale, settings, embedding_oos_stress = NA_real_, candidate
) {
  row <- original_evaluate_candidate(
    representation, x_train, x_eval, x_label, labels_label, types_label,
    reference, gamma_scale, settings, embedding_oos_stress, candidate
  )
  checks <- candidate$score_diagnostics
  for (name in score_diagnostic_columns) {
    row[[name]] <- if (is.null(checks)) NA_real_ else unname(checks[name])
  }
  row
}

evaluate_goix <- function(
    x_train, x_reference, x_eval, x_label, labels_label, types_label,
    seed, settings
) {
  result <- original_evaluate_goix(
    x_train, x_reference, x_eval, x_label, labels_label, types_label,
    seed, settings
  )
  for (name in score_diagnostic_columns) result[[name]] <- NA_real_
  subsets <- goix_subsets(ncol(x_train), seed, settings$goix_subsampling)
  extremes <- range(gamma_scales(settings))

  for (scale in extremes) {
    ensemble <- numeric(nrow(x_label))
    ensemble_eval <- numeric(nrow(x_eval))
    constant_subsets <- 0L
    for (columns in subsets) {
      train <- x_train[, columns, drop = FALSE]
      model <- fit_ocsvm(train, nu = settings$detector$nu,
                         gamma = scale / ncol(train))
      evaluation <- score_ocsvm(model, x_eval[, columns, drop = FALSE])
      label <- score_ocsvm(model, x_label[, columns, drop = FALSE])
      constant_subsets <- constant_subsets + as.integer(score_is_constant(label))
      ensemble <- ensemble + percentile_scores(evaluation, label) / length(subsets)
      ensemble_eval <- ensemble_eval + percentile_scores(evaluation, evaluation) /
        length(subsets)
    }
    diagnostics <- score_diagnostics(ensemble, ensemble_eval)
    diagnostics["goix_constant_subsets"] <- constant_subsets
    i <- which(result$gamma_scale == scale)
    for (name in score_diagnostic_columns) result[[name]][i] <- diagnostics[name]
  }
  result
}

result_dirs <- function(type) {
  root <- file.path("experiments/aim3/widen_result", type, "widen_experiment")
  paths <- list(experiment = file.path(root, "experiment"),
                analysis = file.path(root, "analysis"))
  for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  paths
}

write_report <- function(data, path, digits = 6L) {
  filename <- basename(path)
  if (filename %in% c("syn_bandwidth_exp.csv", "real_bandwidth_exp.csv")) {
    return(invisible(path))
  }
  if (filename %in% c("syn_results.csv", "real_results.csv")) {
    columns <- c(
      "dataset", "run", "truth", "intrinsic_dim", "snr", "original_dim",
      "representation", "gamma_scale", "aumvc_normalized", "aumvc_reliable",
      "zero_occupancy", "min_hit_count", "roc_auc", "pr_auc"
    )
    data <- data[, intersect(columns, names(data)), drop = FALSE]
  }
  original_write_report(data, path, digits)
}

widen_selection <- function(x, groups, tolerance) {
  x <- x[order(x$gamma_scale), , drop = FALSE]
  valid <- which(!is.na(x$aumvc_reliable) & x$aumvc_reliable &
                 is.finite(x$aumvc_normalized))
  raw <- which(is.finite(x$aumvc_normalized))
  raw_chosen <- select_gamma(x, rep(TRUE, nrow(x)), tolerance)
  raw_status <- if (!length(raw)) "unreliable" else if (is.na(raw_chosen)) "tie" else "selected"

  status <- "unreliable"
  chosen <- NA_integer_
  if (length(valid) >= 2L) {
    values <- x$aumvc_normalized[valid]
    minimum <- min(values)
    tied <- valid[abs(values - minimum) <=
                  tolerance * pmax(1, abs(values), abs(minimum))]
    if (length(tied) > 1L) {
      status <- "tie"
    } else {
      chosen <- tied[1L]
      status <- if (length(valid) == nrow(x)) "selected" else "partial"
    }
  }

  pairs <- if (length(valid) >= 2L) utils::combn(valid, 2L) else
    matrix(integer(), nrow = 2L, ncol = 0L)
  matches <- logical(0)
  for (j in seq_len(ncol(pairs))) {
    a <- pairs[1L, j]
    b <- pairs[2L, j]
    roc_order <- sign(x$roc_auc[a] - x$roc_auc[b])
    pr_order <- sign(x$pr_auc[a] - x$pr_auc[b])
    if (!is.finite(roc_order) || !is.finite(pr_order) ||
        roc_order == 0 || roc_order != pr_order) next
    gap <- x$aumvc_normalized[a] - x$aumvc_normalized[b]
    matches <- c(matches, is.finite(gap) &&
                 abs(gap) > tolerance && sign(-gap) == roc_order)
  }

  best_roc <- max(x$roc_auc, na.rm = TRUE)
  roc_winners <- which(is.finite(x$roc_auc) &
                       abs(x$roc_auc - best_roc) <= 1e-12)
  roc_gamma <- if (length(roc_winners) == 1L) x$gamma_scale[roc_winners] else NA_real_
  row <- x[1L, c(groups, "run"), drop = FALSE]
  row$raw_status <- raw_status
  row$raw_chosen_gamma <- raw_chosen
  row$status <- status
  row$chosen_gamma <- if (is.na(chosen)) NA_real_ else x$gamma_scale[chosen]
  row$roc_best_gamma <- roc_gamma
  row$eligible_pairs <- length(matches)
  row$pair_agreement <- if (length(matches)) mean(matches) else NA_real_
  row$roc_regret <- if (is.na(chosen)) NA_real_ else best_roc - x$roc_auc[chosen]
  row$pr_regret <- if (is.na(chosen)) NA_real_ else
    max(x$pr_auc, na.rm = TRUE) - x$pr_auc[chosen]
  row$roc_range <- diff(range(x$roc_auc, na.rm = TRUE))
  row
}

analyze_results <- function(results, bandwidth, type, directory, digits = 6L) {
  groups <- if (type == "real") c("dataset", "representation") else
    c("truth", "intrinsic_dim", "original_dim", "snr", "representation")
  cases <- group_rows(results, c(groups, "run"))
  selection <- do.call(rbind, lapply(cases, function(index) {
    widen_selection(results[index, , drop = FALSE], groups,
                    sqrt(.Machine$double.eps))
  }))
  rownames(selection) <- NULL
  original_write_report(
    selection[, names(selection) != "roc_range", drop = FALSE],
    file.path(directory, "selection.csv"), digits
  )

  summaries <- lapply(group_rows(selection, groups), function(index) {
    x <- selection[index, , drop = FALSE]
    selected <- x$status == "selected"
    row <- x[1L, groups, drop = FALSE]
    row$choice_rate <- mean(selected)
    row$selected_gamma_scale_mode <- mode_unique(x$chosen_gamma[selected])
    row$pair_agreement <- mean_finite(x$pair_agreement)
    row$roc_regret <- mean_finite(x$roc_regret[selected])
    row$pr_regret <- mean_finite(x$pr_regret[selected])
    row
  })
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  original_write_report(summary, file.path(directory, "summary.csv"), digits)

  diagnostic_groups <- if (type == "real") c("dataset", "representation") else
    c("original_dim", "representation")
  boundaries <- range(results$gamma_scale)
  diagnostics <- lapply(group_rows(selection, diagnostic_groups), function(index) {
    x <- selection[index, , drop = FALSE]
    row <- x[1L, diagnostic_groups, drop = FALSE]
    row$n_cases <- nrow(x)
    row$median_roc_range <- stats::median(x$roc_range, na.rm = TRUE)
    row$roc_unique_optima <- sum(is.finite(x$roc_best_gamma))
    row$roc_boundary_optima <- sum(x$roc_best_gamma %in% boundaries, na.rm = TRUE)
    row$raw_aumvc_unique_optima <- sum(x$raw_status == "selected")
    row$raw_aumvc_boundary_optima <- sum(x$raw_chosen_gamma %in% boundaries, na.rm = TRUE)
    row$reliable_aumvc_boundary_optima <- sum(
      x$chosen_gamma[x$status == "selected"] %in% boundaries, na.rm = TRUE
    )
    row$raw_abstentions <- sum(x$raw_status != "selected")
    row$unreliable <- sum(x$status == "unreliable")
    row$partial <- sum(x$status == "partial")
    row$ties <- sum(x$status == "tie")
    row
  })
  diagnostics <- do.call(rbind, diagnostics)
  rownames(diagnostics) <- NULL
  original_write_report(diagnostics, file.path(directory, "grid_diagnostics.csv"), digits)

  extremes <- results$gamma_scale %in% boundaries
  checks <- results[extremes, c(groups, "run", "gamma_scale", score_diagnostic_columns),
                    drop = FALSE]
  utils::write.csv(checks, file.path(directory, "score_checks.csv"),
                   row.names = FALSE, na = "")
  invisible(list(selection = selection, summary = summary,
                 diagnostics = diagnostics, checks = checks))
}

settings <- AIM3_SETTINGS
settings$bandwidth$gamma_scales <- c(2^-6, 2^-2, 1, 2^2, 2^6)
run_synthetic_experiment(settings)
run_real_experiment(settings)
