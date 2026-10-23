source("experiments/aim3/synthetic_experiment.R")

bootstrap_score_grid <- function(prepared, seed, settings) {
  scales <- gamma_scales(settings)
  base <- prepared$prepared
  output <- list()
  methods <- c("ambient", "goix_subsampling", names(prepared$embeddings))

  for (method in methods) {
    if (method %in% c("ambient", "goix_subsampling")) {
      train <- base$train
      evaluation <- base$eval
      label <- base$label
    } else {
      projected <- prepared$embeddings[[method]]
      train <- projected$train
      evaluation <- projected$eval
      label <- projected$label
    }
    predictions <- matrix(0, nrow(label), length(scales))

    if (method == "goix_subsampling") {
      subsets <- goix_subsets(ncol(train), seed + 2000L, settings$goix_subsampling)
      for (columns in subsets) {
        x_train <- train[, columns, drop = FALSE]
        x_eval <- evaluation[, columns, drop = FALSE]
        x_label <- label[, columns, drop = FALSE]
        for (i in seq_along(scales)) {
          model <- fit_ocsvm(x_train, nu = settings$detector$nu,
                             gamma = scales[i] / ncol(x_train))
          cal <- score_ocsvm(model, x_eval)
          predictions[, i] <- predictions[, i] +
            percentile_scores(cal, score_ocsvm(model, x_label)) / length(subsets)
        }
      }
    } else {
      for (i in seq_along(scales)) {
        model <- fit_ocsvm(train, nu = settings$detector$nu,
                           gamma = scales[i] / ncol(train))
        predictions[, i] <- score_ocsvm(model, label)
      }
    }
    output[[method]] <- predictions
  }
  output
}

paired_bootstrap_case <- function(row, run, settings, selection, repetitions) {
  seed <- synthetic_seed(row, run, settings$synthetic, settings$seed)
  case <- generate_synthetic(row$intrinsic_dim, row$original_dim, row$snr,
                             row$truth, seed, settings$synthetic)
  prepared <- prepare_representations(case, seed + 50000L, settings)
  scores <- bootstrap_score_grid(prepared, seed + 50000L, settings)
  labels <- prepared$prepared$labels_label
  normal <- which(labels == 0L)
  anomaly <- which(labels == 1L)
  set.seed(seed + 9876L)
  samples <- replicate(repetitions, c(
    sample(normal, length(normal), replace = TRUE),
    sample(anomaly, length(anomaly), replace = TRUE)
  ), simplify = FALSE)
  scales <- gamma_scales(settings)
  rows <- list()

  for (method in names(scores)) {
    matrix <- scores[[method]]
    observed <- apply(matrix, 2L, function(x) roc_auc_score(labels, x))
    best_set <- which(abs(observed - max(observed)) <= 1e-12)
    best_index <- if (length(best_set) == 1L) best_set else NA_integer_
    bootstrap_roc <- vapply(samples, function(index) {
      apply(matrix[index, , drop = FALSE], 2L, function(x) {
        roc_auc_score(labels[index], x)
      })
    }, numeric(length(scales)))
    chosen <- selection[
      selection$run == run & selection$truth == row$truth &
        selection$intrinsic_dim == row$intrinsic_dim &
        selection$original_dim == row$original_dim &
        selection$snr == row$snr &
        selection$representation == method, , drop = FALSE
    ]
    gamma <- if (nrow(chosen) == 1L && chosen$status == "selected") {
      chosen$chosen_gamma
    } else NA_real_
    selected_index <- match(gamma, scales)
    gap <- c(NA_real_, NA_real_)
    observed_gap <- NA_real_
    regret <- c(NA_real_, NA_real_)
    mean_regret <- NA_real_
    if (!is.na(selected_index)) {
      bootstrap_regret <- apply(bootstrap_roc, 2L, max) -
        bootstrap_roc[selected_index, ]
      regret <- as.numeric(stats::quantile(
        bootstrap_regret, c(0.025, 0.975), names = FALSE
      ))
      mean_regret <- mean(bootstrap_regret)
      if (!is.na(best_index)) {
        paired_gap <- bootstrap_roc[best_index, ] -
          bootstrap_roc[selected_index, ]
        gap <- as.numeric(stats::quantile(
          paired_gap, c(0.025, 0.975), names = FALSE
        ))
        observed_gap <- observed[best_index] - observed[selected_index]
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      run = run, truth = row$truth, intrinsic_dim = row$intrinsic_dim,
      snr = row$snr, original_dim = row$original_dim,
      representation = method, normal_count = length(normal),
      anomaly_count = length(anomaly), selected_gamma = gamma,
      roc_best_gamma = if (is.na(best_index)) NA_real_ else scales[best_index],
      roc_best_tied = is.na(best_index), observed_roc_gap = observed_gap,
      bootstrap_gap_low = gap[1L], bootstrap_gap_high = gap[2L],
      bootstrap_regret_low = regret[1L],
      bootstrap_regret_high = regret[2L],
      bootstrap_mean_regret = mean_regret,
      roc_best_change_rate = mean(vapply(seq_len(ncol(bootstrap_roc)), function(j) {
        winners <- which(abs(bootstrap_roc[, j] - max(bootstrap_roc[, j])) <= 1e-12)
        !identical(winners, best_set)
      }, logical(1)))
    )
  }
  do.call(rbind, rows)
}

bootstrap_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

bootstrap_summary <- function(results) {
  rows <- lapply(split(results, results$representation), function(x) {
    data.frame(
      representation = x$representation[1L],
      cases = nrow(x),
      selected_cases = sum(is.finite(x$selected_gamma)),
      mean_observed_roc_gap = bootstrap_mean(x$observed_roc_gap),
      mean_bootstrap_regret = bootstrap_mean(x$bootstrap_mean_regret),
      mean_roc_best_change_rate = bootstrap_mean(x$roc_best_change_rate),
      intervals_including_zero = sum(
        x$bootstrap_gap_low <= 0 & x$bootstrap_gap_high >= 0,
        na.rm = TRUE
      ),
      paired_intervals_available = sum(is.finite(x$bootstrap_gap_low))
    )
  })
  do.call(rbind, rows)
}

save_bootstrap_pdf <- function(results, summary, path, repetitions) {
  labels <- sprintf(
    "%s q%d d%d SNR%s | %s",
    ifelse(results$truth == "polynomial_interaction", "Poly", "Osc"),
    results$intrinsic_dim, results$original_dim, results$snr,
    results$representation
  )
  colors <- c(ambient = "#BB5B43", goix_subsampling = "#15866D",
              mds = "#3F77AB", isomap = "#9568A5")
  grDevices::pdf(path, width = 11, height = 7.5)
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::plot.new()
  graphics::title("Paired bootstrap | widened gamma", cex.main = 1.45)
  graphics::text(0.03, 0.91,
                 sprintf("%d paired, class-stratified resamples per synthetic case", repetitions),
                 adj = 0, cex = 0.9)
  graphics::text(0.03, 0.86,
                 "OCSVM fits and label scores held fixed; all gamma values share each resampled index set.",
                 adj = 0, cex = 0.82)
  header <- sprintf("%-34s  %8s  %8s  %11s  %23s  %12s",
                    "Case / representation", "Anomalies", "Selected", "ROC best",
                    "Observed gap [95% CI]", "Best changes")
  graphics::text(0.02, 0.78, header, family = "mono", adj = 0, cex = 0.72)
  for (i in seq_len(nrow(results))) {
    z <- results[i, ]
    gap <- if (is.finite(z$observed_roc_gap)) {
      sprintf("%.3f [%.3f, %.3f]", z$observed_roc_gap,
              z$bootstrap_gap_low, z$bootstrap_gap_high)
    } else "not selected"
    line <- sprintf("%-34s  %8d  %8s  %11s  %23s  %11.1f%%",
      labels[i], z$anomaly_count,
      if (is.finite(z$selected_gamma)) format(z$selected_gamma) else "-",
      if (is.finite(z$roc_best_gamma)) format(z$roc_best_gamma) else "tie",
      gap, 100 * z$roc_best_change_rate)
    graphics::text(0.02, 0.76 - i * 0.073, line,
                   family = "mono", adj = 0, cex = 0.68)
  }
  graphics::text(0.03, 0.09,
    "Gap = ROC at observed ROC-best gamma minus ROC at full-grid AUMVC choice.",
    adj = 0, cex = 0.82)
  graphics::text(0.03, 0.05,
    "No AUMVC choice implies missing gap, not zero regret. A tied ROC-best gamma is marked 'tie'.",
    adj = 0, cex = 0.82)

  values <- 100 * results$roc_best_change_rate
  graphics::par(mar = c(4.5, 13, 3.2, 1.2))
  graphics::barplot(rev(values), names.arg = rev(labels),
    col = rev(unname(colors[results$representation])),
    horiz = TRUE, las = 1, xlim = c(0, 100),
    xlab = "Bootstrap samples with a different ROC-optimal gamma (%)",
    main = "Stability of the ROC-optimal bandwidth")
  graphics::mtext("Same resampled observations for all gamma candidates", side = 3,
                  line = 0.25, cex = 0.75)

  eligible <- which(is.finite(results$bootstrap_gap_low) &
                    is.finite(results$bootstrap_gap_high))
  graphics::par(mar = c(4.5, 13, 3.2, 1.2))
  if (length(eligible)) {
    low <- results$bootstrap_gap_low[eligible]
    high <- results$bootstrap_gap_high[eligible]
    mid <- results$observed_roc_gap[eligible]
    y <- seq_along(eligible)
    lim <- range(c(low, high, mid, 0))
    graphics::plot(mid, y, pch = 19, yaxt = "n", xlim = lim,
                   ylim = c(0.5, length(y) + 0.5),
                   ylab = "", xlab = "Paired ROC gap",
                   main = "Observed fixed-gamma ROC gap with 95% bootstrap interval")
    graphics::axis(2, at = y, labels = labels[eligible], las = 1, cex.axis = 0.8)
    graphics::segments(low, y, high, y, lwd = 2)
    graphics::abline(v = 0, lty = 2, col = "grey50")
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No fully selected cases with a unique ROC optimum")
  }
  invisible(path)
}

run_widen_bootstrap <- function(repetitions = 1000L) {
  settings <- AIM3_SETTINGS
  settings$bandwidth$gamma_scales <- c(2^-6, 2^-2, 1, 2^2, 2^6)
  directory <- "experiments/aim3/widen_result/synthetic/widen_experiment/analysis"
  source_file <- file.path(directory, "selection.csv")
  if (!file.exists(source_file)) stop("Run widen_grid_experiment.R first.", call. = FALSE)
  selection <- read.csv(source_file, stringsAsFactors = FALSE)
  cases <- data.frame(
    truth = c("polynomial_interaction", "oscillatory_local"),
    intrinsic_dim = c(3L, 5L), original_dim = c(80L, 400L),
    snr = c(5, 1), stringsAsFactors = FALSE
  )
  results <- lapply(seq_len(nrow(cases)), function(i) {
    paired_bootstrap_case(cases[i, , drop = FALSE], 1L, settings,
                          selection, repetitions)
  })
  results <- do.call(rbind, results)
  rownames(results) <- NULL
  write_report(results, file.path(directory, "paired_bootstrap.csv"),
               settings$output_digits)
  summary <- bootstrap_summary(results)
  write_report(summary, file.path(directory, "paired_bootstrap_summary.csv"),
               settings$output_digits)
  save_bootstrap_pdf(results, summary,
                     file.path(directory, "paired_bootstrap_report.pdf"),
                     repetitions)
  print(summary, row.names = FALSE)
  invisible(results)
}

if (sys.nframe() == 0L) run_widen_bootstrap()
