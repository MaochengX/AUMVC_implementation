representation_levels <- c("ambient", "goix_subsampling", "mds", "isomap")
representation_labels <- c("ambient", "goix", "mds", "isomap")
representation_colors <- c("#D55E00", "#009E73", "#0072B2", "#CC79A7")

results_metrics_real <- c(
  "aumvc_normalized", "aumvc_relative_mc_se", "min_hit_count",
  "zero_occupancy", "low_occupancy", "box_log_volume",
  "roc_auc", "pr_auc", "embedding_oos_stress", "runtime_sec"
)

results_metrics_synthetic <- c(
  results_metrics_real,
  "roc_distributional", "pr_distributional",
  "roc_structural", "pr_structural"
)

bandwidth_metrics <- c(
  "compared_pairs_count", "unreliable_pairs_rate",
  "agreement_pairs_rank_rate", "reliable_agreement_pairs_rank_rate",
  "tie_rate"
)

summary_results <- function(data, type) {
  groups <- if (type == "real") {
    c("dataset", "sample_size", "original_dim", "embedding_dim", "representation", "gamma_scale")
  } else {
    c("truth", "intrinsic_dim", "snr", "sample_size", "original_dim", "embedding_dim", "representation", "gamma_scale")
  }
  metrics <- if (type == "real") results_metrics_real else results_metrics_synthetic
  rows <- lapply(group_rows(data, groups), function(index) {
    part <- data[index, , drop = FALSE]
    row <- part[1L, groups, drop = FALSE]
    row$averaged_detector_dim <- mean_finite(part$detector_dim)
    for (metric in metrics) row[[paste0("averaged_", metric)]] <- mean_finite(part[[metric]])
    row$aumvc_reliable_rate <- mean(part$aumvc_reliable)
    row
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

summary_bandwidth <- function(data, type) {
  groups <- if (type == "real") {
    c("dataset", "sample_size", "original_dim", "embedding_dim", "representation")
  } else {
    c("truth", "intrinsic_dim", "snr", "sample_size", "original_dim", "embedding_dim", "representation")
  }
  rows <- lapply(group_rows(data, groups), function(index) {
    part <- data[index, , drop = FALSE]
    row <- part[1L, groups, drop = FALSE]
    row$averaged_detector_dim <- mean_finite(part$detector_dim)
    for (metric in bandwidth_metrics) row[[paste0("averaged_", metric)]] <- mean_finite(part[[metric]])
    row$selected_gamma_scale_mode <- mode_unique(part$selected_gamma_scale)
    row$selected_gamma_scale_reliable_mode <- mode_unique(part$selected_gamma_scale_reliable)
    row
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

metric_limits <- function(values, probability = FALSE) {
  values <- values[is.finite(values)]
  if (probability) return(c(0, 1))
  if (!length(values)) return(c(0, 1))
  limits <- range(values)
  if (limits[1L] == limits[2L]) {
    pad <- if (limits[1L] == 0) 1 else abs(limits[1L]) * 0.08
    return(limits + c(-pad, pad))
  }
  pad <- diff(limits) * 0.08
  limits + c(-pad, pad)
}

plot_real_points <- function(data, metric, title) {
  methods <- representation_levels[representation_levels %in% data$representation]
  y <- data[[metric]]
  probability <- metric %in% c(
    "zero_occupancy", "low_occupancy",
    "agreement_pairs_rank_rate", "reliable_agreement_pairs_rank_rate"
  )
  plot(
    NA, xlim = c(0.5, length(methods) + 0.5),
    ylim = metric_limits(y, probability), xaxt = "n",
    xlab = "representation", ylab = metric, main = title
  )
  axis(1, at = seq_along(methods), labels = representation_labels[match(methods, representation_levels)])
  grid(nx = NA, ny = NULL, col = "grey90")
  for (i in seq_along(methods)) {
    values <- y[data$representation == methods[i] & is.finite(y)]
    if (!length(values)) next
    offset <- if (length(values) == 1L) 0 else seq(-0.16, 0.16, length.out = length(values))
    points(i + offset, values, pch = 16, col = representation_colors[match(methods[i], representation_levels)])
  }
}

plot_real_analysis <- function(results, bandwidth, directory) {
  for (dataset in unique(results$dataset)) {
    file <- file.path(directory, paste0(dataset, "_figures.pdf"))
    grDevices::pdf(file, width = 8, height = 6)
    on.exit(grDevices::dev.off(), add = TRUE)
    part <- results[results$dataset == dataset, , drop = FALSE]
    for (metric in c("zero_occupancy", "low_occupancy", "runtime_sec")) {
      plot_real_points(part, metric, paste(dataset, metric))
    }
    part <- bandwidth[bandwidth$dataset == dataset, , drop = FALSE]
    for (metric in c("agreement_pairs_rank_rate", "reliable_agreement_pairs_rank_rate")) {
      plot_real_points(part, metric, paste(dataset, metric))
    }
    grDevices::dev.off()
    on.exit(NULL, add = FALSE)
  }
}

plot_synthetic_lines <- function(data, metric, title) {
  truths <- unique(data$truth)
  dimensions <- sort(unique(data$intrinsic_dim))
  snrs <- sort(unique(data$snr))
  probability <- grepl("occupancy|agreement_pairs_rank_rate", metric)
  limits <- metric_limits(data[[metric]], probability)
  old <- par(mfrow = c(length(truths), length(dimensions)), mar = c(4, 4, 3, 1), oma = c(1, 1, 3, 1))
  on.exit(par(old))
  for (truth in truths) {
    for (dimension in dimensions) {
      part <- data[data$truth == truth & data$intrinsic_dim == dimension, , drop = FALSE]
      x <- sort(unique(part$original_dim))
      plot(
        NA, xlim = range(x), ylim = limits, xaxt = "n",
        xlab = "original_dim", ylab = metric,
        main = paste(truth, "q =", dimension)
      )
      axis(1, at = x, labels = x)
      grid(nx = NA, ny = NULL, col = "grey90")
      for (r in seq_along(representation_levels)) {
        for (s in seq_along(snrs)) {
          line <- part[
            part$representation == representation_levels[r] & part$snr == snrs[s],
            , drop = FALSE
          ]
          line <- line[order(line$original_dim), , drop = FALSE]
          if (!nrow(line) || !any(is.finite(line[[metric]]))) next
          lines(
            line$original_dim, line[[metric]], type = "b", pch = 16,
            col = representation_colors[r], lty = s, lwd = 1.5, cex = 0.7
          )
        }
      }
      if (truth == truths[1L] && dimension == dimensions[1L]) {
        legend("topleft", representation_labels, col = representation_colors, lty = 1, bty = "n", cex = 0.65)
        legend("bottomleft", paste0("snr=", snrs), col = "black", lty = seq_along(snrs), bty = "n", cex = 0.65)
      }
    }
  }
  mtext(title, outer = TRUE, side = 3, line = 1, cex = 1.2)
}

plot_synthetic_analysis <- function(results_summary, bandwidth_summary, directory) {
  file <- file.path(directory, "syn_results_figures.pdf")
  grDevices::pdf(file, width = 14, height = 8)
  for (scale in sort(unique(results_summary$gamma_scale))) {
    part <- results_summary[results_summary$gamma_scale == scale, , drop = FALSE]
    for (metric in c("averaged_zero_occupancy", "averaged_low_occupancy", "averaged_runtime_sec")) {
      plot_synthetic_lines(part, metric, paste(metric, "gamma_scale =", scale))
    }
  }
  grDevices::dev.off()
  file <- file.path(directory, "syn_bandwidth_figures.pdf")
  grDevices::pdf(file, width = 14, height = 8)
  for (metric in c("averaged_agreement_pairs_rank_rate", "averaged_reliable_agreement_pairs_rank_rate")) {
    plot_synthetic_lines(bandwidth_summary, metric, metric)
  }
  grDevices::dev.off()
}

analyze_results <- function(results, bandwidth, type, directory, digits = 6L) {
  results_summary <- summary_results(results, type)
  bandwidth_summary <- summary_bandwidth(bandwidth, type)
  prefix <- if (type == "real") "real" else "syn"
  write_report(
    results_summary,
    file.path(directory, paste0(prefix, "_results_summary.csv")),
    digits
  )
  write_report(
    bandwidth_summary,
    file.path(directory, paste0(prefix, "_bandwidth_summary.csv")),
    digits
  )
  if (type == "real") {
    plot_real_analysis(results, bandwidth, directory)
  } else {
    plot_synthetic_analysis(results_summary, bandwidth_summary, directory)
  }
  invisible(list(results_summary = results_summary, bandwidth_summary = bandwidth_summary))
}
