aim3_method_names <- c(
  ambient = "Ambient", goix_subsampling = "Goix subsets",
  mds = "MDS", isomap = "Isomap"
)
aim3_method_colors <- c(
  ambient = "#D55E00", goix_subsampling = "#009E73",
  mds = "#0072B2", isomap = "#CC79A7"
)
aim3_truth_names <- c(
  polynomial_interaction = "Polynomial", oscillatory_local = "Oscillatory"
)
aim3_analysis_metrics <- c(
  "aumvc_normalized", "minimum_hit_count", "zero_occupancy", "low_occupancy",
  "roc_auc", "pr_auc", "pr_delta_ambient", "roc_delta_ambient",
  "embedding_oos_stress", "runtime_seconds"
)
aim3_figure_metrics <- c(
  averaged_zero_occupancy = "Averaged zero occupancy",
  averaged_low_occupancy = "Averaged low occupancy",
  averaged_runtime_seconds = "Averaged runtime (seconds)",
  averaged_embedding_oos_stress = "Averaged embedding OOS stress"
)
aim3_marker_symbols <- c(16, 17, 15, 18, 3, 4, 7, 8, 0, 1, 2, 5, 6)

aim3_group_key <- function(data, columns) {
  do.call(paste, c(data[columns], sep = "\r"))
}

aim3_group_indices <- function(data, columns) {
  split(seq_len(nrow(data)), aim3_group_key(data, columns))
}

aim3_read_analysis_data <- function(input_file, type) {
  if (!file.exists(input_file)) stop("Input file not found: ", input_file)
  data <- read.csv(
    input_file, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("", "NA", "NaN", "N/A")
  )
  if (!nrow(data)) stop("The input file contains no results.")
  if (!"original_dim" %in% names(data) && "ambient_dim" %in% names(data)) {
    data$original_dim <- data$ambient_dim
  }
  case_columns <- if (type == "real") {
    c("original_dim", "representation")
  } else {
    c("intrinsic_dim", "original_dim", "snr", "embedding_dim", "representation", "truth")
  }
  required <- c(case_columns, "run", aim3_analysis_metrics)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }
  numeric_columns <- intersect(c(
    "run", "original_dim", "embedding_dim", "intrinsic_dim", "snr",
    aim3_analysis_metrics
  ), names(data))
  for (column in numeric_columns) {
    converted <- suppressWarnings(as.numeric(data[[column]]))
    if (any(is.na(converted) & !is.na(data[[column]]))) {
      stop("Non-numeric values in column: ", column)
    }
    data[[column]] <- converted
  }
  data$representation <- trimws(data$representation)
  identifiers <- unique(c(
    case_columns, "run", intersect(c("dataset", "embedding_dim"), names(data))
  ))
  if (anyNA(data[identifiers])) stop("Missing experiment identifiers.")
  numeric_identifiers <- intersect(identifiers, numeric_columns)
  if (any(!is.finite(as.matrix(data[numeric_identifiers])))) {
    stop("Non-finite experiment identifiers.")
  }
  if (any(data$run < 1 | data$run != floor(data$run))) {
    stop("Run identifiers must be positive integers.")
  }
  dimensions <- intersect(c("original_dim", "intrinsic_dim", "embedding_dim"), names(data))
  if (any(as.matrix(data[dimensions]) <= 0)) stop("Dimensions must be positive.")
  unknown <- setdiff(unique(data$representation), names(aim3_method_names))
  if (length(unknown)) stop("Unknown representations: ", paste(unknown, collapse = ", "))
  probabilities <- c("aumvc_normalized", "zero_occupancy", "low_occupancy", "roc_auc", "pr_auc")
  for (column in probabilities) {
    values <- data[[column]]
    if (any(is.finite(values) & (values < 0 | values > 1))) {
      stop("Values outside [0, 1] in column: ", column)
    }
  }
  for (column in c("minimum_hit_count", "embedding_oos_stress", "runtime_seconds")) {
    if (any(data[[column]] < 0, na.rm = TRUE)) stop("Negative values in column: ", column)
  }
  if (any(duplicated(aim3_group_key(data, identifiers)))) {
    stop("Duplicate rows for the same setting, run and representation.")
  }
  data
}

aim3_complete_mean <- function(values) {
  if (length(values) && all(is.finite(values))) mean(values) else NA_real_
}

aim3_truth_order <- function(values) {
  c(intersect(names(aim3_truth_names), values), sort(setdiff(values, names(aim3_truth_names))))
}

aim3_summary_table <- function(data, type) {
  group_columns <- if (type == "real") {
    c("original_dim", "representation")
  } else {
    c("intrinsic_dim", "original_dim", "snr", "embedding_dim", "representation", "truth")
  }
  groups <- aim3_group_indices(data, group_columns)
  rows <- lapply(groups, function(indices) {
    part <- data[indices, , drop = FALSE]
    result <- part[1L, group_columns, drop = FALSE]
    for (metric in aim3_analysis_metrics) {
      result[[paste0("averaged_", metric)]] <- aim3_complete_mean(part[[metric]])
    }
    result
  })
  summary <- do.call(rbind, rows)
  method_order <- match(summary$representation, names(aim3_method_names))
  ordering <- if (type == "real") {
    order(summary$original_dim, method_order)
  } else {
    truth_order <- match(summary$truth, aim3_truth_order(unique(summary$truth)))
    order(truth_order, summary$intrinsic_dim, summary$original_dim,
          summary$snr, summary$embedding_dim, method_order)
  }
  summary <- summary[ordering, , drop = FALSE]
  rownames(summary) <- NULL
  summary
}

aim3_write_summary <- function(summary, path) {
  formatted <- summary
  for (column in paste0("averaged_", aim3_analysis_metrics)) {
    values <- summary[[column]]
    scientific <- is.finite(values) & values != 0 & (abs(values) < 0.0001 | abs(values) >= 1000000)
    formatted[[column]] <- ifelse(
      is.finite(values), formatC(values, format = "f", digits = 4), NA_character_
    )
    formatted[[column]][scientific] <- formatC(values[scientific], format = "e", digits = 4)
  }
  write.csv(formatted, path, row.names = FALSE, na = "NA")
}

aim3_dimension_legend <- function(data, type) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  if (!nrow(data)) return(invisible(NULL))
  methods <- names(aim3_method_names)[names(aim3_method_names) %in% data$representation]
  key <- legend(
    0, 0.98, legend = aim3_method_names[methods], col = aim3_method_colors[methods],
    lty = 1, lwd = 2, bty = "n", cex = 0.95, y.intersp = 1.2,
    title = "Representation", title.adj = 0, title.font = 2,
    pch = if (type == "real") aim3_marker_symbols[match(methods, names(aim3_method_names))] else NA
  )
  if (type == "real") return(invisible(NULL))
  snr_values <- sort(unique(data$snr))
  key <- legend(
    0, key$rect$top - key$rect$h - 0.06, legend = paste0("SNR = ", snr_values),
    lty = (seq_along(snr_values) - 1L) %% 6L + 1L, lwd = 2,
    bty = "n", cex = 0.95, y.intersp = 1.2,
    title = "SNR (line style)", title.adj = 0, title.font = 2
  )
  embedding_values <- sort(unique(data$embedding_dim))
  legend(
    0, key$rect$top - key$rect$h - 0.06, legend = paste0("q = ", embedding_values),
    pch = aim3_marker_symbols[(seq_along(embedding_values) - 1L) %% length(aim3_marker_symbols) + 1L],
    bty = "n", cex = 0.95, y.intersp = 1.2,
    title = "Embedding dimension (marker)", title.adj = 0, title.font = 2
  )
}

aim3_dimension_figure <- function(summary, metric, type, title) {
  series_columns <- if (type == "real") {
    "representation"
  } else {
    c("representation", "embedding_dim", "snr")
  }
  groups <- aim3_group_indices(summary, series_columns)
  available <- vapply(groups, function(indices) {
    any(is.finite(summary[[metric]][indices]))
  }, logical(1))
  groups <- groups[available]
  legend_data <- if (length(groups)) {
    summary[unlist(groups, use.names = FALSE), , drop = FALSE]
  } else {
    summary[FALSE, , drop = FALSE]
  }
  dimensions <- sort(unique(summary$original_dim))
  x_limits <- range(dimensions)
  padding <- if (diff(x_limits) > 0) diff(x_limits) * 0.04 else x_limits[1L] * 0.1
  x_limits <- x_limits + c(-padding, padding)
  values <- summary[[metric]]
  y_limits <- if (metric %in% c("averaged_zero_occupancy", "averaged_low_occupancy")) {
    c(0, 1)
  } else if (any(is.finite(values)) && max(values[is.finite(values)]) > 0) {
    c(0, max(values[is.finite(values)]) * 1.08)
  } else {
    c(0, 1)
  }
  layout(matrix(c(1, 2), nrow = 1L), widths = c(0.73, 0.27))
  par(mar = c(4.8, 5, 4.1, 1.1), oma = c(1.8, 0, 0, 0),
      cex.axis = 1, cex.lab = 1.1, las = 1, family = "sans")
  plot(NA_real_, xlim = x_limits, ylim = y_limits, xaxt = "n",
       xlab = "Original ambient dimension", ylab = aim3_figure_metrics[[metric]],
       main = paste(title, aim3_figure_metrics[[metric]], sep = "\n"))
  axis(1, at = dimensions, labels = dimensions)
  grid(nx = NA, ny = NULL, col = "#E5E7EB", lty = 1)
  snr_values <- if (type == "synthetic") sort(unique(legend_data$snr)) else NULL
  embedding_values <- if (type == "synthetic") sort(unique(legend_data$embedding_dim)) else NULL
  for (indices in groups) {
    part <- summary[indices, , drop = FALSE]
    if (anyDuplicated(part$original_dim)) stop("Multiple values at one dimension in a figure line.")
    part <- part[order(part$original_dim), , drop = FALSE]
    method <- part$representation[1L]
    line_type <- if (type == "synthetic") {
      (match(part$snr[1L], snr_values) - 1L) %% 6L + 1L
    } else {
      1L
    }
    marker_index <- if (type == "synthetic") {
      match(part$embedding_dim[1L], embedding_values)
    } else {
      match(method, names(aim3_method_names))
    }
    marker <- aim3_marker_symbols[(marker_index - 1L) %% length(aim3_marker_symbols) + 1L]
    lines(part$original_dim, part[[metric]], type = "b", lwd = 2, lty = line_type,
          pch = marker, col = aim3_method_colors[[method]], cex = 1)
  }
  if (!length(groups)) text(mean(x_limits), mean(y_limits), "No available values")
  aim3_dimension_legend(legend_data, type)
  footer <- if (type == "real") {
    "Each point is a group mean. Real-data sources also change with dimension."
  } else {
    "Each point is a group mean. A line identifies representation, embedding dimension and SNR."
  }
  mtext(footer, side = 1, outer = TRUE, line = 0.5, cex = 0.8)
}

aim3_plot_dimension_results <- function(summary, type) {
  if (type == "real") {
    for (metric in names(aim3_figure_metrics)) {
      aim3_dimension_figure(summary, metric, type, "Aim 3 real results")
    }
    return(invisible(NULL))
  }
  for (truth in aim3_truth_order(unique(summary$truth))) {
    truth_label <- if (truth %in% names(aim3_truth_names)) aim3_truth_names[[truth]] else truth
    part <- summary[summary$truth == truth, , drop = FALSE]
    for (intrinsic in sort(unique(part$intrinsic_dim))) {
      dimension_part <- part[part$intrinsic_dim == intrinsic, , drop = FALSE]
      title <- paste0(truth_label, " | intrinsic dimension = ", intrinsic)
      for (metric in names(aim3_figure_metrics)) {
        aim3_dimension_figure(dimension_part, metric, type, title)
      }
    }
  }
}

aim3_analyze_results <- function(input_file, type = c("real", "synthetic")) {
  type <- match.arg(type)
  data <- aim3_read_analysis_data(input_file, type)
  other_metrics <- setdiff(aim3_analysis_metrics, "embedding_oos_stress")
  if (any(!is.finite(as.matrix(data[other_metrics])))) {
    warning("Missing selected metrics: affected group means remain NA.")
  }
  summary <- aim3_summary_table(data, type)
  output_directory <- file.path(dirname(input_file), "analysis")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_directory)) stop("Cannot create output directory: ", output_directory)
  summary_file <- file.path(output_directory, paste0(type, "_summary.csv"))
  figure_file <- file.path(output_directory, paste0(type, "_figures.pdf"))
  aim3_write_summary(summary, summary_file)
  figure_data <- read.csv(summary_file, stringsAsFactors = FALSE, na.strings = "NA")
  pdf(figure_file, width = 11.7, height = 8.3, onefile = TRUE, family = "Helvetica")
  device <- dev.cur()
  on.exit(if (device %in% dev.list()) dev.off(device), add = TRUE)
  aim3_plot_dimension_results(figure_data, type)
  dev.off(device)
  cat("Summary: ", nrow(summary), " rows, ", ncol(summary), " columns.\n", sep = "")
  cat("Means include all recorded runs; occupancy stays on its original 0-1 scale.\n")
  cat("Figures use the saved summary table; stress lines require available values.\n")
  cat("Saved: ", summary_file, "\nSaved: ", figure_file, "\n", sep = "")
  invisible(summary)
}
