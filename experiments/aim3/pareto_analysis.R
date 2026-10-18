mean_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

sum_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) sum(x) else NA_real_
}

group_rows <- function(data, columns) {
  key <- do.call(
    paste,
    c(lapply(data[columns], function(x) ifelse(is.na(x), "<na>", x)), sep = "\r")
  )
  split(seq_len(nrow(data)), key)
}

read_result <- function(path, required) {
  if (!file.exists(path)) stop("Missing result file: ", path, call. = FALSE)
  data <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "), call. = FALSE)
  data
}

aggregate_runs <- function(results, bandwidth, type) {
  run_groups <- if (type == "real") {
    c("dataset", "run", "representation")
  } else {
    c("truth", "intrinsic_dim", "snr", "original_dim", "run", "representation")
  }
  candidate <- lapply(group_rows(results, run_groups), function(index) {
    part <- results[index, , drop = FALSE]
    row <- part[1L, run_groups, drop = FALSE]
    row$total_runtime_sec <- sum_finite(part$runtime_sec)
    row$zero_occupancy <- mean_finite(part$zero_occupancy)
    row$low_occupancy <- mean_finite(part$low_occupancy)
    row
  })
  candidate <- do.call(rbind, candidate)
  merged <- merge(candidate, bandwidth, by = run_groups, all.x = TRUE, sort = FALSE)
  groups <- setdiff(run_groups, "run")
  rows <- lapply(group_rows(merged, groups), function(index) {
    part <- merged[index, , drop = FALSE]
    row <- part[1L, groups, drop = FALSE]
    for (metric in c(
      "total_runtime_sec", "zero_occupancy", "low_occupancy",
      "unreliable_pairs_rate", "agreement_pairs_rank_rate",
      "reliable_agreement_pairs_rank_rate"
    )) {
      row[[metric]] <- mean_finite(part[[metric]])
    }
    row
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

pareto_flag <- function(x, y, minimize_y = TRUE) {
  valid <- is.finite(x) & is.finite(y)
  result <- rep(NA, length(x))
  if (!any(valid)) return(result)
  xx <- x[valid]
  yy <- if (minimize_y) y[valid] else -y[valid]
  flag <- rep(TRUE, length(xx))
  for (i in seq_along(xx)) {
    dominated <- (xx <= xx[i] & yy <= yy[i]) & (xx < xx[i] | yy < yy[i])
    dominated[i] <- FALSE
    flag[i] <- !any(dominated)
  }
  result[valid] <- flag
  result
}

frontier_table <- function(summary, type) {
  criteria <- list(
    runtime_low_occupancy = c("low_occupancy", "min"),
    runtime_zero_occupancy = c("zero_occupancy", "min"),
    runtime_unreliable_pairs = c("unreliable_pairs_rate", "min"),
    runtime_agreement = c("agreement_pairs_rank_rate", "max"),
    runtime_reliable_agreement = c("reliable_agreement_pairs_rank_rate", "max")
  )
  case_groups <- if (type == "real") {
    "dataset"
  } else {
    c("truth", "intrinsic_dim", "snr", "original_dim")
  }
  rows <- list()
  k <- 0L
  for (index in group_rows(summary, case_groups)) {
    part <- summary[index, , drop = FALSE]
    for (name in names(criteria)) {
      metric <- criteria[[name]][1L]
      minimize <- criteria[[name]][2L] == "min"
      flag <- pareto_flag(part$total_runtime_sec, part[[metric]], minimize)
      for (i in seq_len(nrow(part))) {
        k <- k + 1L
        row <- part[i, c(case_groups, "representation"), drop = FALSE]
        row$criterion <- name
        row$runtime_sec <- part$total_runtime_sec[i]
        row$metric_value <- part[[metric]][i]
        row$metric_direction <- if (minimize) "min" else "max"
        row$pareto_optimal <- flag[i]
        rows[[k]] <- row
      }
    }
  }
  do.call(rbind, rows)
}

frontier_frequency <- function(data) {
  groups <- c("criterion", "representation")
  rows <- lapply(group_rows(data, groups), function(index) {
    part <- data[index, , drop = FALSE]
    valid <- !is.na(part$pareto_optimal)
    data.frame(
      criterion = part$criterion[1L],
      representation = part$representation[1L],
      evaluated_cases = sum(valid),
      pareto_cases = sum(part$pareto_optimal[valid]),
      pareto_rate = if (any(valid)) mean(part$pareto_optimal[valid]) else NA_real_
    )
  })
  do.call(rbind, rows)
}

write_report <- function(data, path, digits = 6L) {
  output <- data
  integer_columns <- intersect(c("intrinsic_dim", "original_dim", "evaluated_cases", "pareto_cases"), names(output))
  for (column in integer_columns) {
    values <- output[[column]]
    output[[column]] <- ifelse(is.na(values), NA_character_, as.character(as.integer(round(values))))
  }
  numeric <- names(output)[vapply(output, is.numeric, logical(1))]
  for (column in numeric) {
    values <- output[[column]]
    output[[column]] <- ifelse(
      is.finite(values),
      formatC(values, format = "f", digits = digits),
      NA_character_
    )
  }
  utils::write.csv(output, path, row.names = FALSE, na = "", quote = FALSE)
}

plot_real_frontiers <- function(data, path) {
  representations <- c("ambient", "goix_subsampling", "mds", "isomap")
  pch <- c(16, 17, 15, 18)
  datasets <- unique(data$dataset)
  criteria <- unique(data$criterion)
  grDevices::pdf(path, width = 12, height = 4)
  for (criterion in criteria) {
    old <- par(mfrow = c(1, length(datasets)), mar = c(4, 4, 3, 1))
    for (dataset in datasets) {
      part <- data[data$dataset == dataset & data$criterion == criterion, , drop = FALSE]
      valid <- is.finite(part$runtime_sec) & is.finite(part$metric_value)
      part <- part[valid, , drop = FALSE]
      if (!nrow(part)) {
        plot.new()
        title(dataset)
        next
      }
      plot(
        part$runtime_sec, part$metric_value,
        xlab = "total runtime (sec)", ylab = criterion,
        main = dataset, pch = pch[match(part$representation, representations)]
      )
      text(part$runtime_sec, part$metric_value, labels = part$representation, pos = 3, cex = 0.7)
      frontier <- part[part$pareto_optimal %in% TRUE, , drop = FALSE]
      if (nrow(frontier) > 1L) {
        frontier <- frontier[order(frontier$runtime_sec), , drop = FALSE]
        lines(frontier$runtime_sec, frontier$metric_value, lwd = 2)
      }
    }
    par(old)
  }
  grDevices::dev.off()
}

plot_frequency <- function(data, path, title) {
  criteria <- unique(data$criterion)
  representations <- c("ambient", "goix_subsampling", "mds", "isomap")
  grDevices::pdf(path, width = 8, height = 5)
  for (criterion in criteria) {
    part <- data[data$criterion == criterion, , drop = FALSE]
    values <- setNames(part$pareto_rate, part$representation)[representations]
    values[is.na(values)] <- 0
    barplot(values, ylim = c(0, 1), ylab = "Pareto-optimal rate", main = paste(title, criterion))
  }
  grDevices::dev.off()
}

run_pareto_analysis <- function() {
  real_root <- "experiments/aim3/result/real"
  syn_root <- "experiments/aim3/result/synthetic"
  required_results <- c("run", "representation", "runtime_sec", "zero_occupancy", "low_occupancy")
  required_bandwidth <- c(
    "run", "representation", "unreliable_pairs_rate",
    "agreement_pairs_rank_rate", "reliable_agreement_pairs_rank_rate"
  )

  real_results <- read_result(file.path(real_root, "experiment/real_results.csv"), c("dataset", required_results))
  real_bandwidth <- read_result(file.path(real_root, "experiment/real_bandwidth_exp.csv"), c("dataset", required_bandwidth))
  syn_results <- read_result(
    file.path(syn_root, "experiment/syn_results.csv"),
    c("truth", "intrinsic_dim", "snr", "original_dim", required_results)
  )
  syn_bandwidth <- read_result(
    file.path(syn_root, "experiment/syn_bandwidth_exp.csv"),
    c("truth", "intrinsic_dim", "snr", "original_dim", required_bandwidth)
  )

  real_summary <- aggregate_runs(real_results, real_bandwidth, "real")
  syn_summary <- aggregate_runs(syn_results, syn_bandwidth, "synthetic")
  real_pareto <- frontier_table(real_summary, "real")
  syn_pareto <- frontier_table(syn_summary, "synthetic")
  real_frequency <- frontier_frequency(real_pareto)
  syn_frequency <- frontier_frequency(syn_pareto)

  real_out <- file.path(real_root, "analysis/pareto")
  syn_out <- file.path(syn_root, "analysis/pareto")
  dir.create(real_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(syn_out, recursive = TRUE, showWarnings = FALSE)

  write_report(real_pareto, file.path(real_out, "real_pareto.csv"))
  write_report(real_frequency, file.path(real_out, "real_pareto_frequency.csv"))
  write_report(syn_pareto, file.path(syn_out, "syn_pareto.csv"))
  write_report(syn_frequency, file.path(syn_out, "syn_pareto_frequency.csv"))

  plot_real_frontiers(real_pareto, file.path(real_out, "real_pareto_figures.pdf"))
  plot_frequency(syn_frequency, file.path(syn_out, "syn_pareto_frequency.pdf"), "synthetic")

  invisible(list(real = real_pareto, synthetic = syn_pareto))
}

if (sys.nframe() == 0L) run_pareto_analysis()
