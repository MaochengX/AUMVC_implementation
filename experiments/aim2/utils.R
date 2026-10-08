aim2_make_splits <- function(labels, counts, split_names, seed) {
  labels <- as.integer(labels)
  counts <- as.integer(counts)

  if (
    length(counts) != length(split_names) ||
    length(counts) < 2L ||
    any(counts < 1L) ||
    sum(counts) > length(labels) ||
    !all(labels %in% c(0L, 1L))
  ) {
    stop("Invalid Aim 2 split request.", call. = FALSE)
  }

  n_label <- tail(counts, 1L)
  positive <- which(labels == 1L)
  negative <- which(labels == 0L)
  n_positive <- round(n_label * length(positive) / length(labels))
  n_positive <- max(1L, min(n_positive, length(positive), n_label - 1L))
  n_negative <- n_label - n_positive

  if (n_negative > length(negative)) {
    n_negative <- length(negative)
    n_positive <- n_label - n_negative
  }
  if (n_positive > length(positive) || n_negative < 1L) {
    stop("The label-evaluation split cannot contain both classes.", call. = FALSE)
  }

  set.seed(seed)
  label_eval <- c(
    sample(positive, n_positive, replace = FALSE),
    sample(negative, n_negative, replace = FALSE)
  )
  label_eval <- sample(label_eval)
  remaining <- setdiff(seq_along(labels), label_eval)
  remaining_counts <- head(counts, -1L)
  selected <- sample(remaining, sum(remaining_counts), replace = FALSE)
  ends <- cumsum(remaining_counts)
  starts <- c(1L, head(ends, -1L) + 1L)
  splits <- Map(function(first, last) selected[first:last], starts, ends)
  splits[[length(split_names)]] <- label_eval
  names(splits) <- split_names
  combined <- unlist(splits, use.names = FALSE)
  if (anyDuplicated(combined)) {
    stop("Aim 2 splits overlap.", call. = FALSE)
  }
  splits
}

aim2_limit_contamination <- function(x, labels, maximum, seed) {
  positive <- which(labels == 1L)
  negative <- which(labels == 0L)
  maximum_positive <- floor(maximum * length(negative) / (1 - maximum))
  if (length(positive) > maximum_positive) {
    set.seed(seed)
    positive <- sample(positive, maximum_positive, replace = FALSE)
  }
  set.seed(seed + 1L)
  keep <- sample(c(negative, positive), replace = FALSE)
  list(
    x = x[keep, , drop = FALSE],
    labels = labels[keep],
    anomaly_fraction = mean(labels[keep] == 1L)
  )
}

aim2_split_counts <- function(n, settings) {
  fractions <- settings$split_fractions
  maximums <- settings$split_maximums[names(fractions)]
  counts <- pmin(floor(n * fractions), maximums)
  counts <- pmax(counts, 1L)
  remaining <- n - sum(counts)
  for (name in c("aumvc", "reference", "label_eval", "detector_train")) {
    addition <- min(remaining, maximums[name] - counts[name])
    counts[name] <- counts[name] + addition
    remaining <- remaining - addition
    if (remaining == 0L) break
  }
  if (sum(counts) > n) {
    counts <- floor(n * fractions / sum(fractions))
    counts <- pmax(counts, 1L)
  }
  counts <- as.integer(counts)
  names(counts) <- names(fractions)
  counts
}

fit_standardizer <- function(x) {
  center <- colMeans(x)
  scale <- apply(x, 2L, sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(center = center, scale = scale)
}

apply_standardizer <- function(x, standardizer) {
  x <- sweep(x, 2L, standardizer$center, "-")
  sweep(x, 2L, standardizer$scale, "/")
}

roc_auc_score <- function(labels, scores) {
  n_positive <- sum(labels == 1L)
  n_negative <- sum(labels == 0L)
  if (n_positive == 0L || n_negative == 0L) {
    stop("ROC-AUC requires both classes.", call. = FALSE)
  }
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1L]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

pr_auc_score <- function(labels, scores) {
  n_positive <- sum(labels == 1L)
  if (n_positive == 0L) stop("PR-AUC requires positive observations.", call. = FALSE)

  index <- order(scores, decreasing = TRUE)
  labels <- labels[index]
  scores <- scores[index]
  true_positive <- cumsum(labels == 1L)
  end_of_tie <- c(which(scores[-length(scores)] != scores[-1L]), length(scores))
  precision <- true_positive[end_of_tie] / end_of_tie
  recall <- true_positive[end_of_tie] / n_positive
  trapezoid_area(c(0, recall), c(1, precision))
}

format_mean_sd <- function(mean_value, sd_value) {
  if (!is.finite(mean_value)) return(as.character(mean_value))
  mean_text <- aim2_format_number(mean_value)
  if (is.na(sd_value)) return(mean_text)
  paste0(mean_text, " +/- ", aim2_format_number(sd_value))
}

aim2_format_number <- function(value) {
  if (is.na(value)) return(NA_character_)
  if (!is.finite(value)) return(as.character(value))
  if (value != 0 && (abs(value) < 1e-4 || abs(value) >= 1e6)) {
    return(sprintf("%.4e", value))
  }
  sprintf("%.4f", value)
}

aim2_next_experiment <- function(root) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  entries <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  numbers <- suppressWarnings(as.integer(sub("^exp", "", entries)))
  next_number <- if (any(!is.na(numbers))) max(numbers, na.rm = TRUE) + 1L else 1L
  path <- file.path(root, paste0("exp", next_number))
  dir.create(path, recursive = FALSE)
  path
}

aim2_save_csv <- function(result, path) {
  numeric_columns <- vapply(result, is.numeric, logical(1))
  numeric_columns[names(numeric_columns) %in% c(
    "minimum_hit_count",
    "reliable_runs",
    "matches",
    "compared"
  )] <- FALSE
  result[numeric_columns] <- lapply(result[numeric_columns], function(column) {
    vapply(column, aim2_format_number, character(1))
  })
  utils::write.csv(result, path, row.names = FALSE, na = "")
  path
}
