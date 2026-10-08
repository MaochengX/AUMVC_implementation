aim3_stratified_indices <- function(strata, size, seed) {
  strata <- as.character(strata)
  groups <- unique(strata)
  available <- vapply(groups, function(group) sum(strata == group), integer(1))
  if (size < length(groups) || any(available < 1L)) {
    stop("The label-evaluation split cannot represent every stratum.", call. = FALSE)
  }

  allocation <- rep(1L, length(groups))
  remaining <- size - sum(allocation)
  capacity <- available - allocation
  if (remaining > 0L) {
    target <- remaining * capacity / sum(capacity)
    addition <- pmin(floor(target), capacity)
    allocation <- allocation + addition
    remaining <- size - sum(allocation)
    priority <- order(target - floor(target), decreasing = TRUE)
    while (remaining > 0L) {
      changed <- FALSE
      for (index in priority) {
        if (allocation[index] < available[index] && remaining > 0L) {
          allocation[index] <- allocation[index] + 1L
          remaining <- remaining - 1L
          changed <- TRUE
        }
      }
      if (!changed) stop("The stratified split is too large.", call. = FALSE)
    }
  }

  set.seed(seed)
  selected <- unlist(lapply(seq_along(groups), function(index) {
    sample(which(strata == groups[index]), allocation[index], replace = FALSE)
  }), use.names = FALSE)
  sample(selected)
}

aim3_make_splits <- function(labels, counts, split_names, seed, strata = labels) {
  labels <- as.integer(labels)
  counts <- as.integer(counts)
  if (
    length(counts) != length(split_names) ||
    length(counts) < 2L ||
    any(counts < 1L) ||
    sum(counts) != length(labels) ||
    !all(labels %in% c(0L, 1L))
  ) {
    stop("Invalid Aim 3 split request.", call. = FALSE)
  }

  if (length(strata) != length(labels) || anyNA(strata)) {
    stop("Invalid Aim 3 stratification values.", call. = FALSE)
  }
  n_label <- tail(counts, 1L)
  label_eval <- aim3_stratified_indices(strata, n_label, seed)
  remaining <- setdiff(seq_along(labels), label_eval)
  remaining_counts <- head(counts, -1L)
  set.seed(seed + 1L)
  selected <- sample(remaining, sum(remaining_counts), replace = FALSE)
  ends <- cumsum(remaining_counts)
  starts <- c(1L, head(ends, -1L) + 1L)
  splits <- Map(function(first, last) selected[first:last], starts, ends)
  splits[[length(split_names)]] <- label_eval
  names(splits) <- split_names
  splits
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

aim3_percentile_scores <- function(calibration_scores, scores) {
  calibration_scores <- sort(validate_scores(calibration_scores))
  scores <- validate_scores(scores)
  findInterval(scores, calibration_scores, rightmost.closed = TRUE) /
    length(calibration_scores)
}

aim3_format_number <- function(value) {
  if (is.na(value)) return(NA_character_)
  if (!is.finite(value)) return(as.character(value))
  if (value != 0 && (abs(value) < 1e-4 || abs(value) >= 1e6)) {
    return(sprintf("%.4e", value))
  }
  sprintf("%.4f", value)
}

aim3_next_experiment <- function(root) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  entries <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  numbers <- suppressWarnings(as.integer(sub("^exp", "", entries)))
  next_number <- if (any(!is.na(numbers))) max(numbers, na.rm = TRUE) + 1L else 1L
  path <- file.path(root, paste0("exp", next_number))
  dir.create(path, recursive = FALSE)
  path
}

aim3_save_csv <- function(result, path) {
  numeric_columns <- vapply(result, is.numeric, logical(1))
  identifiers <- c(
    "run", "intrinsic_dim", "ambient_dim", "original_dim",
    "effective_dim", "sample_size", "embedding_dim", "snr"
  )
  numeric_columns[names(numeric_columns) %in% identifiers] <- FALSE
  result[numeric_columns] <- lapply(result[numeric_columns], function(column) {
    vapply(column, aim3_format_number, character(1))
  })
  utils::write.csv(result, path, row.names = FALSE, na = "")
  path
}
