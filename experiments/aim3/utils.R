stratified_indices <- function(strata, size, seed) {
  if (length(strata) < 1L || anyNA(strata)) {
    stop("Invalid stratification labels.", call. = FALSE)
  }
  size <- as.numeric(size)
  if (length(size) != 1L || !is.finite(size) || size < 1L || size != floor(size)) {
    stop("Invalid stratified split size.", call. = FALSE)
  }
  size <- as.integer(size)
  strata <- as.character(strata)
  groups <- unique(strata)
  available <- vapply(groups, function(x) sum(strata == x), integer(1))
  if (size < length(groups) || size > length(strata)) {
    stop("Invalid stratified split size.", call. = FALSE)
  }
  allocation <- rep(1L, length(groups))
  remaining <- size - length(groups)
  capacity <- available - 1L
  if (remaining > 0L) {
    target <- remaining * capacity / sum(capacity)
    add <- pmin(floor(target), capacity)
    allocation <- allocation + add
    remaining <- size - sum(allocation)
    priority <- order(target - floor(target), decreasing = TRUE)
    while (remaining > 0L) {
      for (i in priority) {
        if (remaining > 0L && allocation[i] < available[i]) {
          allocation[i] <- allocation[i] + 1L
          remaining <- remaining - 1L
        }
      }
    }
  }
  set.seed(seed)
  selected <- unlist(lapply(seq_along(groups), function(i) {
    rows <- which(strata == groups[i])
    rows[sample.int(length(rows), allocation[i])]
  }), use.names = FALSE)
  selected[sample.int(length(selected))]
}

make_splits <- function(labels, counts, seed, strata = labels) {
  names_expected <- c(
    "embedding", "detector_train", "reference", "evaluation", "label_eval"
  )
  if (
    is.null(names(counts)) ||
    anyDuplicated(names(counts)) ||
    !setequal(names(counts), names_expected)
  ) {
    stop("Invalid split counts.", call. = FALSE)
  }
  counts <- counts[names_expected]
  if (
    anyNA(counts) ||
    any(!is.finite(counts)) ||
    any(counts < 1) ||
    any(counts != floor(counts)) ||
    sum(counts) != length(labels) ||
    length(strata) != length(labels)
  ) {
    stop("Invalid split counts.", call. = FALSE)
  }
  counts <- setNames(as.integer(counts), names_expected)
  label_eval <- stratified_indices(strata, counts["label_eval"], seed)
  remaining <- setdiff(seq_along(labels), label_eval)
  set.seed(seed + 1L)
  remaining <- remaining[sample.int(length(remaining))]
  sizes <- counts[names_expected[-5L]]
  ends <- cumsum(sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  splits <- Map(function(a, b) remaining[a:b], starts, ends)
  names(splits) <- names_expected[-5L]
  splits$label_eval <- label_eval
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
  labels <- as.integer(labels)
  positives <- sum(labels == 1L)
  negatives <- sum(labels == 0L)
  if (!positives || !negatives) stop("ROC-AUC requires both classes.", call. = FALSE)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1L]) - positives * (positives + 1) / 2) /
    (positives * negatives)
}

pr_auc_score <- function(labels, scores) {
  labels <- as.integer(labels)
  positives <- sum(labels == 1L)
  if (!positives) stop("PR-AUC requires positive observations.", call. = FALSE)
  index <- order(scores, decreasing = TRUE)
  labels <- labels[index]
  scores <- scores[index]
  true_positive <- cumsum(labels == 1L)
  end_tie <- c(which(diff(scores) != 0), length(scores))
  precision <- true_positive[end_tie] / end_tie
  recall <- true_positive[end_tie] / positives
  trapezoid_area(c(0, recall), c(1, precision))
}

percentile_scores <- function(calibration_scores, scores) {
  calibration_scores <- sort(validate_scores(calibration_scores))
  scores <- validate_scores(scores)
  findInterval(scores, calibration_scores, rightmost.closed = TRUE) /
    length(calibration_scores)
}

attach_metadata <- function(metadata, result) {
  cbind(metadata[rep(1L, nrow(result)), , drop = FALSE], result, row.names = NULL)
}

result_dirs <- function(type) {
  root <- file.path("experiments/aim3/result", type)
  experiment <- file.path(root, "experiment")
  analysis <- file.path(root, "analysis")
  dir.create(experiment, recursive = TRUE, showWarnings = FALSE)
  dir.create(analysis, recursive = TRUE, showWarnings = FALSE)
  list(experiment = experiment, analysis = analysis)
}

write_report <- function(data, path, digits = 6L) {
  output <- data
  integer_columns <- c(
    "run", "sample_size", "original_dim", "embedding_dim", "detector_dim",
    "intrinsic_dim", "min_hit_count", "compared_pairs_count"
  )
  numeric <- names(output)[vapply(output, is.numeric, logical(1))]
  for (column in intersect(integer_columns, numeric)) {
    values <- output[[column]]
    output[[column]] <- ifelse(is.na(values), NA_character_, as.character(as.integer(round(values))))
  }
  for (column in setdiff(numeric, integer_columns)) {
    values <- output[[column]]
    output[[column]] <- ifelse(
      is.finite(values),
      formatC(values, format = "f", digits = digits),
      NA_character_
    )
  }
  utils::write.csv(output, path, row.names = FALSE, na = "", quote = FALSE)
  path
}

mean_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

mode_unique <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  count <- table(x)
  winners <- names(count)[count == max(count)]
  if (length(winners) == 1L) as.numeric(winners) else NA_real_
}

group_rows <- function(data, columns) {
  key <- do.call(paste, c(lapply(data[columns], function(x) ifelse(is.na(x), "<na>", x)), sep = "\r"))
  split(seq_len(nrow(data)), key)
}

validate_runs <- function(n_runs) {
  if (length(n_runs) != 1L || !is.finite(n_runs) || n_runs < 1L || n_runs != floor(n_runs)) {
    stop("n_runs must be a positive integer.", call. = FALSE)
  }
  as.integer(n_runs)
}
