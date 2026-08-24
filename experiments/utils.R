make_splits <- function(n, counts, split_names, seed, require_full = FALSE) {
  counts <- as.integer(counts)

  if (
    length(counts) != length(split_names) ||
    any(counts < 1L) ||
    sum(counts) > n ||
    (require_full && sum(counts) != n)
  ) {
    stop("Invalid split sizes.", call. = FALSE)
  }

  set.seed(seed)
  rows <- sample.int(n, sum(counts))
  ends <- cumsum(counts)
  starts <- c(1L, head(ends, -1L) + 1L)
  splits <- Map(function(a, b) rows[a:b], starts, ends)
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

  if (n_positive == 0L) {
    stop("PR-AUC requires positive observations.", call. = FALSE)
  }

  index <- order(scores, decreasing = TRUE)
  labels <- labels[index]
  scores <- scores[index]
  true_positive <- cumsum(labels == 1L)
  end_of_tie <- c(which(scores[-length(scores)] != scores[-1L]), length(scores))
  precision <- true_positive[end_of_tie] / end_of_tie
  recall <- true_positive[end_of_tie] / n_positive
  trapezoid_area(c(0, recall), c(1, precision))
}

format_mean_sd <- function(mean_value, sd_value, digits = 4L) {
  paste0(
    formatC(mean_value, digits = digits, format = "f"),
    " +/- ",
    formatC(sd_value, digits = digits, format = "f")
  )
}

format_vector_mean_sd <- function(x, digits = 4L) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return("NA")
  if (all(is.infinite(x) & x > 0)) return("Inf")
  if (!all(is.finite(x))) return("non-finite")
  if (length(x) == 1L) return(formatC(x, digits = digits, format = "g"))
  paste0(
    formatC(mean(x), digits = digits, format = "g"),
    " +/- ",
    formatC(sd(x), digits = digits, format = "g")
  )
}
