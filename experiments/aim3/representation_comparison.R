source("experiments/utils.R")
source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")
source("experiments/aim3/embedding.R")

aim3_type_metrics <- function(labels, scores, types) {
  result <- c(roc_distributional = NA_real_, pr_distributional = NA_real_,
              roc_structural = NA_real_, pr_structural = NA_real_)
  if (is.null(types)) return(result)
  for (type in c("distributional", "structural")) {
    keep <- types %in% c("normal", type)
    target <- as.integer(types[keep] == type)
    if (length(unique(target)) != 2L) next
    result[paste0("roc_", type)] <- roc_auc_score(target, scores[keep])
    result[paste0("pr_", type)] <- pr_auc_score(target, scores[keep])
  }
  result
}

aim3_variable_reference_columns <- function(x_reference) {
  lower <- apply(x_reference, 2L, min)
  upper <- apply(x_reference, 2L, max)

  is.finite(lower) &
    is.finite(upper) &
    upper > lower
}

aim3_filter_reference_constants <- function(
    x_train,
    x_reference,
    x_eval,
    representation
) {
  keep <- aim3_variable_reference_columns(x_reference)

  if (!any(keep)) {
    stop(
      representation,
      " has no variable reference coordinates.",
      call. = FALSE
    )
  }

  removed <- sum(!keep)

  if (removed > 0L) {
    message(
      representation,
      ": removed ",
      removed,
      " constant reference coordinates; using ",
      sum(keep),
      "."
    )
  }

  list(
    keep = keep,
    train = x_train[, keep, drop = FALSE],
    reference = x_reference[, keep, drop = FALSE],
    eval = x_eval[, keep, drop = FALSE]
  )
}


aim3_make_reference <- function(x_reference, seed, settings) {
  make_reference(
    x_reference,
    n_reference = settings$n_reference,
    n_mc_repetitions = settings$n_mc_repetitions,
    seed = seed
  )
}

aim3_evaluate_representation <- function(
    representation,
    x_train,
    x_reference,
    x_eval,
    x_label,
    labels_label,
    types_label,
    seed,
    settings
) {
  start <- proc.time()[[3L]]

  filtered <- aim3_filter_reference_constants(
    x_train,
    x_reference,
    x_eval,
    representation
  )

  x_train <- filtered$train
  x_reference <- filtered$reference
  x_eval <- filtered$eval
  x_label <- x_label[, filtered$keep, drop = FALSE]

  model <- fit_ocsvm(
    x_train,
    nu = settings$detector$nu,
    gamma = 1 / ncol(x_train)
  )

  if (!model$converged) {
    warning(representation, " OCSVM solver did not converge.", call. = FALSE)
  }

  reference <- aim3_make_reference(
    x_reference,
    seed + 100L,
    settings
  )

  mv <- aumvc(
    x_eval,
    reference,
    function(x) score_ocsvm(model, x),
    score_direction = "anomaly",
    alpha_grid = settings$aumvc_alpha_grid
  )
  label_scores <- score_ocsvm(model, x_label)
  type_metrics <- aim3_type_metrics(labels_label, label_scores, types_label)

  data.frame(
    representation = representation,
    aumvc = mv$aumvc,
    aumvc_normalized = mv$aumvc_normalized,
    aumvc_mc_se = mv$aumvc_mc_se,
    aumvc_normalized_mc_se = mv$aumvc_normalized_mc_se,
    zero_occupancy = mean(mv$mv_curve$volume_normalized == 0),
    low_occupancy = mean(mv$mv_curve$volume_normalized < 10 / settings$n_reference),
    box_log_volume = reference$box$log_volume,
    roc_auc = roc_auc_score(labels_label, label_scores),
    pr_auc = pr_auc_score(labels_label, label_scores),
    roc_distributional = unname(type_metrics["roc_distributional"]),
    pr_distributional = unname(type_metrics["pr_distributional"]),
    roc_structural = unname(type_metrics["roc_structural"]),
    pr_structural = unname(type_metrics["pr_structural"]),
    runtime_seconds = proc.time()[[3L]] - start
  )
}

source("experiments/aim3/goix_subsampling.R")

aim3_compare_representations <- function(case, seed, settings) {
  if (is.null(case$labels) || length(case$labels) != nrow(case$x)) {
    stop("Aim 3 requires binary labels for held-out evaluation.", call. = FALSE)
  }
  if (anyNA(case$labels) || !all(case$labels %in% c(0L, 1L))) {
    stop("Aim 3 labels must be binary 0/1.", call. = FALSE)
  }
  for (attempt in seq_len(100L)) {
    split <- make_splits(
      nrow(case$x), settings$split_counts,
      names(settings$split_counts), seed + attempt - 1L,
      require_full = TRUE
    )
    if (length(unique(case$labels[split$label_eval])) == 2L) break
  }
  if (length(unique(case$labels[split$label_eval])) != 2L) {
    stop("Could not obtain both classes in the label evaluation split.", call. = FALSE)
  }
  labels_label <- as.integer(case$labels[split$label_eval])
  types_label <- if (is.null(case$outlier_type)) NULL else {
    case$outlier_type[split$label_eval]
  }

  x <- as.matrix(case$x)

  keep <- aim3_variable_reference_columns(
    x[split$reference, , drop = FALSE]
  )

  if (!any(keep)) {
    stop(
      "The reference split has no variable coordinates.",
      call. = FALSE
    )
  }

  removed <- sum(!keep)

  if (removed > 0L) {
    message(
      "Input: removed ",
      removed,
      " constant reference coordinates; using ",
      sum(keep),
      "."
    )
  }

  x <- x[, keep, drop = FALSE]

  minimum_dim <- max(
    as.integer(settings$goix_subsampling$subset_dim),
    as.integer(case$intrinsic_dim)
  )

  if (ncol(x) < minimum_dim) {
    stop(
      "Too few variable coordinates remain for Aim 3.",
      call. = FALSE
    )
  }

  standardizer <- fit_standardizer(
    x[split$embedding, , drop = FALSE]
  )

  x_embedding <- apply_standardizer(
    x[split$embedding, , drop = FALSE],
    standardizer
  )
  x_train <- apply_standardizer(
    x[split$detector_train, , drop = FALSE],
    standardizer
  )
  x_reference <- apply_standardizer(
    x[split$reference, , drop = FALSE],
    standardizer
  )
  x_eval <- apply_standardizer(
    x[split$evaluation, , drop = FALSE],
    standardizer
  )
  x_label <- apply_standardizer(
    x[split$label_eval, , drop = FALSE], standardizer
  )

  ambient <- aim3_evaluate_representation(
    "ambient",
    x_train,
    x_reference,
    x_eval,
    x_label,
    labels_label,
    types_label,
    seed + 1000L,
    settings
  )

  goix <- aim3_evaluate_goix_subsampling(
    x_train,
    x_reference,
    x_eval,
    x_label,
    labels_label,
    types_label,
    seed + 2000L,
    settings
  )

  embedding_start <- proc.time()[[3L]]
  embedding <- aim3_fit_embedding(
    x_embedding,
    case$intrinsic_dim,
    settings$embedding$method
  )
  embedding_standardizer <- fit_standardizer(embedding$points)
  project <- function(x) {
    apply_standardizer(aim3_project_embedding(embedding, x), embedding_standardizer)
  }

  embedding_result <- aim3_evaluate_representation(
    embedding$method,
    project(x_train),
    project(x_reference),
    project(x_eval),
    project(x_label),
    labels_label,
    types_label,
    seed + 3000L,
    settings
  )
  # Include fitting and out-of-sample projections in embedding runtime.
  embedding_result$runtime_seconds <- proc.time()[[3L]] - embedding_start

  rbind(ambient, goix, embedding_result)
}
