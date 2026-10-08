source("experiments/aim3/utils.R")
source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")
source("experiments/aim3/embedding.R")

aim3_type_metrics <- function(labels, scores, types) {
  result <- c(
    roc_distributional = NA_real_,
    pr_distributional = NA_real_,
    roc_structural = NA_real_,
    pr_structural = NA_real_
  )
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
  is.finite(lower) & is.finite(upper) & upper > lower
}

aim3_supported_columns <- function(x, minimum_fraction) {
  if (
    length(minimum_fraction) != 1L ||
    !is.finite(minimum_fraction) ||
    minimum_fraction <= 0 ||
    minimum_fraction >= 0.5
  ) {
    stop("minimum_feature_fraction must be between 0 and 0.5.", call. = FALSE)
  }

  minimum_count <- max(2L, ceiling(nrow(x) * minimum_fraction))
  vapply(seq_len(ncol(x)), function(column) {
    values <- x[, column]
    if (any(!is.finite(values))) return(FALSE)
    counts <- tabulate(match(values, unique(values)))
    length(counts) > 1L && nrow(x) - max(counts) >= minimum_count
  }, logical(1))
}

aim3_filter_reference_constants <- function(
    x_train,
    x_reference,
    x_eval,
    x_label,
    representation
) {
  keep <- aim3_variable_reference_columns(x_reference)
  if (!all(keep)) {
    stop(
      representation,
      " has a constant reference coordinate.",
      call. = FALSE
    )
  }
  list(
    train = x_train,
    reference = x_reference,
    eval = x_eval,
    label = x_label
  )
}

aim3_make_reference <- function(x_reference, seed, settings) {
  make_reference(
    x_reference,
    n_reference = settings$n_reference,
    n_mc_repetitions = settings$n_mc_repetitions,
    seed = seed,
    chunk_size = settings$reference_chunk_size
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
    x_label,
    representation
  )
  model <- fit_ocsvm(
    filtered$train,
    nu = settings$detector$nu,
    gamma = 1 / ncol(filtered$train)
  )
  if (!model$converged) {
    stop(representation, " OCSVM solver did not converge.", call. = FALSE)
  }

  reference <- aim3_make_reference(
    filtered$reference,
    seed + 100L,
    settings
  )
  score_fun <- function(x) score_ocsvm(model, x)
  mv <- aumvc(
    filtered$eval,
    reference,
    score_fun,
    score_direction = "anomaly",
    alpha_grid = settings$aumvc_alpha_grid
  )
  label_scores <- score_fun(filtered$label)
  type_metrics <- aim3_type_metrics(labels_label, label_scores, types_label)

  data.frame(
    representation = representation,
    aumvc = mv$aumvc,
    aumvc_normalized = mv$aumvc_normalized,
    aumvc_mc_se = mv$aumvc_mc_se,
    aumvc_normalized_mc_se = mv$aumvc_normalized_mc_se,
    aumvc_relative_mc_se = mv$aumvc_relative_mc_se,
    aumvc_subset_se = NA_real_,
    minimum_hit_count = mv$minimum_hit_count,
    zero_occupancy = mv$zero_occupancy,
    low_occupancy = mv$low_occupancy,
    reliable = mv$minimum_hit_count >= settings$minimum_hit_count &&
      is.finite(mv$aumvc_relative_mc_se) &&
      mv$aumvc_relative_mc_se <= settings$maximum_relative_mc_se,
    box_log_volume = reference$box$log_volume,
    roc_auc = roc_auc_score(labels_label, label_scores),
    pr_auc = pr_auc_score(labels_label, label_scores),
    roc_distributional = unname(type_metrics["roc_distributional"]),
    pr_distributional = unname(type_metrics["pr_distributional"]),
    roc_structural = unname(type_metrics["roc_structural"]),
    pr_structural = unname(type_metrics["pr_structural"]),
    ocsvm_converged = model$converged,
    ocsvm_kkt_gap = model$kkt_gap,
    embedding_oos_stress = NA_real_,
    runtime_seconds = proc.time()[[3L]] - start
  )
}

source("experiments/aim3/goix_subsampling.R")

aim3_compare_representations <- function(case, seed, settings) {
  if (
    is.null(case$labels) || length(case$labels) != nrow(case$x) ||
    anyNA(case$labels) || !all(case$labels %in% c(0L, 1L))
  ) {
    stop("Aim 3 requires binary labels for held-out evaluation.", call. = FALSE)
  }

  split <- aim3_make_splits(
    case$labels,
    settings$split_counts,
    names(settings$split_counts),
    seed,
    if (is.null(case$outlier_type)) case$labels else case$outlier_type
  )
  labels_label <- as.integer(case$labels[split$label_eval])
  types_label <- if (is.null(case$outlier_type)) {
    NULL
  } else {
    case$outlier_type[split$label_eval]
  }

  x <- validate_matrix(case$x, "case$x")
  keep <- aim3_supported_columns(
    x[split$embedding, , drop = FALSE],
    settings$minimum_feature_fraction
  )
  if (!any(keep)) {
    stop("The embedding split has no supported coordinates.", call. = FALSE)
  }
  x <- x[, keep, drop = FALSE]
  minimum_dim <- max(
    as.integer(settings$goix_subsampling$subset_dim),
    as.integer(case$intrinsic_dim)
  )
  if (ncol(x) < minimum_dim) {
    stop("Too few variable coordinates remain for Aim 3.", call. = FALSE)
  }

  standardizer <- fit_standardizer(x[split$embedding, , drop = FALSE])
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
    x[split$label_eval, , drop = FALSE],
    standardizer
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

  methods <- unique(as.character(settings$embedding$methods))
  if (
    length(methods) < 1L ||
    any(!methods %in% c("mds", "isomap"))
  ) {
    stop("Select at least one valid embedding method.", call. = FALSE)
  }

  embedding_results <- lapply(methods, function(method) {
    embedding_start <- proc.time()[[3L]]
    embedding <- aim3_fit_embedding(
      x_embedding,
      case$intrinsic_dim,
      method,
      settings$embedding
    )
    oos_stress <- aim3_embedding_oos_stress(
      embedding,
      x_embedding,
      x_train
    )
    embedding_standardizer <- fit_standardizer(embedding$points)
    project <- function(value) {
      coordinates <- aim3_project_embedding(embedding, value)
      apply_standardizer(coordinates, embedding_standardizer)
    }
    method_offset <- 2000L + 1000L * match(
      method,
      c("mds", "isomap")
    )
    result <- aim3_evaluate_representation(
      embedding$method,
      project(x_train),
      project(x_reference),
      project(x_eval),
      project(x_label),
      labels_label,
      types_label,
      seed + method_offset,
      settings
    )
    result$embedding_oos_stress <- oos_stress
    result$runtime_seconds <- proc.time()[[3L]] - embedding_start
    result
  })

  results <- do.call(rbind, c(list(ambient, goix), embedding_results))
  attr(results, "effective_dim") <- ncol(x)
  results
}
