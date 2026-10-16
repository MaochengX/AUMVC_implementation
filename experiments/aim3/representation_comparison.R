source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("detectors/ocsvm.R")
source("experiments/aim3/embedding.R")

compute_type_metrics <- function(labels, scores, types) {
  result <- c(
    roc_distributional = NA_real_, pr_distributional = NA_real_,
    roc_structural = NA_real_, pr_structural = NA_real_
  )
  if (is.null(types)) return(result)
  for (type in c("distributional", "structural")) {
    keep <- types %in% c("normal", type)
    target <- as.integer(types[keep] == type)
    if (length(unique(target)) == 2L) {
      result[paste0("roc_", type)] <- roc_auc_score(target, scores[keep])
      result[paste0("pr_", type)] <- pr_auc_score(target, scores[keep])
    }
  }
  result
}

supported_columns <- function(x, minimum_fraction) {
  if (!is.finite(minimum_fraction) || minimum_fraction <= 0 || minimum_fraction >= 0.5) {
    stop("minimum_feature_fraction must be between 0 and 0.5.", call. = FALSE)
  }
  minimum_count <- max(2L, ceiling(nrow(x) * minimum_fraction))
  vapply(seq_len(ncol(x)), function(j) {
    values <- x[, j]
    if (any(!is.finite(values))) return(FALSE)
    counts <- tabulate(match(values, unique(values)))
    length(counts) > 1L && nrow(x) - max(counts) >= minimum_count
  }, logical(1))
}

make_experiment_reference <- function(x_reference, seed, settings) {
  make_reference(
    x_reference,
    n_reference = settings$n_reference,
    n_mc_repetitions = settings$n_mc_repetitions,
    seed = seed,
    chunk_size = settings$reference_chunk_size
  )
}

gamma_scales <- function(settings) {
  scales <- as.numeric(settings$bandwidth$gamma_scales)
  if (length(scales) < 2L || any(!is.finite(scales)) || any(scales <= 0) || anyDuplicated(scales)) {
    stop("gamma_scales must be unique positive values.", call. = FALSE)
  }
  scales
}

aumvc_reliable <- function(aumvc_normalized, relative_se, min_hit_count, settings) {
  is.finite(aumvc_normalized) &&
    is.finite(relative_se) &&
    min_hit_count >= settings$minimum_hit_count &&
    relative_se <= settings$maximum_relative_mc_se
}

prepare_case_data <- function(case, seed, settings) {
  x <- validate_matrix(case$x, "case$x")
  labels <- as.integer(case$labels)
  if (length(labels) != nrow(x) || anyNA(labels) || !all(labels %in% c(0L, 1L))) {
    stop("Aim 3 requires binary labels.", call. = FALSE)
  }
  embedding_dim <- as.integer(case$intrinsic_dim)
  if (length(embedding_dim) != 1L || is.na(embedding_dim) || embedding_dim < 1L) {
    stop("Aim 3 requires a positive embedding dimension.", call. = FALSE)
  }
  strata <- if (is.null(case$outlier_type)) labels else case$outlier_type
  split <- make_splits(labels, settings$split_counts, seed, strata)
  keep <- supported_columns(
    x[split$embedding, , drop = FALSE],
    settings$minimum_feature_fraction
  )
  if (!any(keep)) stop("No supported coordinates remain.", call. = FALSE)
  x <- x[, keep, drop = FALSE]
  minimum_dim <- max(as.integer(settings$goix_subsampling$subset_dim), embedding_dim)
  if (ncol(x) < minimum_dim) stop("Too few coordinates remain for Aim 3.", call. = FALSE)
  standardizer <- fit_standardizer(x[split$embedding, , drop = FALSE])
  transform <- function(index) apply_standardizer(x[index, , drop = FALSE], standardizer)
  list(
    embedding = transform(split$embedding),
    train = transform(split$detector_train),
    reference = transform(split$reference),
    eval = transform(split$evaluation),
    label = transform(split$label_eval),
    labels_label = labels[split$label_eval],
    types_label = if (is.null(case$outlier_type)) NULL else case$outlier_type[split$label_eval],
    embedding_dim = embedding_dim
  )
}

prepare_embedding_data <- function(prepared, method, settings) {
  start <- proc.time()[[3L]]
  model <- fit_embedding(
    prepared$embedding,
    prepared$embedding_dim,
    method,
    settings$embedding
  )
  stress <- embedding_oos_stress(model, prepared$train)
  standardizer <- fit_standardizer(model$points)
  project <- function(x) apply_standardizer(project_embedding(model, x), standardizer)
  list(
    representation = method,
    train = project(prepared$train),
    reference = project(prepared$reference),
    eval = project(prepared$eval),
    label = project(prepared$label),
    embedding_oos_stress = stress,
    setup_runtime = proc.time()[[3L]] - start
  )
}

prepare_representations <- function(case, seed, settings) {
  prepared <- prepare_case_data(case, seed, settings)
  methods <- unique(as.character(settings$embedding$methods))
  if (!length(methods) || any(!methods %in% c("mds", "isomap"))) {
    stop("Invalid embedding methods.", call. = FALSE)
  }
  embeddings <- setNames(
    lapply(methods, function(method) prepare_embedding_data(prepared, method, settings)),
    methods
  )
  list(prepared = prepared, embeddings = embeddings)
}

evaluate_candidate <- function(
    representation, x_train, x_eval, x_label, labels_label, types_label,
    reference, gamma_scale, settings, embedding_oos_stress = NA_real_
) {
  start <- proc.time()[[3L]]
  detector_dim <- ncol(x_train)
  model <- fit_ocsvm(
    x_train,
    nu = settings$detector$nu,
    gamma = gamma_scale / detector_dim
  )
  score_fun <- function(x) score_ocsvm(model, x)
  mv <- aumvc(
    x_eval,
    reference,
    score_fun,
    score_direction = "anomaly",
    alpha_grid = settings$aumvc_alpha_grid
  )
  label_scores <- score_fun(x_label)
  subtype <- compute_type_metrics(labels_label, label_scores, types_label)
  data.frame(
    embedding_dim = if (representation %in% c("mds", "isomap")) detector_dim else NA_integer_,
    detector_dim = detector_dim,
    representation = representation,
    gamma_scale = gamma_scale,
    aumvc_normalized = mv$aumvc_normalized,
    aumvc_relative_mc_se = mv$aumvc_relative_mc_se,
    min_hit_count = mv$minimum_hit_count,
    zero_occupancy = mv$zero_occupancy,
    low_occupancy = mv$low_occupancy,
    box_log_volume = reference$box$log_volume,
    roc_auc = roc_auc_score(labels_label, label_scores),
    pr_auc = pr_auc_score(labels_label, label_scores),
    roc_distributional = unname(subtype["roc_distributional"]),
    pr_distributional = unname(subtype["pr_distributional"]),
    roc_structural = unname(subtype["roc_structural"]),
    pr_structural = unname(subtype["pr_structural"]),
    embedding_oos_stress = embedding_oos_stress,
    candidate_runtime = proc.time()[[3L]] - start,
    aumvc_reliable = aumvc_reliable(
      mv$aumvc_normalized,
      mv$aumvc_relative_mc_se,
      mv$minimum_hit_count,
      settings
    )
  )
}

evaluate_grid <- function(
    representation, x_train, x_reference, x_eval, x_label, labels_label,
    types_label, seed, settings, embedding_oos_stress = NA_real_, setup_runtime = 0
) {
  if (any(apply(x_reference, 2L, function(x) diff(range(x))) <= 0)) {
    stop(representation, " has a constant reference coordinate.", call. = FALSE)
  }
  start <- proc.time()[[3L]]
  reference <- make_experiment_reference(x_reference, seed, settings)
  shared <- setup_runtime + proc.time()[[3L]] - start
  scales <- gamma_scales(settings)
  rows <- lapply(scales, function(scale) {
    evaluate_candidate(
      representation, x_train, x_eval, x_label, labels_label, types_label,
      reference, scale, settings, embedding_oos_stress
    )
  })
  result <- do.call(rbind, rows)
  result$runtime_sec <- result$candidate_runtime + shared / length(scales)
  result$candidate_runtime <- NULL
  result
}

goix_subsets <- function(ambient_dim, seed, settings) {
  subset_dim <- as.integer(settings$subset_dim)
  n_subsets <- as.integer(settings$n_subsets)
  if (subset_dim < 1L || n_subsets < 1L || subset_dim > ambient_dim) {
    stop("Invalid Goix settings.", call. = FALSE)
  }
  set.seed(seed)
  replicate(n_subsets, sample.int(ambient_dim, subset_dim), simplify = FALSE)
}

evaluate_goix <- function(
    x_train, x_reference, x_eval, x_label, labels_label, types_label,
    seed, settings
) {
  scales <- gamma_scales(settings)
  start <- proc.time()[[3L]]
  subsets <- goix_subsets(ncol(x_train), seed, settings$goix_subsampling)
  n_subsets <- length(subsets)
  n_contribution <- settings$n_reference * settings$n_mc_repetitions
  subset_results <- lapply(scales, function(x) vector("list", n_subsets))
  ensemble_scores <- matrix(0, nrow(x_label), length(scales))
  contributions <- matrix(0, n_contribution, length(scales))
  runtime <- numeric(length(scales))
  references <- vector("list", n_subsets)
  for (j in seq_along(subsets)) {
    columns <- subsets[[j]]
    references[[j]] <- make_experiment_reference(
      x_reference[, columns, drop = FALSE],
      seed + 1000L + j,
      settings
    )
  }
  shared <- proc.time()[[3L]] - start
  for (j in seq_along(subsets)) {
    columns <- subsets[[j]]
    train <- x_train[, columns, drop = FALSE]
    eval <- x_eval[, columns, drop = FALSE]
    label <- x_label[, columns, drop = FALSE]
    reference <- references[[j]]
    for (s in seq_along(scales)) {
      candidate_start <- proc.time()[[3L]]
      model <- fit_ocsvm(
        train,
        nu = settings$detector$nu,
        gamma = scales[s] / ncol(train)
      )
      score_fun <- function(x) score_ocsvm(model, x)
      mv <- aumvc(
        eval,
        reference,
        score_fun,
        score_direction = "anomaly",
        alpha_grid = settings$aumvc_alpha_grid
      )
      eval_scores <- score_fun(eval)
      label_scores <- score_fun(label)
      ensemble_scores[, s] <- ensemble_scores[, s] +
        percentile_scores(eval_scores, label_scores) / n_subsets
      contributions[, s] <- contributions[, s] + mv$mc_contribution / n_subsets
      subset_results[[s]][[j]] <- data.frame(
        aumvc_normalized = mv$aumvc_normalized,
        aumvc_relative_mc_se = mv$aumvc_relative_mc_se,
        min_hit_count = mv$minimum_hit_count,
        zero_occupancy = mv$zero_occupancy,
        low_occupancy = mv$low_occupancy,
        box_log_volume = reference$box$log_volume
      )
      runtime[s] <- runtime[s] + proc.time()[[3L]] - candidate_start
    }
  }
  rows <- lapply(seq_along(scales), function(s) {
    subset <- do.call(rbind, subset_results[[s]])
    normalized <- mean(subset$aumvc_normalized)
    mc_se <- sd(contributions[, s]) / sqrt(nrow(contributions))
    relative_se <- if (is.finite(normalized) && normalized > 0) mc_se / normalized else NA_real_
    scores <- ensemble_scores[, s]
    subtype <- compute_type_metrics(labels_label, scores, types_label)
    reliable <- is.finite(normalized) &&
      is.finite(relative_se) &&
      relative_se <= settings$maximum_relative_mc_se &&
      all(is.finite(subset$aumvc_normalized)) &&
      min(subset$min_hit_count) >= settings$minimum_hit_count &&
      all(is.finite(subset$aumvc_relative_mc_se)) &&
      max(subset$aumvc_relative_mc_se) <= settings$maximum_relative_mc_se
    data.frame(
      embedding_dim = NA_integer_,
      detector_dim = as.integer(settings$goix_subsampling$subset_dim),
      representation = "goix_subsampling",
      gamma_scale = scales[s],
      aumvc_normalized = normalized,
      aumvc_relative_mc_se = relative_se,
      min_hit_count = min(subset$min_hit_count),
      zero_occupancy = mean(subset$zero_occupancy),
      low_occupancy = mean(subset$low_occupancy),
      box_log_volume = mean(subset$box_log_volume),
      roc_auc = roc_auc_score(labels_label, scores),
      pr_auc = pr_auc_score(labels_label, scores),
      roc_distributional = unname(subtype["roc_distributional"]),
      pr_distributional = unname(subtype["pr_distributional"]),
      roc_structural = unname(subtype["roc_structural"]),
      pr_structural = unname(subtype["pr_structural"]),
      embedding_oos_stress = NA_real_,
      runtime_sec = runtime[s] + shared / length(scales),
      aumvc_reliable = reliable
    )
  })
  do.call(rbind, rows)
}

near_equal <- function(a, b, tolerance) {
  is.finite(a) && is.finite(b) && abs(a - b) <= tolerance * max(1, abs(a), abs(b))
}

select_gamma <- function(part, eligible, tolerance) {
  index <- which(eligible & is.finite(part$aumvc_normalized))
  if (!length(index)) return(NA_real_)
  values <- part$aumvc_normalized[index]
  best <- min(values)
  tied <- index[vapply(values, function(x) near_equal(x, best, tolerance), logical(1))]
  if (length(tied) == 1L) part$gamma_scale[tied] else NA_real_
}

compare_bandwidths <- function(part, settings) {
  scales <- gamma_scales(settings)
  if (nrow(part) != length(scales) || anyDuplicated(part$gamma_scale) || !setequal(part$gamma_scale, scales)) {
    stop("Incomplete bandwidth grid.", call. = FALSE)
  }
  tolerance <- settings$bandwidth$tie_tolerance
  pairs <- combn(seq_len(nrow(part)), 2L)
  first <- pairs[1L, ]
  second <- pairs[2L, ]
  roc_order <- sign(part$roc_auc[first] - part$roc_auc[second])
  pr_order <- sign(part$pr_auc[first] - part$pr_auc[second])
  consensus <- is.finite(roc_order) & is.finite(pr_order) &
    roc_order != 0 & pr_order != 0 & roc_order == pr_order
  aumvc_finite <- is.finite(part$aumvc_normalized[first]) &
    is.finite(part$aumvc_normalized[second])
  tied <- vapply(seq_along(first), function(i) {
    consensus[i] && near_equal(
      part$aumvc_normalized[first[i]],
      part$aumvc_normalized[second[i]],
      tolerance
    )
  }, logical(1))
  aumvc_order <- sign(part$aumvc_normalized[second] - part$aumvc_normalized[first])
  agree <- consensus & aumvc_finite & !tied & aumvc_order == roc_order
  reliable_pair <- consensus & aumvc_finite &
    part$aumvc_reliable[first] & part$aumvc_reliable[second]
  total <- sum(consensus)
  reliable_total <- sum(reliable_pair)
  data.frame(
    embedding_dim = part$embedding_dim[1L],
    detector_dim = part$detector_dim[1L],
    representation = part$representation[1L],
    compared_pairs_count = total,
    unreliable_pairs_rate = if (total) sum(consensus & !reliable_pair) / total else NA_real_,
    selected_gamma_scale = select_gamma(part, rep(TRUE, nrow(part)), tolerance),
    agreement_pairs_rank_rate = if (total) sum(agree) / total else NA_real_,
    selected_gamma_scale_reliable = select_gamma(part, part$aumvc_reliable, tolerance),
    reliable_agreement_pairs_rank_rate = if (reliable_total) {
      sum(agree & reliable_pair) / reliable_total
    } else {
      NA_real_
    },
    tie_rate = if (total) sum(tied) / total else NA_real_
  )
}

evaluate_case <- function(case, seed, settings) {
  prepared_data <- prepare_representations(case, seed, settings)
  prepared <- prepared_data$prepared
  ambient <- evaluate_grid(
    "ambient", prepared$train, prepared$reference, prepared$eval, prepared$label,
    prepared$labels_label, prepared$types_label, seed + 1000L, settings
  )
  goix <- evaluate_goix(
    prepared$train, prepared$reference, prepared$eval, prepared$label,
    prepared$labels_label, prepared$types_label, seed + 2000L, settings
  )
  embeddings <- lapply(seq_along(prepared_data$embeddings), function(i) {
    item <- prepared_data$embeddings[[i]]
    evaluate_grid(
      item$representation, item$train, item$reference, item$eval, item$label,
      prepared$labels_label, prepared$types_label, seed + 3000L + 1000L * i,
      settings, item$embedding_oos_stress, item$setup_runtime
    )
  })
  results <- do.call(rbind, c(list(ambient, goix), embeddings))
  rownames(results) <- NULL
  bandwidth <- do.call(rbind, lapply(split(results, results$representation), function(part) {
    compare_bandwidths(part, settings)
  }))
  rownames(bandwidth) <- NULL
  list(results = results, bandwidth = bandwidth)
}
