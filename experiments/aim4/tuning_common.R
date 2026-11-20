source("experiments/aim3/real_experiment.R")

AIM4B_SETTINGS <- list(
  seed = 1234L,
  n_runs = as.integer(Sys.getenv("AIM4B_RUNS", "2000")),
  datasets = trimws(strsplit(Sys.getenv("AIM4B_DATASETS", "ecg5000,fashion_mnist"), ",", fixed = TRUE)[[1L]]),
  embedding_dims = c(ecg5000 = 5L, fashion_mnist = 5L),
  n_reference = as.integer(Sys.getenv("AIM4B_MC", "100000")),
  n_mc_repetitions = 1L,
  reference_chunk_size = 5000L,
  alpha_grid = seq(0.9, 0.999, by = 0.0001),
  minimum_hit_count = 10L,
  maximum_relative_mc_se = 0.10,
  minimum_feature_fraction = 0.10,
  tie_tolerance = sqrt(.Machine$double.eps)
)

six_way_splits <- function(labels, seed) {
  names <- c("embedding", "detector_train", "reference", "evaluation", "validation", "test")
  if (length(unique(labels)) != 2L || any(!labels %in% c(0L, 1L))) {
    stop("Expected binary anomaly labels.", call. = FALSE)
  }
  set.seed(seed)
  splits <- setNames(rep(list(integer()), length(names)), names)
  for (class_label in 0:1) {
    rows <- sample(which(labels == class_label))
    if (length(rows) < 6L) stop("Need at least six observations per class.", call. = FALSE)
    allocation <- sample(rep(seq_along(names), length.out = length(rows)))
    for (i in seq_along(names)) {
      splits[[i]] <- c(splits[[i]], rows[allocation == i])
    }
  }
  splits <- lapply(splits, sample)
  if (length(unique(unlist(splits))) != length(labels) ||
      !setequal(unlist(splits), seq_along(labels))) {
    stop("Overlapping or incomplete data splits.", call. = FALSE)
  }
  splits
}

prepare_aim4b_case <- function(dataset, run, loaded, registry, settings) {
  index <- match(dataset, names(settings$embedding_dims))
  seed <- settings$seed + 100000L * run + 10000L * index
  sampled <- registry[[dataset]]$sample(loaded, seed, AIM3_SETTINGS$real[[dataset]])
  labels <- as.integer(sampled$labels)
  split <- six_way_splits(labels, seed + 50000L)
  raw <- as.matrix(sampled$x)
  keep <- supported_columns(raw[split$embedding, , drop = FALSE],
                            settings$minimum_feature_fraction)
  raw <- raw[, keep, drop = FALSE]
  dimension <- as.integer(settings$embedding_dims[[dataset]])
  if (ncol(raw) < dimension) stop("Too few supported input coordinates.", call. = FALSE)

  original_standardizer <- fit_standardizer(raw[split$embedding, , drop = FALSE])
  scaled_embedding <- apply_standardizer(raw[split$embedding, , drop = FALSE], original_standardizer)
  model <- fit_embedding(scaled_embedding, dimension, "mds", AIM3_SETTINGS$embedding)
  embedded_standardizer <- fit_standardizer(model$points)
  project <- function(which_split) {
    scaled <- apply_standardizer(raw[split[[which_split]], , drop = FALSE], original_standardizer)
    apply_standardizer(project_embedding(model, scaled), embedded_standardizer)
  }

  split_names <- names(split)[-1L]
  data <- setNames(lapply(split_names, project), split_names)
  reference <- make_reference(
    data$reference,
    n_reference = settings$n_reference,
    n_mc_repetitions = settings$n_mc_repetitions,
    seed = seed + 80000L,
    chunk_size = settings$reference_chunk_size
  )
  c(list(dataset = dataset, run = run, seed = seed,
         original_dim = ncol(sampled$x), embedding_dim = dimension,
         n_train = length(split$detector_train),
         n_validation = length(split$validation), n_test = length(split$test),
         validation_labels = labels[split$validation],
         test_labels = labels[split$test], reference_mc = reference), data)
}

score_aim4b_candidates <- function(prepared, parameters, evaluation_scores,
                                   validation_scores, test_scores, reference_scores,
                                   settings) {
  parameters <- as.numeric(parameters)
  grids <- list(evaluation_scores, validation_scores, test_scores, reference_scores)
  expected <- c(nrow(prepared$evaluation), nrow(prepared$validation),
                nrow(prepared$test), settings$n_reference * settings$n_mc_repetitions)
  if (any(!vapply(seq_along(grids), function(i) {
    identical(dim(grids[[i]]), c(expected[i], length(parameters))) &&
      all(is.finite(grids[[i]]))
  }, logical(1)))) stop("Invalid candidate score matrices.", call. = FALSE)

  rows <- lapply(seq_along(parameters), function(i) {
    reference_scores_i <- reference_scores[, i]
    mv <- aumvc_scored(
      evaluation_scores[, i], reference_scores_i,
      prepared$reference_mc, score_direction = "anomaly",
      alpha_grid = settings$alpha_grid
    )
    degeneracy <- mv_degeneracy(mv, reference_scores_i)
    reliable <- aumvc_reliable(
      mv$aumvc_normalized, mv$aumvc_relative_mc_se,
      mv$minimum_hit_count, settings, mv, reference_scores_i
    )
    data.frame(
      dataset = prepared$dataset, run = prepared$run,
      original_dim = prepared$original_dim, embedding_dim = prepared$embedding_dim,
      parameter = parameters[i],
      aumvc = mv$aumvc_normalized,
      relative_mc_se = mv$aumvc_relative_mc_se,
      min_hit_count = mv$minimum_hit_count,
      zero_occupancy = mv$zero_occupancy,
      full_box = unname(degeneracy["full_box"]),
      threshold_collapsed = unname(degeneracy["threshold_collapsed"]),
      reference_collapsed = unname(degeneracy["reference_collapsed"]),
      reliable = reliable,
      validation_roc = roc_auc_score(prepared$validation_labels, validation_scores[, i]),
      validation_pr = pr_auc_score(prepared$validation_labels, validation_scores[, i]),
      test_roc = roc_auc_score(prepared$test_labels, test_scores[, i]),
      test_pr = pr_auc_score(prepared$test_labels, test_scores[, i])
    )
  })
  do.call(rbind, rows)
}

choose_label_setting <- function(values, parameters, default) {
  best <- max(values)
  tied <- which(abs(values - best) <= 1e-12)
  if (length(tied) == 1L) return(tied)
  tied[which.min(abs(log(parameters[tied] / default)))]
}

select_aim4b <- function(candidates, default, most_local, elapsed, settings) {
  params <- candidates$parameter
  default_index <- which(params == default)
  if (length(default_index) != 1L) stop("The default must belong to the grid.", call. = FALSE)
  eligible <- which(candidates$reliable & is.finite(candidates$aumvc))
  if (!length(eligible)) {
    status <- "all_unreliable"
    chosen <- default_index
  } else {
    minimum <- min(candidates$aumvc[eligible])
    tied <- eligible[vapply(eligible, function(i) {
      near_equal(candidates$aumvc[i], minimum, settings$tie_tolerance)
    }, logical(1))]
    status <- if (length(tied) == 1L) "selected" else "tie"
    chosen <- if (status == "selected") tied else default_index
  }
  roc_index <- choose_label_setting(candidates$validation_roc, params, default)
  pr_index <- choose_label_setting(candidates$validation_pr, params, default)
  data.frame(
    dataset = candidates$dataset[1L], run = candidates$run[1L],
    selection_status = status,
    n_reliable = length(eligible),
    aumvc_choice = params[chosen],
    roc_label_choice = params[roc_index],
    pr_label_choice = params[pr_index],
    default = default,
    chose_most_local = status == "selected" && params[chosen] == most_local,
    more_local_than_roc = if (status == "selected") {
      if (most_local == max(params)) params[chosen] > params[roc_index] else params[chosen] < params[roc_index]
    } else NA,
    test_roc_aumvc = candidates$test_roc[chosen],
    test_pr_aumvc = candidates$test_pr[chosen],
    roc_loss_vs_labels = candidates$test_roc[roc_index] - candidates$test_roc[chosen],
    pr_loss_vs_labels = candidates$test_pr[pr_index] - candidates$test_pr[chosen],
    roc_gain_vs_default = candidates$test_roc[chosen] - candidates$test_roc[default_index],
    pr_gain_vs_default = candidates$test_pr[chosen] - candidates$test_pr[default_index],
    runtime_sec = elapsed
  )
}

run_aim4b <- function(detector, settings, score_grid, parameters, default, most_local) {
  n_runs <- validate_runs(settings$n_runs)
  if (length(settings$datasets) < 1L ||
      any(!settings$datasets %in% names(settings$embedding_dims)) ||
      anyDuplicated(settings$datasets)) stop("Invalid dataset selection.")
  if (settings$n_reference < 1000L || !default %in% parameters ||
      anyDuplicated(parameters)) stop("Invalid experiment settings.")
  save_every <- max(1L, as.integer(Sys.getenv("AIM4B_SAVE_EVERY", "25")))
  registry <- real_registry()
  output <- file.path("experiments", "aim4", paste0(detector, "_tune"))
  for (dataset in settings$datasets) {
    loaded <- registry[[dataset]]$load(AIM3_SETTINGS$real[[dataset]]$data_dir)
    folder <- file.path(output, dataset)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    candidates <- vector("list", n_runs)
    selections <- vector("list", n_runs)
    for (run in seq_len(n_runs)) {
      start <- proc.time()[[3L]]
      prepared <- prepare_aim4b_case(dataset, run, loaded, registry, settings)
      candidates[[run]] <- score_grid(prepared, parameters, settings)
      selections[[run]] <- select_aim4b(
        candidates[[run]], default, most_local, proc.time()[[3L]] - start, settings
      )
      if (run %% save_every == 0L || run == n_runs) {
        utils::write.csv(do.call(rbind, candidates[seq_len(run)]),
                         file.path(folder, "candidates.csv"), row.names = FALSE)
        utils::write.csv(do.call(rbind, selections[seq_len(run)]),
                         file.path(folder, "selections.csv"), row.names = FALSE)
        cat(sprintf("%s / %s: %d/%d runs\n", detector, dataset, run, n_runs))
      }
    }
  }
  invisible(output)
}
