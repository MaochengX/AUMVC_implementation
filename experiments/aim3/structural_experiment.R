source("experiments/aim3/synthetic_experiment.R")
source("experiments/aim3/structural_generator.R")

structural_basis <- structural_basis_v3

structural_detector_scores <- function(prepared, seed, settings) {
  scales <- gamma_scales(settings)
  methods <- c("ambient", "goix_subsampling", names(prepared$embeddings))
  base <- prepared$prepared
  output <- vector("list", length(methods) * length(scales))
  index <- 0L

  for (method in methods) {
    if (method == "ambient" || method == "goix_subsampling") {
      train <- base$train
      label <- base$label
      evaluation <- base$eval
    } else {
      embedding <- prepared$embeddings[[method]]
      train <- embedding$train
      label <- embedding$label
      evaluation <- embedding$eval
    }

    if (method == "goix_subsampling") {
      subsets <- goix_subsets(ncol(train), seed + 2000L, settings$goix_subsampling)
    }

    for (scale in scales) {
      if (method == "goix_subsampling") {
        scores <- numeric(nrow(label))
        for (columns in subsets) {
          model <- fit_ocsvm(
            train[, columns, drop = FALSE],
            nu = settings$detector$nu,
            gamma = scale / length(columns)
          )
          evaluation_scores <- score_ocsvm(model, evaluation[, columns, drop = FALSE])
          label_scores <- score_ocsvm(model, label[, columns, drop = FALSE])
          scores <- scores + percentile_scores(evaluation_scores, label_scores) / length(subsets)
        }
      } else {
        model <- fit_ocsvm(train, nu = settings$detector$nu,
                           gamma = scale / ncol(train))
        scores <- score_ocsvm(model, label)
      }
      subtype <- compute_type_metrics(base$labels_label, scores, base$types_label)
      index <- index + 1L
      output[[index]] <- data.frame(
        representation = method,
        gamma_scale = scale,
        roc_auc = roc_auc_score(base$labels_label, scores),
        pr_auc = pr_auc_score(base$labels_label, scores),
        roc_distributional = unname(subtype["roc_distributional"]),
        pr_distributional = unname(subtype["pr_distributional"]),
        roc_structural = unname(subtype["roc_structural"]),
        pr_structural = unname(subtype["pr_structural"])
      )
    }
  }
  do.call(rbind, output)
}

structural_geometry <- function(case) {
  reference <- case$x[1:400, , drop = FALSE]
  normal <- case$x[401:450, , drop = FALSE]
  structural <- case$x[case$outlier_type == "structural", , drop = FALSE]
  labels <- c(rep(0L, nrow(normal)), rep(1L, nrow(structural)))
  evaluation <- rbind(normal, structural)
  squared <- outer(rowSums(evaluation^2), rowSums(reference^2), "+") -
    2 * tcrossprod(evaluation, reference)
  distances <- sqrt(pmax(0, apply(squared, 1L, min)))
  magnitude <- sqrt(rowSums(evaluation^2))
  magnitude_auc <- roc_auc_score(labels, magnitude)
  data.frame(
    nearest_normal_roc = roc_auc_score(labels, distances),
    magnitude_separation = max(magnitude_auc, 1 - magnitude_auc)
  )
}

structural_task <- function(task) {
  row <- task$row
  settings <- task$settings
  seed <- synthetic_seed(row, task$run, settings$synthetic, settings$seed)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  case <- generate_synthetic(row$intrinsic_dim, row$original_dim, row$snr,
                             row$truth, seed, settings$synthetic)
  metadata <- data.frame(
    run = task$run, truth = row$truth, intrinsic_dim = row$intrinsic_dim,
    snr = row$snr, original_dim = row$original_dim
  )
  geometry <- cbind(metadata, structural_geometry(case))
  prepared <- prepare_representations(case, seed + 50000L, settings)
  scores <- structural_detector_scores(prepared, seed + 50000L, settings)
  list(geometry = geometry, scores = attach_metadata(metadata, scores))
}

run_structural_experiment <- function(settings = AIM3_SETTINGS) {
  grid <- synthetic_grid(settings$synthetic)
  runs <- validate_runs(settings$n_runs)
  directory <- "experiments/aim3/result/synthetic/structural_final"
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  parallel <- isTRUE(settings$parallel$enabled) && settings$parallel$workers > 1L
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Install future and future.apply.", call. = FALSE)
    }
    previous <- future::plan()
    future::plan(future::multisession, workers = settings$parallel$workers)
    on.exit(future::plan(previous), add = TRUE)
  }

  scored <- vector("list", runs)
  checked <- vector("list", runs)
  for (run in seq_len(runs)) {
    jobs <- lapply(seq_len(nrow(grid)), function(i) {
      list(row = grid[i, , drop = FALSE], run = run, settings = settings)
    })
    cases <- if (parallel) {
      future.apply::future_lapply(jobs, structural_task, future.seed = TRUE)
    } else {
      lapply(jobs, structural_task)
    }
    scored[[run]] <- do.call(rbind, lapply(cases, `[[`, "scores"))
    checked[[run]] <- do.call(rbind, lapply(cases, `[[`, "geometry"))
    write_report(do.call(rbind, scored[seq_len(run)]),
                 file.path(directory, "syn_results.csv"), settings$output_digits)
    write_report(do.call(rbind, checked[seq_len(run)]),
                 file.path(directory, "geometry.csv"), settings$output_digits)
    cat(sprintf("structural runs: %d/%d\n", run, runs))
  }

  results <- do.call(rbind, scored)
  geometry <- do.call(rbind, checked)
  summary <- lapply(group_rows(results, c("truth", "representation")), function(index) {
    x <- results[index, , drop = FALSE]
    row <- x[1L, c("truth", "representation"), drop = FALSE]
    for (metric in c("roc_structural", "pr_structural",
                     "roc_distributional", "pr_distributional")) {
      row[[metric]] <- mean_finite(x[[metric]])
    }
    row
  })
  summary <- do.call(rbind, summary)
  rownames(summary) <- NULL
  write_report(summary, file.path(directory, "summary.csv"), settings$output_digits)
  print(summary, row.names = FALSE)
  invisible(list(results = results, geometry = geometry, summary = summary))
}

if (sys.nframe() == 0L) run_structural_experiment()
