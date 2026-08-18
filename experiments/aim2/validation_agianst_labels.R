source("experiments/settings.R")

source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")
source("aumvc/auemc.R")

source("detectors/ocsvm.R")
source("detectors/lof.R")
source("detectors/isolation_forest.R")

make_aim2_splits <- function(n, counts, seed) {
  counts <- as.integer(counts)

  if (length(counts) != 4L || any(counts < 1L) || sum(counts) > n) {
    stop("Invalid Aim 2 split sizes.", call. = FALSE)
  }

  set.seed(seed)
  rows <- sample.int(n, sum(counts))

  ends <- cumsum(counts)
  starts <- c(1L, head(ends, -1L) + 1L)

  split <- Map(
    function(a, b) rows[a:b],
    starts,
    ends
  )

  names(split) <- c(
    "detector_train",
    "reference",
    "aumvc",
    "label_eval"
  )

  split
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
  n_anomaly <- sum(labels == 1L)
  n_normal <- sum(labels == 0L)

  if (n_anomaly == 0L || n_normal == 0L) {
    stop("ROC-AUC requires both classes.", call. = FALSE)
  }

  ranks <- rank(scores, ties.method = "average")

  (
    sum(ranks[labels == 1L]) -
      n_anomaly * (n_anomaly + 1) / 2
  ) / (n_anomaly * n_normal)
}

pr_auc_score <- function(labels, scores) {
  n_anomaly <- sum(labels == 1L)

  if (n_anomaly == 0L) {
    stop("PR-AUC requires anomalies.", call. = FALSE)
  }

  index <- order(scores, decreasing = TRUE)
  labels <- labels[index]
  scores <- scores[index]
  true_positive <- cumsum(labels == 1L)

  end_of_tie <- c(
    which(scores[-length(scores)] != scores[-1L]),
    length(scores)
  )

  precision <- true_positive[end_of_tie] / end_of_tie
  recall <- true_positive[end_of_tie] / n_anomaly

  trapezoid_area(c(0, recall), c(1, precision))
}

compare_pair <- function(
    aumvc_difference,
    target_difference,
    tolerance = CONCORDANCE_TOLERANCE
) {
  if (abs(target_difference) <= tolerance) return(NA)
  if (abs(aumvc_difference) <= tolerance) return(FALSE)

  sign(aumvc_difference) == -sign(target_difference)
}

pairwise_concordance <- function(detector, aumvc_values, label_values) {
  pairs <- combn(seq_along(detector), 2L)

  do.call(
    rbind,
    lapply(
      seq_len(ncol(pairs)),
      function(k) {
        i <- pairs[1L, k]
        j <- pairs[2L, k]

        aumvc_difference <- aumvc_values[i] - aumvc_values[j]
        label_difference <- label_values[i] - label_values[j]

        data.frame(
          detector_1 = detector[i],
          detector_2 = detector[j],
          agree = compare_pair(
            aumvc_difference,
            label_difference
          )
        )
      }
    )
  )
}

either_concordance <- function(roc_pairs, pr_pairs) {
  if (
    !identical(roc_pairs$detector_1, pr_pairs$detector_1) ||
    !identical(roc_pairs$detector_2, pr_pairs$detector_2)
  ) {
    stop("ROC and PR pair tables do not match.", call. = FALSE)
  }

  agree <- mapply(
    function(roc, pr) {
      values <- c(roc, pr)
      comparable <- !is.na(values)

      if (!any(comparable)) return(NA)
      any(values[comparable])
    },
    roc_pairs$agree,
    pr_pairs$agree
  )

  data.frame(
    detector_1 = roc_pairs$detector_1,
    detector_2 = roc_pairs$detector_2,
    agree = agree
  )
}

fit_aim2_models <- function(x_train, seed) {
  models <- list(
    OCSVM = fit_ocsvm(
      x_train,
      nu = OCSVM_NU,
      gamma = 1 / ncol(x_train)
    ),
    LOF = fit_lof(x_train, k = LOF_K),
    Isolation_Forest = fit_isolation_forest(
      x_train,
      ntrees = IFOREST_NTREES,
      sample_size = min(
        IFOREST_SAMPLE_SIZE,
        nrow(x_train)
      ),
      seed = seed + 200L
    )
  )

  if (!models$OCSVM$converged) {
    warning("OCSVM solver did not converge.", call. = FALSE)
  }

  models
}

make_score_functions <- function(models) {
  list(
    OCSVM = function(x) score_ocsvm(models$OCSVM, x),
    LOF = function(x) score_lof(models$LOF, x),
    Isolation_Forest = function(x) {
      score_isolation_forest(models$Isolation_Forest, x)
    }
  )
}

run_aim2_once <- function(x, labels, counts, seed) {
  split <- make_aim2_splits(nrow(x), counts, seed)

  train_labels <- labels[split$detector_train]
  reference_labels <- labels[split$reference]
  aumvc_labels <- labels[split$aumvc]

  x_train <- x[
    split$detector_train[train_labels == 0L],
    ,
    drop = FALSE
  ]

  x_reference <- x[
    split$reference[reference_labels == 0L],
    ,
    drop = FALSE
  ]

  x_aumvc <- x[
    split$aumvc[aumvc_labels == 0L],
    ,
    drop = FALSE
  ]

  x_label <- x[split$label_eval, , drop = FALSE]
  labels_label <- labels[split$label_eval]

  if (
    nrow(x_train) <= LOF_K ||
    nrow(x_reference) < 2L ||
    nrow(x_aumvc) < 2L ||
    length(unique(labels_label)) != 2L
  ) {
    stop("The Aim 2 split is not usable.", call. = FALSE)
  }

  standardizer <- fit_standardizer(x_train)

  x_train <- apply_standardizer(x_train, standardizer)
  x_reference <- apply_standardizer(x_reference, standardizer)
  x_aumvc <- apply_standardizer(x_aumvc, standardizer)
  x_label <- apply_standardizer(x_label, standardizer)

  reference <- make_reference(
    x_reference,
    n_reference = N_REFERENCE,
    n_mc_repetitions = N_MC_REPETITIONS,
    seed = seed + 100L
  )

  models <- fit_aim2_models(x_train, seed)
  score_functions <- make_score_functions(models)

  results <- do.call(
    rbind,
    lapply(
      names(score_functions),
      function(detector) {
        score_fun <- score_functions[[detector]]
        label_scores <- score_fun(x_label)

        mv <- aumvc(
          x_eval = x_aumvc,
          reference = reference,
          score_fun = score_fun,
          score_direction = "anomaly",
          alpha_grid = AUMVC_ALPHA_GRID
        )

        data.frame(
          detector = detector,
          aumvc = mv$aumvc,
          aumvc_normalized = mv$aumvc_normalized,
          roc_auc = roc_auc_score(labels_label, label_scores),
          pr_auc = pr_auc_score(labels_label, label_scores)
        )
      }
    )
  )

  roc_pairs <- pairwise_concordance(
    results$detector,
    results$aumvc,
    results$roc_auc
  )

  pr_pairs <- pairwise_concordance(
    results$detector,
    results$aumvc,
    results$pr_auc
  )

  either_pairs <- either_concordance(roc_pairs, pr_pairs)

  list(
    results = results,
    roc_pairs = roc_pairs,
    pr_pairs = pr_pairs,
    either_pairs = either_pairs
  )
}

summarize_detector_runs <- function(run_results) {
  combined <- do.call(
    rbind,
    lapply(
      seq_along(run_results),
      function(run) {
        data.frame(
          run = run,
          run_results[[run]]$results
        )
      }
    )
  )

  detectors <- unique(combined$detector)

  summary <- do.call(
    rbind,
    lapply(
      detectors,
      function(detector) {
        x <- combined[
          combined$detector == detector,
          ,
          drop = FALSE
        ]

        data.frame(
          detector = detector,
          aumvc_mean = mean(x$aumvc),
          aumvc_sd = sd(x$aumvc),
          aumvc_normalized_mean = mean(x$aumvc_normalized),
          aumvc_normalized_sd = sd(x$aumvc_normalized),
          roc_mean = mean(x$roc_auc),
          roc_sd = sd(x$roc_auc),
          pr_mean = mean(x$pr_auc),
          pr_sd = sd(x$pr_auc)
        )
      }
    )
  )

  list(
    combined = combined,
    summary = summary
  )
}

summarize_concordance_runs <- function(run_results) {
  metrics <- c(
    "AUMVC vs ROC-AUC",
    "AUMVC vs PR-AUC",
    "AUMVC vs either"
  )

  pair_names <- c(
    "roc_pairs",
    "pr_pairs",
    "either_pairs"
  )

  do.call(
    rbind,
    lapply(
      seq_along(metrics),
      function(i) {
        pairs <- do.call(
          rbind,
          lapply(
            run_results,
            function(x) x[[pair_names[i]]]
          )
        )

        concordant <- sum(pairs$agree, na.rm = TRUE)
        compared <- sum(!is.na(pairs$agree))

        data.frame(
          metric = metrics[i],
          concordant = concordant,
          compared = compared,
          percentage = 100 * concordant / compared
        )
      }
    )
  )
}

format_mean_sd <- function(mean, sd, digits = 4L) {
  paste0(
    formatC(mean, digits = digits, format = "f"),
    " +/- ",
    formatC(sd, digits = digits, format = "f")
  )
}

run_aim2_dataset <- function(x, labels, dataset, counts) {
  x <- validate_matrix(x, "x")
  labels <- as.integer(labels)

  if (
    length(labels) != nrow(x) ||
    anyNA(labels) ||
    !all(labels %in% c(0L, 1L))
  ) {
    stop("labels must be binary 0/1.", call. = FALSE)
  }

  run_results <- lapply(
    seq_len(N_RUNS),
    function(run) {
      run_aim2_once(
        x = x,
        labels = labels,
        counts = counts,
        seed = run_seed(run)
      )
    }
  )

  detector_runs <- summarize_detector_runs(run_results)
  concordance <- summarize_concordance_runs(run_results)

  output <- list(
    dataset = dataset,
    n_runs = N_RUNS,
    base_seed = SEED,
    run_results = run_results,
    detector_results = detector_runs$combined,
    detector_summary = detector_runs$summary,
    concordance = concordance
  )

  dir.create(
    "experiments/aim2/results",
    recursive = TRUE,
    showWarnings = FALSE
  )

  saveRDS(
    output,
    file.path(
      "experiments/aim2/results",
      paste0(tolower(dataset), ".rds")
    )
  )

  display <- data.frame(
    detector = output$detector_summary$detector,
    AUMVC = format_mean_sd(
      output$detector_summary$aumvc_mean,
      output$detector_summary$aumvc_sd
    ),
    ROC_AUC = format_mean_sd(
      output$detector_summary$roc_mean,
      output$detector_summary$roc_sd
    ),
    PR_AUC = format_mean_sd(
      output$detector_summary$pr_mean,
      output$detector_summary$pr_sd
    )
  )

  cat(dataset, " - ", N_RUNS, " runs\n\n", sep = "")
  print(display, row.names = FALSE)

  cat("\n")

  display_concordance <- concordance
  names(display_concordance)[2L] <- "matches"

  print(
    display_concordance[
      ,
      c("metric", "matches", "compared", "percentage")
    ],
    row.names = FALSE
  )

  invisible(output)
}
