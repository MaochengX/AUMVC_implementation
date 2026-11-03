source("experiments/aim3/real_experiment.R")

MV_EXP_RUNS <- 1:3
MV_EXP_DATASETS <- c("ecg5000", "shuttle")
MV_EXP_REPRESENTATIONS <- c("ambient", "goix_subsampling", "mds", "isomap")
MV_EXP_GAMMA <- c(0.25, 1, 2, 4, 8, 16, 32, 64)
MV_EXP_COMPARE <- c(0.25, 4)
MV_EXP_MASSES <- 0.95
MV_EXP_OUTPUT <- "experiments/aim3/mv_curve_exp"

plot_mv_comparison <- function(dataset, run, representation, curves, normal_curves,
                               anomaly_fraction) {
  old <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(old))
  for (normal_only in c(FALSE, TRUE)) {
    selected <- if (normal_only) normal_curves else curves
    vals <- unlist(lapply(selected, function(x) x$volume_normalized))
    plot(NA, xlim = range(selected[[1L]]$alpha),
         ylim = range(c(0, vals), finite = TRUE), xlab = "Target mass",
         ylab = "MC occupancy (volume / box volume)",
         main = if (normal_only) "Normal-only evaluation (diagnostic)" else "Full evaluation")
    if (!normal_only) abline(v = 1 - anomaly_fraction, lty = 3, col = "gray50")
    for (i in seq_along(selected)) {
      lines(selected[[i]]$alpha, selected[[i]]$volume_normalized,
            col = c("#2878B5", "#D55E00")[i], lwd = 2)
    }
    legend("topleft",
           legend = c("gamma_scale 0.25", "gamma_scale 4", if (!normal_only) "normal share"),
           col = c("#2878B5", "#D55E00", if (!normal_only) "gray50"),
           lty = c(1, 1, if (!normal_only) 3), bty = "n", cex = 0.75)
  }
  mtext(sprintf("%s, %s", dataset, representation),
        outer = TRUE, line = 0.4, cex = 0.95)
}

plot_level_slices <- function(case) {
  old <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(old))
  bounds <- rbind(case$reference$box$lower[1:2], case$reference$box$upper[1:2])
  xs <- seq(bounds[1L, 1L], bounds[2L, 1L], length.out = 120L)
  ys <- seq(bounds[1L, 2L], bounds[2L, 2L], length.out = 120L)
  slice <- expand.grid(xs, ys)
  template <- matrix(rep(apply(case$train, 2L, median), each = nrow(slice)),
                     nrow = nrow(slice))
  template[, 1:2] <- as.matrix(slice)
  axis_name <- if (case$representation == "ambient") "Standardized feature" else "Embedding coordinate"
  for (i in seq_along(MV_EXP_COMPARE)) {
    gamma <- as.character(MV_EXP_COMPARE[i])
    scores <- score_ocsvm(case$models[[gamma]], template)
    mv <- case$curves[[gamma]]
    threshold <- mv$threshold[which.min(abs(mv$alpha - 0.95))]
    anomaly_cutoff <- -threshold
    z <- matrix(scores, nrow = length(xs), ncol = length(ys))
    image(xs, ys, (z <= anomaly_cutoff) * 1, breaks = c(-0.5, 0.5, 1.5),
          col = c("white", "#D8EAF3"),
          xlab = paste(axis_name, 1L), ylab = paste(axis_name, 2L),
          main = sprintf("gamma_scale %g, alpha 0.95", MV_EXP_COMPARE[i]))
    if (min(z) < anomaly_cutoff && max(z) > anomaly_cutoff) {
      contour(xs, ys, z, levels = anomaly_cutoff, add = TRUE,
              drawlabels = FALSE, col = "#2878B5", lwd = 2)
    }
    points(case$evaluation[, 1L], case$evaluation[, 2L],
           pch = ifelse(case$eval_labels == 1L, 4L, 1L),
           col = ifelse(case$eval_labels == 1L, "#D55E00", "gray35"), cex = 0.48)
  }
  mtext(sprintf("%s, %s", case$dataset, case$representation),
        outer = TRUE, line = 0.4, cex = 0.95)
}

mv_exp_flags <- function(mv, reference_scores) {
  near_constant <- function(x) diff(range(x)) <= 1e-10 * max(1, max(abs(x)))
  c(
    threshold_collapsed = near_constant(mv$mv_curve$threshold),
    reference_collapsed = near_constant(reference_scores),
    full_box = any(mv$mv_curve$volume_normalized >= 1 - 1e-12)
  )
}

mv_exp_diagnostic <- function(dataset, run, representation, gamma, aumvc,
                                    roc, pr, min_hits, relative_se, flags,
                                    saturation_rate, settings, score_basis,
                                    normal_only = NA_real_, subsets_reliable = TRUE) {
  data.frame(
    dataset = dataset, run = run, representation = representation,
    gamma_scale = gamma, aumvc_normalized = aumvc,
    normal_only_aumvc = normal_only, roc_auc = roc, pr_auc = pr,
    min_hit_count = min_hits, relative_mc_se = relative_se,
    threshold_collapsed = flags[["threshold_collapsed"]],
    reference_collapsed = flags[["reference_collapsed"]],
    full_box = flags[["full_box"]], full_box_mass_rate = saturation_rate,
    reliable = is.finite(aumvc) && is.finite(relative_se) &&
      min_hits >= settings$minimum_hit_count &&
      relative_se <= settings$maximum_relative_mc_se &&
      !any(flags) && subsets_reliable,
    score_basis = score_basis
  )
}

mv_exp_inclusion <- function(dataset, run, representation, gamma, alpha,
                                   eval_labels, n_inside, n_anomalies_inside,
                                   occupancy, normal_occupancy, basis) {
  n_eval <- length(eval_labels)
  n_anomalies <- sum(eval_labels == 1L)
  data.frame(
    dataset = dataset, run = run, representation = representation,
    gamma_scale = gamma, alpha = alpha, n_eval = n_eval,
    n_anomalies_eval = n_anomalies,
    minimum_unavoidable_anomalies = max(0, ceiling(alpha * n_eval - 1e-9) - (n_eval - n_anomalies)),
    n_inside = n_inside, n_anomalies_inside = n_anomalies_inside,
    occupancy_full = occupancy, occupancy_normal_only = normal_occupancy,
    inclusion_basis = basis
  )
}

mv_exp_single <- function(dataset, run, representation, prepared, eval_labels,
                                split_seed, settings) {
  if (representation == "ambient") {
    item <- prepared
    reference_seed <- split_seed + 1000L
  } else {
    item <- prepare_embedding_data(prepared, representation, settings)
    reference_seed <- split_seed + if (representation == "mds") 4000L else 5000L
  }
  reference <- make_experiment_reference(item$reference, reference_seed, settings)
  scored <- prepare_candidate_scores(
    item$train, item$eval, item$label, reference, MV_EXP_GAMMA, settings
  )
  diagnostic <- vector("list", length(MV_EXP_GAMMA))
  inclusion <- list()
  curves <- normal_curves <- models <- list()
  for (i in seq_along(MV_EXP_GAMMA)) {
    gamma <- MV_EXP_GAMMA[i]
    candidate <- scored$candidates[[i]]
    mv <- aumvc_scored(candidate$evaluation_scores, candidate$reference_scores,
                       reference, "anomaly", settings$aumvc_alpha_grid)
    flags <- mv_exp_flags(mv, candidate$reference_scores)
    normal_only <- NA_real_
    if (gamma %in% MV_EXP_COMPARE) {
      normal_mv <- aumvc_scored(candidate$evaluation_scores[eval_labels == 0L],
                                candidate$reference_scores, reference, "anomaly",
                                settings$aumvc_alpha_grid)
      key <- as.character(gamma)
      curves[[key]] <- mv$mv_curve
      normal_curves[[key]] <- normal_mv$mv_curve
      models[[key]] <- candidate$model
      normal_only <- normal_mv$aumvc_normalized
      for (alpha in MV_EXP_MASSES) {
        j <- which.min(abs(mv$mv_curve$alpha - alpha))
        selected <- candidate$evaluation_scores <= -mv$mv_curve$threshold[j]
        inclusion[[length(inclusion) + 1L]] <- mv_exp_inclusion(
          dataset, run, representation, gamma, mv$mv_curve$alpha[j],
          eval_labels, sum(selected), sum(selected & eval_labels == 1L),
          mv$mv_curve$volume_normalized[j], normal_mv$mv_curve$volume_normalized[j],
          "single_model"
        )
      }
    }
    diagnostic[[i]] <- mv_exp_diagnostic(
      dataset, run, representation, gamma, mv$aumvc_normalized,
      roc_auc_score(prepared$labels_label, candidate$label_scores),
      pr_auc_score(prepared$labels_label, candidate$label_scores),
      mv$minimum_hit_count, mv$aumvc_relative_mc_se, flags,
      mean(mv$mv_curve$volume_normalized >= 1 - 1e-12), settings,
      "single_model", normal_only
    )
  }
  list(
    diagnostic = do.call(rbind, diagnostic), inclusion = do.call(rbind, inclusion),
    curves = curves, normal_curves = normal_curves,
    slice = list(dataset = dataset, run = run, representation = representation,
                 reference = reference, train = item$train, evaluation = item$eval,
                 eval_labels = eval_labels, models = models, curves = curves)
  )
}

mv_exp_goix <- function(dataset, run, prepared, eval_labels, split_seed, settings) {
  subsets <- goix_subsets(ncol(prepared$train), split_seed + 2000L,
                          settings$goix_subsampling)
  n_subsets <- length(subsets)
  n_reference <- settings$n_reference * settings$n_mc_repetitions
  gammas <- MV_EXP_COMPARE
  ns <- length(gammas)
  contributions <- matrix(0, n_reference, ns)
  ensemble_scores <- matrix(0, nrow(prepared$label), ns)
  aumvc_values <- min_hits <- saturation_rates <- matrix(NA_real_, n_subsets, ns)
  subset_reliable <- matrix(FALSE, n_subsets, ns)
  subset_flags <- array(FALSE, c(n_subsets, 3L, ns),
                        dimnames = list(NULL, c("threshold_collapsed", "reference_collapsed", "full_box"), NULL))
  averaged_curves <- averaged_normal_curves <- vector("list", ns)
  normal_aumvc <- matrix(NA_real_, n_subsets, ns)
  totals <- matrix(0, length(MV_EXP_MASSES), ns)
  anomalies <- matrix(0, length(MV_EXP_MASSES), ns)
  for (j in seq_along(subsets)) {
    columns <- subsets[[j]]
    reference <- make_experiment_reference(prepared$reference[, columns, drop = FALSE],
                                           split_seed + 3000L + j, settings)
    scored <- prepare_candidate_scores(
      prepared$train[, columns, drop = FALSE],
      prepared$eval[, columns, drop = FALSE],
      prepared$label[, columns, drop = FALSE], reference, gammas, settings
    )
    for (i in seq_along(gammas)) {
      candidate <- scored$candidates[[i]]
      mv <- aumvc_scored(candidate$evaluation_scores, candidate$reference_scores,
                         reference, "anomaly", settings$aumvc_alpha_grid)
      normal_mv <- aumvc_scored(candidate$evaluation_scores[eval_labels == 0L],
                                candidate$reference_scores, reference, "anomaly",
                                settings$aumvc_alpha_grid)
      subset_flags[j, , i] <- mv_exp_flags(mv, candidate$reference_scores)
      subset_reliable[j, i] <- is.finite(mv$aumvc_normalized) &&
        is.finite(mv$aumvc_relative_mc_se) &&
        mv$minimum_hit_count >= settings$minimum_hit_count &&
        mv$aumvc_relative_mc_se <= settings$maximum_relative_mc_se &&
        !any(subset_flags[j, , i])
      aumvc_values[j, i] <- mv$aumvc_normalized
      normal_aumvc[j, i] <- normal_mv$aumvc_normalized
      min_hits[j, i] <- mv$minimum_hit_count
      saturation_rates[j, i] <- mean(mv$mv_curve$volume_normalized >= 1 - 1e-12)
      contributions[, i] <- contributions[, i] + mv$mc_contribution / n_subsets
      ensemble_scores[, i] <- ensemble_scores[, i] +
        percentile_scores(candidate$evaluation_scores, candidate$label_scores) / n_subsets
      if (is.null(averaged_curves[[i]])) {
        averaged_curves[[i]] <- mv$mv_curve[, c("alpha", "volume_normalized")]
        averaged_normal_curves[[i]] <- normal_mv$mv_curve[, c("alpha", "volume_normalized")]
        averaged_curves[[i]]$volume_normalized <- 0
        averaged_normal_curves[[i]]$volume_normalized <- 0
      }
      averaged_curves[[i]]$volume_normalized <- averaged_curves[[i]]$volume_normalized +
        mv$mv_curve$volume_normalized / n_subsets
      averaged_normal_curves[[i]]$volume_normalized <-
        averaged_normal_curves[[i]]$volume_normalized +
        normal_mv$mv_curve$volume_normalized / n_subsets
      for (k in seq_along(MV_EXP_MASSES)) {
        index <- which.min(abs(mv$mv_curve$alpha - MV_EXP_MASSES[k]))
        selected <- candidate$evaluation_scores <= -mv$mv_curve$threshold[index]
        totals[k, i] <- totals[k, i] + sum(selected) / n_subsets
        anomalies[k, i] <- anomalies[k, i] + sum(selected & eval_labels == 1L) / n_subsets
      }
    }
  }
  diagnostic <- inclusion <- list()
  curves <- normal_curves <- list()
  for (i in seq_along(gammas)) {
    gamma <- gammas[i]
    key <- as.character(gamma)
    curves[[key]] <- averaged_curves[[i]]
    normal_curves[[key]] <- averaged_normal_curves[[i]]
    aumvc_value <- mean(aumvc_values[, i])
    relative_se <- if (is.finite(aumvc_value) && aumvc_value > 0) {
      stats::sd(contributions[, i]) / sqrt(n_reference) / aumvc_value
    } else NA_real_
    flags <- setNames(apply(subset_flags[, , i, drop = FALSE], 2L, any),
                      dimnames(subset_flags)[[2L]])
    diagnostic[[i]] <- mv_exp_diagnostic(
      dataset, run, "goix_subsampling", gamma, aumvc_value,
      roc_auc_score(prepared$labels_label, ensemble_scores[, i]),
      pr_auc_score(prepared$labels_label, ensemble_scores[, i]),
      min(min_hits[, i]), relative_se, flags, mean(saturation_rates[, i]), settings,
      "mean_subset_MV__ensemble_ROC", mean(normal_aumvc[, i]),
      all(subset_reliable[, i])
    )
    for (k in seq_along(MV_EXP_MASSES)) {
      alpha <- curves[[key]]$alpha[which.min(abs(curves[[key]]$alpha - MV_EXP_MASSES[k]))]
      inclusion[[length(inclusion) + 1L]] <- mv_exp_inclusion(
        dataset, run, "goix_subsampling", gamma, alpha, eval_labels,
        totals[k, i], anomalies[k, i],
        curves[[key]]$volume_normalized[which.min(abs(curves[[key]]$alpha - alpha))],
        normal_curves[[key]]$volume_normalized[which.min(abs(curves[[key]]$alpha - alpha))],
        paste0("mean_of_", n_subsets, "_subsets")
      )
    }
  }
  list(diagnostic = do.call(rbind, diagnostic), inclusion = do.call(rbind, inclusion),
       curves = curves, normal_curves = normal_curves, slice = NULL)
}

mv_exp_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

mv_exp_summary <- function(diagnostics, inclusions) {
  groups <- split(diagnostics, interaction(diagnostics$dataset,
                                            diagnostics$representation,
                                            diagnostics$gamma_scale, drop = TRUE))
  rows <- lapply(groups, function(group) {
    included <- inclusions[
      inclusions$dataset == group$dataset[1L] &
      inclusions$representation == group$representation[1L] &
      inclusions$gamma_scale == group$gamma_scale[1L], , drop = FALSE
    ]
    data.frame(
      dataset = group$dataset[1L],
      representation = group$representation[1L],
      gamma_scale = group$gamma_scale[1L],
      n_runs = length(unique(group$run)),
      averaged_aumvc_normalized = mv_exp_mean(group$aumvc_normalized),
      averaged_normal_only_aumvc = mv_exp_mean(group$normal_only_aumvc),
      averaged_roc_auc = mv_exp_mean(group$roc_auc),
      averaged_pr_auc = mv_exp_mean(group$pr_auc),
      reliable_rate = mv_exp_mean(group$reliable),
      full_box_rate = mv_exp_mean(group$full_box),
      threshold_collapsed_rate = mv_exp_mean(group$threshold_collapsed),
      reference_collapsed_rate = mv_exp_mean(group$reference_collapsed),
      averaged_anomalies_inside_95 = mv_exp_mean(included$n_anomalies_inside),
      averaged_minimum_unavoidable_95 = mv_exp_mean(included$minimum_unavoidable_anomalies),
      averaged_occupancy_95 = mv_exp_mean(included$occupancy_full),
      score_basis = group$score_basis[1L]
    )
  })
  output <- do.call(rbind, rows)
  output <- output[order(output$dataset, output$representation, output$gamma_scale), ]
  rownames(output) <- NULL
  output
}

run_mv_experiment <- function(runs = MV_EXP_RUNS,
                                        datasets = MV_EXP_DATASETS,
                                        representations = MV_EXP_REPRESENTATIONS,
                                        settings = AIM3_SETTINGS,
                                        output = MV_EXP_OUTPUT) {
  stopifnot(all(representations %in% MV_EXP_REPRESENTATIONS))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  registry <- real_registry()
  stopifnot(all(datasets %in% names(registry)),
            all(datasets %in% settings$real$datasets))
  loaded <- setNames(lapply(datasets, function(name) {
    registry[[name]]$load(settings$real[[name]]$data_dir)
  }), datasets)
  diagnostics <- inclusions <- slices <- list()
  grDevices::pdf(file.path(output, "mv_curves.pdf"), width = 11, height = 5.5)
  tryCatch({
    for (dataset in datasets) for (run in runs) {
      index <- match(dataset, settings$real$datasets)
      seed <- settings$seed + run * 100000L + index * 10000L
      set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
              sample.kind = "Rejection")
      sampled <- registry[[dataset]]$sample(loaded[[dataset]], seed,
                                           settings$real[[dataset]])
      case_settings <- settings
      case_settings$split_counts <- real_split_counts(nrow(sampled$x))
      split_seed <- seed + 50000L
      prepared <- prepare_case_data(
        list(x = sampled$x, labels = sampled$labels,
             intrinsic_dim = as.integer(settings$real$embedding_dim)),
        split_seed, case_settings
      )
      split <- make_splits(sampled$labels, case_settings$split_counts, split_seed)
      eval_labels <- sampled$labels[split$evaluation]
      for (representation in representations) {
        cat(sprintf("MV: %s, run %d, %s\n", dataset, run, representation))
        case <- if (representation == "goix_subsampling") {
          mv_exp_goix(dataset, run, prepared, eval_labels, split_seed, case_settings)
        } else {
          mv_exp_single(dataset, run, representation, prepared,
                              eval_labels, split_seed, case_settings)
        }
        key <- paste(dataset, run, representation, sep = "_")
        diagnostics[[key]] <- case$diagnostic
        inclusions[[key]] <- case$inclusion
        if (!is.null(case$slice)) slices[[key]] <- case$slice
        plot_mv_comparison(dataset, run, representation, case$curves,
                           case$normal_curves, mean(eval_labels == 1L))
      }
    }
  }, finally = grDevices::dev.off())
  if (length(slices)) {
    grDevices::pdf(file.path(output, "level_set_slices.pdf"), width = 11, height = 5.5)
    tryCatch(for (case in slices) plot_level_slices(case),
             finally = grDevices::dev.off())
  }
  summary <- mv_exp_summary(do.call(rbind, diagnostics), do.call(rbind, inclusions))
  utils::write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
  invisible(summary)
}

if (sys.nframe() == 0L) run_mv_experiment()
