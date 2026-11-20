AIM4_TUNING_ROOT <- "experiments/aim4"

mean_se <- function(x) {
  if (!length(x)) return(c(mean = NA_real_, se = NA_real_))
  c(mean = mean(x), se = if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_)
}

check_tuning_results <- function(selections, candidates) {
  needed <- c("run", "selection_status", "n_reliable", "aumvc_choice",
              "roc_label_choice", "pr_label_choice", "default", "test_roc_aumvc",
              "test_pr_aumvc", "roc_loss_vs_labels", "pr_loss_vs_labels",
              "roc_gain_vs_default", "pr_gain_vs_default")
  if (!all(needed %in% names(selections))) stop("Missing selection columns.")
  if (!all(c("run", "parameter", "reliable", "test_roc", "test_pr") %in% names(candidates))) {
    stop("Missing candidate columns.")
  }
  if (anyDuplicated(selections$run) || anyNA(selections[, needed]) ||
      any(!selections$selection_status %in% c("selected", "tie", "all_unreliable"))) {
    stop("Invalid selection rows.")
  }
  by_run <- split(candidates, candidates$run)
  for (i in seq_len(nrow(selections))) {
    s <- selections[i, ]
    group <- by_run[[as.character(s$run)]]
    if (is.null(group) || anyDuplicated(group$parameter) ||
        sum(group$reliable) != s$n_reliable) stop("Candidate mismatch at run ", s$run)
    index <- match(c(s$aumvc_choice, s$roc_label_choice, s$pr_label_choice,
                     s$default), group$parameter)
    if (anyNA(index)) stop("Missing selected parameter at run ", s$run)
    if (s$selection_status != "selected" && s$aumvc_choice != s$default) {
      stop("Fallback does not use default at run ", s$run)
    }
    roc <- group$test_roc[index]
    pr <- group$test_pr[index]
    recorded <- c(s$test_roc_aumvc, s$test_pr_aumvc,
                  s$roc_loss_vs_labels, s$pr_loss_vs_labels,
                  s$roc_gain_vs_default, s$pr_gain_vs_default)
    expected <- c(roc[1L], pr[1L], roc[2L] - roc[1L], pr[3L] - pr[1L],
                  roc[1L] - roc[4L], pr[1L] - pr[4L])
    if (any(!is.finite(expected)) || any(abs(recorded - expected) > 1e-7)) {
      stop("Paired test difference mismatch at run ", s$run)
    }
  }
  invisible(TRUE)
}

summarize_tuning <- function(selections, detector, dataset, scope) {
  selected <- selections$selection_status == "selected"
  rows <- if (scope == "selected_only") selections[selected, , drop = FALSE] else selections
  roc_loss <- mean_se(rows$roc_loss_vs_labels)
  pr_loss <- mean_se(rows$pr_loss_vs_labels)
  roc_gain <- mean_se(rows$roc_gain_vs_default)
  pr_gain <- mean_se(rows$pr_gain_vs_default)
  data.frame(
    detector = detector, dataset = dataset, scope = scope,
    n_runs = nrow(selections), n_used = nrow(rows),
    n_no_choice = sum(!selected),
    no_choice_rate = mean(!selected),
    n_all_unreliable = sum(selections$selection_status == "all_unreliable"),
    n_ties = sum(selections$selection_status == "tie"),
    roc_loss_mean = roc_loss["mean"], roc_loss_se = roc_loss["se"],
    pr_loss_mean = pr_loss["mean"], pr_loss_se = pr_loss["se"],
    roc_gain_default_mean = roc_gain["mean"], roc_gain_default_se = roc_gain["se"],
    pr_gain_default_mean = pr_gain["mean"], pr_gain_default_se = pr_gain["se"],
    row.names = NULL
  )
}

run_tuning_analysis <- function(root = AIM4_TUNING_ROOT) {
  for (detector in c("ocsvm", "lof")) {
    folder <- file.path(root, paste0(detector, "_tune"))
    datasets <- c("ecg5000", "fashion_mnist")
    results <- lapply(datasets, function(dataset) {
      path <- file.path(folder, dataset)
      choices <- utils::read.csv(file.path(path, "selections.csv"))
      candidates <- utils::read.csv(file.path(path, "candidates.csv"))
      check_tuning_results(choices, candidates)
      if (nrow(choices) < 20L) warning("Fewer than 20 runs for ", detector, "/", dataset)
      rbind(summarize_tuning(choices, detector, dataset, "full_rule"),
            summarize_tuning(choices, detector, dataset, "selected_only"))
    })
    report <- do.call(rbind, results)
    utils::write.csv(report, file.path(folder, "evaluation.csv"), row.names = FALSE)
    print(report[, c("dataset", "scope", "n_used", "n_no_choice",
                     "roc_loss_mean", "roc_loss_se", "pr_loss_mean", "pr_loss_se")],
          row.names = FALSE)
  }
}

if (sys.nframe() == 0L) run_tuning_analysis()
