source("experiments/settings.R")

datasets <- c(
  "adult",
  "http",
  "pima",
  "smtp",
  "wilt"
)

paths <- file.path(
  "experiments/aim2/results",
  paste0(datasets, ".rds")
)

if (any(!file.exists(paths))) {
  stop(
    "Run adult.R, http.R, pima.R, smtp.R, and wilt.R first.",
    call. = FALSE
  )
}

runs <- lapply(paths, readRDS)

if (any(vapply(runs, function(x) x$n_runs != AIM2_SETTINGS$n_runs, logical(1)))) {
  stop(
    "Saved results do not use the current Aim 2 run count.",
    call. = FALSE
  )
}

if (any(vapply(runs, function(x) x$base_seed != AIM2_SETTINGS$seed, logical(1)))) {
  stop(
    "Saved results do not use the current Aim 2 seed.",
    call. = FALSE
  )
}

dataset_summary <- do.call(
  rbind,
  lapply(
    runs,
    function(x) {
      data.frame(
        dataset = x$dataset,
        roc = x$concordance$percentage[1L],
        pr = x$concordance$percentage[2L],
        either = x$concordance$percentage[3L]
      )
    }
  )
)

overall <- do.call(
  rbind,
  lapply(
    seq_len(3L),
    function(i) {
      concordant <- sum(
        vapply(
          runs,
          function(x) x$concordance$concordant[i],
          numeric(1)
        )
      )

      compared <- sum(
        vapply(
          runs,
          function(x) x$concordance$compared[i],
          numeric(1)
        )
      )

      data.frame(
        metric = runs[[1L]]$concordance$metric[i],
        concordant = concordant,
        compared = compared,
        percentage = 100 * concordant / compared
      )
    }
  )
)

cat(
  "Aim 2 summary - ",
  AIM2_SETTINGS$n_runs,
  " runs per dataset\n\n",
  sep = ""
)

print(dataset_summary, row.names = FALSE)

cat("\nOverall\n")

display_overall <- overall
names(display_overall)[2L] <- "matches"

print(
  display_overall[
    ,
    c("metric", "matches", "compared", "percentage")
  ],
  row.names = FALSE
)

cat("\nGoix reference for ROC/PR: 73.3%\n")
