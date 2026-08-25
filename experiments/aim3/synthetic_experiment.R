source("experiments/settings.R")
source("experiments/aim3/synthetic.R")
source("experiments/aim3/representation_evaluation.R")

settings <- AIM3_SETTINGS

results_dir <- "experiments/aim3/results"
checkpoint_dir <- file.path(results_dir, "checkpoints")
settings_file <- file.path(checkpoint_dir, "settings.rds")
final_file <- file.path(results_dir, "synthetic.rds")

dir.create(
  checkpoint_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (file.exists(settings_file)) {
  saved_settings <- readRDS(settings_file)

  if (!identical(saved_settings, settings)) {
    stop(
      paste(
        "Aim 3 settings changed since the checkpoints were created.",
        "Delete experiments/aim3/results/checkpoints before rerunning."
      ),
      call. = FALSE
    )
  }
} else {
  saveRDS(settings, settings_file)
}

grid <- expand.grid(
  truth = settings$synthetic$truth_functions,
  intrinsic_dim = settings$synthetic$intrinsic_dims,
  snr = settings$synthetic$snr_levels,
  run = seq_len(settings$n_runs),
  stringsAsFactors = FALSE
)

aim3_case_id <- function(truth, intrinsic_dim, snr, run) {
  sprintf(
    "%s_q%d_snr%s_run%02d",
    truth,
    intrinsic_dim,
    format(snr, scientific = FALSE, trim = TRUE),
    run
  )
}

aim3_format_time <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) {
    return("--:--:--")
  }

  seconds <- as.integer(round(seconds))
  hours <- seconds %/% 3600L
  minutes <- (seconds %% 3600L) %/% 60L
  seconds <- seconds %% 60L

  sprintf("%02d:%02d:%02d", hours, minutes, seconds)
}

aim3_show_progress <- function(done, total, elapsed, mean_case_time) {
  width <- 30L
  fraction <- if (total > 0L) done / total else 1
  filled <- min(width, floor(width * fraction))

  bar <- paste0(
    "[",
    paste(rep("=", filled), collapse = ""),
    paste(rep(" ", width - filled), collapse = ""),
    "]"
  )

  remaining <- total - done
  eta <- if (is.finite(mean_case_time)) {
    remaining * mean_case_time
  } else {
    NA_real_
  }

  cat(
    sprintf(
      "\r%s %d/%d %5.1f%% | elapsed %s | ETA %s",
      bar,
      done,
      total,
      100 * fraction,
      aim3_format_time(elapsed),
      aim3_format_time(eta)
    )
  )
  flush.console()
}

case_ids <- mapply(
  aim3_case_id,
  grid$truth,
  grid$intrinsic_dim,
  grid$snr,
  grid$run,
  USE.NAMES = FALSE
)

checkpoint_files <- file.path(
  checkpoint_dir,
  paste0(case_ids, ".rds")
)

already_done <- file.exists(checkpoint_files)
total_cases <- nrow(grid)
done_cases <- sum(already_done)

cat(
  sprintf(
    "Aim 3 synthetic experiment: %d total cases, %d already complete\n",
    total_cases,
    done_cases
  )
)

start_time <- proc.time()[[3L]]
case_times <- numeric(0)

aim3_show_progress(
  done_cases,
  total_cases,
  elapsed = 0,
  mean_case_time = NA_real_
)

for (i in seq_len(total_cases)) {
  checkpoint_file <- checkpoint_files[i]

  if (file.exists(checkpoint_file)) {
    next
  }

  case_start <- proc.time()[[3L]]

  truth_index <- match(
    grid$truth[i],
    settings$synthetic$truth_functions
  )

  dimension_index <- match(
    grid$intrinsic_dim[i],
    settings$synthetic$intrinsic_dims
  )

  run_seed <- aim_run_seed(settings, grid$run[i])

  case_seed <- run_seed +
    10000L * truth_index +
    1000L * dimension_index

  case <- aim3_generate_synthetic(
    intrinsic_dim = grid$intrinsic_dim[i],
    snr = grid$snr[i],
    truth = grid$truth[i],
    seed = case_seed,
    settings = settings$synthetic
  )

  result <- aim3_evaluate_representations(
    case,
    seed = run_seed + 50000L,
    settings = settings
  )

  result$truth <- grid$truth[i]
  result$intrinsic_dim <- grid$intrinsic_dim[i]
  result$snr <- grid$snr[i]
  result$run <- grid$run[i]
  result$case_seed <- case_seed
  result$analysis_seed <- run_seed + 50000L

  saveRDS(result, checkpoint_file)

  case_times <- c(
    case_times,
    proc.time()[[3L]] - case_start
  )

  done_cases <- done_cases + 1L
  elapsed <- proc.time()[[3L]] - start_time

  aim3_show_progress(
    done_cases,
    total_cases,
    elapsed,
    mean(case_times)
  )
}

cat("\nCombining checkpoints...\n")

results <- lapply(
  checkpoint_files,
  readRDS
)

results <- do.call(rbind, results)

first_columns <- c(
  "truth",
  "intrinsic_dim",
  "snr",
  "run",
  "case_seed",
  "analysis_seed",
  "representation"
)

results <- results[, c(
  first_columns,
  setdiff(names(results), first_columns)
)]

saveRDS(
  list(
    settings = settings,
    results = results
  ),
  final_file
)

cat(
  "Aim 3 synthetic experiment complete\n",
  "Saved final results to: ",
  final_file,
  "\n",
  sep = ""
)
