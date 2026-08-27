path <- "experiments/aim3/results/synthetic.rds"

if (!file.exists(path)) {
  stop(
    "Run experiments/aim3/synthetic_experiment.R first.",
    call. = FALSE
  )
}

saved <- readRDS(path)
results <- saved$results

aim3_mean_sd <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0L) return("NA")
  if (length(x) == 1L) return(sprintf("%.4g", x))
  if (all(is.finite(x))) {
    return(sprintf("%.4g +/- %.4g", mean(x), sd(x)))
  }
  if (all(is.infinite(x) & x > 0)) return("Inf")

  "non-finite"
}

keys <- unique(
  results[, c(
    "truth",
    "intrinsic_dim",
    "snr",
    "representation"
  )]
)

summary_rows <- lapply(
  seq_len(nrow(keys)),
  function(i) {
    keep <-
      results$truth == keys$truth[i] &
      results$intrinsic_dim == keys$intrinsic_dim[i] &
      results$snr == keys$snr[i] &
      results$representation == keys$representation[i]

    x <- results[keep, , drop = FALSE]

    data.frame(
      truth = keys$truth[i],
      intrinsic_dim = keys$intrinsic_dim[i],
      snr = keys$snr[i],
      representation = keys$representation[i],
      AUMVC = aim3_mean_sd(x$aumvc),
      AUMVC_normalized = aim3_mean_sd(x$aumvc_normalized),
      MC_SE = aim3_mean_sd(x$aumvc_mc_se),
      MC_SE_normalized = aim3_mean_sd(x$aumvc_normalized_mc_se),
      mean_occupancy = aim3_mean_sd(x$mean_occupancy),
      zero_occupancy = aim3_mean_sd(x$zero_occupancy_fraction),
      log_box_volume = aim3_mean_sd(x$log_box_volume),
      N_reference = aim3_mean_sd(x$n_reference),
      capped = sum(x$reference_capped),
      log_reference_density = aim3_mean_sd(x$log_reference_density),
      runtime_seconds = aim3_mean_sd(x$runtime_seconds)
    )
  }
)

summary_table <- do.call(rbind, summary_rows)

write.csv(
  summary_table,
  "experiments/aim3/results/summary.csv",
  row.names = FALSE
)

print(
  summary_table,
  row.names = FALSE
)
