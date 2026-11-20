AIM4_REPRODUCTION_ROOT <- "experiments/aim4/reproduction"

analyze_reproduction <- function(root = AIM4_REPRODUCTION_ROOT) {
  studies <- c("boston", "gaussian_mixture")
  results <- lapply(studies, function(study) {
    folder <- file.path(root, study)
    values <- utils::read.csv(file.path(folder, "bandwidth.csv"))
    required <- c("sigma", "gamma", "mv_area", "volume_95", "selected",
                  "valid_calibration", "calibration_max_error")
    if (!all(required %in% names(values)) || any(!is.finite(values$mv_area)) ||
        anyDuplicated(values$sigma) || sum(values$selected) != 1L) {
      stop("Invalid reproduction results for ", study)
    }
    ranked <- values[order(values$mv_area, values$sigma), , drop = FALSE]
    if (!isTRUE(ranked$selected[1L])) stop("Selection does not minimize MV area: ", study)
    top <- head(ranked, 3L)

    if (study == "boston") {
      ordered <- values[order(values$sigma), , drop = FALSE]
      grDevices::pdf(file.path(folder, "selection_curve.pdf"), width = 7, height = 4.8)
      tryCatch({
        plot(ordered$sigma, ordered$mv_area, type = "b", pch = 16, cex = 0.65,
             xlab = expression(sigma), ylab = "Area under MV curve",
             main = "Boston Housing: bandwidth selection")
        abline(v = top$sigma[1L], lty = 2, col = "#C34A31")
        abline(v = 0.42, lty = 3, col = "gray45")
        legend("topright", c(sprintf("selected: %.3f", top$sigma[1L]), "paper: 0.42"),
               lty = c(2, 3), col = c("#C34A31", "gray45"), bty = "n")
      }, finally = grDevices::dev.off())
    }

    result <- data.frame(
      study = study, rank = seq_len(nrow(top)), n_sigma = nrow(values),
      sigma = top$sigma, gamma = top$gamma,
      mv_area = top$mv_area, volume_95 = top$volume_95,
      valid_calibration = top$valid_calibration,
      calibration_max_error = top$calibration_max_error,
      true_mass = rep(NA_real_, nrow(top)),
      symmetric_difference_volume_in_box = rep(NA_real_, nrow(top))
    )
    path <- file.path(folder, "set_comparison.csv")
    if (study == "gaussian_mixture" && file.exists(path)) {
      comparison <- utils::read.csv(path)
      if (nrow(comparison) != 1L || abs(comparison$sigma - top$sigma[1L]) > 1e-9) {
        stop("Gaussian set comparison does not match selected sigma.")
      }
      result$true_mass[1L] <- comparison$true_mass
      result$symmetric_difference_volume_in_box[1L] <-
        comparison$symmetric_difference_volume_in_box
    }
    result
  })
  report <- do.call(rbind, results)
  utils::write.csv(report, file.path(root, "summary.csv"), row.names = FALSE)
  print(report, row.names = FALSE)
  invisible(report)
}

if (sys.nframe() == 0L) analyze_reproduction()
