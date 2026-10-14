source("experiments/aim2/evaluation.R")
source("experiments/aim2/real_data.R")

settings <- AIM2_SETTINGS
display_names <- c(adult = "Adult", http = "HTTP", pima = "Pima",
                   smtp = "SMTP", wilt = "Wilt")
if (length(settings$real_datasets) == 0L) stop("Select at least one real dataset.")
outputs <- lapply(settings$real_datasets, function(dataset) {
  data <- load_real(dataset)
  run_dataset(
    data$x,
    data$labels,
    display_names[[dataset]],
    settings
  )
})

cat("\nTOTAL across ", length(outputs), " datasets and ", settings$n_runs,
    " runs per dataset\n", sep = "")
total <- total_concordance(outputs)
print(total, row.names = FALSE)

rows <- do.call(rbind, lapply(outputs, report_rows))
total_label <- paste0("TOTAL across ", length(outputs), " datasets and ",
                      settings$n_runs, " runs per dataset")
total_rows <- comparison_rows(total_label, total)
directory <- next_experiment("experiments/aim2/result/real")
path <- file.path(directory, "real_results.csv")
save_csv(rbind(rows, total_rows), path)
cat("Saved: ", path, "\n", sep = "")
