source("experiments/aim2/evaluation.R")
source("experiments/aim2/real_data.R")

settings <- AIM2_SETTINGS
display_names <- c(adult = "Adult", http = "HTTP", pima = "Pima",
                   smtp = "SMTP", wilt = "Wilt")
if (length(settings$real_datasets) == 0L) stop("Select at least one real dataset.")
outputs <- lapply(settings$real_datasets, function(dataset) {
  data <- aim2_load_real(dataset)
  aim2_run_dataset(
    data$x, data$labels, display_names[[dataset]],
    settings$split_counts[[dataset]], settings
  )
})

cat("\nTOTAL across ", length(outputs), " datasets and ", settings$n_runs,
    " runs per dataset\n", sep = "")
print(aim2_total_concordance(outputs), row.names = FALSE)
