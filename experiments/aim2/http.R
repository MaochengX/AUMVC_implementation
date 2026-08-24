source("experiments/aim2/validation_agianst_labels.R")

data <- read.csv("dataset/http/http.csv")
features <- setdiff(names(data), "label")

run_aim2_dataset(
  x = as.matrix(data[, features, drop = FALSE]),
  labels = data$label,
  dataset = "HTTP",
  counts = AIM2_SETTINGS$split_counts$http
)
