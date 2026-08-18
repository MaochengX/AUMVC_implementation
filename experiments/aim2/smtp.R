source("experiments/aim2/validation_agianst_labels.R")

data <- read.csv("dataset/smtp/smtp.csv")
features <- setdiff(names(data), "label")

run_aim2_dataset(
  x = as.matrix(data[, features, drop = FALSE]),
  labels = data$label,
  dataset = "SMTP",
  counts = AIM2_SPLIT_COUNTS$smtp
)
