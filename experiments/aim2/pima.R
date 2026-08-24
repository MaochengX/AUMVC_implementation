source("experiments/aim2/validation_agianst_labels.R")
settings <- AIM2_SETTINGS

data <- read.csv("dataset/pima/pima.csv")
features <- setdiff(names(data), "label")

run_aim2_dataset(
  as.matrix(data[, features, drop = FALSE]),
  data$label,
  "Pima",
  settings$split_counts$pima,
  settings
)
