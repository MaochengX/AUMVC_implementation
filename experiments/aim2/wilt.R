source("experiments/aim2/validation_agianst_labels.R")
settings <- AIM2_SETTINGS

data <- rbind(
  read.csv("dataset/wilt/training.csv", stringsAsFactors = FALSE),
  read.csv("dataset/wilt/testing.csv", stringsAsFactors = FALSE)
)
features <- c("GLCM_pan", "Mean_Green", "Mean_Red", "Mean_NIR", "SD_pan")

aim2_run_dataset(
  as.matrix(data[, features, drop = FALSE]),
  as.integer(tolower(trimws(data$class)) == "w"),
  "Wilt",
  settings$split_counts$wilt,
  settings
)
