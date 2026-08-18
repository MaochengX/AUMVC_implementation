source("experiments/aim2/validation_agianst_labels.R")

data <- rbind(
  read.csv(
    "dataset/wilt/training.csv",
    stringsAsFactors = FALSE
  ),
  read.csv(
    "dataset/wilt/testing.csv",
    stringsAsFactors = FALSE
  )
)

features <- c(
  "GLCM_pan",
  "Mean_Green",
  "Mean_Red",
  "Mean_NIR",
  "SD_pan"
)

run_aim2_dataset(
  x = as.matrix(data[, features]),
  labels = as.integer(
    tolower(trimws(data$class)) == "w"
  ),
  dataset = "Wilt",
  counts = AIM2_SPLIT_COUNTS$wilt
)
