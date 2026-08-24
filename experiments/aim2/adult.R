source("experiments/aim2/validation_agianst_labels.R")
settings <- AIM2_SETTINGS

columns <- c(
  "age", "workclass", "fnlwgt", "education", "education_num",
  "marital_status", "occupation", "relationship", "race", "sex",
  "capital_gain", "capital_loss", "hours_per_week", "native_country", "income"
)

read_adult <- function(path) {
  data <- read.table(
    path,
    sep = ",",
    header = FALSE,
    col.names = columns,
    strip.white = TRUE,
    comment.char = "|",
    quote = "",
    stringsAsFactors = FALSE
  )
  data$income <- sub("\\.$", "", trimws(data$income))
  data
}

data <- rbind(
  read_adult("dataset/adult/adult.data"),
  read_adult("dataset/adult/adult.test")
)
features <- c("age", "fnlwgt", "education_num", "capital_gain", "capital_loss", "hours_per_week")

run_aim2_dataset(
  as.matrix(data[, features, drop = FALSE]),
  as.integer(data$income == ">50K"),
  "Adult",
  settings$split_counts$adult,
  settings
)
