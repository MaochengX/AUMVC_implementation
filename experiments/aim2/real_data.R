aim2_load_adult <- function() {
  columns <- c(
    "age", "workclass", "fnlwgt", "education", "education_num",
    "marital_status", "occupation", "relationship", "race", "sex",
    "capital_gain", "capital_loss", "hours_per_week", "native_country", "income"
  )
  read_adult <- function(path) {
    data <- read.table(
      path, sep = ",", header = FALSE, col.names = columns,
      strip.white = TRUE, comment.char = "|", quote = "",
      stringsAsFactors = FALSE
    )
    data$income <- sub("\\.$", "", trimws(data$income))
    data
  }
  data <- rbind(
    read_adult("dataset/adult/adult.data"),
    read_adult("dataset/adult/adult.test")
  )
  features <- c("age", "fnlwgt", "education_num", "capital_gain",
                "capital_loss", "hours_per_week")
  list(x = as.matrix(data[, features, drop = FALSE]),
       labels = as.integer(data$income == ">50K"))
}

aim2_load_labeled_csv <- function(name) {
  data <- read.csv(file.path("dataset", name, paste0(name, ".csv")))
  list(x = as.matrix(data[, names(data) != "label", drop = FALSE]),
       labels = data$label)
}

aim2_load_wilt <- function() {
  data <- rbind(
    read.csv("dataset/wilt/training.csv", stringsAsFactors = FALSE),
    read.csv("dataset/wilt/testing.csv", stringsAsFactors = FALSE)
  )
  features <- c("GLCM_pan", "Mean_Green", "Mean_Red", "Mean_NIR", "SD_pan")
  list(x = as.matrix(data[, features, drop = FALSE]),
       labels = as.integer(tolower(trimws(data$class)) == "w"))
}

aim2_load_real <- function(name) {
  switch(name,
    adult = aim2_load_adult(),
    http = aim2_load_labeled_csv("http"),
    pima = aim2_load_labeled_csv("pima"),
    smtp = aim2_load_labeled_csv("smtp"),
    wilt = aim2_load_wilt(),
    stop("Unknown Aim 2 dataset: ", name, call. = FALSE)
  )
}
