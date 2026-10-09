script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_directory <- if (length(script_argument)) {
  dirname(normalizePath(sub("^--file=", "", script_argument[1L])))
} else {
  normalizePath("experiments/aim3")
}

source(file.path(script_directory, "analysis.R"))
input_file <- file.path(
  script_directory, "result", "synthetic", "exp1", "synthetic_key_results.csv"
)
aim3_analyze_results(input_file, "synthetic")
