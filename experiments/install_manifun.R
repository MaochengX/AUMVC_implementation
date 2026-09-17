if (!requireNamespace("manifun", quietly = TRUE)) {
  if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak", repos = "https://cloud.r-project.org")
  }
  pak::pkg_install("HerrMo/manifun")
}
stopifnot(requireNamespace("manifun", quietly = TRUE))
cat("manifun version:", as.character(utils::packageVersion("manifun")), "\n")
