cran_repo <- "https://cloud.r-project.org"

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = cran_repo)
}

packages <- c(
  "Rtsne",
  "data.table",
  "bioc::diffusionMap",
  "ggplot2",
  "igraph",
  "shiny",
  "umap",
  "vegan",
  "HerrMo/manifun"
)

pak::pkg_install(packages)

if (!requireNamespace("manifun", quietly = TRUE)) {
  stop("manifun installation failed.", call. = FALSE)
}

required <- c("embed", "extract_points")
available <- vapply(
  required,
  exists,
  logical(1),
  where = asNamespace("manifun"),
  inherits = FALSE
)

if (!all(available)) {
  stop("manifun is installed but required functions are missing.", call. = FALSE)
}

cat(
  "manifun installed successfully:",
  as.character(utils::packageVersion("manifun")),
  "\n"
)
