source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")

python <- Sys.which("python3")
if (!nzchar(python)) python <- Sys.which("python")
if (!nzchar(python)) stop("Python was not found.", call. = FALSE)

version_command <- paste(
  "import importlib.metadata;",
  "print(importlib.metadata.version('emmv'))"
)
version <- suppressWarnings(system2(
  python,
  c("-c", shQuote(version_command)),
  stdout = TRUE,
  stderr = TRUE
))
status <- attr(version, "status")
if (is.null(status)) status <- 0L
version <- if (status == 0L) trimws(tail(version, 1L)) else NA_character_

if (is.na(version) || version != "0.0.4") {
  status <- system2(
    python,
    c("-m", "pip", "install", "--upgrade", "emmv==0.0.4")
  )
  if (status != 0L) stop("Could not install emmv==0.0.4.", call. = FALSE)
}

version <- system2(
  python,
  c("-c", shQuote(version_command)),
  stdout = TRUE
)
version <- trimws(tail(version, 1L))
if (version != "0.0.4") stop("emmv 0.0.4 is required.", call. = FALSE)

alpha_grid <- seq(0.9, 0.999, by = 0.0001)
n_reference <- 100000L

validation_box <- function(x) {
  lower <- apply(x, 2L, min)
  upper <- apply(x, 2L, max)
  width <- upper - lower
  list(
    lower = lower,
    width = width,
    volume = prod(width)
  )
}

validation_reference <- function(box, n) {
  x <- matrix(
    runif(n * length(box$lower)),
    nrow = n,
    byrow = TRUE
  )
  x <- sweep(x, 2L, box$width, "*")
  sweep(x, 2L, box$lower, "+")
}

validation_case <- function(name, seed) {
  set.seed(seed)

  if (name == "gaussian_2d") {
    x_box <- matrix(rnorm(600L), ncol = 2L)
    x_eval <- matrix(rnorm(802L), ncol = 2L)
    score <- function(x) -rowSums(x^2)
  } else if (name == "two_clusters_2d") {
    x_box <- rbind(
      cbind(rnorm(150L, -2, 0.7), rnorm(150L, -2, 0.7)),
      cbind(rnorm(150L, 2, 0.7), rnorm(150L, 2, 0.7))
    )
    x_eval <- rbind(
      cbind(rnorm(201L, -2, 0.7), rnorm(201L, -2, 0.7)),
      cbind(rnorm(202L, 2, 0.7), rnorm(202L, 2, 0.7))
    )
    centers <- rbind(c(-2, -2), c(2, 2))
    score <- function(x) {
      d1 <- rowSums(sweep(x, 2L, centers[1L, ], "-")^2)
      d2 <- rowSums(sweep(x, 2L, centers[2L, ], "-")^2)
      -pmin(d1, d2)
    }
  } else if (name == "elliptical_3d") {
    scale <- c(1, 0.6, 1.8)
    x_box <- sweep(matrix(rnorm(1050L), ncol = 3L), 2L, scale, "*")
    x_eval <- sweep(matrix(rnorm(1227L), ncol = 3L), 2L, scale, "*")
    score <- function(x) -rowSums(sweep(x, 2L, scale, "/")^2)
  } else {
    stop("Unknown validation case.", call. = FALSE)
  }

  box <- validation_box(x_box)
  x_reference <- validation_reference(box, n_reference)

  list(
    evaluation_scores = score(x_eval),
    reference_scores = score(x_reference),
    box_volume = box$volume
  )
}

python_script <- tempfile(fileext = ".py")
writeLines(c(
  "import sys",
  "import numpy as np",
  "from emmv import mass_volume",
  "alpha = np.loadtxt(sys.argv[1])",
  "reference_scores = np.loadtxt(sys.argv[2])",
  "evaluation_scores = np.loadtxt(sys.argv[3])",
  "volume = float(sys.argv[4])",
  "auc_value, curve = mass_volume(",
  "    alpha, volume, reference_scores, evaluation_scores, len(reference_scores)",
  ")",
  "np.savetxt(sys.argv[5], curve, fmt='%.17g')",
  "np.savetxt(sys.argv[6], [auc_value], fmt='%.17g')"
), python_script)

write_numbers <- function(values, path) {
  writeLines(sprintf("%.17g", as.numeric(values)), path)
}

run_python_reference <- function(case) {
  files <- replicate(5L, tempfile(fileext = ".txt"))
  on.exit(unlink(files), add = TRUE)
  write_numbers(alpha_grid, files[1L])
  write_numbers(case$reference_scores, files[2L])
  write_numbers(case$evaluation_scores, files[3L])

  output <- suppressWarnings(system2(
    python,
    c(
      shQuote(python_script),
      shQuote(files[1L]),
      shQuote(files[2L]),
      shQuote(files[3L]),
      format(case$box_volume, digits = 17L, scientific = TRUE),
      shQuote(files[4L]),
      shQuote(files[5L])
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L) {
    cat(paste(output, collapse = "\n"), "\n")
    stop("Python emmv calculation failed.", call. = FALSE)
  }

  list(
    curve = scan(files[4L], quiet = TRUE),
    aumvc = scan(files[5L], quiet = TRUE)[1L]
  )
}

validate_case <- function(name, seed) {
  case <- validation_case(name, seed)
  python_reference <- run_python_reference(case)
  r_result <- aumvc_from_scores(
    case$evaluation_scores,
    case$reference_scores,
    case$box_volume,
    "normality",
    alpha_grid
  )
  direction_result <- aumvc_from_scores(
    -case$evaluation_scores,
    -case$reference_scores,
    case$box_volume,
    "anomaly",
    alpha_grid
  )

  curve_difference <- max(abs(r_result$mv_curve$volume - python_reference$curve))
  aumvc_difference <- abs(r_result$aumvc - python_reference$aumvc)
  normalized_difference <- abs(
    r_result$aumvc_normalized - python_reference$aumvc / case$box_volume
  )
  direction_difference <- max(
    abs(r_result$mv_curve$volume - direction_result$mv_curve$volume),
    abs(r_result$aumvc - direction_result$aumvc)
  )

  data.frame(
    case = name,
    r_aumvc = r_result$aumvc,
    python_aumvc = python_reference$aumvc,
    curve_difference = curve_difference,
    aumvc_difference = aumvc_difference,
    normalized_difference = normalized_difference,
    direction_difference = direction_difference,
    passed = max(
      curve_difference,
      aumvc_difference,
      normalized_difference,
      direction_difference
    ) <= 1e-10
  )
}

cases <- c("gaussian_2d", "two_clusters_2d", "elliptical_3d")
results <- do.call(rbind, Map(
  validate_case,
  cases,
  c(1234L, 2234L, 3234L)
))

cat("Python: ", python, "\n", sep = "")
cat("emmv version: ", version, "\n\n", sep = "")
print(results, row.names = FALSE, digits = 10)

if (!all(results$passed)) {
  stop("Python emmv validation failed.", call. = FALSE)
}

cat("\nPython emmv validation passed.\n")
