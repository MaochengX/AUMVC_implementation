source("experiments/settings.R")
source("experiments/aim3/synthetic.R")
source("experiments/aim3/representation_evaluation.R")

settings <- AIM3_SETTINGS
settings$goix_subsampling$n_subsets <- 3L
settings$reference$min_points <- 1000L
settings$reference$max_points <- 1000L
settings$reference$n_mc_repetitions <- 2L

intrinsic_dim <- 3L

case <- aim3_generate_synthetic(
  intrinsic_dim = intrinsic_dim,
  snr = 5,
  truth = "polynomial_interaction",
  seed = settings$seed,
  settings = settings$synthetic
)

stopifnot(nrow(case$x) == settings$synthetic$n)
stopifnot(ncol(case$x) == settings$synthetic$ambient_dim)
stopifnot(all(is.finite(case$x)))

split <- aim3_make_splits(
  nrow(case$x),
  settings$seed,
  settings
)

standardizer <- aim3_fit_standardizer(
  case$x[split$embedding, , drop = FALSE]
)

landmarks <- aim3_apply_standardizer(
  case$x[split$embedding, , drop = FALSE],
  standardizer
)

newdata <- aim3_apply_standardizer(
  case$x[split$detector_train[1:10], , drop = FALSE],
  standardizer
)

mds <- aim3_fit_mds(
  landmarks,
  ndim = intrinsic_dim,
  settings = settings$mds
)

projected <- aim3_project_mds(
  mds,
  newdata
)

stopifnot(all(dim(projected) == c(10L, intrinsic_dim)))
stopifnot(all(is.finite(projected)))

cat("manifun MDS test passed\n")

result <- aim3_evaluate_representations(
  case,
  seed = settings$seed + 500L,
  settings = settings
)

stopifnot(nrow(result) == 3L)
stopifnot(
  setequal(
    result$representation,
    c(
      "ambient",
      "goix_subsampling",
      "mds"
    )
  )
)
stopifnot(all(is.finite(result$aumvc_normalized)))
stopifnot(all(result$aumvc_normalized >= 0))
stopifnot(all(result$n_reference >= 1L))
stopifnot(all(is.finite(result$runtime_seconds)))

cat("Aim 3 synthetic smoke test passed\n")
