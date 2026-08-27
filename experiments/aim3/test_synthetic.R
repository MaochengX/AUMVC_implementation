source("experiments/settings.R")
source("experiments/aim3/synthetic.R")
source("experiments/aim3/representation_comparison.R")

settings <- AIM3_SETTINGS
settings$goix_subsampling$n_subsets <- 3L
settings$n_reference <- 1000L
settings$n_mc_repetitions <- 2L
q <- 3L

z_test <- matrix(
  c(
    -1, 0.5, 1.2,
    0.2, -0.7, 0.4,
    1.1, 0.3, -0.8
  ),
  nrow = 3L,
  byrow = TRUE
)

stopifnot(
  ncol(aim3_basis(z_test, "polynomial_interaction")) == 4L * q,
  ncol(aim3_basis(z_test, "oscillatory_local")) == 5L * q
)

case <- aim3_generate_synthetic(
  q,
  5,
  "polynomial_interaction",
  settings$seed,
  settings$synthetic
)

aim3_show_generated_function(case, output_dim = 1L)

manual <- aim3_apply_basis_scaler(
  aim3_basis(z_test, case$truth),
  case$basis_scaler
) %*% case$coefficients

stopifnot(
  max(abs(manual - aim3_evaluate_function(z_test, case))) < 1e-10,
  nrow(case$x) == settings$synthetic$n,
  ncol(case$x) == settings$synthetic$ambient_dim,
  sum(case$outlier_type == "distributional") == settings$synthetic$n_distributional,
  sum(case$outlier_type == "structural") == settings$synthetic$n_structural,
  all(is.finite(case$x))
)

split <- make_splits(
  nrow(case$x),
  settings$split_counts,
  names(settings$split_counts),
  settings$seed,
  require_full = TRUE
)

standardizer <- fit_standardizer(
  case$x[split$embedding, , drop = FALSE]
)
landmarks <- apply_standardizer(
  case$x[split$embedding, , drop = FALSE],
  standardizer
)

mds <- aim3_fit_mds(landmarks, q)
reprojected <- aim3_project_mds(
  mds,
  landmarks[1:5, , drop = FALSE]
)

stopifnot(
  identical(dim(mds$points), c(length(split$embedding), q)),
  max(abs(reprojected - mds$points[1:5, , drop = FALSE])) < 1e-6
)

comparison <- aim3_compare_representations(
  case,
  settings$seed + 500L,
  settings
)

stopifnot(
  nrow(comparison) == 3L,
  setequal(
    comparison$representation,
    c("ambient", "goix_subsampling", "mds")
  ),
  all(is.finite(comparison$aumvc_normalized)),
  all(is.finite(comparison$aumvc_normalized_mc_se)),
  all(comparison$zero_occupancy >= 0),
  all(comparison$zero_occupancy <= 1),
  all(is.finite(comparison$runtime_seconds))
)

cat("Aim 3 synthetic smoke test passed\n")
