source("experiments/settings.R")

aim3_row_unit <- function(x) {
  x <- as.matrix(x)
  norms <- sqrt(rowSums(x^2))
  norms[norms == 0] <- 1
  x / norms
}

aim3_next_coordinate <- function(q) {
  c(seq.int(2L, q), 1L)
}

aim3_polynomial_interaction_basis <- function(z) {
  z <- as.matrix(z)
  next_index <- aim3_next_coordinate(ncol(z))

  cbind(
    z,
    z^2,
    z^3,
    z * z[, next_index, drop = FALSE]
  )
}

aim3_oscillatory_local_basis <- function(z) {
  z <- as.matrix(z)
  next_index <- aim3_next_coordinate(ncol(z))

  cbind(
    sin(z),
    cos(2 * z),
    tanh(1.5 * z),
    exp(-(z^2) / 2),
    sin(z * z[, next_index, drop = FALSE])
  )
}

aim3_make_basis <- function(z, truth) {
  switch(
    truth,
    polynomial_interaction = aim3_polynomial_interaction_basis(z),
    oscillatory_local = aim3_oscillatory_local_basis(z),
    stop("Unknown Aim 3 truth function.", call. = FALSE)
  )
}

aim3_basis_names <- function(q, truth) {
  z_names <- paste0("z", seq_len(q))
  next_names <- paste0("z", aim3_next_coordinate(q))

  if (truth == "polynomial_interaction") {
    return(c(
      z_names,
      paste0(z_names, "^2"),
      paste0(z_names, "^3"),
      paste0(z_names, "*", next_names)
    ))
  }

  if (truth == "oscillatory_local") {
    return(c(
      paste0("sin(", z_names, ")"),
      paste0("cos(2*", z_names, ")"),
      paste0("tanh(1.5*", z_names, ")"),
      paste0("exp(-", z_names, "^2/2)"),
      paste0("sin(", z_names, "*", next_names, ")")
    ))
  }

  stop("Unknown Aim 3 truth function.", call. = FALSE)
}

aim3_show_basis_function <- function(intrinsic_dim, truth) {
  terms <- aim3_basis_names(as.integer(intrinsic_dim), truth)
  cat("truth =", truth, "\n")
  cat("phi(z) = [", paste(terms, collapse = ", "), "]\n", sep = "")
  invisible(terms)
}

aim3_fit_basis_scaler <- function(phi) {
  center <- colMeans(phi)
  scale <- apply(phi, 2L, sd)
  scale[!is.finite(scale) | scale == 0] <- 1

  list(center = center, scale = scale)
}

aim3_apply_basis_scaler <- function(phi, scaler) {
  phi <- sweep(phi, 2L, scaler$center, "-")
  sweep(phi, 2L, scaler$scale, "/")
}

aim3_make_coefficients <- function(n_basis, ambient_dim, seed) {
  set.seed(seed)

  matrix(
    rnorm(n_basis * ambient_dim),
    nrow = n_basis,
    ncol = ambient_dim
  ) / sqrt(n_basis)
}

aim3_make_distributional_latent <- function(n, q, tail, seed) {
  set.seed(seed)

  direction <- aim3_row_unit(
    matrix(rnorm(n * q), nrow = n, ncol = q)
  )
  probability <- runif(n, tail[1L], tail[2L])
  radius <- sqrt(qchisq(probability, df = q))

  direction * radius
}

aim3_make_structural_directions <- function(n, coefficients, seed) {
  set.seed(seed)

  ambient_dim <- ncol(coefficients)
  basis_space <- qr.Q(qr(t(coefficients)))
  direction <- matrix(rnorm(n * ambient_dim), nrow = n)
  direction <- direction - direction %*% basis_space %*% t(basis_space)

  aim3_row_unit(direction)
}

aim3_evaluate_generated_function <- function(z, generated) {
  z <- as.matrix(z)

  if (ncol(z) != generated$intrinsic_dim) {
    stop("Latent dimension does not match the generated function.", call. = FALSE)
  }

  phi <- aim3_make_basis(z, generated$truth)
  phi <- aim3_apply_basis_scaler(phi, generated$basis_scaler)

  phi %*% generated$coefficients
}

aim3_show_generated_function <- function(
    generated,
    output_dims = 1:3,
    digits = 4L
) {
  output_dims <- as.integer(output_dims)

  if (any(output_dims < 1L | output_dims > generated$ambient_dim)) {
    stop("output_dims is outside the ambient dimension.", call. = FALSE)
  }

  terms <- aim3_basis_names(generated$intrinsic_dim, generated$truth)
  center <- generated$basis_scaler$center
  scale <- generated$basis_scaler$scale
  number_format <- paste0("%.", as.integer(digits), "g")

  cat("truth =", generated$truth, "\n")
  cat("intrinsic dimension =", generated$intrinsic_dim, "\n")
  cat("ambient dimension =", generated$ambient_dim, "\n\n")

  for (j in output_dims) {
    coefficients <- generated$coefficients[, j]
    pieces <- character(length(coefficients))

    for (k in seq_along(coefficients)) {
      sign_text <- if (coefficients[k] < 0) "-" else "+"
      coefficient_text <- sprintf(number_format, abs(coefficients[k]))
      center_text <- sprintf(number_format, center[k])
      scale_text <- sprintf(number_format, scale[k])

      term <- sprintf(
        "((%s - %s) / %s)",
        terms[k],
        center_text,
        scale_text
      )

      pieces[k] <- sprintf(
        "%s %s * %s",
        sign_text,
        coefficient_text,
        term
      )
    }

    pieces[1L] <- sub("^\\+ ", "", pieces[1L])
    cat("f", j, "(z) = ", paste(pieces, collapse = " "), "\n\n", sep = "")
  }

  invisible(NULL)
}

aim3_generate_synthetic <- function(
    intrinsic_dim,
    snr,
    truth,
    seed,
    settings = AIM3_SETTINGS$synthetic
) {
  intrinsic_dim <- as.integer(intrinsic_dim)
  ambient_dim <- as.integer(settings$ambient_dim)
  n_distributional <- as.integer(settings$n_distributional)
  n_structural <- as.integer(settings$n_structural)
  n_normal <- as.integer(settings$n) - n_distributional - n_structural

  if (n_normal < 1L) {
    stop("Aim 3 synthetic sample sizes are invalid.", call. = FALSE)
  }

  if (!intrinsic_dim %in% settings$intrinsic_dims) {
    stop("Unsupported intrinsic dimension.", call. = FALSE)
  }

  if (!snr %in% settings$snr_levels) {
    stop("Unsupported SNR level.", call. = FALSE)
  }

  if (!truth %in% settings$truth_functions) {
    stop("Unsupported truth function.", call. = FALSE)
  }

  set.seed(seed)
  z_normal <- matrix(
    rnorm(n_normal * intrinsic_dim),
    nrow = n_normal,
    ncol = intrinsic_dim
  )

  z_distributional <- aim3_make_distributional_latent(
    n_distributional,
    intrinsic_dim,
    settings$distributional_tail,
    seed + 1L
  )

  set.seed(seed + 2L)
  z_structural <- matrix(
    rnorm(n_structural * intrinsic_dim),
    nrow = n_structural,
    ncol = intrinsic_dim
  )

  phi_normal <- aim3_make_basis(z_normal, truth)
  basis_scaler <- aim3_fit_basis_scaler(phi_normal)

  phi_normal <- aim3_apply_basis_scaler(phi_normal, basis_scaler)
  phi_distributional <- aim3_apply_basis_scaler(
    aim3_make_basis(z_distributional, truth),
    basis_scaler
  )
  phi_structural <- aim3_apply_basis_scaler(
    aim3_make_basis(z_structural, truth),
    basis_scaler
  )

  coefficients <- aim3_make_coefficients(
    ncol(phi_normal),
    ambient_dim,
    seed + 3L
  )

  signal_normal <- phi_normal %*% coefficients
  signal_distributional <- phi_distributional %*% coefficients
  signal_structural <- phi_structural %*% coefficients

  signal_variance <- mean(apply(signal_normal, 2L, var))
  noise_sd <- sqrt(signal_variance / snr)

  structural_direction <- aim3_make_structural_directions(
    n_structural,
    coefficients,
    seed + 4L
  )

  structural_distance <- settings$structural_shift *
    sqrt(ambient_dim) *
    noise_sd

  signal_structural <- signal_structural +
    structural_distance * structural_direction

  set.seed(seed + 5L + as.integer(100 * snr))
  noise <- matrix(
    rnorm(settings$n * ambient_dim, sd = noise_sd),
    nrow = settings$n,
    ncol = ambient_dim
  )

  signal <- rbind(
    signal_normal,
    signal_distributional,
    signal_structural
  )

  outlier_type <- c(
    rep("normal", n_normal),
    rep("distributional", n_distributional),
    rep("structural", n_structural)
  )

  list(
    x = signal + noise,
    labels = as.integer(outlier_type != "normal"),
    outlier_type = outlier_type,
    truth = truth,
    intrinsic_dim = intrinsic_dim,
    ambient_dim = ambient_dim,
    snr = snr,
    noise_sd = noise_sd,
    structural_distance = structural_distance,
    basis_scaler = basis_scaler,
    coefficients = coefficients
  )
}
