aim3_unit_rows <- function(x) {
  x <- as.matrix(x)
  norms <- sqrt(rowSums(x^2))
  norms[norms == 0] <- 1
  x / norms
}

aim3_next_coordinate <- function(q) {
  c(seq.int(2L, q), 1L)
}

aim3_polynomial_basis <- function(z) {
  z <- as.matrix(z)
  next_index <- aim3_next_coordinate(ncol(z))
  cbind(z, z^2, z^3, z * z[, next_index, drop = FALSE])
}

aim3_oscillatory_basis <- function(z) {
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

aim3_basis <- function(z, truth) {
  switch(
    truth,
    polynomial_interaction = aim3_polynomial_basis(z),
    oscillatory_local = aim3_oscillatory_basis(z),
    stop("Unknown Aim 3 truth function.", call. = FALSE)
  )
}

aim3_basis_names <- function(q, truth) {
  z <- paste0("z", seq_len(q))
  next_z <- paste0("z", aim3_next_coordinate(q))

  switch(
    truth,
    polynomial_interaction = c(
      z,
      paste0(z, "^2"),
      paste0(z, "^3"),
      paste0(z, "*", next_z)
    ),
    oscillatory_local = c(
      paste0("sin(", z, ")"),
      paste0("cos(2*", z, ")"),
      paste0("tanh(1.5*", z, ")"),
      paste0("exp(-", z, "^2/2)"),
      paste0("sin(", z, "*", next_z, ")")
    ),
    stop("Unknown Aim 3 truth function.", call. = FALSE)
  )
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
  matrix(rnorm(n_basis * ambient_dim), nrow = n_basis) / sqrt(n_basis)
}

aim3_distributional_latent <- function(n, q, tail, seed) {
  set.seed(seed)
  direction <- aim3_unit_rows(matrix(rnorm(n * q), nrow = n))
  radius <- sqrt(qchisq(runif(n, tail[1L], tail[2L]), df = q))
  direction * radius
}

aim3_structural_directions <- function(n, coefficients, seed) {
  set.seed(seed)
  basis_space <- qr.Q(qr(t(coefficients)))
  direction <- matrix(rnorm(n * ncol(coefficients)), nrow = n)
  direction <- direction - direction %*% basis_space %*% t(basis_space)
  aim3_unit_rows(direction)
}

aim3_evaluate_function <- function(z, generated) {
  z <- as.matrix(z)

  if (ncol(z) != generated$intrinsic_dim) {
    stop("Latent dimension does not match the generated function.", call. = FALSE)
  }

  phi <- aim3_apply_basis_scaler(
    aim3_basis(z, generated$truth),
    generated$basis_scaler
  )

  phi %*% generated$coefficients
}

aim3_show_generated_function <- function(generated, output_dim = 1L, digits = 4L) {
  output_dim <- as.integer(output_dim)

  if (output_dim < 1L || output_dim > ncol(generated$coefficients)) {
    stop("output_dim is outside the ambient dimension.", call. = FALSE)
  }

  terms <- aim3_basis_names(generated$intrinsic_dim, generated$truth)
  coefficients <- generated$coefficients[, output_dim]
  format_number <- function(x) formatC(x, digits = digits, format = "g")

  pieces <- vapply(seq_along(coefficients), function(i) {
    paste0(
      if (coefficients[i] < 0) "- " else "+ ",
      format_number(abs(coefficients[i])),
      " * ((",
      terms[i],
      " - ",
      format_number(generated$basis_scaler$center[i]),
      ") / ",
      format_number(generated$basis_scaler$scale[i]),
      ")"
    )
  }, character(1))

  pieces[1L] <- sub("^\\+ ", "", pieces[1L])
  cat("f", output_dim, "(z) = ", paste(pieces, collapse = " "), "\n", sep = "")
  invisible(NULL)
}

aim3_generate_synthetic <- function(intrinsic_dim, snr, truth, seed, settings) {
  q <- as.integer(intrinsic_dim)
  d <- as.integer(settings$ambient_dim)
  n_distributional <- as.integer(settings$n_distributional)
  n_structural <- as.integer(settings$n_structural)
  n_normal <- as.integer(settings$n) - n_distributional - n_structural

  if (n_normal < 1L) stop("Invalid Aim 3 synthetic sample sizes.", call. = FALSE)
  if (!q %in% settings$intrinsic_dims) stop("Unsupported intrinsic dimension.", call. = FALSE)
  if (!snr %in% settings$snr_levels) stop("Unsupported SNR level.", call. = FALSE)
  if (!truth %in% settings$truth_functions) stop("Unsupported truth function.", call. = FALSE)

  set.seed(seed)
  z_normal <- matrix(rnorm(n_normal * q), nrow = n_normal)
  z_distributional <- aim3_distributional_latent(
    n_distributional,
    q,
    settings$distributional_tail,
    seed + 1L
  )

  set.seed(seed + 2L)
  z_structural <- matrix(rnorm(n_structural * q), nrow = n_structural)

  phi_normal <- aim3_basis(z_normal, truth)
  basis_scaler <- aim3_fit_basis_scaler(phi_normal)
  phi_normal <- aim3_apply_basis_scaler(phi_normal, basis_scaler)
  phi_distributional <- aim3_apply_basis_scaler(
    aim3_basis(z_distributional, truth),
    basis_scaler
  )
  phi_structural <- aim3_apply_basis_scaler(
    aim3_basis(z_structural, truth),
    basis_scaler
  )

  coefficients <- aim3_make_coefficients(ncol(phi_normal), d, seed + 3L)
  signal_normal <- phi_normal %*% coefficients
  signal_distributional <- phi_distributional %*% coefficients
  signal_structural <- phi_structural %*% coefficients

  signal_variance <- mean(apply(signal_normal, 2L, var))
  noise_sd <- sqrt(signal_variance / snr)
  structural_direction <- aim3_structural_directions(
    n_structural,
    coefficients,
    seed + 4L
  )

  signal_structural <- signal_structural +
    settings$structural_shift * sqrt(d) * noise_sd * structural_direction

  set.seed(seed + 5L + as.integer(100 * snr))
  noise <- matrix(
    rnorm(settings$n * d, sd = noise_sd),
    nrow = settings$n,
    ncol = d
  )

  outlier_type <- c(
    rep("normal", n_normal),
    rep("distributional", n_distributional),
    rep("structural", n_structural)
  )

  list(
    x = rbind(signal_normal, signal_distributional, signal_structural) + noise,
    outlier_type = outlier_type,
    truth = truth,
    intrinsic_dim = q,
    basis_scaler = basis_scaler,
    coefficients = coefficients
  )
}
