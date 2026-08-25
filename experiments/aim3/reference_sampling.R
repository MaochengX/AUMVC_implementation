aim3_reference_size <- function(box, settings) {
  if (
    settings$points_per_volume <= 0 ||
    settings$min_points < 1L ||
    settings$max_points < settings$min_points
  ) {
    stop("Invalid Aim 3 reference-density settings.", call. = FALSE)
  }

  log_target <- log(settings$points_per_volume) + box$log_volume
  log_min <- log(settings$min_points)
  log_max <- log(settings$max_points)

  if (log_target <= log_min) {
    return(list(n = as.integer(settings$min_points), capped = FALSE))
  }

  if (log_target >= log_max) {
    return(list(n = as.integer(settings$max_points), capped = TRUE))
  }

  list(
    n = as.integer(ceiling(exp(log_target))),
    capped = FALSE
  )
}

aim3_make_reference <- function(x_reference, seed, settings) {
  box <- fit_reference_box(x_reference)
  size <- aim3_reference_size(box, settings)

  reference <- make_reference(
    x_reference,
    n_reference = size$n,
    n_mc_repetitions = settings$n_mc_repetitions,
    seed = seed
  )

  reference$density_capped <- size$capped
  reference$log_density <- log(reference$n_reference) -
    reference$box$log_volume

  reference
}
