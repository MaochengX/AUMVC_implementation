structural_basis_v3 <- local({
  cache <- new.env(parent = emptyenv())

  perpendicular_unit <- function(vector, x) {
    direction <- matrix(rep(vector, each = nrow(x)), nrow = nrow(x))
    direction <- direction - x * (
      rowSums(direction * x) / pmax(rowSums(x^2), 1e-12)
    )
    direction / pmax(sqrt(rowSums(direction^2)), 1e-12)
  }

  function(z, w, truth, severity) {
    if (!identical(dim(z), dim(w)) || length(severity) != 1L ||
        !is.finite(severity) || severity < 0 || severity > 1) {
      stop("Invalid structural settings.", call. = FALSE)
    }

    key <- paste(truth, ncol(z), sep = "_")
    if (!exists(key, envir = cache, inherits = FALSE)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        previous_seed <- get(".Random.seed", envir = .GlobalEnv)
      } else {
        previous_seed <- NULL
      }
      set.seed(3200L + ncol(z) + if (truth == "oscillatory_local") 100L else 0L)
      ref_z <- matrix(rnorm(2000L * ncol(z)), ncol = ncol(z))
      if (is.null(previous_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      }

      ref <- basis(ref_z, truth)
      center <- colMeans(ref)
      spread <- apply(ref, 2L, sd)
      spread[!is.finite(spread) | spread < 1e-10] <- 1
      ref <- sweep(sweep(ref, 2L, center, "-"), 2L, spread, "/")
      eigenvectors <- eigen(crossprod(ref), symmetric = TRUE)$vectors
      p <- ncol(ref)
      q <- ncol(z)
      assign(key, list(center = center, spread = spread,
                       structural = eigenvectors[, p],
                       visible = eigenvectors[, q]), envir = cache)
    }

    model <- get(key, envir = cache)
    phi <- basis(z, truth)
    standardized <- sweep(sweep(phi, 2L, model$center, "-"),
                          2L, model$spread, "/")
    structural_direction <- perpendicular_unit(model$structural, standardized)
    visible_direction <- perpendicular_unit(model$visible, standardized)
    visible_direction <- visible_direction * ifelse(z[, 1L] >= 0, 1, -1)
    shift <- sqrt(ncol(standardized)) * severity * (
      0.45 * structural_direction + 0.50 * visible_direction
    )
    modified <- standardized + shift
    radius <- sqrt(rowSums(standardized^2))
    modified_radius <- pmax(sqrt(rowSums(modified^2)), 1e-12)
    modified <- modified * (radius / modified_radius)^0.25

    sweep(sweep(modified, 2L, model$spread, "*"),
          2L, model$center, "+")
  }
})
