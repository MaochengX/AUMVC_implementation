source("experiments/aim3/synthetic_experiment.R")
source("experiments/aim3/real_experiment.R")

SCAN_SETTINGS <- AIM3_SETTINGS
SCAN_SETTINGS$bandwidth$gamma_scales <- c(0.25, 1, 4)
SCAN_SETTINGS$synthetic$ambient_dims <- 200L
SCAN_SETTINGS$synthetic$snr_levels <- 5
SCAN_RUNS <- as.integer(Sys.getenv("SCAN_RUNS", "10"))
SCAN_PARTS <- trimws(strsplit(Sys.getenv("SCAN_PARTS", "synthetic,orthogonal,real"), ",", fixed = TRUE)[[1L]])
SCAN_OUTPUT <- Sys.getenv("SCAN_OUTPUT", "experiments/aim3/dim_scan")
SCAN_ORTHOGONAL_SHIFT <- 3

scan_dimensions <- function(source, q = NULL, dataset = NULL) {
  if (source != "real") return(as.integer(q * c(1, 2, 4, 8)))
  if (dataset == "shuttle") return(c(2L, 3L, 5L, 8L))
  c(5L, 10L, 20L, 40L)
}

original_orthogonal_case <- function(row, run, settings) {
  q <- as.integer(row$intrinsic_dim)
  d <- 200L
  snr <- 5
  truth <- as.character(row$truth)
  n_normal <- settings$synthetic$n - settings$synthetic$n_distributional -
    settings$synthetic$n_structural
  n_distributional <- settings$synthetic$n_distributional
  n_structural <- settings$synthetic$n_structural
  seed <- settings$seed + run * 100000L +
    match(truth, settings$synthetic$truth_functions) * 10000L + q * 100L + as.integer(snr)

  set.seed(seed)
  z_normal <- matrix(rnorm(n_normal * q), nrow = n_normal)
  z_distributional <- distributional_latent(
    n_distributional, q, settings$synthetic$distributional_tail, seed + 1L
  )
  set.seed(seed + 2L)
  z_structural <- matrix(rnorm(n_structural * q), nrow = n_structural)

  normal_basis <- basis(z_normal, truth)
  scaler <- fit_basis_scaler(normal_basis)
  normal_basis <- apply_basis_scaler(normal_basis, scaler)
  distributional_basis <- apply_basis_scaler(basis(z_distributional, truth), scaler)
  structural_basis_values <- apply_basis_scaler(basis(z_structural, truth), scaler)

  set.seed(seed + 3L)
  coefficients <- matrix(rnorm(ncol(normal_basis) * d), nrow = ncol(normal_basis)) /
    sqrt(ncol(normal_basis))
  signal_normal <- normal_basis %*% coefficients
  signal_distributional <- distributional_basis %*% coefficients
  signal_structural <- structural_basis_values %*% coefficients
  noise_sd <- sqrt(mean(apply(signal_normal, 2L, var)) / snr)

  set.seed(seed + 4L)
  signal_space <- qr.Q(qr(t(coefficients)))
  directions <- matrix(rnorm(n_structural * d), nrow = n_structural)
  directions <- directions - directions %*% signal_space %*% t(signal_space)
  signal_structural <- signal_structural +
    SCAN_ORTHOGONAL_SHIFT * sqrt(d) * noise_sd * unit_rows(directions)

  set.seed(seed + 5L + as.integer(100 * snr))
  noise <- matrix(rnorm(settings$synthetic$n * d, sd = noise_sd),
                  nrow = settings$synthetic$n, ncol = d)
  types <- c(rep("normal", n_normal), rep("distributional", n_distributional),
             rep("structural", n_structural))
  list(
    case = list(x = rbind(signal_normal, signal_distributional, signal_structural) + noise,
                labels = as.integer(types != "normal"), outlier_type = types,
                intrinsic_dim = q),
    seed = seed
  )
}

scan_one_case <- function(case, seed, dims, settings, metadata) {
  case$intrinsic_dim <- min(dims)
  prepared <- prepare_case_data(case, seed, settings)
  maximum <- min(ncol(prepared$embedding), nrow(prepared$embedding) - 1L)
  if (any(dims > maximum)) {
    stop("Requested dimension exceeds the available embedding rank bound: ", maximum)
  }
  rows <- list()
  ambient <- evaluate_grid(
    "ambient", prepared$train, prepared$reference, prepared$eval, prepared$label,
    prepared$labels_label, prepared$types_label, seed + 1000L, settings
  )
  rows[[1L]] <- attach_metadata(metadata, ambient)
  for (dimension in dims) {
    current <- prepared
    current$embedding_dim <- as.integer(dimension)
    embedded <- prepare_embedding_data(current, "mds", settings)
    measured <- evaluate_grid(
      "mds", embedded$train, embedded$reference, embedded$eval, embedded$label,
      prepared$labels_label, prepared$types_label, seed + 3000L, settings,
      embedded$embedding_oos_stress, embedded$setup_runtime
    )
    rows[[length(rows) + 1L]] <- attach_metadata(metadata, measured)
  }
  do.call(rbind, rows)
}

scan_synthetic_case <- function(row, run, source, settings) {
  if (source == "orthogonal") {
    original <- original_orthogonal_case(row, run, settings)
    case <- original$case
    seed <- original$seed
  } else {
    seed <- synthetic_seed(row, run, settings$synthetic, settings$seed)
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    case <- generate_synthetic(row$intrinsic_dim, 200L, 5, row$truth, seed, settings$synthetic)
  }
  dimensions <- scan_dimensions(source, row$intrinsic_dim)
  metadata <- data.frame(
    source = source, dataset = "synthetic", run = run, truth = row$truth,
    intrinsic_dim = row$intrinsic_dim, snr = 5, sample_size = nrow(case$x),
    original_dim = ncol(case$x), stringsAsFactors = FALSE
  )
  scan_one_case(case, seed + 50000L, dimensions, settings, metadata)
}

scan_real_case <- function(dataset, index, run, loaded, registry, settings) {
  seed <- settings$seed + run * 100000L + index * 10000L
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  sampled <- registry[[dataset]]$sample(loaded[[dataset]], seed, settings$real[[dataset]])
  dimensions <- scan_dimensions("real", dataset = dataset)
  current <- settings
  current$split_counts <- real_split_counts(nrow(sampled$x))
  case <- list(x = sampled$x, labels = sampled$labels, intrinsic_dim = min(dimensions))
  metadata <- data.frame(
    source = "real", dataset = dataset, run = run, truth = NA_character_,
    intrinsic_dim = NA_integer_, snr = NA_real_, sample_size = nrow(sampled$x),
    original_dim = ncol(sampled$x), stringsAsFactors = FALSE
  )
  scan_one_case(case, seed + 50000L, dimensions, current, metadata)
}

scan_average <- function(x, column, gamma, representation, dimension = NULL) {
  keep <- x$gamma_scale == gamma & x$representation == representation
  if (!is.null(dimension)) keep <- keep & !is.na(x$embedding_dim) & x$embedding_dim == dimension
  values <- x[[column]][keep]
  values <- as.numeric(values)
  if (!length(values) || all(is.na(values))) return(NA_real_)
  mean(values, na.rm = TRUE)
}

scan_panel <- function(data, metric, dimensions, heading, x_labels) {
  gamma <- SCAN_SETTINGS$bandwidth$gamma_scales
  colors <- c("#2878B5", "#009E73", "#D55E00")
  plot(NA_real_, xlim = range(dimensions), ylim = c(0, 1), xaxt = "n",
       xlab = "MDS dimensions kept", ylab = metric, main = heading)
  axis(1, at = dimensions, labels = x_labels)
  for (i in seq_along(gamma)) {
    values <- vapply(dimensions, function(d) {
      scan_average(data, metric, gamma[i], "mds", d)
    }, numeric(1))
    lines(dimensions, values, type = "b", pch = 16L, lwd = 2, col = colors[i])
    ambient <- scan_average(data, metric, gamma[i], "ambient")
    if (is.finite(ambient)) abline(h = ambient, col = colors[i], lty = 3L)
  }
  legend("bottomleft", legend = c(paste0("gamma ", gamma), "ambient (dotted)"),
         col = c(colors, "gray50"), lty = c(1, 1, 1, 3),
         pch = c(16, 16, 16, NA), bty = "n", cex = 0.73)
}

scan_figures <- function(results, output) {
  source <- unique(results$source)
  for (kind in c("detection", "reliability")) {
    grDevices::pdf(file.path(output, paste0(kind, "_vs_dimension.pdf")), width = 11, height = 5)
    tryCatch({
      groups <- if (source == "real") unique(results$dataset) else {
        unique(paste(results$truth, results$intrinsic_dim, sep = " | "))
      }
      for (group in groups) {
        part <- if (source == "real") results[results$dataset == group, , drop = FALSE] else {
          selected <- paste(results$truth, results$intrinsic_dim, sep = " | ") == group
          results[selected, , drop = FALSE]
        }
        dimensions <- sort(unique(part$embedding_dim[!is.na(part$embedding_dim)]))
        if (!length(dimensions)) next
        labels <- if (source == "real") as.character(dimensions) else {
          paste0(dimensions / part$intrinsic_dim[1L], "x")
        }
        metrics <- if (kind == "reliability") {
          "aumvc_reliable"
        } else if (source == "real") {
          c("roc_auc", "pr_auc")
        } else if (source == "orthogonal") {
          c("roc_structural", "roc_distributional")
        } else {
          c("roc_distributional", "roc_structural")
        }
        old <- par(mfrow = c(1, length(metrics)), mar = c(4, 4, 3, 1), oma = c(0, 0, 1, 0))
        tryCatch({
          for (metric in metrics) {
            scan_panel(part, metric, dimensions,
                       if (metric == "aumvc_reliable") "AUMVC reliable rate" else metric, labels)
          }
          mtext(paste(source, group, sep = " | "), outer = TRUE, line = -0.2)
        }, finally = par(old))
      }
    }, finally = grDevices::dev.off())
  }
}

run_scan_experiment <- function(settings = SCAN_SETTINGS, runs = SCAN_RUNS,
                                parts = SCAN_PARTS, output = SCAN_OUTPUT) {
  runs <- validate_runs(runs)
  if (!exists("mv_degeneracy", mode = "function")) {
    warning("The active Aim 3 reliability rule does not check full-box/score collapse.")
  }
  if (!length(parts) || any(!parts %in% c("synthetic", "orthogonal", "real"))) {
    stop("SCAN_PARTS must contain synthetic, orthogonal, or real.")
  }
  folders <- c(synthetic = "current_generator", orthogonal = "original_orthogonal_shift",
             real = "real_data")
  start <- proc.time()[[3L]]
  all_results <- list()
  for (source in intersect(names(folders), parts)) {
    folder <- file.path(output, folders[[source]])
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    results <- list()
    if (source != "real") {
      grid <- synthetic_grid(settings$synthetic)
      if (source == "orthogonal") {
        grid <- grid[grid$intrinsic_dim == 5L, , drop = FALSE]
      }
      for (run in seq_len(runs)) for (j in seq_len(nrow(grid))) {
        row <- grid[j, , drop = FALSE]
        cat(sprintf("scan %s | run %d/%d | %s q=%d\n",
                    source, run, runs, row$truth, row$intrinsic_dim))
        results[[length(results) + 1L]] <- scan_synthetic_case(row, run, source, settings)
        write_report(do.call(rbind, results), file.path(folder, "scan_results.csv"), settings$output_digits)
      }
    } else {
      registry <- real_registry()
      datasets <- settings$real$datasets
      loaded <- setNames(lapply(datasets, function(d) {
        registry[[d]]$load(settings$real[[d]]$data_dir)
      }), datasets)
      for (run in seq_len(runs)) for (i in seq_along(datasets)) {
        cat(sprintf("scan real | run %d/%d | %s\n", run, runs, datasets[i]))
        results[[length(results) + 1L]] <- scan_real_case(datasets[i], i, run, loaded, registry, settings)
        write_report(do.call(rbind, results), file.path(folder, "scan_results.csv"), settings$output_digits)
      }
    }
    if (!length(results)) next
    data <- do.call(rbind, results)
    rownames(data) <- NULL
    scan_figures(data, folder)
    all_results[[source]] <- data
  }
  cat(sprintf("Finished in %.1f minutes: %s\n", (proc.time()[[3L]] - start) / 60, output))
  invisible(do.call(rbind, all_results))
}

if (sys.nframe() == 0L) run_scan_experiment()
