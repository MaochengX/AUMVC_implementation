source("experiments/aim3/settings.R")
source("experiments/aim3/utils.R")
source("experiments/aim3/representation_comparison.R")
source("experiments/aim3/analysis.R")

real_split_counts <- function(n) {
  fifth <- n %/% 5L
  c(
    embedding = fifth,
    detector_train = fifth,
    reference = fifth,
    evaluation = fifth,
    label_eval = n - 4L * fifth
  )
}

read_ecg_ts <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  data_line <- which(tolower(lines) == "@data")
  if (length(data_line) != 1L) stop("Invalid ECG200 .ts file.", call. = FALSE)
  rows <- lines[(data_line + 1L):length(lines)]
  rows <- rows[nzchar(rows) & !startsWith(rows, "#")]
  parts <- strsplit(rows, ":", fixed = TRUE)
  if (any(lengths(parts) != 2L)) stop("Invalid ECG200 row.", call. = FALSE)
  x <- do.call(rbind, lapply(parts, function(row) {
    as.numeric(strsplit(row[1L], ",", fixed = TRUE)[[1L]])
  }))
  labels <- vapply(parts, function(row) as.numeric(row[2L]), numeric(1))
  list(x = x, labels = labels)
}

read_ecg_txt <- function(path) {
  data <- as.matrix(utils::read.table(path, header = FALSE, check.names = FALSE))
  list(x = data[, -1L, drop = FALSE], labels = data[, 1L])
}

load_ecg200 <- function(data_dir) {
  train_ts <- file.path(data_dir, "ECG200_TRAIN.ts")
  test_ts <- file.path(data_dir, "ECG200_TEST.ts")
  train_txt <- file.path(data_dir, "ECG200_TRAIN.txt")
  test_txt <- file.path(data_dir, "ECG200_TEST.txt")
  if (file.exists(train_ts) && file.exists(test_ts)) {
    train <- read_ecg_ts(train_ts)
    test <- read_ecg_ts(test_ts)
  } else if (file.exists(train_txt) && file.exists(test_txt)) {
    train <- read_ecg_txt(train_txt)
    test <- read_ecg_txt(test_txt)
  } else {
    stop("ECG200 files not found in ", data_dir, call. = FALSE)
  }
  x <- rbind(train$x, test$x)
  labels <- c(train$labels, test$labels)
  if (nrow(x) != 200L || ncol(x) != 96L || any(!is.finite(x))) {
    stop("Unexpected ECG200 data.", call. = FALSE)
  }
  list(x = x, labels = labels)
}

idx_connection <- function(path) {
  if (grepl("\\.gz$", path)) gzfile(path, "rb") else file(path, "rb")
}

read_idx_labels <- function(path) {
  connection <- idx_connection(path)
  on.exit(close(connection))
  magic <- readBin(connection, integer(), 1L, 4L, endian = "big")
  n <- readBin(connection, integer(), 1L, 4L, endian = "big")
  if (magic != 2049L || n < 1L) stop("Invalid Fashion-MNIST labels.", call. = FALSE)
  readBin(connection, integer(), n, 1L, signed = FALSE)
}

read_idx_images <- function(path) {
  connection <- idx_connection(path)
  on.exit(close(connection))
  header <- readBin(connection, integer(), 4L, 4L, endian = "big")
  if (header[1L] != 2051L) stop("Invalid Fashion-MNIST images.", call. = FALSE)
  pixels <- readBin(
    connection,
    integer(),
    header[2L] * header[3L] * header[4L],
    1L,
    signed = FALSE
  )
  matrix(pixels, nrow = header[2L], byrow = TRUE) / 255
}

load_fashion_mnist <- function(data_dir) {
  images <- file.path(data_dir, c("t10k-images-idx3-ubyte", "t10k-images-idx3-ubyte.gz"))
  labels <- file.path(data_dir, c("t10k-labels-idx1-ubyte", "t10k-labels-idx1-ubyte.gz"))
  image_file <- images[file.exists(images)][1L]
  label_file <- labels[file.exists(labels)][1L]
  if (is.na(image_file) || is.na(label_file)) {
    stop("Fashion-MNIST files not found in ", data_dir, call. = FALSE)
  }
  x <- read_idx_images(image_file)
  y <- read_idx_labels(label_file)
  if (nrow(x) != length(y) || ncol(x) != 784L || any(!is.finite(x))) {
    stop("Unexpected Fashion-MNIST data.", call. = FALSE)
  }
  list(x = x, labels = y)
}

load_shuttle <- function(data_dir) {
  test <- file.path(data_dir, "shuttle.tst")
  train <- file.path(data_dir, "shuttle.trn")
  if (!file.exists(test)) stop("shuttle.tst not found in ", data_dir, call. = FALSE)
  files <- if (file.exists(train)) c(train, test) else test
  data <- do.call(rbind, lapply(files, function(path) {
    as.matrix(utils::read.table(path, header = FALSE, check.names = FALSE))
  }))
  if (ncol(data) != 10L || any(!is.finite(data))) stop("Unexpected Shuttle data.", call. = FALSE)
  list(x = data[, 1:9, drop = FALSE], labels = data[, 10L])
}

sample_ecg200 <- function(data, seed, settings) {
  set.seed(seed)
  rows <- sample.int(nrow(data$x))
  list(x = data$x[rows, , drop = FALSE], labels = as.integer(data$labels[rows] == -1))
}

sample_fashion_mnist <- function(data, seed, settings) {
  size <- as.integer(settings$sample_size)
  normal <- as.integer(settings$normal_class)
  fraction <- settings$anomaly_fraction
  if (size < 5L || fraction <= 0 || fraction >= 1) stop("Invalid Fashion-MNIST settings.", call. = FALSE)
  normal_rows <- which(data$labels == normal)
  anomaly_rows <- which(data$labels != normal)
  n_anomaly <- as.integer(round(size * fraction))
  n_normal <- size - n_anomaly
  if (n_normal > length(normal_rows) || n_anomaly > length(anomaly_rows)) {
    stop("Fashion-MNIST sample is too large.", call. = FALSE)
  }
  set.seed(seed)
  rows <- c(
    sample(normal_rows, n_normal, replace = FALSE),
    sample(anomaly_rows, n_anomaly, replace = FALSE)
  )
  rows <- sample(rows)
  list(x = data$x[rows, , drop = FALSE], labels = as.integer(data$labels[rows] != normal))
}

sample_shuttle <- function(data, seed, settings) {
  keep <- data$labels != as.integer(settings$excluded_class)
  x <- data$x[keep, , drop = FALSE]
  labels <- data$labels[keep]
  size <- as.integer(settings$sample_size)
  if (size > nrow(x)) stop("Shuttle sample is too large.", call. = FALSE)
  set.seed(seed)
  rows <- sample.int(nrow(x), size)
  list(x = x[rows, , drop = FALSE], labels = as.integer(labels[rows] != 1L))
}

real_registry <- function() {
  list(
    ecg200 = list(load = load_ecg200, sample = sample_ecg200),
    fashion_mnist = list(load = load_fashion_mnist, sample = sample_fashion_mnist),
    shuttle = list(load = load_shuttle, sample = sample_shuttle)
  )
}

real_result_columns <- c(
  "dataset", "run", "sample_size", "original_dim", "embedding_dim", "detector_dim",
  "representation", "gamma_scale", "aumvc_normalized", "aumvc_relative_mc_se",
  "min_hit_count", "zero_occupancy", "low_occupancy", "box_log_volume",
  "roc_auc", "pr_auc", "embedding_oos_stress", "runtime_sec", "aumvc_reliable"
)

real_bandwidth_columns <- c(
  "dataset", "run", "sample_size", "original_dim", "embedding_dim", "detector_dim",
  "representation", "compared_pairs_count", "unreliable_pairs_rate",
  "selected_gamma_scale", "agreement_pairs_rank_rate", "selected_gamma_scale_reliable",
  "reliable_agreement_pairs_rank_rate", "tie_rate"
)

run_real_experiment <- function(settings = AIM3_SETTINGS) {
  n_runs <- validate_runs(settings$n_runs)
  registry <- real_registry()
  datasets <- as.character(settings$real$datasets)
  if (!length(datasets) || any(!datasets %in% names(registry))) {
    stop("Invalid real datasets.", call. = FALSE)
  }
  loaded <- setNames(lapply(datasets, function(dataset) {
    registry[[dataset]]$load(settings$real[[dataset]]$data_dir)
  }), datasets)
  total <- n_runs * length(datasets)
  results <- vector("list", total)
  bandwidth <- vector("list", total)
  dirs <- result_dirs("real")
  done <- 0L
  start <- proc.time()[[3L]]
  cat("\naim3 real\n")
  for (run in seq_len(n_runs)) {
    for (dataset in datasets) {
      seed <- settings$seed + run * 100000L + match(dataset, datasets) * 10000L
      sampled <- registry[[dataset]]$sample(loaded[[dataset]], seed, settings$real[[dataset]])
      case_settings <- settings
      case_settings$split_counts <- real_split_counts(nrow(sampled$x))
      case <- list(
        x = sampled$x,
        labels = sampled$labels,
        intrinsic_dim = as.integer(settings$real$embedding_dim)
      )
      evaluated <- evaluate_case(case, seed + 50000L, case_settings)
      metadata <- data.frame(
        dataset = dataset,
        run = run,
        sample_size = nrow(sampled$x),
        original_dim = ncol(sampled$x),
        stringsAsFactors = FALSE
      )
      done <- done + 1L
      results[[done]] <- attach_metadata(metadata, evaluated$results)
      bandwidth[[done]] <- attach_metadata(metadata, evaluated$bandwidth)
      elapsed <- proc.time()[[3L]] - start
      cat(sprintf("\r%d/%d | elapsed %.0fs | eta %.0fs", done, total, elapsed, elapsed * (total - done) / done))
      flush.console()
    }
    current_results <- do.call(rbind, results[seq_len(done)])
    current_bandwidth <- do.call(rbind, bandwidth[seq_len(done)])
    current_results <- current_results[, real_result_columns, drop = FALSE]
    current_bandwidth <- current_bandwidth[, real_bandwidth_columns, drop = FALSE]
    write_report(current_results, file.path(dirs$experiment, "real_results.csv"), settings$output_digits)
    write_report(current_bandwidth, file.path(dirs$experiment, "real_bandwidth_exp.csv"), settings$output_digits)
  }
  results <- do.call(rbind, results)
  bandwidth <- do.call(rbind, bandwidth)
  results <- results[, real_result_columns, drop = FALSE]
  bandwidth <- bandwidth[, real_bandwidth_columns, drop = FALSE]
  write_report(results, file.path(dirs$experiment, "real_results.csv"), settings$output_digits)
  write_report(bandwidth, file.path(dirs$experiment, "real_bandwidth_exp.csv"), settings$output_digits)
  analyze_results(results, bandwidth, "real", dirs$analysis, settings$output_digits)
  cat("\ncomplete: ", dirname(dirs$experiment), "\n", sep = "")
  invisible(list(results = results, bandwidth = bandwidth))
}

if (sys.nframe() == 0L) run_real_experiment()
