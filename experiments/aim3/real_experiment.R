source("experiments/settings.R")
source("experiments/aim3/representation_comparison.R")

settings <- AIM3_SETTINGS
real_settings <- settings$real

results_root <- "experiments/aim3/results/real"
active_file <- file.path(results_root, "active_experiment.txt")
latest_file <- file.path(results_root, "latest_experiment.txt")

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)

aim3_real_experiment_id <- function() {
  base_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  experiment_id <- base_id
  suffix <- 1L

  dataset_roots <- file.path(
    results_root,
    c("ecg200", "fashion_mnist", "shuttle")
  )

  while (any(dir.exists(file.path(dataset_roots, experiment_id)))) {
    suffix <- suffix + 1L
    experiment_id <- paste0(base_id, "_", sprintf("%02d", suffix))
  }

  experiment_id
}

aim3_real_split_counts <- function(n) {
  quarter <- n %/% 4L

  c(
    embedding = quarter,
    detector_train = quarter,
    reference = quarter,
    evaluation = n - 3L * quarter
  )
}

aim3_read_ecg_ts <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  data_line <- which(tolower(lines) == "@data")

  if (length(data_line) != 1L) {
    stop("Invalid ECG200 .ts file.", call. = FALSE)
  }

  rows <- lines[(data_line + 1L):length(lines)]
  rows <- rows[nzchar(rows) & !startsWith(rows, "#")]
  parts <- strsplit(rows, ":", fixed = TRUE)

  if (any(lengths(parts) != 2L)) {
    stop("Unexpected ECG200 .ts row format.", call. = FALSE)
  }

  x <- do.call(
    rbind,
    lapply(
      parts,
      function(row) {
        as.numeric(strsplit(row[1L], ",", fixed = TRUE)[[1L]])
      }
    )
  )

  labels <- vapply(
    parts,
    function(row) as.numeric(row[2L]),
    numeric(1)
  )

  list(x = x, labels = labels)
}

aim3_read_ecg_txt <- function(path) {
  data <- as.matrix(
    utils::read.table(
      path,
      header = FALSE,
      check.names = FALSE
    )
  )

  list(
    x = data[, -1L, drop = FALSE],
    labels = data[, 1L]
  )
}

aim3_load_ecg200 <- function(data_dir) {
  train_ts <- file.path(data_dir, "ECG200_TRAIN.ts")
  test_ts <- file.path(data_dir, "ECG200_TEST.ts")
  train_txt <- file.path(data_dir, "ECG200_TRAIN.txt")
  test_txt <- file.path(data_dir, "ECG200_TEST.txt")

  if (file.exists(train_ts) && file.exists(test_ts)) {
    train <- aim3_read_ecg_ts(train_ts)
    test <- aim3_read_ecg_ts(test_ts)
    files <- c(train_ts, test_ts)
  } else if (file.exists(train_txt) && file.exists(test_txt)) {
    train <- aim3_read_ecg_txt(train_txt)
    test <- aim3_read_ecg_txt(test_txt)
    files <- c(train_txt, test_txt)
  } else {
    stop(
      paste(
        "ECG200 files not found.",
        "Expected ECG200_TRAIN.ts and ECG200_TEST.ts",
        "or the corresponding .txt files in",
        data_dir
      ),
      call. = FALSE
    )
  }

  x <- rbind(train$x, test$x)
  labels <- c(train$labels, test$labels)

  if (nrow(x) != 200L || ncol(x) != 96L || any(!is.finite(x))) {
    stop("Unexpected ECG200 data dimensions or values.", call. = FALSE)
  }

  list(
    x = x,
    labels = labels,
    files = files
  )
}

aim3_idx_connection <- function(path) {
  if (grepl("\\.gz$", path)) {
    gzfile(path, "rb")
  } else {
    file(path, "rb")
  }
}

aim3_read_idx_labels <- function(path) {
  connection <- aim3_idx_connection(path)
  on.exit(close(connection), add = TRUE)

  magic <- readBin(
    connection,
    integer(),
    n = 1L,
    size = 4L,
    endian = "big"
  )

  n <- readBin(
    connection,
    integer(),
    n = 1L,
    size = 4L,
    endian = "big"
  )

  if (magic != 2049L || n < 1L) {
    stop("Invalid Fashion-MNIST label file.", call. = FALSE)
  }

  readBin(
    connection,
    integer(),
    n = n,
    size = 1L,
    signed = FALSE
  )
}

aim3_read_idx_images <- function(path) {
  connection <- aim3_idx_connection(path)
  on.exit(close(connection), add = TRUE)

  header <- readBin(
    connection,
    integer(),
    n = 4L,
    size = 4L,
    endian = "big"
  )

  if (header[1L] != 2051L) {
    stop("Invalid Fashion-MNIST image file.", call. = FALSE)
  }

  n <- header[2L]
  rows <- header[3L]
  columns <- header[4L]

  pixels <- readBin(
    connection,
    integer(),
    n = n * rows * columns,
    size = 1L,
    signed = FALSE
  )

  matrix(
    pixels,
    nrow = n,
    ncol = rows * columns,
    byrow = TRUE
  ) / 255
}

aim3_load_fashion_mnist <- function(data_dir) {
  image_candidates <- file.path(
    data_dir,
    c(
      "t10k-images-idx3-ubyte",
      "t10k-images-idx3-ubyte.gz"
    )
  )

  label_candidates <- file.path(
    data_dir,
    c(
      "t10k-labels-idx1-ubyte",
      "t10k-labels-idx1-ubyte.gz"
    )
  )

  image_file <- image_candidates[file.exists(image_candidates)][1L]
  label_file <- label_candidates[file.exists(label_candidates)][1L]

  if (is.na(image_file) || is.na(label_file)) {
    stop(
      paste(
        "Fashion-MNIST test files not found in",
        data_dir
      ),
      call. = FALSE
    )
  }

  x <- aim3_read_idx_images(image_file)
  labels <- aim3_read_idx_labels(label_file)

  if (
    nrow(x) != length(labels) ||
    ncol(x) != 784L ||
    any(!is.finite(x))
  ) {
    stop("Unexpected Fashion-MNIST data.", call. = FALSE)
  }

  list(
    x = x,
    labels = labels,
    files = c(image_file, label_file)
  )
}

aim3_load_shuttle <- function(data_dir) {
  test_file <- file.path(data_dir, "shuttle.tst")
  train_file <- file.path(data_dir, "shuttle.trn")

  if (!file.exists(test_file)) {
    stop(
      paste(
        "shuttle.tst not found in",
        data_dir
      ),
      call. = FALSE
    )
  }

  files <- test_file

  if (file.exists(train_file)) {
    files <- c(train_file, test_file)
  }

  data <- do.call(
    rbind,
    lapply(
      files,
      function(path) {
        as.matrix(
          utils::read.table(
            path,
            header = FALSE,
            check.names = FALSE
          )
        )
      }
    )
  )

  if (ncol(data) != 10L || any(!is.finite(data))) {
    stop("Unexpected Shuttle data.", call. = FALSE)
  }

  list(
    x = data[, 1:9, drop = FALSE],
    labels = data[, 10L],
    files = files
  )
}

aim3_sample_ecg200 <- function(data, seed, dataset_settings) {
  set.seed(seed)
  rows <- sample.int(nrow(data$x))

  data$x[rows, , drop = FALSE]
}

aim3_sample_fashion_mnist <- function(data, seed, dataset_settings) {
  sample_size <- as.integer(dataset_settings$sample_size)
  normal_class <- as.integer(dataset_settings$normal_class)
  anomaly_fraction <- dataset_settings$anomaly_fraction

  if (
    sample_size < 4L ||
    anomaly_fraction <= 0 ||
    anomaly_fraction >= 1
  ) {
    stop("Invalid Fashion-MNIST settings.", call. = FALSE)
  }

  normal_rows <- which(data$labels == normal_class)
  anomaly_rows <- which(data$labels != normal_class)

  n_anomaly <- as.integer(round(sample_size * anomaly_fraction))
  n_normal <- sample_size - n_anomaly

  if (
    n_normal > length(normal_rows) ||
    n_anomaly > length(anomaly_rows)
  ) {
    stop("Fashion-MNIST sample request exceeds available data.", call. = FALSE)
  }

  set.seed(seed)

  rows <- c(
    sample(normal_rows, n_normal, replace = FALSE),
    sample(anomaly_rows, n_anomaly, replace = FALSE)
  )

  rows <- sample(rows)

  data$x[rows, , drop = FALSE]
}

aim3_sample_shuttle <- function(data, seed, dataset_settings) {
  sample_size <- as.integer(dataset_settings$sample_size)
  excluded_class <- as.integer(dataset_settings$excluded_class)
  keep <- data$labels != excluded_class
  x <- data$x[keep, , drop = FALSE]

  if (sample_size > nrow(x)) {
    stop("Shuttle sample request exceeds available data.", call. = FALSE)
  }

  set.seed(seed)
  rows <- sample.int(nrow(x), sample_size)

  x[rows, , drop = FALSE]
}

aim3_real_data_md5 <- function(files) {
  hashes <- tools::md5sum(files)
  names(hashes) <- basename(files)
  hashes
}

aim3_real_format_time <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) return("--:--:--")

  seconds <- as.integer(round(seconds))
  hours <- seconds %/% 3600L
  minutes <- (seconds %% 3600L) %/% 60L
  seconds <- seconds %% 60L

  sprintf("%02d:%02d:%02d", hours, minutes, seconds)
}

aim3_real_show_progress <- function(
    done,
    total,
    elapsed,
    mean_run_time,
    dataset,
    run
) {
  width <- 30L
  fraction <- done / total
  filled <- floor(width * fraction)

  eta <- if (is.finite(mean_run_time)) {
    (total - done) * mean_run_time
  } else {
    NA_real_
  }

  bar <- paste0(
    "[",
    paste(rep("=", filled), collapse = ""),
    paste(rep(" ", width - filled), collapse = ""),
    "]"
  )

  cat(sprintf(
    "\r%s %d/%d %5.1f%% | %s run %d | elapsed %s | ETA %s",
    bar,
    done,
    total,
    100 * fraction,
    dataset,
    run,
    aim3_real_format_time(elapsed),
    aim3_real_format_time(eta)
  ))

  flush.console()
}

allowed_datasets <- c(
  "ecg200",
  "fashion_mnist",
  "shuttle"
)

if (
  length(real_settings$datasets) < 1L ||
  any(!real_settings$datasets %in% allowed_datasets)
) {
  stop("Invalid Aim 3 real dataset selection.", call. = FALSE)
}

datasets <- list()

if ("ecg200" %in% real_settings$datasets) {
  datasets$ecg200 <- list(
    data = aim3_load_ecg200(real_settings$ecg200$data_dir),
    settings = real_settings$ecg200,
    sample = aim3_sample_ecg200,
    offset = 1000000L
  )
}

if ("fashion_mnist" %in% real_settings$datasets) {
  datasets$fashion_mnist <- list(
    data = aim3_load_fashion_mnist(
      real_settings$fashion_mnist$data_dir
    ),
    settings = real_settings$fashion_mnist,
    sample = aim3_sample_fashion_mnist,
    offset = 2000000L
  )
}

if ("shuttle" %in% real_settings$datasets) {
  datasets$shuttle <- list(
    data = aim3_load_shuttle(real_settings$shuttle$data_dir),
    settings = real_settings$shuttle,
    sample = aim3_sample_shuttle,
    offset = 3000000L
  )
}

if (file.exists(active_file)) {
  experiment_id <- trimws(
    readLines(active_file, warn = FALSE)[1L]
  )

  if (!nzchar(experiment_id)) {
    stop("Invalid active real-data experiment.", call. = FALSE)
  }

  cat("Resuming real-data experiment: ", experiment_id, "\n", sep = "")
} else {
  experiment_id <- aim3_real_experiment_id()
  writeLines(experiment_id, active_file)

  cat("Starting real-data experiment: ", experiment_id, "\n", sep = "")
}

settings_file <- file.path(
  results_root,
  paste0("settings_", experiment_id, ".rds")
)

if (file.exists(settings_file)) {
  if (!identical(readRDS(settings_file), settings)) {
    stop(
      paste(
        "Aim 3 settings changed during an unfinished real-data experiment.",
        "Restore the previous settings or delete",
        active_file,
        "to start a new experiment."
      ),
      call. = FALSE
    )
  }
} else {
  saveRDS(settings, settings_file)
}

for (dataset_name in names(datasets)) {
  dataset_info <- datasets[[dataset_name]]

  experiment_dir <- file.path(
    results_root,
    dataset_name,
    experiment_id
  )

  dir.create(
    experiment_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  md5_file <- file.path(experiment_dir, "data_md5.rds")
  current_md5 <- aim3_real_data_md5(dataset_info$data$files)

  if (file.exists(md5_file)) {
    if (!identical(readRDS(md5_file), current_md5)) {
      stop(
        dataset_name,
        " data files changed during the experiment.",
        call. = FALSE
      )
    }
  } else {
    saveRDS(current_md5, md5_file)
  }
}

total_runs <- length(datasets) * settings$n_runs
done_runs <- 0L
run_times <- numeric(0)
start_time <- proc.time()[[3L]]

for (dataset_name in names(datasets)) {
  dataset_info <- datasets[[dataset_name]]

  experiment_dir <- file.path(
    results_root,
    dataset_name,
    experiment_id
  )

  checkpoint_dir <- file.path(experiment_dir, "checkpoints")

  dir.create(
    checkpoint_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  checkpoint_files <- file.path(
    checkpoint_dir,
    sprintf("run_%02d.rds", seq_len(settings$n_runs))
  )

  done_runs <- done_runs + sum(file.exists(checkpoint_files))

  for (run in seq_len(settings$n_runs)) {
    checkpoint_file <- checkpoint_files[run]

    if (file.exists(checkpoint_file)) next

    run_start <- proc.time()[[3L]]

    run_seed <- settings$seed +
      dataset_info$offset +
      (run - 1L) * 100000L

    x <- dataset_info$sample(
      dataset_info$data,
      run_seed,
      dataset_info$settings
    )

    analysis_settings <- settings
    analysis_settings$split_counts <- aim3_real_split_counts(nrow(x))

    case <- list(
      x = x,
      intrinsic_dim = as.integer(real_settings$mds_dim)
    )

    comparison <- aim3_compare_representations(
      case,
      run_seed + 50000L,
      analysis_settings
    )

    result <- data.frame(
      run = run,
      sample_size = nrow(x),
      ambient_dim = ncol(x),
      mds_dim = as.integer(real_settings$mds_dim),
      comparison
    )

    temporary_file <- paste0(checkpoint_file, ".tmp")
    saveRDS(result, temporary_file)

    if (!file.rename(temporary_file, checkpoint_file)) {
      stop("Could not save real-data checkpoint.", call. = FALSE)
    }

    run_times <- c(
      run_times,
      proc.time()[[3L]] - run_start
    )

    done_runs <- done_runs + 1L

    aim3_real_show_progress(
      done_runs,
      total_runs,
      proc.time()[[3L]] - start_time,
      mean(run_times),
      dataset_name,
      run
    )
  }

  results <- do.call(
    rbind,
    lapply(checkpoint_files, readRDS)
  )

  saveRDS(
    list(
      settings = settings,
      experiment_id = experiment_id,
      dataset = dataset_name,
      results = results
    ),
    file.path(experiment_dir, "real.rds")
  )
}

cat("\n")

writeLines(experiment_id, latest_file)
unlink(active_file)

cat(
  "Completed real-data experiment: ",
  experiment_id,
  "\n",
  sep = ""
)
