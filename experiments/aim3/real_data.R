aim3_real_split_counts <- function(n) {
  fifth <- n %/% 5L
  c(embedding = fifth, detector_train = fifth, reference = fifth,
    evaluation = fifth, label_eval = n - 4L * fifth)
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

  list(x = data$x[rows, , drop = FALSE],
       labels = as.integer(data$labels[rows] == -1))
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

  list(x = data$x[rows, , drop = FALSE],
       labels = as.integer(data$labels[rows] != normal_class))
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

  list(x = x[rows, , drop = FALSE],
       labels = as.integer(data$labels[keep][rows] != 1L))
}

aim3_real_data_md5 <- function(files) {
  hashes <- tools::md5sum(files)
  names(hashes) <- basename(files)
  hashes
}

