source("experiments/aim4/tuning_common.R")

OCSVM_GAMMA_SCALES <- c(1, 2, 4, 6, 8, 10, 12, 14, 16, 64)
OCSVM_DEFAULT_SCALE <- 1
OCSVM_NU <- 0.5

score_ocsvm_grid <- function(prepared, scales, settings) {
  train <- prepared$detector_train
  models <- lapply(scales, function(scale) {
    fit_ocsvm(train, nu = OCSVM_NU, gamma = scale / ncol(train))
  })
  score_all <- function(x) do.call(cbind, lapply(models, score_ocsvm, newdata = x))
  evaluation <- score_all(prepared$evaluation)
  validation <- score_all(prepared$validation)
  test <- score_all(prepared$test)
  score_funs <- lapply(models, function(model) {
    force(model)
    function(x) score_ocsvm(model, x)
  })
  mc <- score_reference_candidates(prepared$reference_mc, score_funs)$scores
  results <- score_aim4b_candidates(
    prepared, scales, evaluation, validation, test, mc, settings
  )
  results$gamma <- scales / ncol(train)
  results$nu <- OCSVM_NU
  results
}

run_ocsvm_tune <- function(settings = AIM4B_SETTINGS) {
  run_aim4b("ocsvm", settings, score_ocsvm_grid,
            OCSVM_GAMMA_SCALES, OCSVM_DEFAULT_SCALE, max(OCSVM_GAMMA_SCALES))
}

if (sys.nframe() == 0L) run_ocsvm_tune()
