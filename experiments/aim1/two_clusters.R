source("aumvc/input_validation.R")
source("aumvc/level_set.R")
source("aumvc/aumvc.R")

set.seed(2027)

n_per_cluster <- 100

x_train <- rbind(
  matrix(
    rnorm(
      n_per_cluster * 2,
      mean = -2,
      sd = 0.7
    ),
    ncol = 2
  ),
  matrix(
    rnorm(
      n_per_cluster * 2,
      mean = 2,
      sd = 0.7
    ),
    ncol = 2
  )
)

x_eval <- rbind(
  matrix(
    rnorm(
      99 * 2,
      mean = -2,
      sd = 0.7
    ),
    ncol = 2
  ),
  matrix(
    rnorm(
      99 * 2,
      mean = 2,
      sd = 0.7
    ),
    ncol = 2
  ),
  matrix(
    c(
      0, 7,
      7, 0
    ),
    ncol = 2,
    byrow = TRUE
  )
)

x_train <- validate_reference_inputs(
  x_train,
  n_reference = 20000L,
  n_mc_repetitions = 5L,
  seed = 2027L
)

centers <- rbind(
  c(-2, -2),
  c(2, 2)
)

squared_distance <- function(x, center) {
  rowSums(
    sweep(
      x,
      2,
      center,
      FUN = "-"
    )^2
  )
}

correct_score <- function(newdata) {
  pmax(
    -squared_distance(
      newdata,
      centers[1L, ]
    ),
    -squared_distance(
      newdata,
      centers[2L, ]
    )
  )
}

reversed_score <- function(newdata) {
  -correct_score(newdata)
}

reference <- make_reference(
  x_train,
  n_reference = 20000L,
  n_mc_repetitions = 5L,
  seed = 2027L
)

x_eval <- validate_aumvc_inputs(
  x_eval,
  reference,
  correct_score,
  "normality",
  seq(0.9, 0.999, length.out = 100L)
)

validate_scores(
  correct_score(x_eval),
  n_expected = nrow(x_eval),
  name = "correct scores"
)

validate_scores(
  reversed_score(x_eval),
  n_expected = nrow(x_eval),
  name = "reversed scores"
)

correct <- aumvc(
  x_eval,
  reference,
  correct_score,
  score_direction = "normality"
)

reversed <- aumvc(
  x_eval,
  reference,
  reversed_score,
  score_direction = "normality"
)

comparison <- data.frame(
  scoring_rule = c(
    "correct",
    "reversed"
  ),
  aumvc = c(
    correct$aumvc,
    reversed$aumvc
  ),
  aumvc_mc_se = c(
    correct$aumvc_mc_se,
    reversed$aumvc_mc_se
  )
)

print(comparison)

stopifnot(
  correct$aumvc < reversed$aumvc
)
