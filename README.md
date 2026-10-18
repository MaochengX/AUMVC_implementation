# AUMVC Implementation

## Overview

This repository implements and validates Mass-Volume (MV) curves and AUMVC in R for label-free evaluation of outlier detectors. Its main focus is how high dimensionality affects Monte Carlo volume estimation and whether computing AUMVC in manifold-embedding spaces improves its reliability.

## Experiments

- **Aim 1 — AUMVC Implementation in R:** Implement AUMVC with sanity checks and validate it against `christian-oleary/emmv`.
- **Aim 2 — Validation against Labels:** Compare AUMVC rankings of OCSVM, LOF, and Isolation Forest with ROC-AUC and PR-AUC.
- **Aim 3 — Dimensionality Robustness:** Study AUMVC's reliability and performance when ambient dimension increases by experimenting in ambient, Goix feature subsampling, MDS, and Isomap spaces.

The AUMVC implementation is in `aumvc/`. Experiments and explanations are in `experiments/aim1/`, `experiments/aim2/`, `experiments/aim3/`, and `experiments/aim4/` (see each aim's `results.md`).

## Dependencies

Requires **R**. External R packages: `e1071`, `dbscan`, `isotree`, `vegan`, `manifun`, `future`, `future.apply`, and `pak`.

Install `manifun` from [HerrMo/manifun](https://github.com/HerrMo/manifun).

From the repository root:

```bash
Rscript -e 'install.packages(c("e1071", "dbscan", "isotree", "vegan", "future", "future.apply", "pak"), repos="https://cloud.r-project.org")'
Rscript experiments/install_manifun.R
```

## Run

Run from the repository root:

```bash
# Aim 1
Rscript experiments/aim1/reference_validation.R
Rscript experiments/aim1/one_cluster.R
Rscript experiments/aim1/two_clusters.R

# Aim 2
Rscript experiments/aim2/synthetic.R
Rscript experiments/aim2/real_experiment.R

# Aim 3
Rscript experiments/test_optimization.R
Rscript experiments/aim3/synthetic_experiment.R
Rscript experiments/aim3/real_experiment.R
Rscript experiments/aim3/pareto_analysis.R
```

For details, see each aim's `results.md`.
