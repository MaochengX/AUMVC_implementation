# Aim 1 — AUMVC validation

Aim 1 checks the R implementation of AUMVC using simple low-dimensional examples and an independent Python reference. Run the following commands from the repository root.

## 1. Reference checks

```bash
Rscript experiments/aim1/reference_validation.R
```

Checks score-direction consistency and confirms that OCSVM, LOF, and Isolation Forest return valid scores. Prints `AUMVC and detector reference validation passed.` when all checks succeed.

## 2. One-cluster example

```bash
Rscript experiments/aim1/one_cluster.R
```

Evaluates AUMVC on a 2D single-cluster dataset with two generated anomalies. Prints a one-row table with `AUMVC` and its Monte Carlo standard error (`MC_SE`).

## 3. Two-cluster example

```bash
Rscript experiments/aim1/two_clusters.R
```

Evaluates AUMVC on a 2D two-cluster dataset. Prints `AUMVC`, `MC_SE`, and `Score-direction test passed` if normality and anomaly score conventions agree.

## 4. Python reference validation

```bash
Rscript experiments/aim1/python_validation.R
```

Compares the R implementation with **`christian-oleary/emmv` (`emmv==0.0.4`)**. The comparison uses three 2D/3D datasets and prints the R and Python AUMVC values, their differences, and a `passed` flag for each case.

**Python is required**, but no separate Python environment or manual Python command is needed. The R script finds `python3` (or `python`) and runs the reference comparison itself; it checks for `emmv==0.0.4` and installs it using that Python installation's `pip` if necessary.

**Outputs:** All four scripts print to the terminal.
