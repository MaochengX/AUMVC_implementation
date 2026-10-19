# Aim 3 — AUMVC in high dimensions

Aim 3 studys how AUMVC's reliability and performance is changed with ambient dimension inceases, under ambient, Goix feature subsampling, MDS, and Isomap spaces. The detector is fixed as One-Class SVM in this experiement. It tests whether AUMVC ranks different OCSVM bandwidths consistently with ROC-AUC and PR-AUC. 

Configurations including parallel workers are in `experiments/aim3/settings.R`.

## 1. Synthetic experiment

```bash
Rscript experiments/aim3/synthetic_experiment.R
```

Generates nonlinear synthetic data with **two generating functions** (`polynomial_interaction`, `oscillatory_local`), intrinsic dimensions **3, 5, 10**, ambient dimensions **80, 200, 400**, and signal-to-noise ratios **1, 5, 10**. Each dataset has 1,000 observations: 900 normal, 50 distributional outliers, and 50 structural outliers.

For each setting, compares four representations with five OCSVM bandwidths. The default is **50 runs**. It produces:

- `result/synthetic/experiment/syn_results.csv` — individual results for every run, setting, representation, and bandwidth.
- `result/synthetic/experiment/syn_bandwidth_exp.csv` — bandwidth selection and ranking agreement.
- `result/synthetic/analysis/syn_results_summary.csv` — averages of the individual metrics across runs.
- `result/synthetic/analysis/syn_bandwidth_summary.csv` — averages of the bandwidth-comparison metrics across runs.
- `result/synthetic/analysis/syn_results_figures.pdf` — occupancy and runtime versus ambient dimension.
- `result/synthetic/analysis/syn_bandwidth_figures.pdf` — bandwidth-ranking agreement versus ambient dimension.

All paths above are relative to `experiments/aim3/`.

## 2. Real-data experiment

```bash
Rscript experiments/aim3/real_experiment.R
```

Evaluates **ECG200** (96 features, 200 observations), **Fashion-MNIST** (784 features, a 1,000-observation sample), and **Shuttle** (9 features, a 1,000-observation sample). The embedding dimension is **5** by default.

In the current setup, ECG200 class `-1` is anomalous; Fashion-MNIST class `0` is normal and other classes are anomalous (sampled at 10% anomalies); Shuttle class `1` is normal, other retained classes are anomalous, and class `4` is excluded. Dataset files must be available under the paths specified in `settings.R`.

The script prints progress and produces:

- `result/real/experiment/real_results.csv` — individual results for every run, dataset, representation, and bandwidth.
- `result/real/experiment/real_bandwidth_exp.csv` — bandwidth selection and ranking agreement.
- `result/real/analysis/real_results_summary.csv` — averages of individual metrics across runs.
- `result/real/analysis/real_bandwidth_summary.csv` — averages of bandwidth-comparison metrics across runs.
- `result/real/analysis/{dataset}_figures.pdf` — occupancy, runtime, and ranking-agreement plots for each dataset.

**Saving behavior:** A new execution overwrites.

## 3. How the experiment works

**Data and evaluation.** Each case uses five disjoint splits: embedding learning, detector training, reference-box construction, AUMVC evaluation, and held-out label evaluation. Standardization is fitted on the embedding split. Labels are used for the label-evaluation split and ROC/PR calculation, **not for AUMVC or bandwidth selection**.

**Representations** :

| `representation` | Description |
|---|---|
| `ambient` | the original feature space. |
| `goix_subsampling` | use 50 subsets, each containing 5 features. Results are aggregated over the subsets. |
| `mds` | OCSVM and AUMVC in classical MDS embedding coordinates. |
| `isomap` | OCSVM and AUMVC in Isomap embedding coordinates. |

MDS and Isomap are fitted using `manifun` on the embedding split and use out-of-sample extensions for all other splits.

**Bandwidth comparison.** For every representation, the code fits five OCSVMs with `nu = 0.5` and `gamma_scale` in **0.25, 0.5, 1, 2, 4**, using `gamma = gamma_scale / detector_dim`. 

**Monte Carlo settings.** The current defaults are 100,000 uniform reference points per reference box, 1 MC repetition, and mass levels from 0.9 to 0.999. Reliability thresholds are a minimum of **10 reference hits** at every mass level and **relative MC SE ≤ 0.10**.

## 4. results CSVs

### A. `syn_results.csv` and `real_results.csv`:

Each row represents **one run × one setting/dataset × one representation × one bandwidth**.

| Column | Description |
|---|---|
| `run` | Repetition number. |
| `dataset` *(real only)* | ECG200, Fashion-MNIST, or Shuttle. |
| `truth` *(synthetic only)* | Generating function: `polynomial_interaction` or `oscillatory_local`. |
| `intrinsic_dim` *(synthetic only)* | Dimension of the latent variables used to generate the data (3, 5, or 10); not the observed feature count. |
| `snr` *(synthetic only)* | Signal-to-noise ratio; larger snr means less added noise. |
| `sample_size` | Number of observations. |
| `original_dim` | Number of observed input features before unsupported-coordinate filtering. |
| `embedding_dim` | Number of embedding coordinates for MDS/Isomap; blank for ambient/Goix. |
| `detector_dim` | Number of coordinates actually supplied to OCSVM (may differ from `original_dim` after filtering). For Goix feature subsampling, it is 5. |
| `representation` | `ambient`, `goix_subsampling`, `mds`, or `isomap`. |
| `gamma_scale` | Bandwidth-grid multiplier; actual OCSVM `gamma = gamma_scale / detector_dim`. |
| `aumvc_normalized` | AUMVC divided by the reference-box volume. |
| `aumvc_relative_mc_se` | MC standard error of normalized AUMVC divided by normalized AUMVC. Blank if this ratio is undefined. |
| `min_hit_count` | Minimum number of MC reference points inside an estimated level set across all mass levels. Small values mean the relevant sets are poorly resolved; for Goix, the minimum is taken across subsets. |
| `zero_occupancy` | Fraction of mass-grid levels with zero MC reference hits. For Goix, averaged over subsets. |
| `low_occupancy` | Fraction of mass-grid levels with fewer than 10 hits, including zero-hit levels; averaged over subsets for Goix. |
| `box_log_volume` | Log of the reference-box volume in the relevant coordinate space; for Goix, mean log-volume across subsets. |
| `roc_auc` | Held-out ROC-AUC. |
| `pr_auc` | Held-out precision–recall area. |
| `roc_distributional`, `pr_distributional` *(synthetic only)* | ROC/PR calculated using **normal observations versus distributional outliers**, excluding structural outliers. |
| `roc_structural`, `pr_structural` *(synthetic only)* | ROC/PR calculated using **normal observations versus structural outliers**, excluding distributional outliers. |
| `embedding_oos_stress` | Out-of-sample embedding-distance discrepancy, measured between held-out detector-training points and embedding landmarks, with optimal distance rescaling. **Lower means better distance preservation**; blank for ambient/Goix. |
| `runtime_sec` | Computational time in seconds. Shared setup and reference-generation time is divided among the five bandwidths. |
| `aumvc_reliable` | `TRUE` when normalized AUMVC and relative MC SE are finite, `min_hit_count >= 10`, and `aumvc_relative_mc_se <= 0.10`; otherwise `FALSE`. It flags estimation quality, not agreement with labels. |


### B. Bandwidth results: `syn_bandwidth_exp.csv` and `real_bandwidth_exp.csv`

Each row combines the **five OCSVM bandwidth candidates** for one run, setting/dataset, and representation. The identifier columns (`run`, dataset or synthetic settings, dimensions, `representation`) have the same meaning as above.

There are **10 possible pairs** among five bandwidths. A pair has **ROC/PR consensus** only if ROC-AUC and PR-AUC both strictly prefer the same candidate. Other pairs are not eligible for the agreement denominator.

| Column | Description |
|---|---|
| `compared_pairs_count` | Number of the bandwidth pairs with ROC/PR consensus (0–10). This is the denominator for the ordinary agreement and unreliable-pair rates. |
| `unreliable_pairs_rate` | Fraction of consensus pairs where **at least one candidate** fails `aumvc_reliable`; blank if there are no consensus pairs. |
| `selected_gamma_scale` | Bandwidth multiplier with the **smallest finite normalized AUMVC**, using all candidates regardless of the reliability flag. Blank if the best value is tied within tolerance. |
| `agreement_pairs_rank_rate` | Fraction of consensus pairs for which the **the candidate with a lower AUMVC also achieves higher ROC-AUC and PR-AUC**. Unreliable candidates remain included; an AUMVC tie does not count as a match. |
| `selected_gamma_scale_reliable` | Smallest-AUMVC bandwidth **among candidates passing** `aumvc_reliable`. Blank if none passes or the best reliable candidate is tied. |
| `reliable_agreement_pairs_rank_rate` | Agreement rate **only among consensus pairs where both candidates are reliable**. Its denominator may be smaller than `compared_pairs_count`; blank if no pair qualifies. |
| `tie_rate` | Fraction of consensus pairs whose two AUMVC values are equal within the tie tolerance; such pairs count as nonmatches in ordinary agreement. |

**Example:** If 8 bandwidth pairs have ROC/PR consensus and AUMVC agrees on 6 of them, then `compared_pairs_count = 8` and `agreement_pairs_rank_rate = 0.75`. If only 4 of those 8 pairs's candidates are reliable, then `reliable_agreement_pairs_rank_rate` is calculated from **those 4 pairs**, not 8.

### C. summary files: `*_results_summary.csv` and `*_bandwidth_summary.csv`

These files group identical dataset/synthetic settings, representations, and (for candidate results) `gamma_scale`, and summarize them **across experimental runs**:

- **`averaged_<metric>`**: arithmetic average across runs of the corresponding per-run metric, using finite values. For example, `averaged_zero_occupancy` averages the per-run `zero_occupancy`.
- **`aumvc_reliable_rate`**: proportion of runs whose `aumvc_reliable` flag is `TRUE`. For example, `0.80` means 80% of runs passed the estimation thresholds.
- **`selected_gamma_scale_mode`**, **`selected_gamma_scale_reliable_mode`**: the most frequently selected bandwidth multiplier across runs (ignoring missing selections). Blank if there is no unique most-frequent value.


## 5. Pareto analysis

Run after both the synthetic and real experiments finishes, because the script reads both sets of detailed CSVs:

```bash
Rscript experiments/aim3/pareto_analysis.R
```

This compares each representation's total measured computational time across its five bandwidths with five outcomes: low occupancy, zero occupancy, fraction of unreliable candidate pairs, AUMVC–ROC/PR agreement, and agreement restricted to reliable pairs. Runtime and the first three outcomes are **minimized**; the two agreement rates are **maximized**.



### Generated files

All paths are relative to `experiments/aim3/`.

| File | Description |
|---|---|
| `result/real/analysis/pareto/real_pareto.csv` | Pareto status for each real-data case and representation. |
| `result/synthetic/analysis/pareto/syn_pareto.csv` | Pareto status for each synthetic case and representation. |
| `result/real/analysis/pareto/real_pareto_frequency.csv` | Pareto-optimal frequency across real-data cases. |
| `result/synthetic/analysis/pareto/syn_pareto_frequency.csv` | Pareto-optimal frequency across synthetic cases. |
| `result/real/analysis/pareto/real_pareto_figures.pdf` | Runtime–outcome plots and Pareto frontiers for real datasets. |
| `result/synthetic/analysis/pareto/syn_pareto_frequency.pdf` | Pareto-frequency plots for synthetic settings. |

**Key columns:**

- `criterion`: Outcome metric used in the Pareto comparison.
- `runtime_sec`: Execution time.
- `metric_value`: Value of the outcome metric.
- `metric_direction`: Whether lower (`min`) or higher (`max`) is better.
- `pareto_optimal`: Whether the representation is Pareto-optimal (`TRUE`/`FALSE`).
- `evaluated_cases`: Number of cases with valid comparisons.
- `pareto_cases`: Number of cases where the representation is Pareto-optimal.
- `pareto_rate`: `pareto_cases / evaluated_cases`.



