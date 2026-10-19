# Aim 2 — Validation against labels

Aim 2 checks whether AUMVC agrees with ROC-AUC and PR-AUC when comparing One-Class SVM, Local Outlier Factor, and Isolation Forest. Run the commands from the repository root.

## 1. Synthetic experiment

```bash
Rscript experiments/aim2/synthetic.R
```

Generates one 2D dataset with a central normal cluster and anomalies in a surrounding ring. After evaluation, prints results for each dataset and an overall agreement summary, and saves:

`experiments/aim2/result/synthetic/expN/synthetic_results.csv`

## 2. Real-data experiment

```bash
Rscript experiments/aim2/real_experiment.R
```

Evaluates Adult, HTTP, Pima, SMTP, and Wilt (as selected in `experiments/aim2/settings.R`). Prints results for each dataset and an overall agreement summary, and saves:

`experiments/aim2/result/real/expN/real_results.csv`

Each execution creates a new `expN` folder, preserving previous results. The default is **50 runs per dataset**.

## 3. Evaluation and Comparison 

Each run uses separate data splits for detector training, construction of the Monte Carlo reference box, AUMVC evaluation, and label-based evaluation. Features are standardized. Labels are used in ROC-AUC/PR-AUC; **labels are not used to calculate AUMVC**.

By default, each run uses 100,000 Monte Carlo reference points, one MC repetition, mass levels from 0.9 to 0.999, and a maximum anomaly fraction of 10% (see `settings.R`).

For each detector, lower AUMVC is better, while higher ROC-AUC and PR-AUC are better. The agreement calculation compares **pairs of detectors** only when ROC-AUC and PR-AUC agree with each other.

## 4. The CSV results

The CSV contains five rows each dataset: **detector rows** (three per dataset) followed by **comparison rows** (two per dataset). Empty cells are intentional: the two row types report different information.

### Detector rows

| Column | Meaning |
|---|---|
| `dataset` | Dataset name and total number of runs, e.g. `Synthetic - 5 runs`. |
| `detector` | `OCSVM`, `LOF`, or `Isolation_Forest`. |
| `AUMVC` | Mean raw AUMVC ± standard deviation; smaller is better. |
| `AUMVC_normalized` | Mean AUMVC divided by the reference-box volume ± standard deviation. Removes the overall box-volume scale; smaller is better. |
| `AUMVC_normalized_MC_SE` | Mean of the Monte Carlo standard errors of normalized AUMVC. Smaller means less MC sampling uncertainty. |
| `AUMVC_relative_MC_SE` | Mean of `normalized MC SE / normalized AUMVC` across runs; measures MC uncertainty relative to the estimate. |
| `minimum_hit_count` | The smallest number of sampled reference points inside an estimated level set. Low counts indicate poorly resolved volumes. |
| `reliable_runs` | Number of runs meeting both thresholds: `minimum_hit_count >= 10` and `AUMVC_relative_MC_SE <= 0.10`. It doesn't affect agreement calculations. |
| `ROC_AUC` | Mean ROC-AUC ± standard deviation; higher is better. |
| `PR_AUC` | Mean precision–recall AUC ± standard deviation; higher is better. |


### Comparison rows

| Column | Meaning |
|---|---|
| `comparison` | `strict` or `uncertainty_aware` (defined below). |
| `metric` | Always `ROC/PR consensus`: ROC-AUC and PR-AUC must prefer the same detector in a pair. |
| `matching_pairs` | Number of eligible detector pairs where AUMVC chooses the same winner as ROC/PR consensus. |
| `compared_pairs` | Number of eligible pairs actually compared. A pair is excluded if ROC-AUC and PR-AUC disagree or tie; `uncertainty_aware` also excludes unresolved AUMVC differences. |
| `total_pairs` | All possible detector-pair comparisons: **3 pairs × number of runs** per dataset. With 5 runs, this is 15. |
| `percentage` | `100 × matching_pairs / compared_pairs`. Blank when no pairs can be compared. |

**Strict:** Compares the AUMVC values directly. An exact tie produces no AUMVC winner.

**Uncertainty-aware:** Compares AUMVC using a **paired 95% Monte Carlo confidence interval** for the difference between two detectors (based on their shared reference-point contributions). If the interval includes zero, that pair is not counted. Thus, `compared_pairs` can be smaller than under `strict`.

**Example:** `matching_pairs = 11`, `compared_pairs = 13`, `total_pairs = 15` gives **84.62% agreement**. Then, it means 2 pairs have no ROC/PR consensus; among the 13 eligible pairs, AUMVC matches 11 pairs comparison.

For real data, the final `TOTAL ...` rows **sum the pair counts across datasets**, then recompute the percentage. They do not average dataset percentages.

**Outputs:** Both scripts print to the terminal and save one CSV each.
