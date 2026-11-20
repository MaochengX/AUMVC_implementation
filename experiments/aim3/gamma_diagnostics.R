# Aim 3: targeted gamma diagnostics for the existing interpretation.
# The output consists only of CSV measurements read by aim3_figures.py.
# NO level-set or MV-curve images are generated here: the original thesis
# already contains representative examples.
#
# Pilot (two real datasets, two embeddings, one run):
#   Rscript experiments/aim3/gamma_diagnostics.R
#
# Complete widened-grid check, only when needed for the exclusion analysis:
#   AIM3_DIAG_RUNS=10 AIM3_DIAG_GRID=widened \
#   AIM3_DIAG_DATASETS=ecg5000,shuttle,fashion_mnist,synthetic \
#   AIM3_DIAG_REPS=ambient,goix_subsampling,mds,isomap \
#   Rscript experiments/aim3/gamma_diagnostics.R
#
# The baseline original-grid and widened-grid results are read directly by
# aim3_figures.py and are NOT rerun by default. Additional grids can be run
# explicitly if their candidate-level exclusion diagnostics are needed.

source("experiments/aim3/real_experiment.R")
source("experiments/aim3/synthetic_experiment.R")
root <- "experiments/aim3"
grid_name <- Sys.getenv("AIM3_DIAG_GRID", "widened")
grids <- list(
  original = c(.25, .5, 1, 2, 4),
  widened = c(2^-6, 2^-2, 1, 4, 2^6),
  intermediate = c(.25, 1, 2, 4, 8, 16, 32, 64)
)
if (!grid_name %in% names(grids)) stop("Unknown gamma grid: ", grid_name)
settings <- AIM3_SETTINGS
settings$parallel$enabled <- FALSE
settings$bandwidth$gamma_scales <- grids[[grid_name]]
scales <- settings$bandwidth$gamma_scales
runs <- as.integer(Sys.getenv("AIM3_DIAG_RUNS", "1"))
if (length(runs) != 1L || !is.finite(runs) || runs < 1L) stop("Invalid run count")
parse_list <- function(key, default) trimws(strsplit(Sys.getenv(key, default), ",", fixed=TRUE)[[1L]])
datasets <- parse_list("AIM3_DIAG_DATASETS", "shuttle,ecg5000")
reps <- parse_list("AIM3_DIAG_REPS", "mds,isomap")
if (any(!datasets %in% c("synthetic", "ecg5000", "fashion_mnist", "shuttle"))) stop("Invalid dataset")
if (any(!reps %in% c("ambient", "goix_subsampling", "mds", "isomap"))) stop("Invalid representation")

out <- file.path(root, "interpretation", "diagnostics", grid_name)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
# Avoid showing stale gamma choices if a new diagnostic run is interrupted.
previous_choices <- file.path(out, "gamma_choices.csv")
if (file.exists(previous_choices)) file.remove(previous_choices)

# A: direct tests of the EXACT original R pair-ranking function.
check_sign_convention <- function() {
  test_settings <- AIM3_SETTINGS
  g <- test_settings$bandwidth$gamma_scales
  part <- data.frame(embedding_dim=5L, detector_dim=5L,
      representation="mds", gamma_scale=g,
      aumvc_normalized=seq_along(g), roc_auc=seq(.9,.5,length.out=length(g)),
      pr_auc=seq(.9,.5,length.out=length(g)), aumvc_reliable=TRUE)
  ans <- compare_bandwidths(part,test_settings)
  stopifnot(ans$compared_pairs_count==10L,
            isTRUE(all.equal(ans$agreement_pairs_rank_rate,1)))
  part$aumvc_normalized <- rev(part$aumvc_normalized)
  stopifnot(isTRUE(all.equal(compare_bandwidths(part,test_settings)$agreement_pairs_rank_rate,0)))
  part$aumvc_normalized <- c(1,1,3,4,5)
  stopifnot(isTRUE(all.equal(compare_bandwidths(part,test_settings)$tie_rate,.1)))
  part$roc_auc[1] <- part$roc_auc[2]
  stopifnot(compare_bandwidths(part,test_settings)$compared_pairs_count==9L)
  part$roc_auc <- seq(.9,.5,length.out=5)
  part$aumvc_reliable <- c(TRUE,FALSE,FALSE,FALSE,FALSE)
  stopifnot(is.na(compare_bandwidths(part,test_settings)$reliable_agreement_pairs_rank_rate))
  cat("Sign/tie/reliability tests: PASS; AUMVC down, ROC/PR up.\n")
}
check_sign_convention()

is_constant <- function(x) {
  if (!length(x) || any(!is.finite(x))) return(TRUE)
  diff(range(x)) <= 1e-10 * max(1, max(abs(x)))
}
is_score_good <- function(mv, ref_scores, eval_scores, label_scores) {
  flags <- mv_degeneracy(mv, ref_scores)
  numerical <- is.finite(mv$aumvc_normalized) &&
    is.finite(mv$aumvc_relative_mc_se) &&
    mv$minimum_hit_count >= settings$minimum_hit_count &&
    mv$aumvc_relative_mc_se <= settings$maximum_relative_mc_se
  extras <- c(reference_constant=is_constant(ref_scores),
              evaluation_constant=is_constant(eval_scores),
              label_constant=is_constant(label_scores))
  c(numerical=as.integer(numerical),
    threshold_collapsed=as.integer(flags[["threshold_collapsed"]]),
    reference_collapsed=as.integer(flags[["reference_collapsed"]]),
    full_box=as.integer(flags[["full_box"]]),
    reference_constant=as.integer(extras[["reference_constant"]]),
    evaluation_constant=as.integer(extras[["evaluation_constant"]]),
    label_constant=as.integer(extras[["label_constant"]]),
    valid=as.integer(numerical && !any(flags) &&
      !extras[["reference_constant"]] && !extras[["evaluation_constant"]]))
}

# Independent evaluation points set the normality threshold, NOT detector training points.
# The train-based MV curve below is a diagnostic counterfactual, not AUMVC tuning.
all_candidate <- list(); all_mass <- list(); all_include <- list()
all_rank <- list(); all_box <- list(); all_dups <- list()
append_record <- function(collection, row) {
  collection[[length(collection)+1L]] <- row
  collection
}

write_rows <- function(x, filename) {
  if (!length(x)) return(invisible(NULL))
  data <- do.call(rbind,x)
  write.csv(data,file.path(out,filename),row.names=FALSE,na="")
}

add_box <- function(dataset, case_id, run, rep, reference, x_reference) {
  lo <- reference$box$lower; hi <- reference$box$upper
  q1 <- apply(x_reference,2L,quantile,probs=.25)
  q3 <- apply(x_reference,2L,quantile,probs=.75)
  iqr <- q3-q1
  all_box[[length(all_box)+1L]] <<- data.frame(dataset=dataset,case_id=case_id,
     run=run,representation=rep,coordinate=seq_along(lo),
     lower=lo,upper=hi,box_width=hi-lo,reference_iqr=iqr,
     box_to_iqr=ifelse(iqr>0,(hi-lo)/iqr,NA_real_),gamma_grid=grid_name)
}

training_neighbours <- function(item, dataset,case_id,run,representation) {
  x <- item$eval; y <- item$train
  d2 <- outer(rowSums(x*x),rowSums(y*y),"+")-2*tcrossprod(x,y)
  distances <- sqrt(pmax(0,apply(d2,1L,min))/ncol(x))
  all_dups[[length(all_dups)+1L]] <<- data.frame(
     dataset=dataset,case_id=case_id,run=run,representation=representation,gamma_grid=grid_name,
     exact_pct=100*mean(distances < 1e-9),
     near_001_pct=100*mean(distances < .001),
     near_01_pct=100*mean(distances < .01),
     near_05_pct=100*mean(distances < .05),
     near_10_pct=100*mean(distances < .10),
     median_rms_distance=median(distances))
}

# We deliberately do not render any new level-set slices here. The existing
# matching Shuttle-Isomap MV and level-set figures remain the thesis examples.
# Reference-box bounds are recorded below to annotate those figures.

run_single <- function(ds, case_id, run, rep, item, p, eval_labels, label_labels, seed) {
  reference <- make_experiment_reference(item$reference,seed,settings)
  add_box(ds,case_id,run,rep,reference,item$reference)
  scored <- prepare_candidate_scores(item$train,item$eval,item$label,reference,scales,settings)
  for(i in seq_along(scales)) {
    gamma <- scales[i]; cand <- scored$candidates[[i]]
    mv <- aumvc_scored(cand$evaluation_scores,cand$reference_scores,reference,
                       score_direction="anomaly",alpha_grid=settings$aumvc_alpha_grid)
    checks <- is_score_good(mv,cand$reference_scores,cand$evaluation_scores,cand$label_scores)
    all_candidate[[length(all_candidate)+1L]] <<- data.frame(
       dataset=ds,case_id=case_id,run=run,representation=rep,gamma_grid=grid_name,
       gamma_scale=gamma,aumvc=mv$aumvc_normalized,
       min_hits=mv$minimum_hit_count,relative_se=mv$aumvc_relative_mc_se,
       as.list(checks),stringsAsFactors=FALSE)
    if(gamma %in% c(.25,4)) {
      training_scores <- score_ocsvm(cand$model,item$train)
      mv_train <- aumvc_scored(training_scores,cand$reference_scores,reference,
                               score_direction="anomaly",alpha_grid=settings$aumvc_alpha_grid)
      for(src in c("heldout","training")) {
        curve <- if(src=="heldout") mv$mv_curve else mv_train$mv_curve
        all_mass[[length(all_mass)+1L]] <<- data.frame(
          dataset=ds,case_id=case_id,run=run,representation=rep,gamma_grid=grid_name,
          gamma_scale=gamma,mass=curve$alpha,threshold_source=src,
          threshold=curve$threshold,occupancy=curve$volume_normalized,
          hit_count=curve$hit_count)
      }
      cutoff <- -mv$mv_curve$threshold[which.min(abs(mv$mv_curve$alpha-.95))]
      all_include[[length(all_include)+1L]] <<- data.frame(
         dataset=ds,case_id=case_id,run=run,representation=rep,gamma_grid=grid_name,
         gamma_scale=gamma,observation_id=seq_along(eval_labels),label=eval_labels,
         inside_mass95=as.integer(cand$evaluation_scores<=cutoff))
      all_rank[[length(all_rank)+1L]] <<- data.frame(
         dataset=ds,case_id=case_id,run=run,representation=rep,gamma_grid=grid_name,
         gamma_scale=gamma,observation_id=seq_along(label_labels),label=label_labels,
         anomaly_score_percentile=rank(cand$label_scores,ties.method="average")/length(cand$label_scores))
    }
  }
}

# Goix uses 50 different boxes; record its aggregate eligibility and individual
# subset failures. No single inverted/embedded Goix level set exists.
run_goix <- function(ds,case_id,run,p,seed) {
  subsets <- goix_subsets(ncol(p$train),seed,settings$goix_subsampling)
  values <- matrix(NA_real_,nrow=length(subsets),ncol=length(scales))
  flags <- matrix(FALSE,nrow=length(subsets),ncol=length(scales))
  numeric_pass <- matrix(FALSE,nrow=length(subsets),ncol=length(scales))
  min_hits <- matrix(NA_real_,nrow=length(subsets),ncol=length(scales))
  rel_se <- matrix(NA_real_,nrow=length(subsets),ncol=length(scales))
  reference_constant <- matrix(FALSE,nrow=length(subsets),ncol=length(scales))
  evaluation_constant <- matrix(FALSE,nrow=length(subsets),ncol=length(scales))
  label_constant <- matrix(FALSE,nrow=length(subsets),ncol=length(scales))
  full <- matrix(FALSE,nrow=length(subsets),ncol=length(scales))
  for(j in seq_along(subsets)) {
    columns <- subsets[[j]]
    ref <- make_experiment_reference(p$reference[,columns,drop=FALSE],seed+1000L+j,settings)
    scored <- prepare_candidate_scores(p$train[,columns,drop=FALSE],
               p$eval[,columns,drop=FALSE],p$label[,columns,drop=FALSE],ref,scales,settings)
    for(i in seq_along(scales)) {
      ca <- scored$candidates[[i]]
      mv <- aumvc_scored(ca$evaluation_scores,ca$reference_scores,ref,
                         score_direction="anomaly",alpha_grid=settings$aumvc_alpha_grid)
      ck <- is_score_good(mv,ca$reference_scores,ca$evaluation_scores,ca$label_scores)
      values[j,i] <- mv$aumvc_normalized
      flags[j,i] <- as.logical(ck[["valid"]])
      numeric_pass[j,i] <- as.logical(ck[["numerical"]])
      min_hits[j,i] <- mv$minimum_hit_count
      rel_se[j,i] <- mv$aumvc_relative_mc_se
      reference_constant[j,i] <- as.logical(ck[["reference_constant"]])
      evaluation_constant[j,i] <- as.logical(ck[["evaluation_constant"]])
      label_constant[j,i] <- as.logical(ck[["label_constant"]])
      full[j,i] <- as.logical(ck[["full_box"]])
    }
  }
  for(i in seq_along(scales)) {
    all_candidate[[length(all_candidate)+1L]] <<- data.frame(
      dataset=ds,case_id=case_id,run=run,representation="goix_subsampling",
      gamma_grid=grid_name,gamma_scale=scales[i],aumvc=mean(values[,i]),
      min_hits=min(min_hits[,i]),relative_se=max(rel_se[,i]),
      numerical=as.integer(all(numeric_pass[,i])),threshold_collapsed=NA_integer_,
      reference_collapsed=NA_integer_,full_box=as.integer(any(full[,i])),
      reference_constant=as.integer(any(reference_constant[,i])),
      evaluation_constant=as.integer(any(evaluation_constant[,i])),
      label_constant=as.integer(any(label_constant[,i])),
      valid=as.integer(all(flags[,i])))
  }
}

registry <- real_registry()
real_names <- intersect(datasets,names(registry))
real_data <- lapply(real_names,function(ds) registry[[ds]]$load(settings$real[[ds]]$data_dir))
names(real_data) <- real_names

for(run in seq_len(runs)) {
  for(ds in datasets) {
    tasks <- list()
    if(ds=="synthetic") {
      for(truth in settings$synthetic$truth_functions) {
        row <- data.frame(truth=truth,intrinsic_dim=5L,original_dim=200L,snr=5)
        seed <- synthetic_seed(row,run,settings$synthetic,settings$seed)
        tasks[[truth]] <- list(case=generate_synthetic(5,200,5,truth,seed,settings$synthetic),
                              seed=seed+50000L,split_counts=AIM3_SETTINGS$split_counts)
      }
    } else {
      idx <- match(ds,settings$real$datasets)
      seed <- settings$seed+run*100000L+idx*10000L
      set.seed(seed,kind="Mersenne-Twister",normal.kind="Inversion",sample.kind="Rejection")
      sampled <- registry[[ds]]$sample(real_data[[ds]],seed,settings$real[[ds]])
      tasks[[ds]] <- list(case=list(x=sampled$x,labels=sampled$labels,
        intrinsic_dim=as.integer(settings$real$embedding_dim)),
        seed=seed+50000L,split_counts=real_split_counts(nrow(sampled$x)))
    }
    for(case_id in names(tasks)) {
      item <- tasks[[case_id]]
      settings$split_counts <- item$split_counts
      strata <- if(is.null(item$case$outlier_type)) item$case$labels else item$case$outlier_type
      ix <- make_splits(item$case$labels,settings$split_counts,item$seed,strata)
      labels_eval <- item$case$labels[ix$evaluation]
      labels_test <- item$case$labels[ix$label_eval]
      p <- prepare_case_data(item$case,item$seed,settings)
      stopifnot(identical(as.integer(p$labels_label),as.integer(labels_test)))
      for(rep in reps) {
        cat(sprintf("Gamma grid=%s, run %d/%d, dataset=%s, case=%s, space=%s\n",
                    grid_name,run,runs,ds,case_id,rep))
        if(rep=="goix_subsampling") {
          run_goix(ds,case_id,run,p,item$seed+2000L)
        } else {
          emb <- if(rep=="ambient") p else prepare_embedding_data(p,rep,settings)
          training_neighbours(emb,ds,case_id,run,rep)
          refseed <- item$seed+switch(rep,ambient=1000L,mds=4000L,isomap=5000L)
          run_single(ds,case_id,run,rep,emb,p,labels_eval,labels_test,refseed)
        }
      }
    }
  }
  # Checkpoint after every run so a long experiment does not lose all results.
  write_rows(all_candidate,"gamma_candidate_diagnostics.csv")
  write_rows(all_mass,"mass_profiles.csv")
  write_rows(all_include,"anomaly_inclusion.csv")
  write_rows(all_rank,"score_ranks.csv")
  write_rows(all_box,"reference_box_bounds.csv")
  write_rows(all_dups,"duplicate_checks.csv")
}

# Critical: keep separate synthetic generator cases within each run.
all <- do.call(rbind,all_candidate)
bycase <- split(all,interaction(all$dataset,all$case_id,all$run,all$representation,drop=TRUE))
selected <- lapply(bycase,function(h) {
  # Candidates excluded for MC failures, full box, collapsed thresholds, and
  # near-constant scores. No anomaly labels enter the exclusion decision.
  ok <- !is.na(h$valid) & as.logical(h$valid) & is.finite(h$aumvc)
  n <- sum(ok)
  choice <- NA_real_
  status <- if (n == 0L) "no_eligible" else if (n == 1L) "one_eligible" else "tie"
  if (n >= 2L) {
    z <- h[ok, , drop=FALSE]
    best <- min(z$aumvc)
    tolerance <- settings$bandwidth$tie_tolerance
    tied <- abs(z$aumvc-best) <= tolerance * pmax(1, abs(best), abs(z$aumvc))
    if (sum(tied) == 1L) {
      choice <- z$gamma_scale[which(tied)]
      status <- if (n == nrow(h)) "selected_complete" else "selected_partial"
    }
  }
  data.frame(dataset=h$dataset[1],case_id=h$case_id[1],run=h$run[1],
      representation=h$representation[1],gamma_grid=grid_name,
      n_eligible=n,n_excluded=nrow(h)-n,
      chosen_gamma=choice,status=status)
})
write_rows(selected,"gamma_choices.csv")
cat("Done. Diagnostic CSVs are in: ",out,"\n",sep="")
cat("No figures were generated by this diagnostic R script.\n")
