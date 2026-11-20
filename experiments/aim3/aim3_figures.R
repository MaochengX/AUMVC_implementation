draw_bars <- function(labels, series, main='',ylab='',ylim=NULL,fmt='%.1f',
                      xlabel='',zero_label=TRUE, legend_here=FALSE) {
  n <- length(labels);m <- length(series)
  values <- do.call(rbind,lapply(series,function(s) as.numeric(s$values)))
  if(ncol(values)!=n) stop('Plot matrix does not match axis labels: ',main)
  maxfinite <- suppressWarnings(max(values[is.finite(values)],na.rm=TRUE))
  if(!is.finite(maxfinite)) maxfinite <- 1
  ymax <- if(is.null(ylim)) max(1,maxfinite*1.17) else ylim
  bar_w <- min(.78/m,.24)
  graphics::plot(NA,NA,xlim=c(.4,n+.6),ylim=c(0,ymax),type='n',xaxt='n',yaxt='n',
                 xlab=xlabel,ylab=ylab,bty='n',xaxs='i',yaxs='i',las=1)
  gridby <- if(ymax>=100) seq(0, floor(ymax/20)*20, by=20) else pretty(c(0,ymax),n=5)
  abline(h=gridby,col='#eeeeee',lwd=.8)
  axis(1,at=seq_len(n),labels=labels,lwd=1,cex.axis=.83)
  at <- gridby[gridby>=0 & gridby<=ymax];axis(2,at=at,labels=format(at,trim=TRUE),lwd=1,cex.axis=.82)
  segments(.4,0,n+.6,0,col='#111111',lwd=.95)
  segments(.4,0,.4,ymax,col='#111111',lwd=.95)
  for(j in seq_len(m)) for(i in seq_len(n)) {
    x <- i + (j-(m+1)/2)*bar_w
    v <- values[j,i]; if(!is.finite(v)) {
      text(x,.025*ymax,'NA',srt=90,adj=c(0,.5),cex=.65,col='#64727d');next
    }
    if(v!=0) rect(x-.47*bar_w,0,x+.47*bar_w,v,col=series[[j]]$color,border=NA)
    if(v==0 && !zero_label) next
    lab <- sprintf(fmt,v)
    if(v>.79*ymax) text(x,v-.027*ymax,lab,srt=90,adj=c(1,.5),col='white',cex=.70,font=2)
    else text(x,max(v,.008*ymax)+.013*ymax,lab,srt=90,adj=c(0,.5),col='#3c4b59',cex=.68,xpd=NA)
  }
  if(nzchar(main)) title(main=main,font.main=1,cex.main=.91,line=.6)
  if(legend_here) legend('topright',legend=vapply(series,`[[`,character(1),'name'),
                          fill=vapply(series,`[[`,character(1),'color'),bty='n',cex=.72)
  invisible(values)
}

spec <- function(name,color,values) list(name=name,color=color,values=values)

header <- function(title,leg=NULL,colors=NULL,cex=1.0) {
  par(mar=c(0,0,0,0));plot.new();plot.window(c(0,1),c(0,1))
  text(.5,.74,title,cex=cex)
  if(!is.null(leg)) legend('center',legend=leg,fill=colors,horiz=TRUE,
                           bty='n',cex=.8,x.intersp=.65,y.intersp=.65)
}

layout_panels <- function(nrow,ncol,title,leg=NULL,colors=NULL,height=0.52) {
  head <- matrix(1L,nrow=1,ncol=ncol)
  content <- matrix(seq_len(nrow*ncol)+1L,nrow=nrow,ncol=ncol,byrow=TRUE)
  layout(rbind(head,content),heights=c(height,rep(4,nrow)))
  header(title,leg,colors)
  par(mar=c(4.3,4.5,2.2,.8),las=1)
}

draw_line <- function(x,y,cols,leg,main='',ylim=c(0,1.07),xticks=NULL,xticklabels=NULL) {
  plot(NA,NA,xlim=range(x,finite=TRUE),ylim=ylim,xaxt='n',xlab='',ylab='ROC-AUC / reliable fraction',bty='n')
  abline(h=pretty(ylim),col='#eeeeee')
  if(is.null(xticks)) axis(1) else axis(1,at=xticks,labels=xticklabels)
  for(j in seq_along(y)) lines(x,y[[j]],type='b',pch=15+j,col=cols[j],lwd=2,lty=if(j==3) 2 else 1)
  title(main=main,font.main=1,cex.main=.91)
}

args <- commandArgs(trailingOnly = TRUE)

arg <- function(key, default) {
  p <- args[startsWith(args, paste0(key, '='))]
  if (length(p)) substring(tail(p, 1), nchar(key) + 2L) else default
}

ROOT <- arg('--aim3', 'experiments/aim3')

OUT <- arg('--out', file.path(ROOT, 'interpretation', 'visulization'))

FIGS <- file.path(OUT, 'figures'); TABLES <- file.path(OUT, 'tables')

dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)

dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)

REPS <- c('ambient','goix_subsampling','mds','isomap')

NAMES <- c('Ambient','Goix','MDS','Isomap')

DS <- c('synthetic','ecg5000','fashion_mnist','shuttle')

DSN <- setNames(c('Synthetic','ECG5000','Fashion-MNIST','Shuttle'),DS)

GEN <- c(polynomial_interaction='Polynomial',oscillatory_local='Oscillatory')

COLORS <- c('#bd6047','#168676','#326eaf','#8b59a7')

METH <- c('AUMVC','ROC','PR')

METHOD_COLORS <- c(AUMVC='#AD6E45',ROC='#316FA5',PR='#168676')

GAMMAS <- c(.25,.5,1,2,4)

WIDE_GAMMAS <- c(1/64,.25,1,4,64)

INPUTS <- c(
 syn='with ECG5000/50run/synthetic/experiment/syn_results.csv',
 sb='with ECG5000/50run/synthetic/experiment/syn_bandwidth_exp.csv',
 real='with ECG5000/50run/real/experiment/real_results.csv',
 rb='with ECG5000/50run/real/experiment/real_bandwidth_exp.csv',
 ws='with ECG5000/widen_result/synthetic/widen_experiment/experiment/syn_results.csv',
 wss='with ECG5000/widen_result/synthetic/widen_experiment/analysis/selection.csv',
 wsc='with ECG5000/widen_result/synthetic/widen_experiment/analysis/score_checks.csv',
 wsg='with ECG5000/widen_result/synthetic/widen_experiment/analysis/grid_diagnostics.csv',
 wr='with ECG5000/widen_result/real/widen_experiment/experiment/real_results.csv',
 wrs='with ECG5000/widen_result/real/widen_experiment/analysis/selection.csv',
 wrc='with ECG5000/widen_result/real/widen_experiment/analysis/score_checks.csv',
 wrg='with ECG5000/widen_result/real/widen_experiment/analysis/grid_diagnostics.csv',
 scan_current='dim_scan/current_generator/scan_results.csv',
 scan_original='dim_scan/original_orthogonal_shift/scan_results.csv',
 scan_real='dim_scan/real_data/scan_results.csv',
 mv='mv_curve_exp/summary.csv',
 bootstrap='with ECG200/widen_experiment/synthetic/widen_experiment/analysis/paired_bootstrap.csv'
)

KEEP_FIGURES <- c(
  'synthetic_low_hit',
  'synthetic_pair_agreement',
  'intrinsic_dimension_agreement',
  'real_reliability_agreement',
  'real_runtime',
  'synthetic_outlier_types',
  'current_generator_dimension_scan',
  'orthogonal_shift_dimension_scan',
  'gamma64_diagnostics',
  'gamma_selection_synthetic',
  'gamma_selection_ecg5000',
  'widened_roc_variation',
  'widened_roc_grid_edge',
  'paired_bootstrap',
  'gamma_profile_ecg5000',
  'gamma_profile_fashion_mnist',
  'gamma_profile_shuttle',
  'gamma_thresholds_volumes',
  'gamma4_choice_only',
  'gamma_choice_original_grid',
  '1g_regret_default_random',
  '1f_real_mds_scan'
)

KEEP_TABLES <- c('synthetic_candidate_reliability','widened_selection_archived_screen','anomaly_score_rank_comparison')

FIGURE_TITLES <- c(
 synthetic_low_hit='Low-Hit Mass Levels (%)',
 synthetic_pair_agreement='Synthetic Pairwise Ranking Agreement (%)',
 intrinsic_dimension_agreement='Pairwise Agreement by Intrinsic Dimension (%)',
 real_reliability_agreement='Real-Data Reliability and Pairwise Agreement',
 real_runtime='Real Runtime',
 synthetic_outlier_types='Distributional and Structural Outlier Detection',
 current_generator_dimension_scan='Current Generator: MDS Dimension Scan',
 orthogonal_shift_dimension_scan='Orthogonal Shift: MDS Dimension Scan',
 widened_roc_variation='ROC-AUC Variation across the Gamma Grid',
 widened_roc_grid_edge='ROC-Best Gamma at Grid Boundary (%)',
 gamma64_diagnostics='Gamma 64: Score Degeneracy and Full-Box Rate',
 paired_bootstrap='Paired-Bootstrap ROC Gaps',
 gamma_selection_synthetic='Synthetic Gamma Selection Frequencies (%)',
 gamma_selection_ecg5000='ECG5000 Gamma Selection Frequencies (%)',
 gamma_profile_ecg5000='ECG5000: Original-Grid AUMVC and ROC by Gamma',
 gamma_profile_fashion_mnist='Fashion-MNIST: Original-Grid AUMVC and ROC by Gamma',
 gamma_profile_shuttle='Shuttle: Original-Grid AUMVC and ROC by Gamma',
 gamma_thresholds_volumes='Mass Thresholds and Reference-Box Volumes',
 gamma4_choice_only='Gamma choice {0.015625, 0.25, 1, 4, 64}',
 gamma_choice_original_grid='Gamma choice {0.25, 0.5, 1, 2, 4}'
)

read_csv <- function(path, optional=FALSE) {
  if (!file.exists(path)) {
    if (optional) return(NULL)
    stop('Missing original experiment CSV: ',path,call.=FALSE)
  }
  read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)
}

read_input <- function(key) read_csv(file.path(ROOT, INPUTS[[key]]))

num <- function(x) {if (is.logical(x)) return(as.numeric(x)); suppressWarnings(as.numeric(as.character(x)))}

average <- function(x) {x <- num(x);x <- x[is.finite(x)];if(length(x)) mean(x) else NA_real_}

se <- function(x) {x <- num(x);x <- x[is.finite(x)];if(length(x)>1) sd(x)/sqrt(length(x)) else NA_real_}

as_pct <- function(x) 100*average(x)

count_finite <- function(x) sum(is.finite(num(x)))

median_finite <- function(x) {x <- num(x);x <- x[is.finite(x)];if(length(x)) median(x) else NA_real_}

require_cols <- function(x,cols,label) {
  miss <- setdiff(cols,names(x))
  if(length(miss)) stop(label,' is missing columns: ',paste(miss,collapse=', '),call.=FALSE)
}

rbind_fill <- function(...) {
  parts <- list(...)
  if(length(parts)==1L && is.list(parts[[1L]]) && !is.data.frame(parts[[1L]])) parts <- parts[[1L]]
  parts <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x)>0, parts)
  if(!length(parts)) return(data.frame())
  allcols <- unique(unlist(lapply(parts,names)))
  parts <- lapply(parts,function(p) {for(col in setdiff(allcols,names(p))) p[[col]] <- NA;p[,allcols,drop=FALSE]})
  z <- do.call(rbind,parts);rownames(z) <- NULL;z
}

case_id <- function(d,cols) {
  require_cols(d,cols,'case ID')
  if(!nrow(d)) return(character())
  do.call(paste,c(lapply(cols,function(c) ifelse(is.na(d[[c]]),'[NA]',as.character(d[[c]]))),sep='\037'))
}

group_rows <- function(d,cols,fun) {
  if(!nrow(d)) return(data.frame())
  splits <- split(seq_len(nrow(d)),case_id(d,cols),drop=TRUE)
  z <- lapply(splits,function(i) {
    g <- d[i,,drop=FALSE]
    result <- as.data.frame(fun(g),check.names=FALSE,stringsAsFactors=FALSE)
    data.frame(g[1,cols,drop=FALSE],result,check.names=FALSE)
  })
  rbind_fill(z)
}

lookup <- function(d, criteria, col) {

  if (!is.data.frame(d))
    stop('lookup expected a data.frame for column ', col, '; received ',
         paste(class(d), collapse='/'), call.=FALSE)
  if (!col %in% names(d)) stop('lookup column is missing: ', col, call.=FALSE)
  missing <- setdiff(names(criteria), names(d))
  if (length(missing)) stop('lookup grouping columns are missing: ', paste(missing, collapse=', '), call.=FALSE)
  z <- rep.int(TRUE, NROW(d))
  for (k in names(criteria)) {
    target <- criteria[[k]]
    z <- z & !is.na(d[[k]]) & as.character(d[[k]]) == as.character(target)
  }
  a <- d[z, col, drop=TRUE]
  if (!length(a)) return(NA_real_)
  if (length(a)>1L) stop('Non-unique lookup for ', col, ' / ',
    paste(names(criteria), unlist(criteria), sep='=', collapse=', '), call.=FALSE)
  num(a[1])
}

table_out <- function(key,d,source='') {
  if(!(key %in% KEEP_TABLES)) return(invisible(d))
  write.csv(d,file.path(TABLES,paste0(key,'.csv')),row.names=FALSE,na='NA',fileEncoding='UTF-8')
  TABLE_WRITTEN <<- unique(c(TABLE_WRITTEN,key))
  invisible(d)
}

TABLE_WRITTEN <- character(); FIGURE_WRITTEN <- character()

with_png <- function(key, width, height, code, legend=NULL, panels=NULL) {
  if (!(key %in% KEEP_FIGURES)) return(invisible(NULL))
  png(file.path(FIGS,paste0(key,'.png')),width=width,height=height,units='in',res=190,bg='white')
  on.exit(dev.off(),add=TRUE)
  par(family='sans',fg='black',col.axis='black',col.lab='black',las=1)
  eval(substitute(code),envir=parent.frame())
  FIGURE_WRITTEN <<- unique(c(FIGURE_WRITTEN,key))
  invisible(key)
}

run_mean_se <- function(d,metric,by) {
  temp <- group_rows(d,c(by,'run'),function(g) data.frame(value=average(g[[metric]])))
  group_rows(temp,by,function(g) data.frame(mean=average(g$value),sd=sd(g$value),n_runs=count_finite(g$value),se=se(g$value)))
}

S <- read_input('syn');SB <- read_input('sb'); R <- read_input('real');RB <- read_input('rb')

WS <- read_input('ws');WR <- read_input('wr');WSS <- read_input('wss');WRS <- read_input('wrs')

WSC <- read_input('wsc');WRC <- read_input('wrc');WSG <- read_input('wsg');WRG <- read_input('wrg')

SC <- read_input('scan_current');SO <- read_input('scan_original');SR <- read_input('scan_real')

MV <- read_input('mv');BOOT <- read_input('bootstrap')

pair_counts <- function(archived,keys) {
  require_cols(archived,c(keys,'compared_pairs_count','unreliable_pairs_rate',
    'agreement_pairs_rank_rate','reliable_agreement_pairs_rank_rate','tie_rate'),'Saved R bandwidth rates')
  z <- archived[,keys,drop=FALSE]
  z$compared_pairs <- as.integer(na_zero(num(archived$compared_pairs_count)))
  uv <- num(archived$unreliable_pairs_rate);uv[!is.finite(uv)] <- 1
  z$reliable_compared_pairs <- as.integer(round(z$compared_pairs*(1-uv)))
  ag <- num(archived$agreement_pairs_rank_rate);ag[!is.finite(ag)] <- 0
  rag <- num(archived$reliable_agreement_pairs_rank_rate);rag[!is.finite(rag)] <- 0
  tie <- num(archived$tie_rate);tie[!is.finite(tie)] <- 0
  z$matching_pairs <- as.integer(round(z$compared_pairs*ag))
  z$reliable_matching_pairs <- as.integer(round(z$reliable_compared_pairs*rag))
  z$tied_pairs <- as.integer(round(z$compared_pairs*tie))
  z$raw_pair_agreement <- num(archived$agreement_pairs_rank_rate)
  z$reliable_pair_agreement <- num(archived$reliable_agreement_pairs_rank_rate)
  if(any(z$matching_pairs>z$compared_pairs) || any(z$reliable_matching_pairs>z$reliable_compared_pairs))
    stop('Invalid archived pair-count reversal')
  good <- z$compared_pairs>0 & is.finite(z$raw_pair_agreement)
  if(any(abs(z$matching_pairs[good]/z$compared_pairs[good]-z$raw_pair_agreement[good])>6e-7))
    stop('Raw pair count inversion inconsistent with archived rates')
  good <- z$reliable_compared_pairs>0 & is.finite(z$reliable_pair_agreement)
  if(any(abs(z$reliable_matching_pairs[good]/z$reliable_compared_pairs[good]-z$reliable_pair_agreement[good])>6e-7))
    stop('Reliable pair count inversion inconsistent with archived rates')
  z
}

na_zero <- function(x) {x[!is.finite(x)] <- 0;x}

PS0 <- pair_counts(SB,c('truth','intrinsic_dim','original_dim','snr','run','representation'))

PR0 <- pair_counts(RB,c('dataset','run','representation'))

summarize_pairs <- function(d,by) group_rows(d,by,function(g) {
  a <- sum(g$compared_pairs);b <- sum(g$reliable_compared_pairs)
  data.frame(n_cases=nrow(g),n_consensus=a,n_matching=sum(g$matching_pairs),
    n_valid_pairs=b,n_valid_matching=sum(g$reliable_matching_pairs),
    n_tied=sum(g$tied_pairs),n_supported_cases=count_finite(g$reliable_pair_agreement),
    mean_raw_agreement=average(g$raw_pair_agreement),
    mean_valid_agreement=average(g$reliable_pair_agreement),
    reliable_pair_coverage_pct=if(a>0) 100*b/a else NA_real_,
    weighted_raw_agreement_pct=if(a>0) 100*sum(g$matching_pairs)/a else NA_real_,
    weighted_valid_agreement_pct=if(b>0) 100*sum(g$reliable_matching_pairs)/b else NA_real_,
    supported_case_pct=100*count_finite(g$reliable_pair_agreement)/nrow(g))
})

ps <- summarize_pairs(PS0,c('original_dim','representation'))

pr <- summarize_pairs(PR0,c('dataset','representation'))

summarize_candidates <- function(d,by) group_rows(d,by,function(g)
  data.frame(n=nrow(g),low_hit_pct=as_pct(g$low_occupancy),
   zero_hit_pct=as_pct(g$zero_occupancy),reliable_pct=as_pct(g$aumvc_reliable),
   roc=average(g$roc_auc),pr=average(g$pr_auc),runtime=average(g$runtime_sec),
   mc_se=average(g$aumvc_relative_mc_se),min_hit=average(g$min_hit_count),
   oos_stress=average(g$embedding_oos_stress)))

cs <- summarize_candidates(S,c('original_dim','representation'))

cr <- summarize_candidates(R,c('dataset','representation'))

subtype <- group_rows(S,c('truth','representation'),function(g)
  data.frame(n=nrow(g),roc_dist=average(g$roc_distributional),
    roc_struct=average(g$roc_structural),pr_dist=average(g$pr_distributional),
    pr_struct=average(g$pr_structural)))

scan_synthetic <- function(d) {
 d$dimension_ratio <- ifelse(d$representation=='ambient',0,num(d$embedding_dim)/num(d$intrinsic_dim))
 group_rows(d[num(d$gamma_scale)==1,,drop=FALSE],c('truth','dimension_ratio'),function(g)
  data.frame(n=nrow(g),distributional_roc=average(g$roc_distributional),
    structural_roc=average(g$roc_structural),reliable_pct=as_pct(g$aumvc_reliable)))
}

curr <- scan_synthetic(SC);ortho <- scan_synthetic(SO)

WSG$dataset <- 'synthetic'

diag <- rbind_fill(WSG,WRG)

extremes <- rbind_fill(lapply(c('synthetic','real'),function(kind){
 candidates <- if(kind=='synthetic') WS else WR
 checks <- if(kind=='synthetic') WSC else WRC
 by <- if(kind=='synthetic') 'representation' else c('dataset','representation')
 f <- group_rows(candidates[num(candidates$gamma_scale)==64,,drop=FALSE],by,function(g)
  data.frame(old_reliable_pct=as_pct(g$aumvc_reliable),n=nrow(g),mean_aumvc=average(g$aumvc_normalized)))
 z <- group_rows(checks[num(checks$gamma_scale)==64,,drop=FALSE],by,function(g)
  data.frame(reference_constant_pct=as_pct(g$reference_near_constant),
   eval_constant_pct=as_pct(g$eval_near_constant),label_constant_pct=as_pct(g$label_near_constant)))
 f <- merge(f,z,by=by,all.x=TRUE,sort=FALSE)
 if(kind=='synthetic') f$dataset <- 'synthetic'
 f
}))

mv64 <- MV[num(MV$gamma_scale)==64 & MV$score_basis=='single_model',,drop=FALSE]

ps200 <- pair_counts(SB[num(SB$original_dim)==200,,drop=FALSE],
 c('truth','intrinsic_dim','original_dim','snr','run','representation'))

q200 <- group_rows(ps200,c('intrinsic_dim','representation'),function(g)
 data.frame(n_cases=nrow(g),n_consensus=sum(g$compared_pairs),
   n_reliable_pairs=sum(g$reliable_compared_pairs),
   n_supported=count_finite(g$reliable_pair_agreement),
   raw_agreement=average(g$raw_pair_agreement),
   reliable_agreement=average(g$reliable_pair_agreement)))

pivot_reliability <- function(df,group,values) {
  keys <- c(group,c('ambient','goix_subsampling','mds','isomap'))
  rows <- lapply(values,function(v) {
    x <- c(as.character(v),vapply(REPS,function(rep) {
      z <- df[num(df[[group]])==v & df$representation==rep,,drop=FALSE]
      if(nrow(z)!=1L) stop('Missing/duplicate candidate reliability for ',group,'=',v,', ',rep)
      sprintf('%.1f%%',num(z$reliable_pct[1]))
    },character(1)))
    as.list(setNames(x,c(if(group=='original_dim') 'd' else 'q',NAMES)))
  })
  out <- rbind_fill(lapply(rows,as.data.frame,stringsAsFactors=FALSE))
  names(out) <- c(if(group=='original_dim') 'd' else 'q',NAMES)
  out
}

table_out('synthetic_candidate_reliability',pivot_reliability(cs,'original_dim',c(80,200,400)))

regret_rows <- list()

for(ds in DS) {
 choices <- if(ds=='synthetic') WSS else WRS[WRS$dataset==ds,,drop=FALSE]
 candidates <- if(ds=='synthetic') WS else WR[WR$dataset==ds,,drop=FALSE]
 keys <- c('representation','run',if(ds=='synthetic') c('truth','intrinsic_dim','original_dim','snr'))
 cgroups <- split(seq_len(nrow(candidates)),case_id(candidates,keys),drop=TRUE)
 for(rep in REPS) {
   h <- choices[choices$representation==rep,,drop=FALSE]
   if(!nrow(h)) next
   rows <- list()
   good <- h[h$status=='selected' & is.finite(num(h$chosen_gamma)),,drop=FALSE]
   for(i in seq_len(nrow(good))) {
     rec <- good[i,,drop=FALSE];id <- case_id(rec,keys)
     idx <- cgroups[[id]]
     if(is.null(idx)) next
     x <- candidates[idx,,drop=FALSE]
     if(nrow(x)!=5L) next
     chosen <- which(abs(num(x$gamma_scale)-num(rec$chosen_gamma))<1e-12)
     base <- which(abs(num(x$gamma_scale)-1)<1e-12)
     if(length(chosen)!=1L || length(base)!=1L) next
     roc_scores <- num(x$roc_auc);pr_scores <- num(x$pr_auc)
     if(any(!is.finite(roc_scores)) || any(!is.finite(pr_scores))) next
     ran <- max(roc_scores)-min(roc_scores)
     rows[[length(rows)+1L]] <- data.frame(run=rec$run,
       roc_regret=max(roc_scores)-roc_scores[chosen],pr_regret=max(pr_scores)-pr_scores[chosen],
       default_regret=max(roc_scores)-roc_scores[base],random_regret=max(roc_scores)-mean(roc_scores),
       relative_regret=if(ran>0) (max(roc_scores)-roc_scores[chosen])/ran else NA_real_)
   }
   cases <- rbind_fill(rows);n_ok <- nrow(cases);n_total <- nrow(h)
   n_runs <- length(unique(h$run));eligible_runs <- if(n_ok) length(unique(cases$run)) else 0L
   result <- data.frame(dataset=ds,representation=rep,n_cases=n_total,
     n_eligible=n_ok,fail_rate=100*(n_total-n_ok)/n_total,
     n_runs=n_runs,eligible_runs=eligible_runs)
   for(metric in c('roc_regret','pr_regret','default_regret','random_regret','relative_regret')) {
     if(n_ok) {
       byrun <- group_rows(cases,'run',function(g) data.frame(v=average(g[[metric]])))
       result[[paste0(metric,'_mean')]] <- average(byrun$v)
       result[[paste0(metric,'_se')]] <- se(byrun$v)
     } else {
       result[[paste0(metric,'_mean')]] <- NA_real_
       result[[paste0(metric,'_se')]] <- NA_real_
     }
   }
   regret_rows[[length(regret_rows)+1L]] <- result
 }
}

regret <- rbind_fill(regret_rows)

widened_screen <- rbind_fill(lapply(DS,function(ds){
 choices <- if(ds=='synthetic') WSS else WRS[WRS$dataset==ds,,drop=FALSE]
 rbind_fill(lapply(REPS,function(rep){
  g <- choices[choices$representation==rep,,drop=FALSE]
  selected <- g[g$status=='selected' & is.finite(num(g$chosen_gamma)),,drop=FALSE]
  out <- data.frame(Dataset=DSN[[ds]],Representation=NAMES[match(rep,REPS)],
    Cases=nrow(g),Selected=nrow(selected),No_choice=nrow(g)-nrow(selected),check.names=FALSE)
  for(gamma in WIDE_GAMMAS) out[[format(gamma,scientific=FALSE,trim=TRUE)]] <-
    if(nrow(g)) 100*sum(abs(num(selected$chosen_gamma)-gamma)<1e-12)/nrow(g) else NA_real_
  out
 }))
}))

table_out('widened_selection_archived_screen',widened_screen)

choose_max_unique <- function(x,gammas) {
 z <- num(x);valid <- which(is.finite(z))
 if(!length(valid)) return(NA_real_)
 winners <- valid[z[valid]==max(z[valid])]
 if(length(winners)!=1L) NA_real_ else num(gammas[winners])
}

selection_all <- rbind_fill(lapply(DS,function(ds) {
 syn <- ds=='synthetic'
 candidates <- if(syn) S else R[R$dataset==ds,,drop=FALSE]
 bw <- if(syn) SB else RB[RB$dataset==ds,,drop=FALSE]
 candidates$dataset <- ds;bw$dataset <- ds
 keys <- c('dataset','representation','run',if(syn) c('truth','intrinsic_dim','original_dim','snr'))
 decision_id <- case_id(bw,keys)
 if(anyDuplicated(decision_id)) stop('Duplicated saved AUMVC decisions for ',ds)
 choices <- setNames(num(bw$selected_gamma_scale_reliable),decision_id)
 group_rows(candidates,keys,function(g) {

  g <- g[!duplicated(g$gamma_scale),,drop=FALSE]
  g <- g[order(num(g$gamma_scale)),,drop=FALSE]
  id <- case_id(g[1,,drop=FALSE],keys)
  data.frame(ROC=choose_max_unique(g$roc_auc,g$gamma_scale),
    PR=choose_max_unique(g$pr_auc,g$gamma_scale),
    AUMVC=if(id %in% names(choices)) unname(choices[[id]]) else NA_real_)
 })
}))

frequency <- rbind_fill(lapply(DS,function(ds) rbind_fill(lapply(REPS,function(rep){
 h <- selection_all[selection_all$dataset==ds & selection_all$representation==rep,,drop=FALSE]
 rbind_fill(lapply(METH,function(method) rbind_fill(lapply(GAMMAS,function(gamma){
  n <- sum(num(h[[method]])==gamma,na.rm=TRUE)
  data.frame(dataset=ds,representation=rep,method=method,gamma_scale=gamma,
    n=n,total=nrow(h),pct=if(nrow(h)) 100*n/nrow(h) else NA_real_,
    n_missing=sum(is.na(h[[method]])))
 }))))
}))))

with_png('gamma_choice_original_grid',11.2,7.4,{
 layout_panels(2,2,FIGURE_TITLES[['gamma_choice_original_grid']],
   c('0.25','0.5','1','2','4'),
   c('#C9A27B','#B57743','#8F654D','#A66C3D','#C0946D'))
 gamma_cols <- c(0.25,0.5,1,2,4)
 gamma_fill <- c('#C9A27B','#B57743','#8F654D','#A66C3D','#C0946D')
 for(ds in DS) {
  bars <- lapply(seq_along(gamma_cols), function(k) {
   vals <- vapply(REPS,function(rep) {
    row <- frequency[frequency$dataset==ds & frequency$representation==rep &
                     frequency$method=='AUMVC' & abs(num(frequency$gamma_scale)-gamma_cols[k])<1e-12,,drop=FALSE]
    if(!nrow(row)) return(NA_real_)
    num(row$pct[1])
   },numeric(1))
   spec(format(gamma_cols[k],trim=TRUE,scientific=FALSE), gamma_fill[k], vals)
  })
  draw_bars(NAMES,bars,DSN[[ds]],'Selected cases (%)',104)
 }
})

with_png('synthetic_low_hit',8.9,3.9,{
 par(mar=c(5,5,3.0,1))
 series <- lapply(seq_along(REPS),function(j) spec(NAMES[j],COLORS[j],
   vapply(c(80,200,400),function(d) lookup(cs,list(original_dim=d,representation=REPS[j]),'low_hit_pct'),numeric(1))))
 draw_bars(c('80','200','400'),series,main=FIGURE_TITLES[['synthetic_low_hit']],
           ylab='Low-hit mass levels (%)',xlabel='Ambient dimension d',ylim=106)
 par(xpd=NA)
 legend('top',inset=c(0,-.10),legend=NAMES,fill=COLORS,bty='n',horiz=TRUE,cex=.7)
 par(xpd=FALSE)
})

with_png('synthetic_pair_agreement',13.5,4.3,{
 layout_panels(1,3,FIGURE_TITLES[['synthetic_pair_agreement']],
  c('Raw','Reliable-only'),c('#8db0ca','#207BA2'))
 for(d in c(80,200,400)) {
  raw <- vapply(REPS,function(rep) 100*lookup(ps,list(original_dim=d,representation=rep),'mean_raw_agreement'),numeric(1))
  rel <- vapply(REPS,function(rep) 100*lookup(ps,list(original_dim=d,representation=rep),'mean_valid_agreement'),numeric(1))
  draw_bars(NAMES,list(spec('Raw','#8db0ca',raw),spec('Reliable-only','#207BA2',rel)),
    paste('d =',d),'Pair agreement (%)',104)
 }
})

with_png('intrinsic_dimension_agreement',13.5,4.3,{
 layout_panels(1,3,FIGURE_TITLES[['intrinsic_dimension_agreement']],
  c('Raw','Reliable-only'),c('#8db0ca','#207BA2'))
 for(q in c(3,5,10)) {
  raw <- vapply(REPS,function(rep) 100*lookup(q200,list(intrinsic_dim=q,representation=rep),'raw_agreement'),numeric(1))
  rel <- vapply(REPS,function(rep) 100*lookup(q200,list(intrinsic_dim=q,representation=rep),'reliable_agreement'),numeric(1))
  draw_bars(NAMES,list(spec('Raw','#8db0ca',raw),spec('Reliable-only','#207BA2',rel)),
    paste('intrinsic_dim =',q),'Pair agreement (%)',104)
 }
})

with_png('real_reliability_agreement',11.5,11.0,{
 layout_panels(3,2,FIGURE_TITLES[['real_reliability_agreement']])
 for(ds in DS[-1]) {
  reliable <- vapply(REPS,function(rep) lookup(cr,list(dataset=ds,representation=rep),'reliable_pct'),numeric(1))
  raw <- vapply(REPS,function(rep) 100*lookup(pr,list(dataset=ds,representation=rep),'mean_raw_agreement'),numeric(1))
  good <- vapply(REPS,function(rep) 100*lookup(pr,list(dataset=ds,representation=rep),'mean_valid_agreement'),numeric(1))
  draw_bars(NAMES,list(spec('Candidate reliability','#307A8A',reliable)),
    paste(DSN[[ds]],'- reliability'),'Reliable fits (%)',104)
  draw_bars(NAMES,list(spec('Raw','#8db0ca',raw),spec('Reliable-only','#207BA2',good)),
    paste(DSN[[ds]],'- pair agreement'),'Pair agreement (%)',118,legend_here=TRUE)
 }
})

with_png('real_runtime',10.5,4.0,{
 par(mar=c(5,5,2.8,1))
 series <- lapply(seq_along(REPS),function(j) spec(NAMES[j],COLORS[j],
    vapply(DS[-1],function(ds) lookup(cr,list(dataset=ds,representation=REPS[j]),'runtime'),numeric(1))))
 draw_bars(unname(DSN[DS[-1]]),series,main=FIGURE_TITLES[['real_runtime']],
           ylab='Logged seconds per candidate',fmt='%.2f')
 par(xpd=NA);legend('top',inset=c(0,-.08),legend=NAMES,fill=COLORS,horiz=TRUE,bty='n',cex=.7);par(xpd=FALSE)
})

with_png('synthetic_outlier_types',10.5,4.1,{
 layout_panels(1,2,FIGURE_TITLES[['synthetic_outlier_types']],
  c('Distributional','Structural'),c('#316FA5','#C98345'))
 for(truth in names(GEN)) {
  x <- vapply(REPS,function(rep) lookup(subtype,list(truth=truth,representation=rep),'roc_dist'),numeric(1))
  y <- vapply(REPS,function(rep) lookup(subtype,list(truth=truth,representation=rep),'roc_struct'),numeric(1))
  draw_bars(NAMES,list(spec('Distributional','#316FA5',x),spec('Structural','#C98345',y)),
    GEN[[truth]],'ROC-AUC',1.07,fmt='%.3f')
  abline(h=.5,lty=3,col='#777777')
 }
})

for(k in c('current_generator','orthogonal_shift')) {
 key <- if(k=='current_generator') 'curr' else 'ortho'
 scan <- if(key=='curr') curr else ortho
 fig_name <- paste0(k,'_dimension_scan')
 with_png(fig_name,10.3,3.9,{
  layout_panels(1,2,FIGURE_TITLES[[fig_name]],
    c('Distributional ROC','Structural ROC','Reliable fits'),c('#316FA5','#C98345','#8158AA'))
  for(truth in names(GEN)) {
    g <- scan[scan$truth==truth & num(scan$dimension_ratio)>0,,drop=FALSE]
    g <- g[order(num(g$dimension_ratio)),,drop=FALSE]
    if(!nrow(g)) {plot.new();title(main=GEN[[truth]]);next}
    x <- num(g$dimension_ratio)
    ys <- list(num(g$distributional_roc),num(g$structural_roc),num(g$reliable_pct)/100)
    draw_line(x,ys,c('#316FA5','#C98345','#8158AA'),
              c('Distributional ROC','Structural ROC','Reliable fits'),GEN[[truth]],
              xticks=c(1,2,4,8),xticklabels=c('q','2q','4q','8q'))
    title(xlab='MDS dimensions retained',cex.lab=.85)
  }
 })
}

mean_finite <- function(x) {x <- x[is.finite(x)];if(length(x)) mean(x) else NA_real_}

real_scan <- read_input('scan_real')
real_scan <- real_scan[real_scan$representation=='mds' & real_scan$gamma_scale==1,,drop=FALSE]
with_png('1f_real_mds_scan',12.4,4.1,{
  layout_panels(1,3,'Real-data MDS dimension scan',
    c('ROC-AUC','PR-AUC','Reliable fits'),c('#316FA5','#168676','#8158AA'))
  par(mar=c(4.4,4.6,2.3,.9))
  for(ds in DS[-1]) {
    z <- real_scan[real_scan$dataset==ds,,drop=FALSE];q <- sort(unique(z$embedding_dim))
    ys <- lapply(c('roc_auc','pr_auc','aumvc_reliable'),function(metric)
      vapply(q,function(d) mean_finite(z[z$embedding_dim==d,metric]),numeric(1)))
    draw_line(q,ys,c('#316FA5','#168676','#8158AA'),
      c('ROC-AUC','PR-AUC','Reliable fits'),DSN[[ds]],
      ylim=c(0,1.07),xticks=q,xticklabels=q)
    title(xlab='MDS dimension q',cex.lab=.85)

  }
})

with_png('gamma4_choice_only',11.2,7.4,{
 layout_panels(2,2,FIGURE_TITLES[['gamma4_choice_only']],
   c('0.015625','0.25','1','4','64'),
   c('#C9A27B','#B57743','#8F654D','#A66C3D','#C0946D'))
 gamma_cols <- c('0.015625','0.25','1','4','64')
 gamma_fill <- c('#C9A27B','#B57743','#8F654D','#A66C3D','#C0946D')
 for(ds in DS) {
  bars <- lapply(seq_along(gamma_cols), function(k) {
   vals <- vapply(REPS,function(rep) {
    row <- widened_screen[widened_screen$Dataset==DSN[[ds]] &
                          widened_screen$Representation==NAMES[match(rep,REPS)],,drop=FALSE]
    if(!nrow(row) || !(gamma_cols[k] %in% names(row))) return(NA_real_)
    num(row[[gamma_cols[k]]][1])
   },numeric(1))
   spec(gamma_cols[k], gamma_fill[k], vals)
  })
  draw_bars(NAMES,bars,DSN[[ds]],'Selected cases (%)',104)
 }
})

range_syn <- group_rows(WS,
  c('truth','intrinsic_dim','original_dim','snr','run','representation'),
  function(g) data.frame(dataset='synthetic',roc_range=diff(range(num(g$roc_auc),na.rm=TRUE))))

range_real <- group_rows(WR,c('dataset','run','representation'),
  function(g) data.frame(roc_range=diff(range(num(g$roc_auc),na.rm=TRUE))))

all_ranges <- rbind_fill(range_syn[,c('dataset','representation','roc_range')],
                          range_real[,c('dataset','representation','roc_range')])

range_summary <- group_rows(all_ranges,c('dataset','representation'),function(g)
  data.frame(n=nrow(g),median_roc_range=median_finite(g$roc_range),mean_roc_range=average(g$roc_range)))

with_png('widened_roc_variation',11.0,7.5,{
 layout_panels(2,2,FIGURE_TITLES[['widened_roc_variation']])
 for(ds in DS) {
  v <- vapply(REPS,function(rep) lookup(range_summary,list(dataset=ds,representation=rep),'median_roc_range'),numeric(1))
  good <- v[is.finite(v)]
  top <- max(.08,if(length(good)) max(good)*1.3 else .08)
  draw_bars(NAMES,list(spec('ROC range','#5187A7',v)),DSN[[ds]],
      'Within-run ROC-AUC range',top,fmt='%.3f')
 }
})

with_png('widened_roc_grid_edge',11.0,7.3,{
 layout_panels(2,2,FIGURE_TITLES[['widened_roc_grid_edge']])
 for(ds in DS) {
  v <- vapply(REPS,function(rep) {
    g <- diag[diag$dataset==ds & diag$representation==rep,,drop=FALSE]
    if(!nrow(g)) return(NA_real_)
    100*sum(num(g$roc_boundary_optima),na.rm=TRUE)/sum(num(g$n_cases),na.rm=TRUE)
  },numeric(1))
  draw_bars(NAMES,list(spec('ROC best at edge','#B7803F',v)),DSN[[ds]],
      'ROC optimum at grid endpoint (%)',105)
 }
})

with_png('gamma64_diagnostics',12.4,12.6,{
 layout_panels(4,2,FIGURE_TITLES[['gamma64_diagnostics']])
 for(ds in DS) {
  vals <- function(key) vapply(REPS,function(rep) lookup(extremes,list(dataset=ds,representation=rep),key),numeric(1))
  draw_bars(NAMES,list(spec('Numerically eligible','#8158AA',vals('old_reliable_pct')),
    spec('Reference nearly constant','#C76252',vals('reference_constant_pct')),
    spec('Evaluation nearly constant','#DAA34F',vals('eval_constant_pct'))),
    paste(DSN[[ds]],': score checks'),'Fits (%)',104)
  v1 <- vapply(REPS,function(rep) 100*lookup(mv64,list(dataset=ds,representation=rep),'full_box_rate'),numeric(1))
  v2 <- vapply(REPS,function(rep) 100*lookup(mv64,list(dataset=ds,representation=rep),'reliable_rate'),numeric(1))
  draw_bars(NAMES,list(spec('Full-box','#C76252',v1),spec('Valid after all checks','#258579',v2)),
    paste(DSN[[ds]],': full-box'),'Fits (%)',104)
 }
})

with_png('1g_regret_default_random',11.8,7.5,{
  layout_panels(2,2,'ROC regret: AUMVC, default gamma = 1, random gamma (mean +/- run-level SE)',
    c('AUMVC','Default gamma=1','Random gamma'),c('#2871a5','#e58b2a','#188675'))
  for(ds in DS){

    r <- regret[FALSE,,drop=FALSE]
    for(rep in REPS){
      p <- regret[regret$dataset==ds & regret$representation==rep,,drop=FALSE]
      if(nrow(p)!=1L) stop('Missing regret row for ',ds,'/',rep)
      r <- rbind(r,p)
    }
    means <- rbind(r$roc_regret_mean,r$default_regret_mean,r$random_regret_mean)
    errors <- rbind(r$roc_regret_se,r$default_regret_se,r$random_regret_se)
    top <- suppressWarnings(max(means + ifelse(is.finite(errors),errors,0),na.rm=TRUE))
    if(!is.finite(top)||top<=0) top <- .02
    bars <- list(spec('AUMVC','#2871a5',means[1,]),
                 spec('Default gamma=1','#e58b2a',means[2,]),
                 spec('Random gamma','#188675',means[3,]))
    draw_bars(NAMES,bars,DSN[[ds]],'ROC regret',top*1.35,fmt='%.3f')

    width <- min(.78/3,.24)
    for(j in seq_len(3)) for(i in seq_len(4)){
      if(!is.finite(means[j,i])||!is.finite(errors[j,i])) next
      x <- i+(j-2)*width
      arrows(x,max(0,means[j,i]-errors[j,i]),x,means[j,i]+errors[j,i],
             angle=90,code=3,length=.035,col='#333333',lwd=.8)
    }
  }
})

with_png('paired_bootstrap',10.4,4.1,{
 b <- BOOT[BOOT$representation!='ambient' & is.finite(num(BOOT$selected_gamma)),,drop=FALSE]
 b <- b[order(b$truth,b$representation),,drop=FALSE]
 if(!nrow(b)) {
   plot.new();text(.5,.5,'No eligible bootstrap cases')
 } else {
   labels <- vapply(seq_len(nrow(b)),function(i) paste0(
     GEN[[as.character(b$truth[i])]],' d=',as.integer(b$original_dim[i]),
     ' / ',NAMES[match(b$representation[i],REPS)]),character(1))
   obs <- num(b$observed_roc_gap); lo <- num(b$bootstrap_gap_low); hi <- num(b$bootstrap_gap_high)
   lim <- range(c(0,obs,lo,hi),na.rm=TRUE);xpad <- diff(lim)*.15 + .003
   par(mar=c(5.0,max(8.0,max(nchar(labels))*.15),1,1))
   plot(obs,rev(seq_len(nrow(b))),xlim=lim+c(-xpad,xpad),ylim=c(.5,nrow(b)+.5),
     yaxt='n',pch=19,col='#326FA6',xlab='ROC-best gamma AUC minus AUMVC-selected gamma AUC',
     ylab='',bty='n')
   axis(2,at=rev(seq_len(nrow(b))),labels=labels,las=1,cex.axis=.8)
   abline(v=0,col='#777777',lty=2)
   arrows(lo,rev(seq_len(nrow(b))),hi,rev(seq_len(nrow(b))),angle=90,code=3,
     length=.06,col='#326FA6',lwd=1.6)
   abline(v=pretty(lim),col='#eeeeee')
 }
})

for(ds in c('synthetic','ecg5000')) {
 fig_key <- paste0('gamma_selection_',ds)
 with_png(fig_key,10.9,7.4,{
   layout_panels(2,2,FIGURE_TITLES[[fig_key]],METH,unname(METHOD_COLORS[METH]),height=.60)
   for(rep in REPS) {
     specs <- lapply(METH,function(method){
       values <- vapply(GAMMAS,function(gamma)
         lookup(frequency,list(dataset=ds,representation=rep,
                               method=method,gamma_scale=gamma),'pct'),numeric(1))
       spec(method,METHOD_COLORS[[method]],values)
     })
     draw_bars(c('0.25','0.5','1','2','4'),specs,NAMES[match(rep,REPS)],
       'Selected cases (%)',105,xlabel='gamma_scale')
   }
 })
}

for(ds in DS[-1]) {
 dat <- R[R$dataset==ds,,drop=FALSE]
 gp_a <- run_mean_se(dat,'aumvc_normalized',c('representation','gamma_scale'))
 gp_r <- run_mean_se(dat,'roc_auc',c('representation','gamma_scale'))
 key <- paste0('gamma_profile_',ds)
 with_png(key,10.6,10.5,{
  layout_panels(4,2,FIGURE_TITLES[[key]])
  for(rep in REPS) for(which in c('AUMVC','ROC')) {
   g <- if(which=='AUMVC') gp_a else gp_r
   ytitle <- if(which=='AUMVC') 'AUMVC (minimize)' else 'ROC-AUC (maximize)'
   g <- g[g$representation==rep,,drop=FALSE]
   g <- g[order(num(g$gamma_scale)),,drop=FALSE]
   y <- num(g$mean); er <- num(g$se);er[!is.finite(er)] <- 0
   xr <- seq_len(nrow(g));ylo <- min(y-er,na.rm=TRUE);yhi <- max(y+er,na.rm=TRUE)
   padding <- max(.025,(yhi-ylo)*.12)
   plot(xr,y,type='b',pch=19,lwd=1.5,col='#26799b',xaxt='n',
     ylim=c(max(0,ylo-padding),yhi+padding),xlab='gamma_scale',
     ylab=if(which=='AUMVC') 'AUMVC' else 'ROC-AUC',bty='n')
   abline(h=pretty(c(ylo-padding,yhi+padding)),col='#eeeeee')
   axis(1,at=xr,labels=format(num(g$gamma_scale),trim=TRUE))
   arrows(xr,y-er,xr,y+er,angle=90,code=3,length=.03,col='#26799b')
   win <- if(which=='AUMVC') which.min(y) else which.max(y)
   points(xr[win],y[win],pch=1,cex=2,col='#c76b41',lwd=2)
   text(xr[win],y[win]+padding*.55,'best',cex=.64)
   title(main=paste(NAMES[match(rep,REPS)],':',ytitle),font.main=1,cex.main=.83)
  }
 })
}

diag_base <- file.path(ROOT,'interpretation','diagnostics')

mass_folder <- NULL

for(grid in c('widened','original')) {
 p <- file.path(diag_base,grid,'mass_profiles.csv')
 if(file.exists(p)) {mass_folder <- file.path(diag_base,grid);break}
}

if(!is.null(mass_folder)) {
 profile_path <- file.path(mass_folder,'mass_profiles.csv')
 profiles <- read_csv(profile_path)
 require_cols(profiles,c('dataset','representation','gamma_scale','mass',
  'threshold_source','threshold','occupancy','run'),profile_path)
 selected_cases <- list()
 for(case in list(c('shuttle','isomap'),c('ecg5000','mds'))) {
  ds <- case[1];rep <- case[2]
  p <- profiles[profiles$dataset==ds & profiles$representation==rep &
                  num(profiles$gamma_scale) %in% c(.25,4),,drop=FALSE]
  if(nrow(p)) selected_cases[[length(selected_cases)+1L]] <- list(ds=ds,rep=rep,data=p)
 }
 if(length(selected_cases)) {
  with_png('gamma_thresholds_volumes',10,3.45*length(selected_cases),{
    layout_panels(length(selected_cases),2,FIGURE_TITLES[['gamma_thresholds_volumes']])
    for(case in selected_cases) for(metric in c('threshold','occupancy')) {
      p <- case$data
      yy <- num(p[[metric]]);xx <- num(p$mass)
      range_y <- range(yy[is.finite(yy)],na.rm=TRUE)
      if(!all(is.finite(range_y)) || diff(range_y)==0) range_y <- c(0,1)
      plot(NA,NA,xlim=range(xx,na.rm=TRUE),ylim=range_y,xlab='Target mass',
        ylab=if(metric=='threshold') 'Normality threshold' else 'Volume / reference box',bty='n')
      abline(h=pretty(range_y),col='#eeeeee')
      for(gamma in c(.25,4)) for(src in c('heldout','training')) {
        z <- p[num(p$gamma_scale)==gamma & p$threshold_source==src,,drop=FALSE]
        if(!nrow(z)) next
        avg <- group_rows(z,'mass',function(g) data.frame(y=average(g[[metric]])))
        avg <- avg[order(num(avg$mass)),,drop=FALSE]
        lines(num(avg$mass),num(avg$y),col=if(gamma==.25) '#287c9f' else '#bd744d',
          lty=if(src=='heldout') 1 else 2,lwd=1.5)
      }
      title(main=paste(DSN[[case$ds]],NAMES[match(case$rep,REPS)],':',
                 if(metric=='threshold') 'Normality threshold' else 'Volume / reference box'),
            font.main=1,cex.main=.81)
      if(case$ds=='shuttle' && metric=='threshold') legend('topright',
        legend=c('gamma 0.25, heldout','gamma 0.25, training','gamma 4, heldout','gamma 4, training'),
        col=c('#287c9f','#287c9f','#bd744d','#bd744d'),lty=c(1,2,1,2),
        bty='n',cex=.65)
    }
  })
 }

 rank_path <- file.path(mass_folder,'score_ranks.csv')
 if(file.exists(rank_path)) {
  z <- read_csv(rank_path)
  ids <- c('dataset','case_id','representation','run','observation_id','label')
  needed <- c(ids,'gamma_scale','anomaly_score_percentile')
  if(all(needed %in% names(z))) {
   z <- z[z$dataset %in% c('shuttle','ecg5000') & z$representation %in% c('isomap','mds') &
        num(z$gamma_scale) %in% c(.25,4) & num(z$label)==1,,drop=FALSE]
   z <- group_rows(z,c(ids,'gamma_scale'),function(g)
         data.frame(rank=average(g$anomaly_score_percentile)))
   by_obs <- group_rows(z,ids,function(g){
     a <- g$rank[num(g$gamma_scale)==.25];b <- g$rank[num(g$gamma_scale)==4]
     data.frame(rank025=if(length(a)) average(a) else NA_real_,
       rank4=if(length(b)) average(b) else NA_real_)
   })
   by_obs <- by_obs[is.finite(by_obs$rank025) & is.finite(by_obs$rank4),,drop=FALSE]
   if(nrow(by_obs)) {
    byrun <- group_rows(by_obs,c('dataset','representation','run'),function(g)
      data.frame(rank025=average(g$rank025),rank4=average(g$rank4),
                 delta=average(g$rank4-g$rank025)))
    bygroup <- group_rows(byrun,c('dataset','representation'),function(g)
      data.frame(n=nrow(g),rank025=average(g$rank025),rank4=average(g$rank4),
        delta=average(g$delta),delta_se=se(g$delta)))
    rank_table <- rbind_fill(lapply(seq_len(nrow(bygroup)),function(i){
      g <- bygroup[i,,drop=FALSE]
      data.frame('Dataset'=DSN[[g$dataset]],'Representation'=NAMES[match(g$representation,REPS)],
        'Runs'=g$n,'Mean rank gamma 0.25'=sprintf('%.3f',g$rank025),
        'Mean rank gamma 4'=sprintf('%.3f',g$rank4),
        'Rank change ± SE'=if(is.finite(g$delta_se))
           sprintf('%+.3f ± %.3f',g$delta,g$delta_se) else sprintf('%+.3f (SE NA)',g$delta),
        check.names=FALSE)
    }))
    table_out('anomaly_score_rank_comparison',rank_table)
   }
  }
 }

}

required <- KEEP_FIGURES

missing <- setdiff(required,FIGURE_WRITTEN)

if(length(missing)) warning('Missing historical Python figures: ',paste(missing,collapse=', '))

missing_tables <- setdiff(KEEP_TABLES,TABLE_WRITTEN)

if(length(missing_tables)) warning('Missing historical tables (check optional diagnostics): ',
                                   paste(missing_tables,collapse=', '))
