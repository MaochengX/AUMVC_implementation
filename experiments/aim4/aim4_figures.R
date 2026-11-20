root <- if (dir.exists('experiments/aim4/reproduction')) 'experiments/aim4' else

        if (dir.exists('aim4/reproduction')) 'aim4' else

        stop('Run from the repository root containing experiments/aim4.',call.=FALSE)

fig_dir <- file.path(root,'interpretation','figures')

tab_dir <- file.path(root,'interpretation','tables')

dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)

dir.create(tab_dir,recursive=TRUE,showWarnings=FALSE)

DET <- c('ocsvm','lof'); DATASETS <- c('ecg5000','fashion_mnist')

DN <- c(ocsvm='OCSVM',lof='LOF'); SN <- c(ecg5000='ECG5000',fashion_mnist='Fashion-MNIST')

PAL <- c(AUMVC='#A66C3D',ROC='#2871A5',PR='#188675',Default='#E58B2A',
         Selected='#2871A5',Full_rule='#E58B2A')

read_file <- function(path) {
  if (!file.exists(path)) stop('Missing saved result: ',path,call.=FALSE)
  read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
}

need <- function(d,cols,where) {
  missing <- setdiff(cols,names(d))
  if(length(missing)) stop(where,': missing ',paste(missing,collapse=', '),call.=FALSE)
}

write_table <- function(name,d) {
  write.csv(d,file.path(tab_dir,paste0(name,'.csv')),row.names=FALSE,na='NA')
  invisible(d)
}

mean_ok <- function(x) {x<-x[is.finite(x)];if(length(x)) mean(x) else NA_real_}

se_ok <- function(x) {x<-x[is.finite(x)];if(length(x)>1L) sd(x)/sqrt(length(x)) else NA_real_}

check_close <- function(actual,expected,label,tol=2e-7) {
  if(length(actual)!=length(expected) || anyNA(actual) || anyNA(expected) ||
     any(!is.finite(actual)) || any(!is.finite(expected)) ||
     any(abs(actual-expected)>tol)) stop('Saved-results check failed: ',label,call.=FALSE)
}

validate_saved <- function(s,c,det,ds) {
  need(s,c('dataset','run','selection_status','n_reliable','aumvc_choice',
           'roc_label_choice','pr_label_choice','default','test_roc_aumvc',
           'test_pr_aumvc','roc_loss_vs_labels','pr_loss_vs_labels',
           'roc_gain_vs_default','pr_gain_vs_default'),paste(det,ds,'selections'))
  need(c,c('dataset','run','parameter','reliable','aumvc','validation_roc',
           'validation_pr','test_roc','test_pr'),paste(det,ds,'candidates'))
  if(nrow(s)<20L || anyDuplicated(s$run) || anyNA(s$run) ||
     any(s$dataset!=ds) || any(c$dataset!=ds) ||
     any(!s$selection_status %in% c('selected','tie','all_unreliable')))
    stop('Invalid or insufficient runs: ',det,'/',ds,call.=FALSE)
  params <- sort(unique(c$parameter))
  if(anyDuplicated(paste(c$run,c$parameter)) ||
     nrow(c)!=nrow(s)*length(params) ||
     any(as.integer(table(factor(c$run,levels=s$run)))!=length(params)) ||
     any(!c$run %in% s$run)) stop('Incomplete candidate grid: ',det,'/',ds)
  index <- function(par) match(paste(s$run,par),paste(c$run,c$parameter))
  ai <- index(s$aumvc_choice);ri <- index(s$roc_label_choice)
  pi <- index(s$pr_label_choice);di <- index(s$default)
  if(anyNA(c(ai,ri,pi,di))) stop('Missing chosen candidate: ',det,'/',ds)
  if(any(s$selection_status!='selected' & s$aumvc_choice!=s$default))
    stop('Incorrect no-choice fallback: ',det,'/',ds)
  if(any(s$n_reliable!=as.integer(tapply(c$reliable,c$run,sum)[as.character(s$run)])))
    stop('Reliability candidate counts disagree: ',det,'/',ds)
  check_close(s$test_roc_aumvc,c$test_roc[ai],'AUMVC test ROC')
  check_close(s$test_pr_aumvc,c$test_pr[ai],'AUMVC test PR')
  check_close(s$roc_loss_vs_labels,c$test_roc[ri]-c$test_roc[ai],'ROC paired gap')
  check_close(s$pr_loss_vs_labels,c$test_pr[pi]-c$test_pr[ai],'PR paired gap')
  check_close(s$roc_gain_vs_default,c$test_roc[ai]-c$test_roc[di],'ROC default gap')
  check_close(s$pr_gain_vs_default,c$test_pr[ai]-c$test_pr[di],'PR default gap')
  for(run in s$run) {
    g <- c[c$run==run,,drop=FALSE]
    z <- s[s$run==run,,drop=FALSE]

    if(g$validation_roc[match(z$roc_label_choice,g$parameter)] <
         max(g$validation_roc)-1e-12 ||
       g$validation_pr[match(z$pr_label_choice,g$parameter)] <
         max(g$validation_pr)-1e-12)
      stop('Label-selected parameter not validation-optimal: ',det,'/',ds,' run ',run)
    eligible <- which(g$reliable & is.finite(g$aumvc))
    if(z$selection_status=='selected' &&
       (!length(eligible) || !(match(z$aumvc_choice,g$parameter) %in% eligible) ||
        g$aumvc[match(z$aumvc_choice,g$parameter)] > min(g$aumvc[eligible])+1e-10))
      stop('AUMVC selection not eligible/optimal: ',det,'/',ds,' run ',run)
  }
  invisible(TRUE)
}

runs <- list(); candidates <- list()

for(det in DET) for(ds in DATASETS) {
  name <- paste(det,ds,sep='_')
  folder <- file.path(root,paste0(det,'_tune'),ds)
  s <- read_file(file.path(folder,'selections.csv'))
  c <- read_file(file.path(folder,'candidates.csv'))
  validate_saved(s,c,det,ds)
  s$detector <- det;s$dataset <- ds

  s$roc_delta_vs_roc <- -s$roc_loss_vs_labels
  s$pr_delta_vs_pr <- -s$pr_loss_vs_labels
  s$roc_delta_vs_default <- s$roc_gain_vs_default
  s$pr_delta_vs_default <- s$pr_gain_vs_default
  s$valid_choice <- s$selection_status=='selected'
  s$raw_aumvc_choice <- ifelse(s$valid_choice,s$aumvc_choice,NA_real_)
  s$more_local_than_pr <- if(det=='ocsvm') {
    ifelse(s$valid_choice,s$aumvc_choice>s$pr_label_choice,NA)
  } else {
    ifelse(s$valid_choice,s$aumvc_choice<s$pr_label_choice,NA)
  }
  c$detector <- det;c$dataset <- ds
  runs[[name]] <- s;candidates[[name]] <- c
}

summary_rows <- list();locality_rows <- list()

for(det in DET) for(ds in DATASETS) {
  s <- runs[[paste(det,ds,sep='_')]];v <- s[s$valid_choice,,drop=FALSE]
  locality_rows[[length(locality_rows)+1L]] <- data.frame(
    detector=det,dataset=ds,n_runs=nrow(s),n_valid=nrow(v),
    n_no_choice=sum(!s$valid_choice),no_choice_pct=100*mean(!s$valid_choice),
    n_all_unreliable=sum(s$selection_status=='all_unreliable'),
    n_ties=sum(s$selection_status=='tie'),
    aumvc_more_local_than_roc_pct=100*mean(v$more_local_than_roc),
    aumvc_more_local_than_pr_pct=100*mean(v$more_local_than_pr),
    aumvc_equals_roc_pct=100*mean(v$aumvc_choice==v$roc_label_choice),
    aumvc_equals_pr_pct=100*mean(v$aumvc_choice==v$pr_label_choice),
    mean_runtime_sec=mean_ok(s$runtime_sec))
  for(scope in c('selected_only','full_rule')) {
    z <- if(scope=='selected_only') v else s
    summary_rows[[length(summary_rows)+1L]] <- data.frame(
      detector=det,dataset=ds,scope=scope,n_runs=nrow(s),n_used=nrow(z),
      no_choice_n=sum(!s$valid_choice),no_choice_pct=100*mean(!s$valid_choice),
      roc_delta_vs_roc_mean=mean_ok(z$roc_delta_vs_roc),
      roc_delta_vs_roc_se=se_ok(z$roc_delta_vs_roc),
      pr_delta_vs_pr_mean=mean_ok(z$pr_delta_vs_pr),
      pr_delta_vs_pr_se=se_ok(z$pr_delta_vs_pr),
      roc_delta_vs_default_mean=mean_ok(z$roc_delta_vs_default),
      roc_delta_vs_default_se=se_ok(z$roc_delta_vs_default),
      pr_delta_vs_default_mean=mean_ok(z$pr_delta_vs_default),
      pr_delta_vs_default_se=se_ok(z$pr_delta_vs_default))
  }
}

sumtab <- do.call(rbind,summary_rows)

loc <- do.call(rbind,locality_rows)

for(det in DET) {
  reference <- read_file(file.path(root,paste0(det,'_tune'),'evaluation.csv'))
  for(ds in DATASETS) for(scope in c('selected_only','full_rule')) {
    saved <- reference[reference$dataset==ds & reference$scope==scope,,drop=FALSE]
    calculated <- sumtab[sumtab$detector==det & sumtab$dataset==ds &
                         sumtab$scope==scope,,drop=FALSE]
    if(nrow(saved)!=1L || nrow(calculated)!=1L ||
       saved$n_used!=calculated$n_used) stop('Evaluation summary scope mismatch')
    for(metric in c('roc','pr')) {
      check_close(calculated[[paste0(metric,'_delta_vs_',metric,'_mean')]],
                  -saved[[paste0(metric,'_loss_mean')]],'Signed summary mean')
      check_close(calculated[[paste0(metric,'_delta_vs_',metric,'_se')]],
                  saved[[paste0(metric,'_loss_se')]],'Signed summary SE')
      check_close(calculated[[paste0(metric,'_delta_vs_default_mean')]],
                  saved[[paste0(metric,'_gain_default_mean')]],'Default summary mean')
      check_close(calculated[[paste0(metric,'_delta_vs_default_se')]],
                  saved[[paste0(metric,'_gain_default_se')]],'Default summary SE')
    }
  }
}

locality_thesis <- data.frame(
  'Detector / dataset'=paste(DN[loc$detector],SN[loc$dataset],sep=' / '),
  'Valid choices'=sprintf('%s/%s',format(loc$n_valid,big.mark=',',scientific=FALSE,trim=TRUE),
                          format(loc$n_runs,big.mark=',',scientific=FALSE,trim=TRUE)),
  'No choice'=sprintf('%s (%.2f%%)',format(loc$n_no_choice,big.mark=',',scientific=FALSE,trim=TRUE),
                       loc$no_choice_pct),
  'More local vs ROC'=sprintf('%.1f%%',loc$aumvc_more_local_than_roc_pct),
  'More local vs PR'=sprintf('%.1f%%',loc$aumvc_more_local_than_pr_pct),
  check.names=FALSE,stringsAsFactors=FALSE)

stopifnot(nrow(locality_thesis)==length(DET)*length(DATASETS),
          all(loc$n_valid+loc$n_no_choice==loc$n_runs),
          all(is.finite(loc$aumvc_more_local_than_roc_pct)),
          all(is.finite(loc$aumvc_more_local_than_pr_pct)))

write_table('tuning_locality_summary',locality_thesis)

freq_list <- list()

for(det in DET) for(ds in DATASETS) {
  name <- paste(det,ds,sep='_')
  s <- runs[[name]];c <- candidates[[name]]
  params <- sort(unique(c$parameter))
  for(scope in c('selected_only')) {
    z <- s[s$valid_choice,,drop=FALSE]
    for(method in c('AUMVC','ROC','PR')) {
      choices <- z[[c(AUMVC='aumvc_choice',ROC='roc_label_choice',PR='pr_label_choice')[[method]]]]
      for(p in params) freq_list[[length(freq_list)+1L]] <- data.frame(
        detector=det,dataset=ds,scope=scope,method=method,parameter=p,
        n=sum(choices==p),denominator=nrow(z),pct=100*mean(choices==p))
    }
  }

}

freq <- do.call(rbind,freq_list)

boston <- read_file(file.path(root,'reproduction','boston','bandwidth.csv'))

gauss <- read_file(file.path(root,'reproduction','gaussian_mixture','bandwidth.csv'))

compare <- read_file(file.path(root,'reproduction','gaussian_mixture','set_comparison.csv'))

for(nm in c('sigma','mv_area','volume_95','valid_calibration','selected')) {
  need(boston,nm,'Boston');need(gauss,nm,'Gaussian mixture')
}

if(sum(boston$selected)!=1L || sum(gauss$selected)!=1L || nrow(compare)!=1L ||
   abs(gauss$sigma[gauss$selected]-compare$sigma)>1e-8 ||
   !isTRUE(boston$valid_calibration[boston$selected]) ||
   !isTRUE(gauss$valid_calibration[gauss$selected]) ||
   which.min(boston$mv_area)!=which(boston$selected) ||
   which.min(gauss$mv_area)!=which(gauss$selected))

  stop('Reproduction selected bandwidth/set comparison mismatch',call.=FALSE)

need(compare,c('target_mass','true_mass','true_mass_se',
 'symmetric_difference_volume_in_box','symmetric_difference_mc_se',
 'relative_symmetric_difference','reference_box_mass'),'Gaussian set comparison')

write_table('gaussian_set_thesis_table',data.frame(
  measure=c('Selected bandwidth (sigma)','Target probability mass',
            'True mass of selected set','Symmetric-difference volume (within box)',
            'Relative symmetric-difference volume'),
  estimate=c(sprintf('%.4f',compare$sigma),
             sprintf('%.2f%%',100*compare$target_mass),
             sprintf('%.2f%% +/- %.2f percentage points',100*compare$true_mass,
                     100*compare$true_mass_se),
             sprintf('%.3f +/- %.3f',compare$symmetric_difference_volume_in_box,
                     compare$symmetric_difference_mc_se),
             sprintf('%.2f%%',100*compare$relative_symmetric_difference)),
  check.names=FALSE
))

with_png <- function(name,w,h,code) {
  png(file.path(fig_dir,paste0(name,'.png')),width=w,height=h,units='in',res=190,bg='white')
  on.exit(dev.off(),add=TRUE)
  par(family='sans',fg='#151515',col.axis='#151515',col.lab='#151515',las=1,
      cex.axis=.84,cex.lab=.89)
  eval(substitute(code),envir=parent.frame())
  invisible(name)
}

header <- function(title,leg=NULL,cols=NULL) {
  par(mar=c(0,0,0,0));plot.new();plot.window(c(0,1),c(0,1))
  text(.5,.75,title,cex=1.04)
  if(!is.null(leg)) legend('center',legend=leg,fill=cols,horiz=TRUE,
     bty='n',cex=.81,x.intersp=.66,y.intersp=.65)
}

panels <- function(nrow,ncol,title,leg=NULL,cols=NULL,head_height=.65) {
  top <- matrix(1L,nrow=1,ncol=ncol)
  main <- matrix(seq_len(nrow*ncol)+1L,nrow=nrow,ncol=ncol,byrow=TRUE)
  layout(rbind(top,main),heights=c(head_height,rep(4,nrow)))
  header(title,leg,cols)
  par(mar=c(4.1,4.6,2.25,.9),las=1)
}

axis_grid <- function(ymin,ymax,at=NULL) {
  if(is.null(at)) at <- pretty(c(ymin,ymax),n=5)
  at<-at[at>=ymin & at<=ymax]
  abline(h=at,col='#ececec',lwd=.8)
  axis(2,at=at,labels=format(round(at,3),trim=TRUE),lwd=.7)
  abline(h=0,col='#707b83',lwd=.9)
}

frequency_panel <- function(det,ds) {
  g <- freq[freq$detector==det & freq$dataset==ds & freq$scope=='selected_only',,drop=FALSE]
  params <- sort(unique(g$parameter));m <- length(params)
  z <- sapply(c('AUMVC','ROC','PR'),function(method)
      g$pct[match(paste(method,params),paste(g$method,g$parameter))])
  ylim <- c(0,max(105, max(z,na.rm=TRUE)*1.11))
  plot(NA,NA,type='n',xlim=c(.35,m+.65),ylim=ylim,xaxt='n',yaxt='n',
       xlab=if(det=='ocsvm') 'Gamma scale (larger = more local)' else 'LOF k (smaller = more local)',
       ylab='Selection frequency (%)',bty='n',xaxs='i',yaxs='i')
  axis_grid(0,ylim[2],at=seq(0,100,20));axis(1,at=seq_len(m),labels=params,cex.axis=.84)
  for(j in seq_len(3L)) for(i in seq_len(m)) {
    val <- z[i,j];xx <- i+(j-2)*.245
    rect(xx-.11,0,xx+.11,val,col=PAL[[colnames(z)[j]]],border=NA)
  }
  n <- unique(g$denominator)
  title(main=paste0(SN[[ds]],'  (',n,' successful AUMVC runs)'),font.main=1,cex.main=.91)
}

for(det in DET) {
 with_png(paste0('parameter_choices_',det),13.1,7.6,{
   panels(2,1,paste(DN[[det]],': gamma / k choices on independent tuning and validation splits'),
     c('AUMVC','ROC validation','PR validation'),unname(PAL[c('AUMVC','ROC','PR')]))
   for(ds in DATASETS) frequency_panel(det,ds)
 })
}

delta_cols <- c('roc_delta_vs_roc','pr_delta_vs_pr',
                'roc_delta_vs_default','pr_delta_vs_default')

delta_labs <- c('ROC vs ROC-tuned','PR vs PR-tuned',
                'ROC vs default','PR vs default')

signed_panel <- function(det,ds) {
  z <- sumtab[sumtab$detector==det & sumtab$dataset==ds,,drop=FALSE]
  z <- z[match(c('selected_only','full_rule'),z$scope),,drop=FALSE]
  if(nrow(z)!=2L || anyNA(z$scope)) stop('Missing signed summary scope')
  v <- sapply(delta_cols,function(key) z[[paste0(key,'_mean')]])
  err <- sapply(delta_cols,function(key) z[[paste0(key,'_se')]])
  if(!all(dim(v)==c(2L,4L)) || any(!is.finite(v)))
    stop('Missing paired difference summary for ',det,'/',ds)
  span <- max(.02,abs(c(v,err)),na.rm=TRUE)
  low <- min(0,v-ifelse(is.finite(err),err,0))-.12*span
  high <- max(0,v+ifelse(is.finite(err),err,0))

  right <- max(high,0)+.90*span
  plot(NA,NA,type='n',xlim=c(low,right),ylim=c(.45,4.65),
       xaxt='n',yaxt='n',xlab='Test AUC difference (AUMVC - comparator)',
       ylab='',bty='n',xaxs='i',yaxs='i')
  at <- pretty(c(low,max(0,high)),n=4)
  at <- at[at>=low & at<=max(0,high)]
  abline(v=at,col='#ededed',lwd=.8)
  abline(v=0,col='#76818b',lty=2,lwd=1.1)
  abline(h=c(1.5,2.5,3.5),col='#f0f0f0',lwd=.8)
  axis(1,at=at,labels=format(round(at,3),trim=TRUE),lwd=0,cex.axis=.78)
  axis(2,at=4:1,labels=delta_labs,las=1,tick=FALSE,cex.axis=.83)
  text(.055*span,4.54,'Mean +/- SE',adj=c(0,0.5),
       col='#56616c',font=2,cex=.76)
  for(i in seq_len(4)) for(j in seq_len(2)) {
    y <- 5-i + c(.16,-.16)[j]
    value <- v[j,i];e <- err[j,i]
    color <- PAL[[c('Selected','Full_rule')[j]]]
    if(is.finite(e)) {
      segments(value-e,y,value+e,y,col=color,lwd=2.2)
      segments(value-e,y-.060,value-e,y+.060,col=color,lwd=1.2)
      segments(value+e,y-.060,value+e,y+.060,col=color,lwd=1.2)
    }
    points(value,y,pch=16,col=color,cex=1.15)
    label <- if(is.finite(e)) sprintf('%+.3f +/- %.3f',value,e) else
                             sprintf('%+.3f (SE NA)',value)
    text(.055*span,y,label,adj=c(0,.5),col=color,cex=.74)
  }
  title(main=sprintf('%s / %s  (%d/%d valid AUMVC choices)',
                     DN[[det]],SN[[ds]],z$n_used[1],z$n_runs[1]),
        cex.main=.91,font.main=1,line=.8)
}

with_png('tuning_paired_test_differences',14.2,9.0,{
  panels(2,2,'Test performance: AUMVC minus comparator (mean +/- run-level SE)',
         c('Valid AUMVC choices','All runs, including default fallback'),
         unname(PAL[c('Selected','Full_rule')]),head_height=.73)
  for(det in DET) for(ds in DATASETS) {
    par(mar=c(4.3,7.25,2.4,1.0),las=1)
    signed_panel(det,ds)
  }
})

with_png('tuning_locality_summary',13.4,2.5,{
  par(mar=c(0,0,0,0));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
  xs <- c(.025,.25,.40,.57,.79,.975)
  ys <- seq(.93,.07,length.out=6L)
  rect(xs[1],ys[2],xs[6],ys[1],col='#eef2f5',border=NA)
  for(i in c(2L,4L)) rect(xs[1],ys[i+1L],xs[6],ys[i],col='#f8fafb',border=NA)
  for(x in xs) segments(x,ys[6],x,ys[1],col='#b8c2ca',lwd=.9)
  for(y in ys) segments(xs[1],y,xs[6],y,col='#b8c2ca',lwd=.9)
  for(j in seq_len(ncol(locality_thesis))) {
    cx <- (xs[j]+xs[j+1L])/2
    text(if(j==1) xs[j]+.01 else cx,(ys[1]+ys[2])/2,
         names(locality_thesis)[j],adj=if(j==1) c(0,.5) else c(.5,.5),
         cex=.87,font=2,col='#26323b')
    for(i in seq_len(nrow(locality_thesis))) {
      text(if(j==1) xs[j]+.01 else cx,(ys[i+1L]+ys[i+2L])/2,
           locality_thesis[[j]][i],adj=if(j==1) c(0,.5) else c(.5,.5),
           cex=.84,col='#26323b')
    }
  }
})

with_png('reproduction_mv_bandwidth',12.5,5.4,{
 panels(1,2,'Paper reproduction: calibrated OCSVM mass-volume area by bandwidth')
 for(study in c('boston','gaussian_mixture')) {
   d <- if(study=='boston') boston else gauss
   d <- d[order(d$sigma),,drop=FALSE]
   top3 <- d[which(as.logical(d$valid_calibration)),,drop=FALSE]
   top3 <- head(top3[order(top3$mv_area,top3$sigma),,drop=FALSE],3L)
   stopifnot(nrow(top3)==3L, all(is.finite(top3$mv_area)))
   rr <- range(d$mv_area);pad <- max(.04,diff(rr)*.15)
   plot(d$sigma,d$mv_area,type='n',xlab=expression('Kernel bandwidth ' * sigma),
        ylab='Area under MV curve',ylim=c(rr[1]-pad*.36,rr[2]+pad*.95),bty='n')
   abline(h=pretty(c(rr[1]-pad,rr[2]+pad)),col='#ececec')
   lines(d$sigma,d$mv_area,col=PAL[['ROC']],lwd=1.8)
   points(d$sigma,d$mv_area,pch=16,col=ifelse(d$valid_calibration,PAL[['ROC']],'#a9adb2'),cex=.60)
   points(d$sigma[d$selected],d$mv_area[d$selected],pch=21,bg=PAL[['Default']],
          col='#353f49',lwd=1,cex=1.55)
   abline(v=d$sigma[d$selected],lty=2,col=PAL[['Default']],lwd=1.25)
   if(study=='boston') abline(v=.42,lty=3,col='#555b66',lwd=1.2)
   mark_cols <- c('#C98542','#366D99','#168676')
   for(k in seq_len(3L)) {
     points(top3$sigma[k],top3$mv_area[k],pch=21,bg='white',
            col=mark_cols[k],lwd=2,cex=1.55)
     text(top3$sigma[k],top3$mv_area[k],labels=k,
          col=mark_cols[k],font=2,cex=.67)
   }
   title(main=paste0(if(study=='boston') 'Boston Housing' else '2D Gaussian mixture',
                     '  (chosen ',sprintf('%.4f',d$sigma[d$selected]),')'),
         font.main=1,cex.main=.91)
   text_rows <- sprintf('%d.  sigma %.4f ; MV %.4f',seq_len(3L),top3$sigma,top3$mv_area)
   legend('topright',legend=text_rows,title='Three lowest valid MV areas',
          text.col=mark_cols,bty='o',box.col='#d8e0e5',bg='white',cex=.78,
          x.intersp=.25,inset=.014)
   if(study=='boston') {
     u <- par('usr')
     text(u[1]+.04*diff(u[1:2]),u[3]+.83*diff(u[3:4]),
          'Paper reference: sigma ~ 0.42 (dotted)',adj=c(0,1),
          cex=.68,col='#677781')
   }
 }
})

if(Sys.getenv('AIM4_SET_OVERLAY','1')!='0') {
  model_source <- file.path(root,'reproduction.R')
  dependencies <- c('detectors/ocsvm.R','aumvc/input_validation.R','aumvc/level_set.R')
  if(file.exists(model_source) && all(file.exists(dependencies))) {
    source(model_source,local=TRUE)
    setup <- REPRODUCTION_SETTINGS
    set.seed(setup$seed+10L)
    x <- reproduction_gaussian_mixture(setup$gaussian_mixture_n)
    splits <- reproduction_splits(nrow(x),setup$gaussian_mixture_models,
                                  setup$train_fraction,setup$seed+11L)

    masses <- sort(unique(c(setup$mass_grid,.90,.95)))
    fit <- reproduction_fit(x,compare$sigma,splits,
       matrix(colMeans(x),nrow=1L),masses,setup$nu,keep_models=TRUE,
       solver_tolerance=setup$solver_tolerance)
    set.seed(setup$seed+13L)
    truth_draw <- reproduction_gaussian_mixture(setup$truth_n)
    density_cutoff <- unname(quantile(reproduction_gaussian_mixture_density(truth_draw),
                                      probs=1-setup$alpha))
    box <- fit_reference_box(x)
    box_low <- box$lower;box_high<-box$upper
    grid_n<-140L
    xs<-seq(box_low[1],box_high[1],length.out=grid_n)
    ys<-seq(box_low[2],box_high[2],length.out=grid_n)
    xy<-as.matrix(expand.grid(xs,ys))
    true_density<-reproduction_gaussian_mixture_density(xy)
    score<-numeric(nrow(xy))
    for(start in seq.int(1L,nrow(xy),by=1000L)) {
      idx<-start:min(start+999L,nrow(xy))
      score[idx]<-reproduction_predict(fit,xy[idx,,drop=FALSE])
    }
    cutoff<-fit$cutoff[match(.95,fit$masses)]
    true<-true_density>=density_cutoff;estimated<-score>=cutoff
    type<-matrix(as.integer(true)+2L*as.integer(estimated),nrow=grid_n,ncol=grid_n)
    with_png('gaussian_true_vs_selected_set',8.6,7.5,{
      par(mar=c(4.5,4.7,3.7,1.0))
      image(xs,ys,type,col=c('#ffffff','#EEC2AF','#BDD6EA','#C9E4D5'),
            breaks=c(-.5,.5,1.5,2.5,3.5),xlab='x1',ylab='x2',
            main='Gaussian mixture: true versus selected 95% MV set',asp=1,
            useRaster=TRUE)
      contour(xs,ys,matrix(as.numeric(true_density),grid_n,grid_n),
              levels=density_cutoff,add=TRUE,drawlabels=FALSE,
              col='#B96048',lwd=1.6)
      contour(xs,ys,matrix(score,grid_n,grid_n),levels=cutoff,
              add=TRUE,drawlabels=FALSE,col=PAL[['ROC']],lwd=1.6,lty=2)
      legend('topleft',c('True only','Estimated only','Overlap'),
             fill=c('#EEC2AF','#BDD6EA','#C9E4D5'),
             border=NA,bty='n',cex=.80)
      legend('bottomright',c('True boundary','Selected boundary'),
             lty=c(1,2),lwd=1.6,col=c('#B96048',PAL[['ROC']]),
             bty='n',cex=.80)

    })
  } else {
    warning('Gaussian set overlay skipped: original detector/AUMVC code is missing. ',
            'Tables still report the archived true mass and symmetric difference. ',
            'Run from repository root to reconstruct the selected set.')
  }
}
