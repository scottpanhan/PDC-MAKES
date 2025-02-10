library("energy")
library("MFSIS")
library("doParallel")
library("knockoff")
##########################################################################################################
##
#       function for PDC-Screen
#
# Input: - y: The response vector or a data matrix.
#        - x: The design matrix of dimensions n * p. 
#        - subset: if subset=FALSE, then all covariates except for X_i will be used as conditional 
#variables, if subset=TRUE, then a bootstrap based method will be used for subset selection for
#conditional variables.
#        - quant_boot: a quantile level used in bootstrap based subset selection method.
#        - dn: number of predictors recruited by PDC-Screen.  
#
# Output: - select: the labels of first dn largest active set of all predictors
#         - Time: the computatioanl times of PDC-Screen.
#
##########################################################################################################
PDC_SIS<-function(y,x,subset=FALSE,quant_boot,dn){
  if(subset==FALSE){
    p=ncol(x)
    n=nrow(x)
    start_time<-Sys.time()
    PDC=c()
    for (i in 1:p){
      PDC[i]=(pdcor(y,x[,i],x[,-i]))^2
    }
    PDCC=order(PDC,decreasing = T)[1:dn]
    end_time<-Sys.time()
    times=end_time-start_time
  }
  if(subset==TRUE){
  y=as.matrix(y)
  n=nrow(y)
  p_y=ncol(y)
  p=ncol(x)
  y_boot=matrix(NA,nrow = n,ncol = p_y)
  for (i in 1:p_y) {
    y_boot[,i]=sample(y[,i],n,replace = T)
  }
  PD_sel1=c()
  PD_sel2=c()
  for (i in 1:p){
    PD_sel1[i]=(energy::dcor(y_boot,x[,i]))^2
    PD_sel2[i]=(energy::dcor(y,x[,i]))^2
  }
  PD_sel=which(PD_sel2>=quantile(PD_sel1,quant_boot))
  
  start_time<-Sys.time()
  ynew=y
  xnew=x
  
  PDC=c()
  for (i in 1:p){
    index_set=ifelse(length(which(PD_sel==i))==0,PD_sel,PD_sel[-which(PD_sel==i)])
    PDC[i]=(pdcor(ynew,xnew[,i],xnew[,index_set]))^2
  }
  PDCC=order(PDC,decreasing = T)[1:dn]
  end_time<-Sys.time()
  times=end_time-start_time
  }
  return(list(select=PDCC,Time=times))
}

##########################################################################################################
##
#       function for rPDC-Screen
#
# Input: - y: The response vector or a data matrix.
#        - x: The design matrix of dimensions n * p. 
#        - subset: if subset=FALSE, then all covariates except for X_i will be used as conditional 
#variables, if subset=TRUE, then a bootstrap based method will be used for subset selection for
#conditional variables.
#        - quant_boot: a quantile level used in bootstrap based subset selection method.
#        - dn: number of predictors recruited by PDC-Screen.  
#
# Output: - select: the labels of first dn largest active set of all predictors
#         - Time: the computatioanl times of rPDC-Screen.
#
##########################################################################################################
rPDC_SIS<-function(y,x,subset=FALSE,quant_boot,dn){
  if(subset==FALSE){
    y=as.matrix(y)
    p=ncol(x)
    n=nrow(x)
    start_time<-Sys.time()
    ynew=matrix(NA,nrow = n,ncol = ncol(y))
    for (i in 1:ncol(y)) {
      ynew[,i]=sapply(1:n,function(k) (length(which(y[,i]<=y[k,i]))/n))
    }
    xnew=matrix(NA,nrow = n,ncol = p)
    for (i in 1:p) {
      xnew[,i]=sapply(1:n,function(k) (length(which(x[,i]<=x[k,i]))/n))
    }
    PDC=c()
    for (i in 1:p){
      PDC[i]=(pdcor(ynew,xnew[,i],xnew[,-i]))^2
    }
    PDCC=order(PDC,decreasing = T)[1:dn]
    end_time<-Sys.time()
    times=end_time-start_time
  }
  if(subset==TRUE){
    y=as.matrix(y)
    n=nrow(y)
    p_y=ncol(y)
    p=ncol(x)
    y_boot=matrix(NA,nrow = n,ncol = p_y)
    for (i in 1:p_y) {
      y_boot[,i]=sample(y[,i],n,replace = T)
    }
    PD_sel1=c()
    PD_sel2=c()
    for (i in 1:p){
      PD_sel1[i]=(energy::dcor(y_boot,x[,i]))^2
      PD_sel2[i]=(energy::dcor(y,x[,i]))^2
    }
    PD_sel=which(PD_sel2>=quantile(PD_sel1,quant_boot))
    
    start_time<-Sys.time()
    ynew=matrix(NA,nrow = n,ncol = p_y)
    for (i in 1:p_y) {
      ynew[,i]=sapply(1:n,function(k) (length(which(y[,i]<=y[k,i]))/n))
    }
    xnew=matrix(NA,nrow = n,ncol = p)
    for (i in 1:p) {
      xnew[,i]=sapply(1:n,function(k) (length(which(x[,i]<=x[k,i]))/n))
    }
    
    PDC=c()
    for (i in 1:p){
      index_set=ifelse(length(which(PD_sel==i))==0,PD_sel,PD_sel[-which(PD_sel==i)])
      PDC[i]=(pdcor(ynew,xnew[,i],xnew[,index_set]))^2
    }
    PDCC=order(PDC,decreasing = T)[1:dn]
    end_time<-Sys.time()
    times=end_time-start_time
  }
  return(list(select=PDCC,Time=times))
}

##########################################################################################################
##
#       function for PC-Screen
#
# Input: - x: The design matrix of dimensions n * p. 
#        - y: The response vector or a data matrix.
#        - multi: under multi-repsonse settings, we choose multi=TRUE, otherwise, multi=FALSE.
#        - dn: number of predictors recruited by PDC-Screen.  
#
# Output: - select: the labels of first dn largest active set of all predictors
#         - Time: the computatioanl times of rPDC-Screen.
#
##########################################################################################################
PC_SIS<-function(x,y,multi=FALSE,dn){
  if (multi==FALSE) {
    start_time<-Sys.time()
    PCC=PCSIS(x,y,nsis = dn)
    end_time<-Sys.time()
    times=end_time-start_time
  }
  if(multi==TRUE){
    start_time<-Sys.time()
    Arccos<-function (x) 
    {
      py_path = system.file("python", "PCSIS.py", package = "MFSIS")
      reticulate::source_python(py_path, envir = globalenv())
      ax = get_arccos_1d(x)
      return(ax)
    }
    Arccos_md<-function (x) 
    {
      py_path = system.file("python", "PCSIS.py", package = "MFSIS")
      reticulate::source_python(py_path, envir = globalenv())
      ax = get_arccos_md(x)
      return(ax)
    }
    pcsis<-function (X, Y, nsis = (dim(X)[1])/log(dim(X)[1])) 
    {
      if (dim(X)[1] != dim(Y)[1]) {
        stop("X and Y should have same number of rows!")
      }
      if (missing(X) | missing(Y)) {
        stop("The data is missing!")
      }
      # if (TRUE %in% (is.na(X) | is.na(Y) | is.na(nsis))) {
      #   stop("The input vector or matrix cannot have NA!")
      # }
      if (inherits(Y, "Surv")) {
        stop("PCSIS can not implemented with object  of Surv")
      }
      req_py()
      reticulate::py_config()
      py_path = system.file("python", "PCSIS.py", package = "MFSIS")
      reticulate::source_python(py_path, envir = globalenv())
      n = dim(X)[1]
      p = dim(X)[2]
      A_y = get_arccos_md(Y)
      Cor<-function (Xj, A_y, n) 
      {
        A_x = get_arccos(Xj)
        corr = projection_corr(A_x, A_y, n)
        return(corr)
      }
      if (n * p <= 1e+05) {
        result = vector(mode = "numeric", length = p)
        for (j in 1:p) {
          result[j] = Cor(X[, j], A_y, n)
        }
      }
      else {
        cores = detectCores(logical = FALSE)
        cl = makeCluster(cores)
        registerDoParallel(cl, cores = cores)
        j = NULL
        result = foreach::foreach(j = 1:p, .combine = "c", .export = c("get_arccos", 
                                                                       "projection_corr"), .packages = c("reticulate")) %dopar% 
          Cor(X[, j], A_y, n)
        stopImplicitCluster()
        stopCluster(cl)
      }
      A = order(result, decreasing = TRUE)
      return(A[1:nsis])
    }
    PCC=pcsis(x,y,nsis = p)
    end_time<-Sys.time()
    times=end_time-start_time
  }
  return(list(select=PCC,Time=times))
}
##########################################################################################################
##
#       function for the determination of the knockoff threshold
#
# Input: - W: the knockoff statistics. 
#        - fdr: The predetermined fdr level.
#        - offset: 1 for knockoff+ procedure and 0 for knockoff procedure.
#
# Output: - The threshold of koockoff method.
#
##########################################################################################################
alphakn_threshold <- function(W, fdr, offset) {
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}

##########################################################################################################
##
#       function for the determination of the threshold (9) in our main text.
#
# Input: - W: the knockoff statistics. 
#        - gamma: The predetermined fdr level.
#        - offset: 1 for knockoff+ procedure and 0 for knockoff procedure.
#
# Output: - The threshold (9) in our main text.
#
##########################################################################################################
stop_early <- function(W, gamma, offset){
  
  tau <- alphakn_threshold(W, fdr =  gamma, offset = offset) 
  ord_W <- order(abs(W), decreasing = TRUE)
  sorted_W <- W[ord_W]
  
  if(sum(W>0) >= 1 / gamma){
    pos_ind <- which(sorted_W > 0)
    tau1 <- sorted_W[pos_ind[ceiling(1/gamma)-1]]
  }else{
    tau1 <- 0
  }
  tau <- min(tau,tau1) 
  
  return(tau)
}

##########################################################################################################
##
#       function for the determination of the rejection set of e-BH procedure.
#
# Input: - E: the e-values. 
#        - alpha: The predetermined fdr level.
#
# Output: - The rejection set of e-BH procedure.
#
##########################################################################################################
ebh <- function(E, alpha){
  
  p <- length(E)
  E_ord <- order(E, decreasing = TRUE)
  E <- sort(E, decreasing = TRUE)
  comp <- E >= (p / alpha / (1:p))
  id <- suppressWarnings(max(which(comp>0)))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }
  return(list(rej = rej))
}

##########################################################################################################
##
#       function for PDC-MAKES
#
# Input: - y: The response vector or a data matrix.
#        - x: The design matrix of dimensions n * p. 
#        - subset: if subset=FALSE, then all covariates except for X_i will be used as conditional 
#variables, if subset=TRUE, then a bootstrap based method will be used for subset selection for
#conditional variables.
#        - rep: the repeatition times of the PDC-MAKES procedure.
#        - quant_boot: a quantile level used in bootstrap based subset selection method.
#        - knock_method: either "equi", "sdp" or "asdp" (default: "asdp"). This determines the method 
#that will be used to minimize the correlation between the original variables and the knockoffs.
#        - gam: The predetermined fdr level for knockoff procedure.
#        - alpha: The final predetermined fdr level.
#        - 1 for knockoff+ procedure and 0 for knockoff procedure.
#        - signal: The true active set.
#
# Output: - select: the selected features.
#         - fdr: the estimated fdr.
#         - power: the estimated power.
#         - Time: the computatioanl times of PDC-MAKES.
#
##########################################################################################################

PDC_MAKES<-function(y,x,subset,rep,quant_boot,knock_method,gam,alpha,offset,signal){
  if(subset==FALSE){
    n=nrow(x)
    p=ncol(x)
    start_time<-Sys.time()
    PDC=c()
    for (i in 1:p){
      PDC[i]=(pdcor(y,x[,i],x[,-i]))^2
    }
    scr_con1=PDC
    
    times=rep
    E=matrix(NA,nrow = times,ncol=p)
    for (tim in 1:times) {
      x_knn=x
      x_kn=create.second_order(x_knn,method = knock_method) 
      
      ynew=y
      xnew=x_kn
      
      PDC=c()
      for (i in 1:p){
        PDC[i]=(pdcor(ynew,xnew[,i],xnew[,-i]))^2
      }
      scr_con2=PDC
      
      
      w_kn=scr_con1-scr_con2 
      thr=stop_early(W=w_kn,gamma = gam,offset)
      
      
      E[tim,] <- (w_kn >= thr) / (1 + sum(w_kn <= -thr))
    }
    eBH <- p * colMeans(E)
    select=ebh(eBH,alpha)$rej
    end_time<-Sys.time()
    runtime=end_time-start_time

    if(length(select)==0){
      fdr=0
      pow=0
    }
    if(length(select)!=0){
      fdr=length(setdiff(select,signal))/length(select)
      pow=length(intersect(signal,select))/length(signal)
    }
  }
  
  if(subset==TRUE){
    y=as.matrix(y)
    n=nrow(x)
    p=ncol(x)
    p_y=ncol(y)
    start_time<-Sys.time()
    y_boot=matrix(NA,nrow = n,ncol = p_y)
    for (i in 1:p_y) {
      y_boot[,i]=sample(y[,i],n,replace = T)
    }
    PD_sel1=c()
    PD_sel2=c()
    for (i in 1:p){
      PD_sel1[i]=(energy::dcor(y_boot,x[,i]))^2
      PD_sel2[i]=(energy::dcor(y,x[,i]))^2
    }
    
    PD_sel=which(PD_sel2>=quantile(PD_sel1,quant_boot))
  
    
    ynew=y
    xnew=x
    
    PDC=c()
    for (i in 1:p){
      index_set=ifelse(length(which(PD_sel==i))==0,PD_sel,PD_sel[-which(PD_sel==i)])
      if(is.na(index_set)){
        PDC[i]=(pdcor(ynew,xnew[,i],xnew[,-i]))^2}
      if(!is.na(index_set)){
        PDC[i]=(pdcor(ynew,xnew[,i],xnew[,index_set]))^2}
    }
    scr_con1=PDC
    
    times=rep
    E=matrix(NA,nrow = times,ncol=p)
    for (tim in 1:times) {
      x_knn=x
      x_kn=create.second_order(x_knn,method = knock_method) 
      
      y_boot=matrix(NA,nrow = n,ncol = p_y)
      for (i in 1:p_y) {
        y_boot[,i]=sample(y[,i],n,replace = T)
      }
      PD_sel1=c()
      PD_sel2=c()
      for (i in 1:p){
        PD_sel1[i]=(energy::dcor(y_boot,x_kn[,i]))^2
        PD_sel2[i]=(energy::dcor(y,x_kn[,i]))^2
      }
      
      PD_sel=which(PD_sel2>=quantile(PD_sel1,quant_boot))
      
      ynew=y
      xnew=x_kn
      
      PDC=c()
      for (i in 1:p){
        index_set=ifelse(length(which(PD_sel==i))==0,PD_sel,PD_sel[-which(PD_sel==i)])
        if(is.na(index_set)){
          PDC[i]=(pdcor(ynew,xnew[,i],xnew[,-i]))^2}
        if(!is.na(index_set)){
          PDC[i]=(pdcor(ynew,xnew[,i],xnew[,index_set]))^2}
      }
      scr_con2=PDC
      
      
      w_kn=scr_con1-scr_con2 
      thr=stop_early(W=w_kn,gamma = gam,offset)
      
      
      E[tim,] <- (w_kn >= thr) / (1 + sum(w_kn <= -thr))
    }
    eBH <- p * colMeans(E)
    select=ebh(eBH,alpha)$rej
    end_time<-Sys.time()
    runtime=end_time-start_time
    
    if(length(select)==0){
      fdr=0
      pow=0
    }
    if(length(select)!=0){
      fdr=length(setdiff(select,signal))/length(select)
      pow=length(intersect(signal,select))/length(signal)
    }
  }
  return(list(select=select,fdr=fdr,power=pow,Time=runtime))
}

##########################################################################################################
##
#       function for rPDC-MAKES
#
# Input: - y: The response vector or a data matrix.
#        - x: The design matrix of dimensions n * p. 
#        - subset: if subset=FALSE, then all covariates except for X_i will be used as conditional 
#variables, if subset=TRUE, then a bootstrap based method will be used for subset selection for
#conditional variables.
#        - rep: the repeatition times of the rPDC-MAKES procedure.
#        - quant_boot: a quantile level used in bootstrap based subset selection method.
#        - knock_method: either "equi", "sdp" or "asdp" (default: "asdp"). This determines the method 
#that will be used to minimize the correlation between the original variables and the knockoffs.
#        - gam: The predetermined fdr level for knockoff procedure.
#        - alpha: The final predetermined fdr level.
#        - 1 for knockoff+ procedure and 0 for knockoff procedure.
#        - signal: The true active set.
#
# Output: - select: the selected features.
#         - fdr: the estimated fdr.
#         - power: the estimated power.
#         - Time: the computatioanl times of rPDC-MAKES.
#
##########################################################################################################

rPDC_MAKES<-function(y,x,subset,rep,quant_boot,knock_method,gam,alpha,offset,signal){
  if(subset==FALSE){
    y=as.matrix(y)
    n=nrow(x)
    p=ncol(x)
    start_time<-Sys.time()
    ynew=matrix(NA,nrow = n,ncol = ncol(y))
    for (i in 1:ncol(y)) {
      ynew[,i]=sapply(1:n,function(k) (length(which(y[,i]<=y[k,i]))/n))
    }
    xnew=matrix(NA,nrow = n,ncol = p)
    for (i in 1:p) {
      xnew[,i]=sapply(1:n,function(k) (length(which(x[,i]<=x[k,i]))/n))
    }
    PDC=c()
    for (i in 1:p){
      PDC[i]=(pdcor(ynew,xnew[,i],xnew[,-i]))^2
    }
    scr_con1=PDC
    
    times=rep
    E=matrix(NA,nrow = times,ncol=p)
    for (tim in 1:times) {
      x_knn=x
      x_kn=create.second_order(x_knn,method = knock_method) 
      
      ynew=matrix(NA,nrow = n,ncol = ncol(y))
      for (i in 1:ncol(y)) {
        ynew[,i]=sapply(1:n,function(k) (length(which(y[,i]<=y[k,i]))/n))
      }
      xnew=matrix(NA,nrow = n,ncol = p)
      for (i in 1:p) {
        xnew[,i]=sapply(1:n,function(k) (length(which(x_kn[,i]<=x_kn[k,i]))/n))
      }
      
      PDC=c()
      for (i in 1:p){
        PDC[i]=(pdcor(ynew,xnew[,i],xnew[,-i]))^2
      }
      scr_con2=PDC
      
      
      w_kn=scr_con1-scr_con2 
      thr=stop_early(W=w_kn,gamma = gam,offset)
      
      
      E[tim,] <- (w_kn >= thr) / (1 + sum(w_kn <= -thr))
    }
    eBH <- p * colMeans(E)
    select=ebh(eBH,alpha)$rej
    end_time<-Sys.time()
    runtime=end_time-start_time
    
    if(length(select)==0){
      fdr=0
      pow=0
    }
    if(length(select)!=0){
      fdr=length(setdiff(select,signal))/length(select)
      pow=length(intersect(signal,select))/length(signal)
    }
  }
  
  if(subset==TRUE){
    y=as.matrix(y)
    n=nrow(x)
    p=ncol(x)
    p_y=ncol(y)
    start_time<-Sys.time()
    y_boot=matrix(NA,nrow = n,ncol = p_y)
    for (i in 1:p_y) {
      y_boot[,i]=sample(y[,i],n,replace = T)
    }
    PD_sel1=c()
    PD_sel2=c()
    for (i in 1:p){
      PD_sel1[i]=(energy::dcor(y_boot,x[,i]))^2
      PD_sel2[i]=(energy::dcor(y,x[,i]))^2
    }
    
    PD_sel=which(PD_sel2>=quantile(PD_sel1,quant_boot))
    
    
    ynew=matrix(NA,nrow = n,ncol = ncol(y))
    for (i in 1:ncol(y)) {
      ynew[,i]=sapply(1:n,function(k) (length(which(y[,i]<=y[k,i]))/n))
    }
    xnew=matrix(NA,nrow = n,ncol = p)
    for (i in 1:p) {
      xnew[,i]=sapply(1:n,function(k) (length(which(x[,i]<=x[k,i]))/n))
    }
    
    PDC=c()
    for (i in 1:p){
      index_set=ifelse(length(which(PD_sel==i))==0,PD_sel,PD_sel[-which(PD_sel==i)])
      if(is.na(index_set)){
        PDC[i]=(pdcor(ynew,xnew[,i],xnew[,-i]))^2}
      if(!is.na(index_set)){
        PDC[i]=(pdcor(ynew,xnew[,i],xnew[,index_set]))^2}
    }
    scr_con1=PDC
    
    times=rep
    E=matrix(NA,nrow = times,ncol=p)
    for (tim in 1:times) {
      x_knn=x
      x_kn=create.second_order(x_knn,method = knock_method) 
      
      y_boot=matrix(NA,nrow = n,ncol = p_y)
      for (i in 1:p_y) {
        y_boot[,i]=sample(y[,i],n,replace = T)
      }
      PD_sel1=c()
      PD_sel2=c()
      for (i in 1:p){
        PD_sel1[i]=(energy::dcor(y_boot,x_kn[,i]))^2
        PD_sel2[i]=(energy::dcor(y,x_kn[,i]))^2
      }
      
      PD_sel=which(PD_sel2>=quantile(PD_sel1,quant_boot))
      
      ynew=matrix(NA,nrow = n,ncol = ncol(y))
      for (i in 1:ncol(y)) {
        ynew[,i]=sapply(1:n,function(k) (length(which(y[,i]<=y[k,i]))/n))
      }
      xnew=matrix(NA,nrow = n,ncol = p)
      for (i in 1:p) {
        xnew[,i]=sapply(1:n,function(k) (length(which(x_kn[,i]<=x_kn[k,i]))/n))
      }
      
      PDC=c()
      for (i in 1:p){
        index_set=ifelse(length(which(PD_sel==i))==0,PD_sel,PD_sel[-which(PD_sel==i)])
        if(is.na(index_set)){
          PDC[i]=(pdcor(ynew,xnew[,i],xnew[,-i]))^2}
        if(!is.na(index_set)){
          PDC[i]=(pdcor(ynew,xnew[,i],xnew[,index_set]))^2}
      }
      scr_con2=PDC
      
      
      w_kn=scr_con1-scr_con2 
      thr=stop_early(W=w_kn,gamma = gam,offset)
      
      
      E[tim,] <- (w_kn >= thr) / (1 + sum(w_kn <= -thr))
    }
    eBH <- p * colMeans(E)
    select=ebh(eBH,alpha)$rej
    end_time<-Sys.time()
    runtime=end_time-start_time
    
    if(length(select)==0){
      fdr=0
      pow=0
    }
    if(length(select)!=0){
      fdr=length(setdiff(select,signal))/length(select)
      pow=length(intersect(signal,select))/length(signal)
    }
  }
  return(list(select=select,fdr=fdr,power=pow,Time=runtime))
}

##########################################################################################################
##
#       function for PC-Knockoff
#
# Input: - y: The response vector or a data matrix.
#        - x: The design matrix of dimensions n * p. 
#        - multi: under multi-repsonse settings, we choose multi=TRUE, otherwise, multi=FALSE.
#        - knock_method: either "equi", "sdp" or "asdp" (default: "asdp"). This determines the method 
#that will be used to minimize the correlation between the original variables and the knockoffs.
#        - alpha: The predetermined fdr level.
#        - signal: The true active set.
#
# Output: - select: the selected features.
#         - fdr: the estimated fdr.
#         - power: the estimated power.
#         - Time: the computatioanl times of PC-Knockoff.
#
##########################################################################################################
PC_knockoff<-function(y,x,multi,knock_method,alpha,signal){
  if (multi==FALSE) {
    n=nrow(x)
    p=ncol(x)
    Arccos<-function (x) 
    {
      py_path = system.file("python", "PCSIS.py", package = "MFSIS")
      reticulate::source_python(py_path, envir = globalenv())
      ax = get_arccos_1d(x)
      return(ax)
    }
    Arccos_md<-function (x) 
    {
      py_path = system.file("python", "PCSIS.py", package = "MFSIS")
      reticulate::source_python(py_path, envir = globalenv())
      ax = get_arccos_md(x)
      return(ax)
    }
    pccoff<-function(x,xk,y)
    {
      w=c()
      yn=Arccos_md(y)
      for(i in 1:p){
        xnew=Arccos(x[,i])
        xknew=Arccos(xk[,i])
        w[i]=projection_corr(yn,xnew,n)^2-projection_corr(yn,xknew,n)^2
      }
      return(w)
    }
    knockoffs = function(X) create.second_order(X, method = knock_method)
    start_time<-Sys.time()
    result = knockoff.filter(x, y, knockoffs=knockoffs, statistic=pccoff,fdr=alpha)
    end_time<-Sys.time()
    runtime=end_time-start_time
    select=result$selected
    
    if(length(select)==0){
      fdr=0
      pow=0
    }
    if(length(select)!=0){
      fdr=length(setdiff(select,signal))/length(select)
      pow=length(intersect(signal,select))/length(signal)
    }
  }
  if (multi==TRUE) {
    n=nrow(x)
    p=ncol(x)
    Arccos<-function (x) 
    {
      py_path = system.file("python", "PCSIS.py", package = "MFSIS")
      reticulate::source_python(py_path, envir = globalenv())
      ax = get_arccos_1d(x)
      return(ax)
    }
    Arccos_md<-function (x) 
    {
      py_path = system.file("python", "PCSIS.py", package = "MFSIS")
      reticulate::source_python(py_path, envir = globalenv())
      ax = get_arccos_md(x)
      return(ax)
    }
    pccoff<-function(x,xk,y)
    {
      w=c()
      yn=Arccos_md(y)
      for(i in 1:p){
        xnew=Arccos(x[,i])
        xknew=Arccos(xk[,i])
        w[i]=projection_corr(yn,xnew,n)^2-projection_corr(yn,xknew,n)^2
      }
      return(w)
    }
    knockoff.filter<-function (X, y, knockoffs = create.second_order, statistic = stat.glmnet_coefdiff, 
                               fdr = 0.1, offset = 1) 
    {
      if (is.data.frame(X)) {
        X.names = names(X)
        X = as.matrix(X, rownames.force = F)
      }
      else if (is.matrix(X)) {
        X.names = colnames(X)
      }
      else {
        stop("Input X must be a numeric matrix or data frame")
      }
      if (!is.numeric(X)) 
        stop("Input X must be a numeric matrix or data frame")
      if (!is.factor(y) && !is.numeric(y)) {
        stop("Input y must be either of numeric or factor type")
      }
      # if (is.numeric(y)) 
      #   y = as.vector(y)
      if (offset != 1 && offset != 0) {
        stop("Input offset must be either 0 or 1")
      }
      if (!is.function(knockoffs)) 
        stop("Input knockoffs must be a function")
      if (!is.function(statistic)) 
        stop("Input statistic must be a function")
      n = nrow(X)
      p = ncol(X)
      # stopifnot(length(y) == n)
      if (identical(knockoffs, create.fixed)) 
        knockoffs = function(x) create.fixed(x, y = y)
      knock_variables = knockoffs(X)
      if (is(knock_variables, "knockoff.variables")) {
        X = knock_variables$X
        Xk = knock_variables$Xk
        if (!is.null(knock_variables$y)) 
          y = knock_variables$y
        rm(knock_variables)
      }
      else if (is(knock_variables, "matrix")) {
        Xk = knock_variables
        rm(knock_variables)
      }
      else {
        stop("Knockoff variables of incorrect type")
      }
      W = statistic(X, Xk, y)
      t = knockoff.threshold(W, fdr = fdr, offset = offset)
      selected = sort(which(W >= t))
      if (!is.null(X.names)) 
        names(selected) = X.names[selected]
      structure(list(call = match.call(), X = X, Xk = Xk, y = y, 
                     statistic = W, threshold = t, selected = selected), 
                class = "knockoff.result")
    }
    knockoffs = function(X) create.second_order(X, method = knock_method)
    start_time<-Sys.time()
    result = knockoff.filter(x, y, knockoffs=knockoffs, statistic=pccoff,fdr=alpha)
    end_time<-Sys.time()
    runtime=end_time-start_time
    select=result$selected
    
    if(length(select)==0){
      fdr=0
      pow=0
    }
    if(length(select)!=0){
      fdr=length(setdiff(select,signal))/length(select)
      pow=length(intersect(signal,select))/length(signal)
    }
  }
  return(list(select=select,fdr=fdr,power=pow,Time=runtime))
  }
