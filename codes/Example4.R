library("mvtnorm")
source("./functions.R")
sim=200
sel_pdc=matrix(NA,nrow = sim,ncol = 7)
sel_rpdc=matrix(NA,nrow = sim,ncol = 7)
sel_pc=matrix(NA,nrow = sim,ncol = 7)
fdr_pdc=c()
fdr_rpdc=c()
fdr_pc=c()
pow_pdc=c()
pow_rpdc=c()
pow_pc=c()
Msize_pdc=c()
Msize_rpdc=c()
Msize_pc=c()
for (m in 1:sim) {
  print(m)
  n=200
  p=2000
  mu=rep(0,p)
  pho1=0.5
  sigma_x=matrix(NA,nrow = p,ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma_x[i,j]=pho1^(abs(i-j))
    }
  }
  x=rmvnorm(n, mean = mu, sigma = sigma_x,method = "chol")
  y=matrix(NA,nrow = n,ncol = 3)
  for (i in 1:3) {
    y[,i]=x[,(2*i-1)]^3+x[,(2*i)]^3+rnorm(n)
  }
  
  signal=rep(1:6)
  
  PDCMAKES=PDC_MAKES(y,x,subset=TRUE,rep=50,quant_boot=0.95,knock_method="equi",gam=0.15,alpha=0.2,offset=1,signal)
  for (i in 1:6) {
    sel_pdc[m,i]=ifelse(i%in%PDCMAKES$select,1,0)
  }
  sel_pdc[m,7]=ifelse(sum(sel_pdc[m,1:6])==6,1,0)
  fdr_pdc[m]=PDCMAKES$fdr
  Msize_pdc[m]=length(PDCMAKES$select)
  pow_pdc[m]=PDCMAKES$power
  
  rPDCMAKES=rPDC_MAKES(y,x,subset=TRUE,rep=50,quant_boot=0.95,knock_method="equi",gam=0.15,alpha=0.2,offset=1,signal)
  for (i in 1:6) {
    sel_rpdc[m,i]=ifelse(i%in%rPDCMAKES$select,1,0)
  }
  sel_rpdc[m,7]=ifelse(sum(sel_rpdc[m,1:6])==6,1,0)
  fdr_rpdc[m]=rPDCMAKES$fdr
  Msize_rpdc[m]=length(rPDCMAKES$select)
  pow_rpdc[m]=rPDCMAKES$power
  
  PCK=PC_knockoff(y,x,multi=TRUE,knock_method="asdp",alpha=0.2,signal)
  for (i in 1:6) {
    sel_pc[m,i]=ifelse(i%in%PCK$select,1,0)
  }
  sel_pc[m,7]=ifelse(sum(sel_pc[m,1:6])==6,1,0)
  fdr_pc[m]=PCK$fdr
  Msize_pc[m]=length(PCK$select)
  pow_pc[m]=PCK$power
}
colSums(sel_pdc)/sim
mean(fdr_pdc)
mean(pow_pdc)
mean(Msize_pdc)
colSums(sel_rpdc)/sim
mean(fdr_rpdc)
mean(pow_rpdc)
mean(Msize_rpdc)
colSums(sel_pc)/sim
mean(fdr_pc)
mean(pow_pc)
mean(Msize_pc)