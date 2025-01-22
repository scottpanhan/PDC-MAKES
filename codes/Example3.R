## code for Example 2 in main text###
library("mvtnorm")
library("Ball")
source("./functions.R")
sim=200
MMS_rPDC=c()
MMS_PDC=c()
MMS_DC=c()
MMS_BC=c()
MMS_PC=c()
p_all1=c()
pow1=c()
p_all2=c()
pow2=c()
p_all3=c()
pow3=c()
p_all4=c()
pow4=c()
p_all5=c()
pow5=c()
rPDC_time=c()
PDC_time=c()
DC_time=c()
BC_time=c()
PC_time=c()
for (m in 1:sim) {
  print(m)
  n=100
  p=2000
  #p=5000
  #p=10000
  mu=rep(0,p)
  pho1=0.25
  sigma_x=matrix(pho1,nrow = p,ncol = p)
  sigma_x[4,]=rep(sqrt(pho1),p)
  sigma_x[,4]=rep(sqrt(pho1),p)
  sigma_x[5,]=rep(0,p)
  sigma_x[,5]=rep(0,p)
  diag(sigma_x)=rep(1,p)
  #sigma_x<- matrix(0.5,nrow=p,ncol=p)
  #diag(sigma_x)=1
  x=rmvnorm(n, mean = mu, sigma = sigma_x,method = "chol")
  x1=x[,1]
  x2=x[,2]
  x3=x[,3]
  x4=x[,4]
  x5=x[,5]
  beta=2.5
  y=beta*x1+beta*x2+beta*x3-3*beta*sqrt(pho1)*x4+0.4*beta*x5+rnorm(n,mean = 0,sd=sqrt((beta^2)*pho1))
  
  
  signal=c(1,2,3,4,5)
  
  ###rPDC-Screen###
  rPDC_select=rPDC_SIS(y,x,subset=TRUE,quant_boot=0.7,dn=p)
  PDCC=rPDC_select$select
  rec_PDC=c()
  for (i in 1:length(signal)) {
    rec_PDC[i]=which(PDCC==signal[i])
  }
  MMS_rPDC[m]=max(rec_PDC)
  dn=round(n/log(n))
  p_all1[m]=ifelse(MMS_rPDC[m]<=dn,1,0)
  pow1[m]=length(intersect(signal,PDCC[1:dn]))/length(signal)
  rPDC_time[m]=rPDC_select$Time
  
  ###PDC-Screen###
  PDC_select=PDC_SIS(y,x,subset=TRUE,quant_boot=0.7,dn=p)
  PDCC2=PDC_select$select
  rec_PDC2=c()
  for (i in 1:length(signal)) {
    rec_PDC2[i]=which(PDCC2==signal[i])
  }
  MMS_PDC[m]=max(rec_PDC2)
  p_all2[m]=ifelse(MMS_PDC[m]<=dn,1,0)
  pow2[m]=length(intersect(signal,PDCC2[1:dn]))/length(signal)
  PDC_time[m]=PDC_select$Time
  
  ###DC-SIS###
  start_time<-Sys.time()
  DC=c()
  for (i in 1:p) {
    DC[i]=(energy::dcor(y,x[,i]))^2
  }
  DCC=order(DC,decreasing = T)
  end_time<-Sys.time()
  rec_DC=c()
  for (i in 1:length(signal)) {
    rec_DC[i]=which(DCC==signal[i])
  }
  MMS_DC[m]=max(rec_DC)
  p_all3[m]=ifelse(MMS_DC[m]<=dn,1,0)
  pow3[m]=length(intersect(signal,DCC[1:dn]))/length(signal)
  DC_time[m]=end_time-start_time
  
  ###Bcor-SIS###
  start_time<-Sys.time()
  BCC=bcorsis(x = x, y = y, method = "standard", d = p)$ix
  end_time<-Sys.time()
  rec_BC=c()
  for (i in 1:length(signal)) {
    rec_BC[i]=which(BCC==signal[i])
  }
  MMS_BC[m]=max(rec_BC)
  p_all4[m]=ifelse(MMS_BC[m]<=dn,1,0)
  pow4[m]=length(intersect(signal,BCC[1:dn]))/length(signal)
  BC_time[m]=end_time-start_time
  
  ###PC-Screen
  PC_select=PC_SIS(x,y,multi=FALSE,dn=p)
  PCC=PC_select$select
  rec_PC=c()
  for (i in 1:length(signal)) {
    rec_PC[i]=which(PCC==signal[i])
  }
  MMS_PC[m]=max(rec_PC)
  p_all5[m]=ifelse(MMS_PC[m]<=dn,1,0)
  pow5[m]=length(intersect(signal,PCC[1:dn]))/length(signal)
  PC_time[m]=PC_select$Time
}
quantile(MMS_rPDC,c(0.05,0.25,0.5,0.75,0.95))
quantile(MMS_PDC,c(0.05,0.25,0.5,0.75,0.95))
quantile(MMS_DC,c(0.05,0.25,0.5,0.75,0.95))
quantile(MMS_BC,c(0.05,0.25,0.5,0.75,0.95))
quantile(MMS_PC,c(0.05,0.25,0.5,0.75,0.95))
mean(p_all1)
mean(pow1)
mean(p_all2)
mean(pow2)
mean(p_all3)
mean(pow3)
mean(p_all4)
mean(pow4)
mean(p_all5)
mean(pow5)
mean(rPDC_time)
mean(PDC_time)
mean(DC_time)
mean(BC_time)
mean(PC_time)