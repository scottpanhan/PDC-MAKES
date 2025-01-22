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
  pho1=0.5
  sigma_x=matrix(NA,nrow = p,ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma_x[i,j]=pho1^(abs(i-j))
    }
  }
  x_ori=rmvnorm(n, mean = mu, sigma = sigma_x,method = "chol")
  eve=seq(2,p,by=2)
  x=x_ori
  for (i in eve) {
    for (j in 1:n) {
      if(x_ori[j,i]<0){x[j,i]=0}
      if(x_ori[j,i]>=0&x_ori[j,i]<=1.5){x[j,i]=1}
      if(x_ori[j,i]>1.5){x[j,i]=2}
    }
    x[,i]=x[,i]-mean(x[,i])
  }
  
  y_o1=matrix(NA,nrow = n,ncol = 10)
  beta_1=(-1)^rbinom(1,size=1,prob=0.5)*runif(1,min = 1,max = 2)
  for (k in 1:10) {
    y_o1[,k]=beta_1*(1+x[,(2*k-1)]+x[,(2*k)])^3*exp(-x[,(2*k)])+rnorm(n)
  }
  
  y_o2=matrix(NA,nrow = n,ncol = 10)
  beta_2=(-1)^rbinom(1,size=1,prob=0.5)*runif(1,min = 1,max = 2)
  for (k in 1:10) {
    y_o2[,k]=beta_2*(1+x[,(2*k-1)]+x[,(2*k)])^3*exp(-x[,(2*k)])+rnorm(n)
  }
  
  y_o3=matrix(NA,nrow = n,ncol = 10)
  beta_3=(-1)^rbinom(1,size=1,prob=0.5)*runif(1,min = 1,max = 2)
  for (k in 1:10) {
    y_o3[,k]=beta_3*(1+x[,(2*k-1)]+x[,(2*k)])^3*exp(-x[,(2*k)])+rnorm(n)
  }
  
  y_o4=matrix(NA,nrow = n,ncol = 10)
  beta_4=(-1)^rbinom(1,size=1,prob=0.5)*runif(1,min = 1,max = 2)
  for (k in 1:10) {
    y_o4[,k]=beta_4*(1+x[,(2*k-1)]+x[,(2*k)])^3*exp(-x[,(2*k)])+rnorm(n)
  }
  
  y_o5=matrix(NA,nrow = n,ncol = 10)
  beta_5=(-1)^rbinom(1,size=1,prob=0.5)*runif(1,min = 1,max = 2)
  for (k in 1:10) {
    y_o5[,k]=beta_5*(1+x[,(2*k-1)]+x[,(2*k)])^3*exp(-x[,(2*k)])+rnorm(n)
  }
  
  y=cbind(y_o1,y_o2,y_o3,y_o4,y_o5)
  
  signal=rep(1:20)
  
  ###rPDC-Screen###
  rPDC_select=rPDC_SIS(y,x,subset=TRUE,quant_boot=0.99,dn=p)
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
  PDC_select=PDC_SIS(y,x,subset=TRUE,quant_boot=0.99,dn=p)
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
  PC_select=PC_SIS(x,y,multi=TRUE,dn=p)
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