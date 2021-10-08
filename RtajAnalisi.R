setwd("/home/irpi/ITC/Tajikistan/Rv2")
source("/home/irpi/Dropbox/CNR/ITC/Rv2/TajFunctionsCV.R")
library(raster)
library(rgdal)
library(GISTools)
library(matrixStats)
library(boot)
library(parallel)
library(ROCR)
library(ggcorrplot)
library(bcaboot)
library(caret)
library(Hmisc)
library(INLA)
library(rlist)
load("data_scaled.Rdata")
Slope.CLASS = inla.group(data_scaled$S_mean, n = 20)
Rlf.CLASS = inla.group(data_scaled$Rlf_mean, n = 20)
Rain.CLASS = inla.group(data_scaled$Rn_mean, n = 20)
nonlinear=cbind(Slope.CLASS=Slope.CLASS, Rain.CLASS=Rain.CLASS,Rlf.CLASS=Rlf.CLASS)

formulaINLA = y~ -1 + intercept + PGA_stdev+PGA_mean+TnC_stdev+TnC_mean+PrC_mean+PrC_stdev+
  Rlf_stdev+S_stdev+Rn_stdev+surfM+
  f(Slope.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4) +
  f(Rain.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4)+
  f(Rlf.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4) +
  f(Litho, model="iid",hyper=hyper.iid, constr=T)+
  f(LU, model="iid",hyper=hyper.iid, constr=T)
#############################################

ROC1000=c()
library(Hmisc)
vp=10#cross-validation
iters=10#numero di ripetizioni dell'analisi
sampleP=5
perc.test=sampleP#%
giro=floor(100/sampleP)#riduzione numero di frane
summ.fixed=list()
inla.tot=list()
summ.categorical=list()
summ.mat=list()
summ.index=list()
#counta=1
for (k in 1:iters){
  data.mat5=c()
  data.mat5=data_scaled
  pp=length(data.mat5$LSD[data.mat5$LSD==1])
  ROCa=c()
  fixed=c()
  TOT=c()
  categorical=c()
  frane=c()
  index=c()
  mat=c()
  #summ.fixed=list()
  countb=1
  for (i in 1:giro){
    if (i>1){
      uni = data.mat5[data.mat5$LSD == 1,]
      campione=sample(uni$ID, size=floor(sampleP/100*pp), replace =F)
      data.mat5$LSD[data.mat5$ID %in% campione]=0
    }
    
    
    #####################################################################
    
    people=crossValidation(data.mat5,vp,data.mat5$LSD)
    idxlist=people$idxlist
    idxlist.p=people$idxlist.p
    idxlist.n=people$idxlist.n
    
    #To verify the ID frequency in idxlist
    n_occur = data.frame(table(unlist(idxlist,recursive=F)))
    n_occur[n_occur$Freq > 1,]
    # the result should be this: <0 rows> (or 0-length row.names)
    
    # 10-fold cross validation
    name=toString(i)
    
    person=kfold(data.mat5,vp,idxlist,data.mat5$LSD,nonlinear,name,formulaINLA)
    Northridge.cv=person$matt
    tot=person$tot
    fix=person$cvruns
    random=person$random
    
    ROCanalysis=froc(vp,idxlist,data.mat5$LSD,Northridge.cv,'CVruns_')
    
    
    
    #rownames(fix)=paste(rownames,toString(i),sep='')
    
    ######################################################################
    
    ROCa=cbind(ROCa,cv=ROCanalysis)#as.numeric(ROCanalysis$auc)
    TOT=append(TOT,list(tot))
    fixed=append(fixed,fix)
    categorical=append(categorical,random)
    mat=append(mat,list(data.mat5))
    index=append(index,idxlist)
  }
  inla.tot=append(inla.tot,list(TOT))
  summ.fixed=append(summ.fixed,list(fixed))
  summ.categorical=append(summ.categorical,list(categorical))
  ROC1000=append(ROC1000,list(ROCa))
  summ.mat=append(summ.mat,list(data.frame(mat)))
  summ.index=append(summ.index,list(index))
  #counta=counta+1
}

#df5=list.rbind(ROC1000)
#df2=data.frame(ROC1000)
#df5=data.frame(matrix((ROC1000), 200, 20, byrow = T))
#matplot(t(df5),type="l")
#df3=rowMeans(df2)
#df4=data.frame(matrix(df3, 200, 20, byrow = T))
#matplot(t(df4),type="l")
save(ROC1000,file="AUClsd1000.Rdata")
#save(df4,file="AUClsd1000media.Rdata")
save(inla.tot,file='inla1000.Rdata')
save(summ.fixed,file='summ.fixed_lsd1000.Rdata')
save(summ.categorical,file='summ.categorical_lsd1000.Rdata')
save(summ.index,file='index1000.Rdata')
save(summ.mat,file='mat1000.Rdata')
