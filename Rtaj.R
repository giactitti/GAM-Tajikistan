setwd("/work/folder")
source("/function/folder/TajFunctionsCV.R")
load("data_scaled.Rdata")
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

#fitting
Slope.CLASS = inla.group(data_scaled$S_mean, n = 20)



Rlf.CLASS = inla.group(data_scaled$Rlf_mean, n = 20)
Rain.CLASS = inla.group(data_scaled$Rn_mean, n = 20)

# load('data_df.Rdata')
Rain.deCLASS = inla.group(df$Rn_mean, n = 20)
sort(unique(Rain.deCLASS))
Slope.deCLASS = inla.group(df$S_mean, n = 20)
sort(unique(Slope.deCLASS))
Rlf.deCLASS = inla.group(df$Rlf_mean, n = 20)
sort(unique(Rlf.deCLASS))

Dataset1 = data.frame(intercept = 1, Ntrials = 1, data_scaled, Slope.CLASS=Slope.CLASS,Rain.CLASS=Rain.CLASS,Rlf.CLASS=Rlf.CLASS)
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))
hyper.iid = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))

attach(Dataset1)

Formula1.INLA = LSD~ -1 + intercept+
  PGA_stdev+PGA_mean+TnC_stdev+TnC_mean+PrC_mean+PrC_stdev+
  Rlf_stdev+S_stdev+Rn_stdev+surfM+
  f(Slope.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4) +
  f(Rain.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4)+
  f(Rlf.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4) +
  f(Litho, model="iid",hyper=hyper.iid, constr=T)+
  f(LU, model="iid",hyper=hyper.iid, constr=T)

Mod.INLA1 = inla(formula = Formula1.INLA, family = "binomial", data=Dataset1,
                 control.fixed=list(prec=.1),
                 num.threads = 2,
                 Ntrials = Ntrials,
                 control.family = list(control.link = list(model = "logit")),
                 control.mode = list(restart = T),
                 control.inla = list(int.strategy = "eb"),
                 control.predictor = list(compute = T, link = 1), verbose = T)
save(Mod.INLA1,file='susc100.Rdata')



plot(Mod.INLA1)



table(Dataset1$LSD)
par(mfrow=c(2,1))
plot(density(Mod.INLA1$summary.fitted.values$mean))
m1=Mod.INLA1$summary.fitted.values$mean
ci=Mod.INLA1$summary.fitted.values$'0.975quant'-Mod.INLA1$summary.fitted.values$'0.025quant'
plot(m1,ci)
dev.off()

#############
extractSignificantEffects=function(results.df){
  print(results.df[sign(results.df$`0.025quant`)==sign(results.df$`0.975quant`),])
}
extractSignificantEffects(Mod.INLA1$summary.fixed)
########

library(pROC)
ROC=roc(data_scaled$LSD~Mod.INLA1$summary.fitted.values$mean)
save(ROC,file='ROC100.Rdata')
plot.roc(data_scaled$LSD~Mod.INLA1$summary.fitted.values$mean)
dev.off()



#########################################
susc4map=data.frame(ID=data_scaled$ID,suscM=m1,suscCI=ci)
Q1=quantile(susc4map$suscM,probs=0.25)
Q3=quantile(susc4map$suscM,probs=0.75)
IQR=Q3-Q1
med=median(susc4map$suscM)
classes=matrix(c(Q1-(1.5*IQR),Q1,med,Q3,Q3+(1.5*IQR)),ncol=5, byrow=TRUE)
colnames(classes) = c('Q1-1.5*IQR','Q1','Median','Q3','Q3+1.5*IQR')
write.table(classes,'boxplots.txt',col.names = T,row.names = F,sep = '\t')


# ###########################################
# plot(density(Mod.INLA1$summary.fitted.values$mean))
# abline(v=quantile(Mod.INLA1$summary.fitted.values$mean,probs=0.05),col='red')
# abline(v=quantile(Mod.INLA1$summary.fitted.values$mean,probs=0.63),col='red')
# abline(v=quantile(Mod.INLA1$summary.fitted.values$mean,probs=0.95),col='red')
classes=matrix(c(0,quantile(susc4map$suscM,probs=0.05),
                 quantile(susc4map$suscM,probs=0.63),
                 quantile(susc4map$suscM,probs=0.95),1),ncol=5, byrow=TRUE)
colnames(classes) = c('0%','5%','63%','95%','100%')
write.table(classes,'classes.txt',col.names = T,row.names = F,sep = '\t')
###########################################################

CM=readOGR("/home/irpi/Dropbox/CNR/ITC/catchements/intersection_medium_valid_last_v2.shp")
susc4map=merge(CM, susc4map, by = "ID")
writeOGR(obj=susc4map, dsn="FitSusc", layer="susc4map", driver="ESRI Shapefile", overwrite_layer = T)

########################################
#cross validation
vp=10
formulaINLA = y~ -1 + intercept + PGA_stdev+PGA_mean+TnC_stdev+TnC_mean+PrC_mean+PrC_stdev+
  Rlf_stdev+S_stdev+Rn_stdev+surfM+
  f(Slope.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4) +
  f(Rain.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4)+
  f(Rlf.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4) +
  f(Litho, model="iid",hyper=hyper.iid, constr=T)+
  f(LU, model="iid",hyper=hyper.iid, constr=T)


people=crossValidation(data_scaled,vp,data_scaled$LSD)
idxlist=people$idxlist
#save(idxlist,file='idx100.Rdata')

#To verify the ID frequency in idxlist
n_occur = data.frame(table(unlist(idxlist,recursive=F)))
n_occur[n_occur$Freq > 1,]
# the result should be this: <0 rows> (or 0-length row.names)

# 10-fold cross validation
nonlinear=cbind(Slope.CLASS=Slope.CLASS, Rain.CLASS=Rain.CLASS,Rlf.CLASS=Rlf.CLASS)

person=kfold(data_scaled,vp,idxlist,data_scaled$LSD,nonlinear,'cv',formulaINLA)
Northridge.cv=person$matt
tot=person$tot
fix=person$cvruns
random=person$random

ROCanalysis=froc(vp,idxlist,data_scaled$LSD,Northridge.cv,'CV')

save(Northridge.cv,file='summ_lsd.Rdata')
save(fix,file='summ.fixed_lsd.Rdata')
save(random,file='summ.categorical_lsd.Rdata')
save(idxlist,file='idx.Rdata')
save(tot,file='inla_lsd.Rdata')
save(ROCanalysis,file="AUClsd.Rdata")
