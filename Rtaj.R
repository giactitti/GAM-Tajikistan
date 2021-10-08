setwd("/home/irpi/ITC/Tajikistan/Rv2")
source("/home/irpi/Dropbox/CNR/ITC/Rv2/TajFunctionsCV.R")
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
###############################àà






#####################################
# #save and plot CV
# color=rainbow(n=10)
# plot(roc(data_scaled$LSD~Mod.INLA1$summary.fitted.values$mean),lwd=1)
# #ROCtest=c()
# for(i in 1:vp){
#   idx.holdout=idxlist[[i]]
#   tmp.y=data_scaled$LSD
#   #tmp.y[idx.holdout]=NA
#   #tmp.y=na.omit(tmp.y)
#   ROCtest=roc(tmp.y[-idx.holdout]~Northridge.cv[-idx.holdout])
#   #save(ROCtest,file=paste('CV_ROCtest',toString(i),'.Rdata'))
#   #plot(roc(tmp.y[-idx.holdout]~Northridge.cv[-idx.holdout]),col=color[i],add = TRUE,lwd=0.8)
#   ROCvalid=roc(tmp.y[idx.holdout]~Northridge.cv[idx.holdout])
#   #save(ROCvalid,file=paste('CV_ROCvalid',toString(i),'.Rdata'))
#   #plot(roc(tmp.y[idx.holdout]~Northridge.cv[idx.holdout]),col=color[i],add = TRUE,lwd=0.8)
#   print(ggroc(ROCvalid))
# }




#ROCv=froc(vp,idxlist,data_scaled$LSD,Northridge.cv,'CV')
####################################################################################################


####################################################################################################

# corIn=Dataset1[,c('PGA_stdev','PGA_mean','PlC_mean','PlC_stdev','PrC_mean','PrC_stdev',
# 'Rlf_mean','Rlf_stdev','S_stdev','Flw_mean','Flw_stdev','topo_surf','Slope.CLASS',
# 'Rain.CLASS',"Rn_surf","Flw_surf",'PGA_surf')]
# 
# require(corrplot)
# 
# corAn=cor(corIn)
# 
# corrplot(corAn,order='AOE',method='color',addCoef.col='gray')
# 
# 
# 
# extractSignificantEffects=function(results.df){
#   print(results.df[sign(results.df$`0.025quant`)==sign(results.df$`0.975quant`),])
# }
# 
# extractSignificantEffects(Mod.INLA1$summary.fixed)
# 
# slope.info2 = Mod.INLA1$summary.random$Slope.CLASS
# 
# par(mfrow=c(1,2))
# plot(Mod.INLA1$summary.random$Slope.CLASS$ID, Mod.INLA1$summary.random$Slope.CLASS$mean, 
#      type = "l", xlab = "Slope ID", ylab = "Regression Coefficients", ylim = c(-4,4))
# lines(Mod.INLA1$summary.random$Slope.CLASS$ID, Mod.INLA1$summary.random$Slope.CLASS$`0.025quant`, pch = 19,
#       col = "red")
# lines(Mod.INLA1$summary.random$Slope.CLASS$ID, Mod.INLA1$summary.random$Slope.CLASS$`0.975quant`, pch = 19,
#       col = "red")
# abline(h = 0, col = "grey")
# 
# plot(Mod.INLA1$summary.random$Slope.CLASS$ID, Mod.INLA1$summary.random$Slope.CLASS$mean, 
#      pch = 19, xlab = "Slope ID", ylab = "Regression Coefficients", ylim = c(-4,4))
# points(Mod.INLA1$summary.random$Slope.CLASS$ID, Mod.INLA1$summary.random$Slope.CLASS$`0.025quant`, pch = 19,
#        col = "red")
# points(Mod.INLA1$summary.random$Slope.CLASS$ID, Mod.INLA1$summary.random$Slope.CLASS$`0.975quant`, pch = 19,
#        col = "red")
# abline(h = 0, col = "grey")
# 
# 
# pga.info2 = Mod.INLA1$summary.random$PGA.CLASS
# 
# par(mfrow=c(1,2))
# plot(Mod.INLA1$summary.random$PGA.CLASS$ID, Mod.INLA1$summary.random$PGA.CLASS$mean, 
#      type = "l", xlab = "PGA ID", ylab = "Regression Coefficients", ylim = c(-4,4))
# lines(Mod.INLA1$summary.random$PGA.CLASS$ID, Mod.INLA1$summary.random$PGA.CLASS$`0.025quant`, pch = 19,
#       col = "red")
# lines(Mod.INLA1$summary.random$PGA.CLASS$ID, Mod.INLA1$summary.random$PGA.CLASS$`0.975quant`, pch = 19,
#       col = "red")
# abline(h = 0, col = "grey")
# 
# plot(Mod.INLA1$summary.random$PGA.CLASS$ID, Mod.INLA1$summary.random$PGA.CLASS$mean, 
#      pch = 19, xlab = "PGA ID", ylab = "Regression Coefficients", ylim = c(-4,4))
# points(Mod.INLA1$summary.random$PGA.CLASS$ID, Mod.INLA1$summary.random$PGA.CLASS$`0.025quant`, pch = 19,
#        col = "red")
# points(Mod.INLA1$summary.random$PGA.CLASS$ID, Mod.INLA1$summary.random$PGA.CLASS$`0.975quant`, pch = 19,
#        col = "red")
# abline(h = 0, col = "grey")
# 
# 
# Mod.INLA1$summary.fitted.values
# 
# ## or to see better
# 
# plot(density(Mod.INLA1$summary.fitted.values$mean)) ## What do you see?
# 
# plot(Mod.INLA1$summary.fixed$mean[2],ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
# text(1,0.6,"surf", font = 2, cex = 1.4)
# points(Mod.INLA1$summary.fixed$`0.025quant`[2],col="blue",pch=18, cex = 2)
# points(Mod.INLA1$summary.fixed$`0.975quant`[2],col="blue",pch=18, cex = 2)
# segments(1,Mod.INLA1$summary.fixed$`0.025quant`[2],1,Lushan.MIonly$summary.fixed$`0.975quant`[2],lwd=2, lty = 3, col = "blue")
# abline(h=0,lty=2,lwd=2,col="gray50")
# 
# ## However, if we want a comprehensive view of all fixed effect, we can plot them in a subpaneled structure as follows:
# par(mfrow=c(3,3))
# 
# ## Panel 
# plot(Mod.INLA1$summary.fixed$mean[2], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
# text(1,0.6,"surf", font = 2, cex = 1.4)
# points(Lushan.MIonly$summary.fixed$`0.025quant`[2],col="blue",pch=18, cex = 2)
# points(Lushan.MIonly$summary.fixed$`0.975quant`[2],col="blue",pch=18, cex = 2)
# segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[2],1,Lushan.MIonly$summary.fixed$`0.975quant`[2],lwd=2, lty = 3, col = "blue")
# abline(h=0,lty=2,lwd=2,col="gray50")
# ## Panel 
# plot(Mod.INLA1$summary.fixed$mean[3], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
# text(1,0.6,"PGA_stdev", font = 2, cex = 1.4)
# points(Lushan.MIonly$summary.fixed$`0.025quant`[3],col="blue",pch=18, cex = 2)
# points(Lushan.MIonly$summary.fixed$`0.975quant`[3],col="blue",pch=18, cex = 2)
# segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[3],1,Lushan.MIonly$summary.fixed$`0.975quant`[3],lwd=2, lty = 3, col = "blue")
# abline(h=0,lty=2,lwd=2,col="gray50")
# ## Panel 
# plot(Mod.INLA1$summary.fixed$mean[4], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
# text(1,0.6,"PlC_mean", font = 2, cex = 1.4)
# points(Lushan.MIonly$summary.fixed$`0.025quant`[4],col="blue",pch=18, cex = 2)
# points(Lushan.MIonly$summary.fixed$`0.975quant`[4],col="blue",pch=18, cex = 2)
# segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[4],1,Lushan.MIonly$summary.fixed$`0.975quant`[4],lwd=2, lty = 3, col = "blue")
# abline(h=0,lty=2,lwd=2,col="gray50")
# ## Panel 
# plot(Mod.INLA1$summary.fixed$mean[5], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
# text(1,0.6,"PlC_stdev", font = 2, cex = 1.4)
# points(Lushan.MIonly$summary.fixed$`0.025quant`[5],col="blue",pch=18, cex = 2)
# points(Lushan.MIonly$summary.fixed$`0.975quant`[5],col="blue",pch=18, cex = 2)
# segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[5],1,Lushan.MIonly$summary.fixed$`0.975quant`[5],lwd=2, lty = 3, col = "blue")
# abline(h=0,lty=2,lwd=2,col="gray50")
# ## Panel 
# plot(Mod.INLA1$summary.fixed$mean[6], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
# text(1,0.6,"PrC_mean", font = 2, cex = 1.4)
# points(Lushan.MIonly$summary.fixed$`0.025quant`[6],col="blue",pch=18, cex = 2)
# points(Lushan.MIonly$summary.fixed$`0.975quant`[6],col="blue",pch=18, cex = 2)
# segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[6],1,Lushan.MIonly$summary.fixed$`0.975quant`[6],lwd=2, lty = 3, col = "blue")
# abline(h=0,lty=2,lwd=2,col="gray50")
# ## Panel 
# plot(Mod.INLA1$summary.fixed$mean[7], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
# text(1,0.6,"PrC_stdev", font = 2, cex = 1.4)
# points(Lushan.MIonly$summary.fixed$`0.025quant`[7],col="blue",pch=18, cex = 2)
# points(Lushan.MIonly$summary.fixed$`0.975quant`[7],col="blue",pch=18, cex = 2)
# segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[7],1,Lushan.MIonly$summary.fixed$`0.975quant`[7],lwd=2, lty = 3, col = "blue")
# abline(h=0,lty=2,lwd=2,col="gray50")
# ## Panel 
# plot(Mod.INLA1$summary.fixed$mean[8], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
# text(1,0.6,"Rlf_mean", font = 2, cex = 1.4)
# points(Mod.INLA1$summary.fixed$`0.025quant`[8],col="blue",pch=18, cex = 2)
# points(Mod.INLA1$summary.fixed$`0.975quant`[8],col="blue",pch=18, cex = 2)
# segments(1,Mod.INLA1$summary.fixed$`0.025quant`[8],1,Mod.INLA1$summary.fixed$`0.975quant`[8],lwd=2, lty = 3, col = "blue")
# abline(h=0,lty=2,lwd=2,col="gray50")
# ## Panel 
# plot(Mod.INLA1$summary.fixed$mean[9], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
# text(1,0.6,"Rlf_stdev", font = 2, cex = 1.4)
# points(Mod.INLA1$summary.fixed$`0.025quant`[9],col="blue",pch=18, cex = 2)
# points(Mod.INLA1$summary.fixed$`0.975quant`[9],col="blue",pch=18, cex = 2)
# segments(1,Mod.INLA1$summary.fixed$`0.025quant`[9],1,Mod.INLA1$summary.fixed$`0.975quant`[9],lwd=2, lty = 3, col = "blue")
# abline(h=0,lty=2,lwd=2,col="gray50")
# ## Panel 
# plot(Mod.INLA1$summary.fixed$mean[10], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
# text(1,0.6,"S_stdev", font = 2, cex = 1.4)
# points(Mod.INLA1$summary.fixed$`0.025quant`[10],col="blue",pch=18, cex = 2)
# points(Mod.INLA1$summary.fixed$`0.975quant`[10],col="blue",pch=18, cex = 2)
# segments(1,Mod.INLA1$summary.fixed$`0.025quant`[10],1,Mod.INLA1$summary.fixed$`0.975quant`[10],lwd=2, lty = 3, col = "blue")
# abline(h=0,lty=2,lwd=2,col="gray50")
# 
# 
# 
# ## With a similar concept in mind, we can asses the effect of non linear covariates, by plotting each regression coefficient
# ## First of all, we need to load a file where each Lithotype name recorded against the numeric ID we used in the calculations
# # Lushan.MIonly$summary.random$Lithology ## Let's see directly what I mean
# # LitoNames = read.delim("Lithologies.txt", sep = "\t", header = T)
# # View(LitoNames)
# # plot(1:5+0.22,Lushan.MIonly$summary.random$Lithology$mean,xlab="Lithology",ylab="Linear predictor",lwd=2,ylim=c(-2.7,2.7),pch=19,xaxt="n", col="black", cex.lab=1.5, cex = 2)
# # axis(1,at=1:5+0.22,labels=c("D","P","Pt","Q","T"), cex.axis = 1.5)
# # points(1:5+0.22,Lushan.MIonly$summary.random$Lithology$"0.025quant",lwd=2,col="blue",pch=18, cex = 2)
# # points(1:5+0.22,Lushan.MIonly$summary.random$Lithology$"0.975quant",lwd=2,col="blue",pch=18, cex = 2)
# # segments(1:5+0.22,Lushan.MIonly$summary.random$Lithology$"0.025quant",1:5+0.22,Lushan.MIonly$summary.random$Lithology$"0.975quant",lwd=2, lty = 3, col = "blue")
# # abline(h=0,lty=2,lwd=2,col="gray50")
# # text(1:5+0.15,y=Lushan.MIonly$summary.random$Lithology$mean,as.character(LitoNames$Class),srt=90, font = 2, cex = 1.75)
# 
# save(Mod.INLA1, file = "ModINLA1.Rdata")
# 
# library(pROC)
# roc(data_scaled$LSD~Mod.INLA1$summary.fitted.values$mean)

