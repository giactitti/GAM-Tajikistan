setwd("/home/irpi/ITC/Tajikistan/Rv2")
library(raster)
library(rgdal)

#carico land use e linearizzo
CMlu=readOGR("/home/irpi/Dropbox/CNR/ITC/catchements/intersection_medium_LUmasked_v2.shp")

dflu = as(CMlu, "data.frame")
dflu$ID=as.integer(dflu$ID)

dflu$CI=dflu[,c('HISTO_1')]
dflu$CR=dflu[,c('HISTO_2')]
dflu$Cv=rowSums(dflu[,c('HISTO_3','HISTO_4')])
dflu$Fo=rowSums(dflu[,c('HISTO_5','HISTO_7','HISTO_8','HISTO_9','HISTO_10')])
dflu$FSG=dflu[,c('HISTO_11')]
dflu$GFS=dflu[,c('HISTO_12')]
#dflu$SC=dflu[,c('HISTO_13')]
dflu$GS=dflu[,c('HISTO_14')]
dflu$Ba=rowSums(dflu[,c('HISTO_15','HISTO_19')])
#dflu$FRF=dflu[,c('HISTO_16')]
dflu$GRF=dflu[,c('HISTO_17')]
dflu$U=rowSums(dflu[,c('HISTO_18','HISTO_23','HISTO_24','HISTO_25','HISTO_27')])
dflu$WB=dflu[,c('HISTO_20')]
dflu$SI=dflu[,c('HISTO_21')]
dflu$ND=dflu[,c('HISTO_22')]
dflu$AOV=dflu[,c('HISTO_26')]
dflu$QU=dflu[,c('HISTO_28')]

da=dflu[,c('CI','CR','Cv','Fo','FSG','GFS','GS','Ba','GRF','U','WB','SI','AOV','QU')]
#da=da[,-13]#tolgo no data
dm=data.matrix(da)

#mode(dm) = "integer"
#round(colSums(dm)/sum(colSums(dm))*100,2)#calcolo la % di area di ciascuna classe

# dm = t(apply(dm, 1, function(i) i/sum(i)*100))
#sum1=rowSums(dm)#clacolo la somma iniziale 100%

# dm[is.na(dm)]=0

#tolgo le colonne di cui non ho bisogno
#ix=which(colnames(dm) %in% c('ND','AOV','QU','GRF','GFS','FSG','SI','WB','U','UAG'))
#dm = dm[,-ix] 

# sum2=rowSums(dm)#calcolo di nuovo la somma
# r1=sum2/sum1
# r1=round(r1,0)
# row=which(r1 == 100, arr.ind=TRUE)#questi sono le feature ancolra al 100%
# dm[dm==0]=NA
# 
# #scelgo i minimi da escludere per tutti quelli che raggiungono ancora il 100%
# column=apply(dm,1,which.min)
# column1=as.numeric(column)
# dm1=dm
# for (i in 1:length(row)) {dm1[row[i],column1[row[i]]]=0}
# dm[is.na(dm)]=0
# dm1[is.na(dm1)]=0
# sum3=rowSums(dm1)
# r2=sum3/sum1
# #############verifiche finali
# which(r2 == 100)#verifico che sia vettore nullo
# diff=setdiff(names(sum3[sum3==0]),names(sum2[sum2==0]))#####
# diff=as.integer(diff)
# for (ii in 1:length(diff)) {dm1[diff[ii]+1,column1[diff[ii]+1]]=dm[diff[ii]+1,column1[1+diff[ii]]]}
# sum4=rowSums(dm1)
# setdiff(names(sum4[sum4==0]),names(sum2[sum2==0]))#####deve avere lunghezza nulla 

#df.LU=as.data.frame(dm1)

M=cbind(1:nrow(dm), max.col(dm, 'first'))
df.LU=as.data.frame(dm)
somma=rowSums(df.LU)
idx=which(somma==0)
df.LU$LU=M[,2]
df.LU$LU[idx]=NaN

df.LU=data.frame(ID=dflu$ID,df.LU)
df.LU=df.LU[,c('ID','LU')]


save(df.LU,file="dfLUpred_v2.Rda")
