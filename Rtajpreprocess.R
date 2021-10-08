setwd("/home/irpi/ITC/Tajikistan/Rv2")
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

load("dfLUpred_v2.Rda")
lsd=readOGR("/home/irpi/Dropbox/CNR/ITC/Test Landslide area/LSPoints.shp")
CM=readOGR("/home/irpi/Dropbox/CNR/ITC/catchements/intersection_medium_valid_last_v2.shp")

library(GISTools)
lsdcount=poly.counts(lsd,CM)

CM.df = as(CM, "data.frame")
CM.df$ID = as.numeric(CM.df$ID)
dflsd = data.frame(LSD=lsdcount)

wh = matrix(c(23.49322153090000143,-23.49322153090000143#,
              # 23.49322153090000143,-23.49322153090000143,
              # 23.49322153090000143,-23.49322153090000143,
              # 83.69415908030080686,-113.3896980265187295,
              # 23.69018716713384265,-31.10176606334843186,
              # 89.98329621380000276,-89.98329621380000276,
              # 23.6901871670999995,-23.6901871670999995,
              # 23.69018716713384265,-23.69018716632081478
              ),ncol=2,byrow=TRUE)
colnames(wh) = c("width","height")
rownames(wh) = c("S")
u = as.table(wh)

df=data.frame(LSD=dflsd$LSD,surf=CM.df$area, CM.df, Litho=CM.df$Rt4)
ix=which(colnames(df) %in% c("Rt4","area"))
df=df[,-ix]
df$surfM=df$S_count*abs(u['S','width']*u['S','height'])#/df$surf
df$area=df$surfM/df$surf
df=merge(df, df.LU, by = "ID")
df=na.omit(df)

save(df,file='data_df.Rdata')

vars2scale=c("PGA_mean","PGA_stdev","TnC_mean","TnC_stdev","PrC_mean","PrC_stdev",
             "Rlf_mean","Rlf_stdev","S_mean","S_stdev","Rn_mean",
             "Rn_stdev","surfM","area")
#lm(originale~scalato,data=scalare)
#y=27.165+8.236*x#smaean
data_scaled=df
data_scaled[,vars2scale]=apply(data_scaled[,vars2scale],2,scale)
data_scaled=na.omit(data_scaled)
data_scaled$LSD[data_scaled$LSD>0]=1
save(data_scaled,file='data_scaled.Rdata')
