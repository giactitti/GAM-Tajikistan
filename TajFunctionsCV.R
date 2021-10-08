#########################################################################################
crossValidation = function(data_scaled,Vpercentuale,bin) {
  library(Hmisc)
  
  data.mat4=data_scaled
  perc.test=Vpercentuale#%
  pos = data.mat4[bin == 1, ]
  neg = data.mat4[bin == 0, ]
  myset = 1:nrow(neg)
  
  p=length(pos[,1])
  p.partition=c(rep(ceiling(p/perc.test),perc.test-1),(p-(ceiling(p/perc.test))*(perc.test-1)))
  idxlist.p=partition.vector(sample(1:p,size=p),p.partition)
  count=1
  for(i in idxlist.p){
    for(ii in 1:length(i)){
      i[ii]=which(data.mat4$ID==pos[i[ii],c('ID')])
    }
    idxlist.p[[count]]=i
    count=count+1
  }
  n=length(neg[,1])
  n.partition=c(rep(ceiling(n/perc.test),perc.test-1),(n-(ceiling(n/perc.test))*(perc.test-1)))
  idxlist.n=partition.vector(sample(1:n,size=n),n.partition)
  count=1
  for(i in idxlist.n){
    for(ii in 1:length(i)){
      i[ii]=which(data.mat4$ID==neg[i[ii],c('ID')])
    }
    idxlist.n[[count]]=i
    count=count+1
  }
  
  idxlist=Map(c, idxlist.p, idxlist.n)
  
  
  return(list("idxlist"=idxlist, "idxlist.p"=idxlist.p,"idxlist.n"=idxlist.n))
  #return(idxlist)
}

###################################################################################
kfold = function(data_scaled,Vpercentuale,idx,bin,nonlinearclasses,name,formulaINLA){
  Northridge.cv=c()
  #browser()
  #ROC10=c()
  library(pROC)
  fix.cv=c()
  random.cv=c()
  for(i in 1:Vpercentuale){
    idx.holdout=idx[[i]]
    tmp.y=bin
    tmp.y[idx.holdout]=NA
    tmp.dataset = data.frame(cbind(intercept = 1, Ntrials = 1, y = tmp.y, data_scaled, nonlinearclasses))
    #tmp.dataset=na.omit(tmp.dataset)
    #browser()
    hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))
    hyper.iid = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))
    
    formula=formulaINLA
    
    CV.runs = inla(formula = formula, family = "binomial", data=tmp.dataset,
                   #c(as.list(tmp.dataset), list(Dataset1=Dataset1)),####controllare???????
                   control.fixed=list(prec=.1),
                   num.threads = 2, 
                   Ntrials = Ntrials,
                   control.family = list(control.link = list(model = "logit")), 
                   control.mode = list(restart = T),
                   control.inla = list(int.strategy = "eb"), 
                   control.predictor = list(compute = T, link = 1), verbose = T)
    
    Northridge.cv[idx.holdout]=CV.runs$summary.fitted.values$mean[idx.holdout]
    fix.cv=append(fix.cv,list(CV.runs$summary.fixed))
    random.cv=append(random.cv,list(CV.runs$summary.random))
    #random.cv=append(random.cv,tmp.dataset)
    #print(CV.runs$summary.fixed)
    #browser()
  }
  save(Northridge.cv,file=paste("CV.runs",toString(name),".Rdata",sep=''))
  return(list("matt"=Northridge.cv, "cvruns"=fix.cv,"random"=random.cv,"tot"=CV.runs))
  #return(Northridge.cv,CV.runs$summary.fixed$mean)
  #browser()
}

############################################################################save ROC
froc=function(Vpercentuale,idx,bin,summary,name){
  ROC=list()
  for(i in 1:Vpercentuale){
    idx.holdout=idxlist[[i]]
    tmp.y=bin
    #ROCtest=roc(tmp.y[-idx.holdout]~Northridge.cv[-idx.holdout])
    #save(ROCtest,file=paste('CV_ROCtest',toString(i),'.Rdata'))
    ROCvalid=roc(tmp.y[idx.holdout]~summary[idx.holdout])
    ROC=append(ROC,list(ROCvalid))
    #ROC=cbind(ROC,ROCvalid)
  }
  save(ROC,file=paste(name,'ROCvalid.Rdata'))
  return(ROC)
}