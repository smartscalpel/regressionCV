library(drake)
library(ggplot2)
library(ggpmisc)
library(data.table)
library(dplyr)
library(xtable)
library(MALDIquant)
library(xgboost)
library(iml)
library(caret)
library(randomForest)
library(doParallel)

#path<-'~/Downloads/peak2019.full/'
dpath<-'/Users/lptolik/Dropbox/Скальпель/DBData/regression/'
#dpath<-'~/regression/'

#' Get peak file names from peak files directory.
#'
#' @return list of file names for peaks
get_peaks<-function(dpath){
  fl<-dir(path = dpath,pattern = '*.peak.rds')
  if(grepl('/$',dpath)){
    p<-dpath
  }else{
    p<-paste0(dpath,'/')
  }
  return(paste0(p,fl))
  # return(fl)
}

#' Prepare feature matrix.
#' Method load files named in the peaks parameter, and convert them into feature matrix.
#' Add all annotation required for the further data analysis including resolution, mode, 
#' exp_setup, diagnosis, percentage of norm tissue etc.
#'
#' @param peaks -- list of peak files to be converted into feature matrix
#'
prepare_feature_matrix<-function(peaks){
  # n<-peaks[[1]]
  # cat('prepare_feature_matrix',n)
  # d<-data.frame(name=n,MZ_1=rnorm(10),MZ_2=2*rnorm(10))
  # #write.csv(d,paste0(n,'.csv'))
  # return(d)
    getMD<-function(p){
      as.data.frame(metaData(p))
    }
    peaksL<-readRDS(peaks)
    dl<-lapply(peaksL, getMD)
    md<-do.call(rbind,dl)
    md$norm.p<-as.numeric(as.character(md$norm.p))
    md$tumor.p<-as.numeric(as.character(md$tumor.p))
    md$necro.p<-as.numeric(as.character(md$necro.p))
    md$fname<-basename(peaks)
    wf<-determineWarpingFunctions(peaksL,
                                  method="lowess",
                                  plot=FALSE,minFrequency=0.05)
    aPeaks<-warpMassPeaks(peaksL,wf)
    bPeaks <- binPeaks(aPeaks, tolerance=2e-3)
    fpeaks <- filterPeaks(bPeaks,
                          labels=md$diag,
                          minFrequency=0.25, mergeWhitelists=TRUE)
    featureMatrix <- intensityMatrix(fpeaks)
    idNA<-which(is.na(featureMatrix),arr.ind =TRUE)
    featureMatrix[idNA]<-0
    colnames(featureMatrix)<-paste0('MZ_',round(as.numeric(colnames(featureMatrix)),3))
    return(cbind(md,featureMatrix))
}

normtypes<-factor(c('None','Autoscaling','Pareto'))

#' Title
#'
#' @param fm feature matrix
#' @param normtype type of norm to apply
#'
#' @return normalized fm
normalize<-function(fm,normtype){
  cat(length(fm),'\n')
  #fm<-fml[[1]]
  #cat(names(fm),'\n')
  cidx<-grep('MZ_.*',names(fm))
  #cat(cidx,'\n')
  mdt<-fm[,-cidx]
  mz<-fm[,cidx]
  #cat(normtype,dim(mdt),dim(mz),'\n')
  mdt$Norm<-normtype
  cat(unique(as.character(mdt$Norm)),unique(mdt$fname),'\n')
  mz<-switch(as.character(normtype),
  None=mz,
  Autoscaling=scale(mz),
  Pareto=paretoscale(mz))
  return(cbind(mdt,mz))
}

paretoscale<-function(mz){
  s<-apply(mz,2,sd)
  a<-apply(mz,2,mean)
  return(scale(mz,center = a,scale = sqrt(s)))
}

get_mdt<-function(fm){
  mdt<-fm %>% dplyr::select(spectrumid,patientid,diagnosis,t.id,smpl.id) %>% unique
  return(mdt)
}

groups<-factor(c('train','test'))
smpl_split_fm<-function(fm){
  mdt<-get_mdt(fm)
  trainIndexSmpl <- createDataPartition(mdt$smpl.id, p = .6,
                                    list = FALSE,
                                    times = 1)
  test_smpl<-mdt$smpl.id[-trainIndexSmpl]
  fm$grp<-groups[1]
  fm$grp[fm$smpl.id %in% test_smpl]<-groups[2]
  return(fm)
}

train_model<-function(fm,modeltype){
  idx<-grep("(MZ_.*|norm.p)",names(fm))
  train<-fm[fm$grp==groups[1],idx]
  train<-train[sample.int(dim(train)[1],size = smpl),]
  res<-switch (modeltype,
    rf=train_rf(train),
    xgb=train_xgb(train)
  )
  return(res)
}

smpl<-100

train_rf<-function(train){
  fitCV10<-trainControl(method = "repeatedcv",
                        number = 10,
                        repeats = 10)
  # cl <- makePSOCKcluster(ncores)
  # registerDoParallel(cl)
  rfFitCVpat <- train(norm.p ~ ., data = train,
                      method = "rf",
                      trControl = fitCV10,
                      verbose = FALSE)
  #stopCluster(cl)
  return(rfFitCVpat)
}

train_xgb<-function(train){
  fitCV10<-trainControl(method = "repeatedcv",
                        number = 10,
                        repeats = 10)
  # cl <- makePSOCKcluster(ncores)
  # registerDoParallel(cl)
  xboostFitCVspec <- train(norm.p ~ ., data = train,
                           method = "xgbDART",
                           trControl = fitCV10,
                           verbose = FALSE)
  #stopCluster(cl)
  return(xboostFitCVspec)
}

train_trigger<-function(fm){
  return(length(unique(fm$norm.p))>2)
}
#' Prepare panel of three PCA plots: 1-2, 2-3, 1-3
#'
#' @param fm feature matrix to plot
#' @param color name of the column to color plot with
#'
plot_pca<-function(fm, color){
  data=fm[[1]]
  stop('Not implemented yet')
}



get_model_fname<-function(ms_setup,diag,expt,method,idx){
  res=ms_setup[1]
  mode=ms_setup[2]
  mz=ms_setup[3]
  dev=ifelse(expt==1,2,4)
  fname<-sprintf('peak2019.diag_%d.expt_%d.res_%d.mode_%d.dev_%d.mz_%d.%s.cv10.%s.fmodel.rds',
                 diag,expt,res,mode,dev,mz,method,idx)
  return(fname)
}

get_fm_fname<-function(ms_setup,ddiag,dexpt){
  res=ms_setup[1]
  mode=ms_setup[2]
  mz=ms_setup[3]
  dev=ifelse(dexpt==1,2,4)
  fname<-sprintf('peak2019.diag_%d.expt_%d.res_%d.mode_%d.dev_%d.mz_%d.fm.rds',
                 ddiag,dexpt,res,mode,dev,mz)
  return(fname)
}

get_model<-function(ms_setup,diag,expt,method,idx){
  fname<-get_model_fname(ms_setup,diag,expt,method,idx)
  fpath<-paste0(path,fname)
  if(file.exists(fpath)){
    m<- readRDS(fpath)
    return(list(model=m,ms_setup=ms_setup,diag=diag,expt=expt,method=method,idx=idx))
  }else{
    return(list())
  }
}

mapMZ<-function(model,fm){
  modMZ<-as.numeric(sub('^MZ_','',model$xNames))
  dataMZ<-as.numeric(sub('^MZ_','',names(fm)[grep('^MZ',names(fm))]))
}

predict_dataset<-function(model,ms_setup,ddiag,dexpt){
  if(length(mode)==0) return(data.frame())
  m<-model$model
  fm<-load_dataset(model,ms_setup,ddiag,dexpt)
  if(dim(fm)[1]==0) return(data.frame())
  fmapped<-mapMZ(m,fm)
}

load_dataset<-function(ms_setup,ddiag,dexpt){
  fmfname<-get_fm_fname(ms_setup,ddiag,dexpt)
  fpath<-paste0(path,fmfname)
  if(file.exists(fpath)){
    fm<- readRDS(fpath)
  }else{
    return(data_frame())
  }
}

get_pat_df<-function(fm){
  if(dim(fm)[2]>2){
  return(unique(fm[,c("patientid","diagnosis")]))
  }else{
    return(data.frame())
  }
}

get_spec_df<-function(fm){
  if(dim(fm)[2]>2){
    return(unique(fm[,c("spectrumid","patientid","diagnosis")]))
  }else{
    return(data.frame())
  }
}

get_dim<-function(fm){
  return(data.frame(nrow=dim(fm)[1],ncol=dim(fm)[2]))
}
