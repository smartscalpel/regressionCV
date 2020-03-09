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

path<-'~/Downloads/peak2019.full/'

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
  dev=ifelse(expt==1,2,4)
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
  
}

predict_dataset<-function(model,ms_setup,ddiag,dexpt){
  if(length(mode)==0) return(data.frame())
  m<-model$model
  fm<-load_dataset(model,ms_setup,ddiag,dexpt)
  if(dim(fm)[1]==0) return(data.frame())
  fmapped<-mapMZ(m,fm)
}
load_dataset<-function(model,ms_setup,ddiag,dexpt){
  fmfname<-get_fm_fname(ms_setup,ddiag,dexpt)
  fpath<-paste0(path,fmfname)
  if(file.exists(fpath)){
    fm<- readRDS(fpath)
  }else{
    return(data_frame())
  }
}