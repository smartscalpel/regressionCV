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
path<-'/Users/lptolik/Dropbox/Скальпель/DBData/regression/RegReports/'

#' Get peak file names from peak files directory.
#'
#' @return
get_peaks<-function(){
  stop('Not implemented yet')
}

#' Prepare feature matrix.
#' Method load files named in the peaks parameter, and convert them into feature matrix.
#' Add all annotation required for the further data analysis including resolution, mode, 
#' exp_setup, diagnosis, percentage of norm tissue etc.
#'
#' @param peaks -- list of peak files to be converted into feature matrix
#'
prepare_feature_matrix<-function(peaks){
  stop('Not implemented yet')
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