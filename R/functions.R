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
  cat(format(Sys.time(), "%b %d %X"),'Function: get_peaks("',dpath,'") starts.\n')
  fl<-dir(path = dpath,pattern = '*.peak.rds')
  fl<-fl[grep('diag_32',fl,inver=TRUE)]
  idx32<-sapply(fl,function(.x)file.exists(paste0(dpath,sub('diag_[0-9]+\\.','diag_32.',.x))))
  fl<-fl[idx32]
  if(grepl('/$',dpath)){
    p<-dpath
  }else{
    p<-paste0(dpath,'/')
  }
  cat(format(Sys.time(), "%b %d %X"),'Function: get_peaks("',dpath,'") finish.\n')
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
prepare_feature_matrix<-function(peaks,norm_shift=0){
  cat(format(Sys.time(), "%b %d %X"),'Function: prepare_feature_matrix("',peaks,'",',norm_shift,') starts.\n')
  # n<-peaks[[1]]
  # cat('prepare_feature_matrix',n)
  # d<-data.frame(name=n,MZ_1=rnorm(10),MZ_2=2*rnorm(10))
  # #write.csv(d,paste0(n,'.csv'))
  # return(d)
    getMD<-function(p){
      as.data.frame(metaData(p))
    }
    norm<-sub('diag_[0-9]+\\.','diag_32.',peaks)
    l<-lapply(list(peaks,norm),readRDS)
    peaksL<-do.call(c,l)
    dl<-lapply(peaksL, getMD)
    md<-do.call(rbind,dl)
    md$norm.p<-as.numeric(as.character(md$norm.p))
    md$norm.p[md$diagnosis==32]<- md$norm.p[md$diagnosis==32]-norm_shift
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
    cat(format(Sys.time(), "%b %d %X"),'Function: prepare_feature_matrix("',peaks,'",',norm_shift,') finish.\n')
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
  cat(format(Sys.time(), "%b %d %X"),'Function: normalize(fm,"',normtype,'") starts.\n')
  cat('dim(fm)=',dim(fm),'\n')
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
  cat(format(Sys.time(), "%b %d %X"),'Function: normalize(fm,"',normtype,'") finish.\n')
  return(cbind(mdt,mz))
}

paretoscale<-function(mz){
  cat(format(Sys.time(), "%b %d %X"),'Function: paretoscale',' starts.\n')
  s<-apply(mz,2,sd)
  a<-apply(mz,2,mean)
  return(scale(mz,center = a,scale = sqrt(s)))
}

get_mdt<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_mdt("',fm$fname[1],'") starts.\n')
  mdt<-fm %>% dplyr::select(spectrumid,patientid,diagnosis,t.id,smpl.id) %>% unique
  cat(format(Sys.time(), "%b %d %X"),'Function: get_mdt("',fm$fname[1],'") finish.\n')
  return(mdt)
}

groups<-factor(c('train','test'))
smpl_split_fm<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: smpl_split_fm("',fm$fname[1],'") starts.\n')
  mdt<-get_mdt(fm)
  trainIndexSmpl <- createDataPartition(mdt$smpl.id, p = .6,
                                    list = FALSE,
                                    times = 1)
  test_smpl<-mdt$smpl.id[-trainIndexSmpl]
  fm$grp<-groups[1]
  fm$grp[fm$smpl.id %in% test_smpl]<-groups[2]
  cat(format(Sys.time(), "%b %d %X"),'Function: smpl_split_fm("',fm$fname[1],'") finish.\n')
  return(fm)
}

train_model<-function(fm,modeltype){
  cat(format(Sys.time(), "%b %d %X"),'Function: train_model("',fm$fname[1],'","',modeltype,'") starts.\n')
  idx<-grep("(MZ_.*|norm.p)",names(fm))
  train<-fm[fm$grp==groups[1],idx]
  train<-train[sample.int(dim(train)[1],size = smpl),]
  res<-switch (modeltype,
    rf=train_rf(train),
    xgb=train_xgb(train)
  )
  cat(format(Sys.time(), "%b %d %X"),'Function: train_model("',fm$fname[1],'","',modeltype,'") finish.\n')
  return(res)
}

smpl<-10

test_model<-function(fm,model){
  cat(format(Sys.time(), "%b %d %X"),'Function: test_model("',fm$fname[1],'","',model$method,'") starts.\n')
  idx<-grep("(MZ_.*|norm.p)",names(fm))
  test<-fm[,idx]
  res<-predict(model,newdata=test)
  fm$predict<-res
  fm$method<-model$method
  cat(format(Sys.time(), "%b %d %X"),'Function: test_model("',fm$fname[1],'","',model$method,'") finish.\n')
  return(fm)
}

plot_test<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_test("',fm$fname[1],'","',fm$method[1],'") starts.\n')
  test<-fm[fm$grp==groups[1],]
  
  my.formula <- y ~ x
  p <- ggplot(data = test, aes(x = norm.p, y = predict)) +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")),
                 parse = TRUE) +
    geom_point()
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_test("',fm$fname[1],'","',fm$method[1],'") finish\n')
  return(p)
}
train_rf<-function(train){
  cat(format(Sys.time(), "%b %d %X"),'Function: train_rf',' starts.\n')
  fitCV10<-trainControl(method = "repeatedcv",
                        number = 10,
                        repeats = 10)
  if(!exists('ncores')){
    ncores<- detectCores()
  }
  cat(format(Sys.time(), "%b %d %X"),'ncores=',ncores,'.\n')
  cl <- makePSOCKcluster(ncores)
  registerDoParallel(cl)
  rfFitCVpat <- train(norm.p ~ ., data = train,
                      method = "rf",
                      trControl = fitCV10,
                      verbose = FALSE)
  stopCluster(cl)
  cat(format(Sys.time(), "%b %d %X"),'Function: train_rf',' finish.\n')
  return(rfFitCVpat)
}

train_xgb<-function(train){
  cat(format(Sys.time(), "%b %d %X"),'Function: train_xgb',' starts.\n')
  fitCV10<-trainControl(method = "repeatedcv",
                        number = 10,
                        repeats = 10)
  if(!exists('ncores')){
    ncores<- detectCores()
  }
  cat(format(Sys.time(), "%b %d %X"),'ncores=',ncores,'.\n')
  cl <- makePSOCKcluster(ncores)
  registerDoParallel(cl)
  xboostFitCVspec <- train(norm.p ~ ., data = train,
                           method = "xgbDART",
                           trControl = fitCV10,
                           verbose = FALSE)
  stopCluster(cl)
  cat(format(Sys.time(), "%b %d %X"),'Function: train_xgb',' finish.\n')
  return(xboostFitCVspec)
}

train_trigger<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: train_trigger',' starts.\n')
  return(length(unique(fm$norm.p))>2)
}
#' Prepare panel of three PCA plots: 1-2, 2-3, 1-3
#'
#' @param fm feature matrix to plot
#' @param color name of the column to color plot with
#'
plot_pca<-function(fm, color){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_pca','\n')
  data=fm[[1]]
  stop('Not implemented yet')
}



get_model_fname<-function(ms_setup,diag,expt,method,idx){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_model_fname','\n')
  res=ms_setup[1]
  mode=ms_setup[2]
  mz=ms_setup[3]
  dev=ifelse(expt==1,2,4)
  fname<-sprintf('peak2019.diag_%d.expt_%d.res_%d.mode_%d.dev_%d.mz_%d.%s.cv10.%s.fmodel.rds',
                 diag,expt,res,mode,dev,mz,method,idx)
  return(fname)
}

get_fm_fname<-function(ms_setup,ddiag,dexpt){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_fm_fname','\n')
  res=ms_setup[1]
  mode=ms_setup[2]
  mz=ms_setup[3]
  dev=ifelse(dexpt==1,2,4)
  fname<-sprintf('peak2019.diag_%d.expt_%d.res_%d.mode_%d.dev_%d.mz_%d.fm.rds',
                 ddiag,dexpt,res,mode,dev,mz)
  return(fname)
}

get_model<-function(ms_setup,diag,expt,method,idx){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_model','\n')
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
  cat(format(Sys.time(), "%b %d %X"),'Function: mapMZ','\n')
  modMZ<-as.numeric(sub('^MZ_','',model$xNames))
  dataMZ<-as.numeric(sub('^MZ_','',names(fm)[grep('^MZ',names(fm))]))
}

predict_dataset<-function(model,ms_setup,ddiag,dexpt){
  cat(format(Sys.time(), "%b %d %X"),'Function: predict_dataset','\n')
  if(length(mode)==0) return(data.frame())
  m<-model$model
  fm<-load_dataset(model,ms_setup,ddiag,dexpt)
  if(dim(fm)[1]==0) return(data.frame())
  fmapped<-mapMZ(m,fm)
}

load_dataset<-function(ms_setup,ddiag,dexpt){
  cat(format(Sys.time(), "%b %d %X"),'Function: load_dataset','\n')
  fmfname<-get_fm_fname(ms_setup,ddiag,dexpt)
  fpath<-paste0(path,fmfname)
  if(file.exists(fpath)){
    fm<- readRDS(fpath)
  }else{
    return(data_frame())
  }
}

get_pat_df<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_pat_df','\n')
  if(dim(fm)[2]>2){
  return(unique(fm[,c("patientid","diagnosis")]))
  }else{
    return(data.frame())
  }
}

get_spec_df<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_spec_df','\n')
  if(dim(fm)[2]>2){
    return(unique(fm[,c("spectrumid","patientid","diagnosis")]))
  }else{
    return(data.frame())
  }
}

get_dim<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_dim','\n')
  return(data.frame(nrow=dim(fm)[1],ncol=dim(fm)[2]))
}
