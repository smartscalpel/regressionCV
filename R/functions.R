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
source('./myRF.R')

#path<-'~/Downloads/peak2019.full/'
dpath<-'/Users/lptolik/Documents/Projects/MSpeaks/data/regression/'
#dpath<-'/Users/lptolik/Dropbox/Скальпель/DBData/regression/'
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
  cat(format(Sys.time(), "%b %d %X"),'Function: prepare_feature_matrix("',peaks,'",',norm_shift,') starts. Mem:',getFreeMem(),'GB\n')
  # n<-peaks[[1]]
  # cat('prepare_feature_matrix',n)
  # d<-data.frame(name=n,MZ_1=rnorm(10),MZ_2=2*rnorm(10))
  # #write.csv(d,paste0(n,'.csv'))
  # return(d)
    getMD<-function(p){
      as.data.frame(metaData(p))
    }
    getRDS<-function(f){
      cat(format(Sys.time(), "%b %d %X"),'Function: getRDS("',f,'") starts. Mem:',getFreeMem(),'GB\n')
      res<-try(readRDS(f))
      cat(format(Sys.time(), "%b %d %X"),class(res),'.\n')
      if(inherits(res, "try-error")){
        return(list())
      }else{
        return(res)
      }
    }
    norm<-sub('diag_[0-9]+\\.','diag_32.',peaks)
    idx<-sapply(norm,file.exists)
    l<-lapply(c(peaks,norm[idx]),getRDS)
    peaksL<-do.call(c,l)
    dl<-lapply(peaksL, getMD)
    md<-do.call(rbind,dl)
    md$norm.p<-as.numeric(as.character(md$norm.p))
    md$tumor.p<-as.numeric(as.character(md$tumor.p))
    md$necro.p<-as.numeric(as.character(md$necro.p))
    if(grepl('othr.p',names(md))){
      md$othr.p<-as.numeric(as.character(md$othr.p))
      md$target<-md$norm.p+md$othr.p
    }else{
      md$target<-md$norm.p
    }
    md$target[md$diagnosis==32]<- md$target[md$diagnosis==32]+norm_shift
    md$fname<-basename(peaks[1])
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

normtypes<-factor(c('None'))#,'Autoscaling','Pareto'))
filtertypes<-c('None','ZVar','Corr')

feature_filter<-function(fm,ftype){
  cat(format(Sys.time(), "%b %d %X"),'Function: feature_filter("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',ftype,'") starts.\n')
  idx<-grep("MZ_.*",names(fm))
  features<-fm[,idx]
  mdt<-fm[,-idx]
  mdt$Filter<-ftype
  res<-switch (ftype,
               None=features,
               ZVar=filter_nzv(features),
               Corr=filter_corr(features)
  )
  res<-cbind(mdt,res)
  cat(format(Sys.time(), "%b %d %X"),'Function: feature_filter("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',ftype,'") finish.\n')
  return(res)
  
}

filter_nzv<-function(fm){
  nzv <- nearZeroVar(fm)
  res<-fm[,-nzv]
  return(res)
}
filter_corr<-function(fm,cutoff = .8){
  fm1<-filter_nzv(fm)
  descrCor <- cor(fm1)
  highlyCorDescr <- findCorrelation(descrCor, cutoff = cutoff)
  res<-fm1[,-highlyCorDescr]
  return(res)
}
#' Title
#'
#' @param fm feature matrix
#' @param normtype type of norm to apply
#'
#' @return normalized fm
normalize<-function(fm,normtype){
  cat(format(Sys.time(), "%b %d %X"),'Function: normalize(fm,"',normtype,'") starts. Mem:',getFreeMem(),'GB\n')
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
  cat(format(Sys.time(), "%b %d %X"),'Function: paretoscale',' starts. Mem:',getFreeMem(),'GB\n')
  s<-apply(mz,2,sd)
  a<-apply(mz,2,mean)
  return(scale(mz,center = a,scale = sqrt(s)))
}

get_mdt<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_mdt("',fm$fname[1],'") starts. Mem:',getFreeMem(),'GB\n')
  mdt<-fm %>% dplyr::select(spectrumid,patientid,diagnosis,t.id,smpl.id,target) %>% unique
  cat(format(Sys.time(), "%b %d %X"),'Function: get_mdt("',fm$fname[1],'") finish.\n')
  return(mdt)
}

groups<-factor(c('train','test'))
smpl_split_fm<-function(fm,split=0.6){
  cat(format(Sys.time(), "%b %d %X"),'Function: smpl_split_fm("',fm$fname[1],'","',as.character(fm$Norm[1]),'") starts. Mem:',getFreeMem(),'GB\n')
  mdt<-get_mdt(fm)
  smpl<-mdt %>% dplyr::select(smpl.id,target) %>% unique
  trainIndexSmpl <- createDataPartition(smpl$target, p = split,
                                    list = FALSE,
                                    times = 1)
  test_smpl<-smpl$smpl.id[-trainIndexSmpl]
  fm$grp<-groups[1]
  fm$grp[fm$smpl.id %in% test_smpl]<-groups[2]
  cat(format(Sys.time(), "%b %d %X"),'Function: smpl_split_fm("',fm$fname[1],'","',as.character(fm$Norm[1]),'") finish.\n')
  return(fm)
}

train_model<-function(fm,modeltype){
  cat(format(Sys.time(), "%b %d %X"),'Function: train_model("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$Filter[1],'","',modeltype,'") starts. Mem:',getFreeMem(),'GB\n')
  fm$was.trained<-0
  idx<-grep("(MZ_.*|target)",names(fm))
  trdx<-which(fm$grp==groups[1])
  if(smpl<length(trdx)){
  jdx<-trdx[sample.int(length(trdx),size = smpl)]
  }else{
    jdx<-trdx
  }
  fm$was.trained[jdx]<-1
  train<-fm[jdx,idx]
  cat(format(Sys.time(), "%b %d %X"),'train dataset',dim(train),'\n')
  res<-switch (modeltype,
    rf=train_rf(train),
    xgb=train_xgb(train)
  )
  cat(format(Sys.time(), "%b %d %X"),'Function: train_model("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$Filter[1],'","',modeltype,'") finish.\n')
  return(list(model=res,data=fm))
}

smpl<-5e6

test_model<-function(mod){
  fm<-mod$data
  model<-mod$model
  cat(format(Sys.time(), "%b %d %X"),'Function: test_model("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$Filter[1],'","',model$method,'") starts. Mem:',getFreeMem(),'GB\n')
  idx<-grep("(MZ_.*|target)",names(fm))
  test<-fm[,idx]
  res<-predict(model,newdata=test)
  fm$predict<-res
  fm$method<-model$method
  cat(format(Sys.time(), "%b %d %X"),'Function: test_model("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$Filter[1],'","',model$method,'") finish.\n')
  return(fm)
}

make_point_plot<-function(test){
  my.formula <- y ~ x
  p <- ggplot(data = test, aes(x = target, y = predict)) +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")),
                 parse = TRUE)
  return(p)
}
plot_test_point<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_test("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  test<-fm[fm$grp==groups[2],]
  p<-make_point_plot(test)+geom_point() +
    geom_jitter()
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_test("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}
plot_train_point<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_train("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  test<-fm[fm$was.trained==0,]
  p<-make_point_plot(test)+geom_point()+
    geom_jitter(aes(x = target, y = predict,color='blue'),data=fm[fm$was.trained==1,])
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_train("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}
plot_test_box<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_test("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  test<-fm[fm$grp==groups[2],]
  p<-make_point_plot(test)+geom_boxplot(aes(group=target))
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_test("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}
plot_train_box<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_train("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  test<-fm[fm$was.trained==0,]
  p<-make_point_plot(test)+geom_boxplot(aes(group=target))+
    geom_jitter(aes(x = target, y = predict,color='blue'),data=fm[fm$was.trained==1,])
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_train("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}

train_rf<-function(train){
  cat(format(Sys.time(), "%b %d %X"),'Function: train_rf',' starts. Mem:',getFreeMem(),'GB\n')
  fitCV10<-trainControl(method = "repeatedcv",
                        number = 10,
                        repeats = 3)
  N<-dim(train)[1]
  p<-dim(train)[2]-1
  #tunegrid <- expand.grid(.mtry=c(1:(p/3)),.ntree=c(500,1000,1500))
  tunegrid <- expand.grid(mtry=sample(c(1:(p/3)),size=10))#,.ntree=c(500,1000,1500))
  
  if(!exists('ncores')){
    ncores<- detectCores()
  }
  cat(format(Sys.time(), "%b %d %X"),'ncores=',ncores,'.\n')
  cl <- makePSOCKcluster(ncores)
  registerDoParallel(cl)
  rfFitCVpat <- train(target ~ ., data = train,
                      method = "rf",#customRF,
                      trControl = fitCV10,
                      tuneGrid=tunegrid, 
                      verbose = FALSE)
  stopCluster(cl)
  cat(format(Sys.time(), "%b %d %X"),'Function: train_rf',' finish.\n')
  return(rfFitCVpat)
}

train_xgb<-function(train){
  cat(format(Sys.time(), "%b %d %X"),'Function: train_xgb',' starts. Mem:',getFreeMem(),'GB\n')
  fitCV10<-trainControl(method = "repeatedcv",
                        number = 10,
                        repeats = 3)
  if(!exists('ncores')){
    ncores<- detectCores()
  }
  cat(format(Sys.time(), "%b %d %X"),'ncores=',ncores,'.\n')
  cl <- makePSOCKcluster(ncores)
  registerDoParallel(cl)
  xboostFitCVspec <- train(target ~ ., data = train,
                           method = "xgbDART",
                           trControl = fitCV10,
                           verbose = FALSE)
  stopCluster(cl)
  cat(format(Sys.time(), "%b %d %X"),'Function: train_xgb',' finish.\n')
  return(xboostFitCVspec)
}

train_trigger<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: train_trigger',' starts. Mem:',getFreeMem(),'GB\n')
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

get_fm_fname<-function(res,mode,mz,ddiag,path){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_fm_fname','\n')
  fpatt<-sprintf('peak2019.diag_%d.expt_.*.res_%d.mode_%d.dev_.*.mz_%d.peak.rds',
                 ddiag,res,mode,mz)
  fname<-dir(path = path,pattern = fpatt)
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

load_dataset<-function(res,mode,mz,diag){
  cat(format(Sys.time(), "%b %d %X"),'Function: load_dataset','\n')
  path<-dpath
  fmfname<-get_fm_fname(res,mode,mz,diag,path)
  fpath<-paste0(path,fmfname)
  return(fpath)
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
