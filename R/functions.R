library(drake)
library(ggplot2)
library(plyr)
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
library(SHAPforxgboost)
#source('./myRF.R')

#path<-'~/Downloads/peak2019.full/'
dpath<-'/Users/lptolik/Documents/Projects/MSpeaks/dataN2TIC/regression/'
#dpath<-'/Users/lptolik/Dropbox/Скальпель/DBData/regression/'
#dpath<-'~/regression/'

getFreeMem<-function(){
  #as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE))/1e6
  return(0)
}
#' Get peak file names from peak files directory.
#'
#' @return list of file names for peaks
get_peaks<-function(dpath){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_peaks("',dpath,'") starts. Mem:',getFreeMem(),'GB\n')
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
prepare_feature_matrix<-function(peaks,norm_shift=0,monoisotopic=FALSE,size=3L:10L){
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
  if(all(grepl('res_2',peaks))){
    tol=5e-4
  }else{
    tol=5e-5
  }
  if(monoisotopic){
    if(all(grepl('mode_2',peaks))){
      K=TRUE
      Cl=FALSE
    }else{
      K=FALSE
      Cl=TRUE
    }
    peaksL<-myPeakList(peaksL, minCor=0.95, tolerance=tol, size=size,Cl=Cl,K=K) 
  }
  dl<-lapply(peaksL, getMD)
  md<-do.call(rbind,dl)
  md$norm.p<-as.numeric(as.character(md$norm.p))
  md$tumor.p<-as.numeric(as.character(md$tumor.p))
  md$necro.p<-as.numeric(as.character(md$necro.p))
  if(any(grepl('othr.p',names(md)))){
    md$othr.p<-as.numeric(as.character(md$othr.p))
    md$target<-md$norm.p+md$othr.p
  }else{
    md$othr.p<-0
    md$target<-md$norm.p
  }
  md$target[md$diagnosis==32]<- md$target[md$diagnosis==32]+norm_shift
  md$fname<-basename(peaks[1])
  wf<-determineWarpingFunctions(peaksL,
                                method="lowess",
                                plot=FALSE,minFrequency=0.05)
  aPeaks<-warpMassPeaks(peaksL,wf)
  bPeaks <- binPeaks(aPeaks,  method="strict",tolerance=tol)
  fpeaks <- filterPeaks(bPeaks,
                        labels=md$diag,
                        minFrequency=0.25, mergeWhitelists=TRUE)
  featureMatrix <- intensityMatrix(fpeaks)
  idNA<-which(is.na(featureMatrix),arr.ind =TRUE)
  featureMatrix[idNA]<-0
  colnames(featureMatrix)<-paste0('MZ_',round(as.numeric(colnames(featureMatrix)),3))
  fm<-cbind(md,featureMatrix)
  tot<-fm$norm.p+fm$tumor.p+fm$necro.p+fm$othr.p
  idx<-which(tot>=100)
  fm<-fm[idx,]
  cat(format(Sys.time(), "%b %d %X"),'feature matrix',dim(featureMatrix),'\n')
  cat(format(Sys.time(), "%b %d %X"),'filtered matrix',dim(fm),'\n')
  cat(format(Sys.time(), "%b %d %X"),'Function: prepare_feature_matrix("',peaks,'",',norm_shift,') finish.\n')
  return(fm)
}

#normtypes<-factor(c('None'))#,'Pareto','Autoscaling'))
#filtertypes<-c('None')#,'ZVar','Corr')

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
  cat(format(Sys.time(), "%b %d %X"),'Function: filter_nzv; dim(fm)=',dim(fm),'\n')
  nzv <- nearZeroVar(fm)
  if(length(nzv)>0){
    res<-fm[,-nzv]
  }else{
    res<-fm
  }
  cat(format(Sys.time(), "%b %d %X"),'Function: filter_nzv; length(nzv)=',length(nzv),'; dim(res)=',dim(res),'\n')
  return(res)
}
filter_corr<-function(fm,cutoff = .8){
  cat(format(Sys.time(), "%b %d %X"),'Function: filter_corr; dim(fm)=',dim(fm),'\n')
  fm1<-filter_nzv(fm)
  descrCor <- cor(fm1)
  highlyCorDescr <- findCorrelation(descrCor, cutoff = cutoff)
  if(length(highlyCorDescr)>0){
    res<-fm1[,-highlyCorDescr]
  }else{
    res<-fm1
  }
  cat(format(Sys.time(), "%b %d %X"),'Function: filter_corr; length(highlyCorDescr)=',length(highlyCorDescr),'; dim(res)=',dim(res),'\n')
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

eval_model<-function(tst){
  return(postResample(tst[!tst$was.trained, "predict"], tst[!tst$was.trained, "target"]))
}

apply_model<-function(mod,newdata){
  fm<-newdata
  model<-mod$model
  cat(format(Sys.time(), "%b %d %X"),'Function: test_model("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$Filter[1],'","',model$method,'") starts. Mem:',getFreeMem(),'GB\n')
  idx<-match(model$finalModel$xNames,names(fm))
  if(any(is.na(idx))){
    i<-which(is.na(idx))
    stop('Model parameters [',model$finalModel$xNames[i],'] are missing from the dataset.\n')
  }
  test<-fm[,c('target',names(fm)[idx])]
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

#### TCP plots ####
plot_tcp_box<-function(fm,theme=theme_grey()){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_tcp_train("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  fmtcp<-fm
  fmtcp$tcp<-100-fm$target
  fmtcp$predict<-100-fm$predict
  test<-fmtcp[fmtcp$was.trained==0,]
  p<-make_point_plot(test)+geom_boxplot(aes(group=target))+
    geom_jitter(aes(x = target, y = predict,color='blue'),data=fmtcp[fmtcp$was.trained==1,])+
    xlab('Assigned TCP')+ylab('Predicted TCP')+theme
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_tcp_train("',fmtcp$fname[1],'","',
      as.character(fmtcp$Norm[1]),'","',fmtcp$method[1],'") finish\n')
  return(p)
}

smpl_box_theme<-theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size=20),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = c(0.99, 0.5),
    legend.position = c(0.99, 0.5))

plot_tcp_smpl_box<-function(fm,theme=theme_grey(base_size=14),palette=NULL){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_tcp_smpl_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  fmtcp<-fm
  fmtcp$target<-100-fm$target
  fmtcp$predict<-100-fm$predict
  sum_tst<-ddply(fmtcp,.(smpl.id,target,was.trained),
                 summarise,min.pred=min(predict),mean.pred=mean(predict),
                 max.pred=max(predict))
  fmtcp$samples<-factor(fmtcp$smpl.id,levels = sum_tst$smpl.id[order(sum_tst$target)])
  fmtcp$Set<-factor(fmtcp$was.trained,labels=c('Validation','Train'))
  sum_tst$samples<-factor(sum_tst$smpl.id,levels = sum_tst$smpl.id[order(sum_tst$target)])
  p <- ggplot(fmtcp, aes(x=samples, y=predict,color=Set)) +
    geom_boxplot()+
    geom_point(data=sum_tst,
               aes(x=samples, y=target),
               color='black',
               shape='+',size=6)+
    coord_flip()+ylab('Predicted TCP')+xlab('Sample ID')
  p<-p+theme
  if(!is.null(palette)){
    p<-p+scale_colour_brewer(type = "seq", palette =palette)
  }
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_tcp_smpl_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}

plot_median_spectrum_tcp_box<-function(fm,theme=theme_grey(base_size=26)){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_median_spectrum_tcp_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  fmtcp<-fm
  fmtcp$target<-100-fm$target
  fmtcp$predict<-100-fm$predict
  spec_tst<-make_mean_spectrum(fmtcp)
  test<-spec_tst[spec_tst$was.trained==0,]
  p<-make_median_spectrum_point_plot(test)+geom_boxplot(aes(group=target))+
    geom_jitter(aes(x = target, y = median.pred),color='blue', size=2,
                data=spec_tst[spec_tst$was.trained==1,])+
    xlab('Assigned TCP')+ylab('Median of predicted TCP')+theme
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_median_spectrum_tcp_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}

plot_mean_spectrum_tcp_box<-function(fm,theme=theme_grey(base_size=26)){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_mean_spectrum_tcp_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  fmtcp<-fm
  fmtcp$target<-100-fm$target
  fmtcp$predict<-100-fm$predict
  spec_tst<-make_mean_spectrum(fmtcp)
  test<-spec_tst[spec_tst$was.trained==0,]
  p<-make_mean_spectrum_point_plot(test)+geom_boxplot(aes(group=target))+
    geom_jitter(aes(x = target, y = mean.pred),color='blue',size=2,
                data=spec_tst[spec_tst$was.trained==1,])+
    xlab('Assigned TCP')+ylab('Mean of predicted TCP')+theme
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_mean_spectrum_tcp_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}

#### Spectrum aggregation #####
make_mean_spectrum<-function(tst){
  spec_tst<-ddply(tst,.(spectrumid,patientid,diagnosis,smpl.id,t.id,norm.p,tumor.p,necro.p,norm.type,othr.p,target,was.trained),summarise,min.pred=min(predict),mean.pred=mean(predict),median.pred=median(predict),max.pred=max(predict))
  return(spec_tst)
}
#### Spectrum median plot #####
make_median_spectrum_point_plot<-function(test){
  my.formula <- y ~ x
  p <- ggplot(data = test, aes(x = target, y = median.pred)) +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")),
                 parse = TRUE)
  return(p)
}
plot_median_spectrum_test_box<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_median_spectrum_test_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  spec_tst<-make_mean_spectrum(fm)
  test<-spec_tst[spec_tst$was.trained==0,]
  p<-make_median_spectrum_point_plot(test)+geom_boxplot(aes(group=target))
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_median_spectrum_test_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}
plot_median_spectrum_train_box<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_median_spectrum_train_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  spec_tst<-make_mean_spectrum(fm)
  test<-spec_tst[spec_tst$was.trained==0,]
  p<-make_median_spectrum_point_plot(test)+geom_boxplot(aes(group=target))+
    geom_jitter(aes(x = target, y = median.pred),color='blue',data=spec_tst[spec_tst$was.trained==1,])
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_median_spectrum_train_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}

#### Spectrum mean plot #####
make_mean_spectrum_point_plot<-function(test){
  my.formula <- y ~ x
  p <- ggplot(data = test, aes(x = target, y = mean.pred)) +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")),
                 parse = TRUE)
  return(p)
}
plot_mean_spectrum_test_box<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_mean_spectrum_test_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  spec_tst<-make_mean_spectrum(fm)
  test<-spec_tst[spec_tst$was.trained==0,]
  p<-make_mean_spectrum_point_plot(test)+geom_boxplot(aes(group=target))
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_mean_spectrum_test_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}
plot_mean_spectrum_train_box<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_mean_spectrum_train_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  spec_tst<-make_mean_spectrum(fm)
  test<-spec_tst[spec_tst$was.trained==0,]
  p<-make_mean_spectrum_point_plot(test)+geom_boxplot(aes(group=target))+
    geom_jitter(aes(x = target, y = mean.pred),color='blue',data=spec_tst[spec_tst$was.trained==1,])
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_mean_spectrum_train_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
  return(p)
}


plot_train_smpl_box<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_train_smpl_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") starts. Mem:',getFreeMem(),'GB\n')
  sum_tst<-ddply(fm,.(smpl.id,target,was.trained),summarise,min.pred=min(predict),mean.pred=mean(predict),max.pred=max(predict))
  fm$samples<-factor(fm$smpl.id,levels = sum_tst$smpl.id[order(sum_tst$target)])
  sum_tst$samples<-factor(sum_tst$smpl.id,levels = sum_tst$smpl.id[order(sum_tst$target)])
  p <- ggplot(fm, aes(x=samples, y=predict,color=factor(was.trained))) +
    geom_boxplot()+
    geom_point(data=sum_tst,
               aes(x=samples, y=target),
               color='black',
               shape='+',size=6)+
    coord_flip()
  cat(format(Sys.time(), "%b %d %X"),'Function: plot_train_smpl_box("',fm$fname[1],'","',as.character(fm$Norm[1]),'","',fm$method[1],'") finish\n')
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

xgb_importance<-function(mod){
  cat(format(Sys.time(), "%b %d %X"),'Function: xgb_importance starts. \n')
  fm<-mod$data
  model<-mod$model
  xgb.importance(model =model$finalModel)->imp_mat
  cat(format(Sys.time(), "%b %d %X"),'Function: xgb_importance  finish. \n')
  return(list(data=fm,model=model,importance=imp_mat))
}

xgb_plot_importance<-function(imp){
  cat(format(Sys.time(), "%b %d %X"),'Function: xgb_plot_importance starts. \n')
  fm<-imp$data
  model<-imp$model$finalModel
  imp_mat<-imp$importance
  p<-xgb.ggplot.importance(imp_mat, rel_to_first = TRUE, xlab = "Relative importance")->p
  cat(format(Sys.time(), "%b %d %X"),'Function: xgb_plot_importance  finish. \n')
  return(p)
}
get_shap_values<-function(imp){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_values starts. \n')
  fm<-imp$data
  #cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_values dim(fm)',dim(fm),' \n')
  model<-imp$model$finalModel
  #cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_values length(model)',length(model),' \n')
  imp<-imp$importance
  #cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_values dim(imp)',dim(imp),' \n')
  idx<-grep("(MZ_.*)",names(fm))
  tm<-as.matrix(fm[,idx])
  #cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_values dim(tm)',dim(tm),' \n')
  shap_values <- SHAPforxgboost::shap.values(xgb_model = model, X_train = tm)
  shap_values$tm<-tm
  return(shap_values)
}
get_shap_dep_plot<-function(shap_values){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_dep_plot starts. \n')
  shap_values$tm->tm
  shap_long <- SHAPforxgboost::shap.prep(shap_contrib = shap_values$shap_score, X_train = tm)
  fig_list <- lapply(names(shap_values$mean_shap_score)[1:16],
                     shap.plot.dependence, data_long = shap_long)
  p<-gridExtra::grid.arrange(grobs = fig_list, ncol = 4)
  return(p)
}
get_shap_plot_data<-function(shap_values){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_plot_data starts. \n')
  plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score, top_n = 10, n_groups = 10)
  return(plot_data)
}
get_shap_force_plot<-function(plot_data){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_force_plot starts. \n')
  p<-shap.plot.force_plot(plot_data,zoom_in = FALSE)
  return(p)
}

get_shap_force_group_plot<-function(plot_data){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_force_group_plot starts. \n')
  p<-shap.plot.force_plot_bygroup(plot_data)
  return(p)
}

get_xgb.shap_plot<-function(imp){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_xgb.shap_plot starts. \n')
  fm<-imp$data
  model<-imp$model$finalModel
  imp<-imp$importance
  idx<-grep("(MZ_.*)",names(fm))
  tm<-as.matrix(fm[,idx])
  contr<-predict(model,newdata=tm, predcontrib = TRUE)
  shap<-xgb.plot.shap(tm,contr,model=model$finalModel,
                      features = imp$Feature[1:10],
                      plot = FALSE)
  cat(format(Sys.time(), "%b %d %X"),'Function: get_xgb.shap_plot  finish. \n')
  return(shap)
}

get_shap_plot<-function(imp){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_plot starts. \n')
  fm<-imp$data
  model<-imp$model$finalModel
  imp<-imp$importance
  idx<-grep("(MZ_.*)",names(fm))
  tm<-as.matrix(fm[,idx])
  sh_res<-shap.score.rank(xgb_model = model, 
                          X_train =tm,
                          shap_approx = F
  )
  sh_long<-shap.prep(shap = sh_res,
                     X_train = tm , 
                     top_n = 10
  )
  p<-plot.shap.summary(data_long = sh_long)
  cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_plot  finish. \n')
  return(p)
}

train_trigger<-function(fm){
  cat(format(Sys.time(), "%b %d %X"),'Function: train_trigger',' starts. Mem:',getFreeMem(),'GB\n')
  return(length(unique(fm$norm.p))>2)
}

#### Train reduced model ####
get_reduced_fm<-function(imp,threshold=5){
  cat(format(Sys.time(), "%b %d %X"),'Function: get_shap_values starts. \n')
  fm<-imp$data
  shap_values <- get_shap_values(imp)
  idxMZ<-grep('MZ_',names(fm))
  shval<-cbind(fm[,-idxMZ],shap_values$shap_score)
  idxMZ<-grep('MZ_',names(shval))
  sh_mean<-ddply(shval,.(target),function(.x){apply(.x[,idxMZ],2,mean)})
  sh_mean_long<-melt(sh_mean,id='target')
  idxMZ<-grep('MZ_',names(fm))
  idxOpt<-match(as.character(unique(
    sh_mean_long$variable[abs(sh_mean_long$value)>threshold])),
    names(fm))
  fm_opt<-cbind(fm[,-idxMZ],fm[,idxOpt])
  return(fm_opt)
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
