getCorrStruct<-function(X,cutoff=0.85){
  cm<-cor(X)
  h<-hclust(as.dist(1-abs(cm)),method = 'ward.D2')
  c<-cutree(h,h=1)
  ncl<-max(c);ncl
  mmm<-min(sapply(1:ncl,function(k){j<-which(c==k);return(min(min(abs(cm[j,j]))))}))
  while(mmm<cutoff){
    ccl<-ncl+1
    c<-cutree(h,k=ccl)
    ncl<-max(c);ncl
    mmm<-min(sapply(1:ncl,function(k){j<-which(c==k);return(min(min(abs(cm[j,j]))))}))
  }
  return(c)
}

groupCATscore<-function(X,class,cutoff=0.85){
  tstat = catscore(as.matrix(X), class)
  c<-getCorrStruct(X,cutoff)
  laply(c,function(jj){j<-which(c==jj);tmp<-sign(tstat[j,])*sqrt(sum(tstat[j,]^2));tmp[j==j[1]]})->tmp
colnames(tmp)<-paste0('g',colnames(tmp))
res<-cbind(as.data.frame(tstat[TRUE,]),tmp)
return(res[order(abs(tmp[,1]),decreasing = TRUE),])

}