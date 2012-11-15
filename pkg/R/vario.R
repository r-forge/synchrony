vario <- function (nbins=20, extent=0.5, data, data2=NULL, is.latlon=TRUE, centered=FALSE,
                       nrands=0, 
                       type=c("semivar", "cov", 
                              "pearson", "spearman", "kendall", "moran", "geary"),
                       mult.test.corr=c(FALSE, "holm", "hochberg", "sidak", 
                                        "bonferroni")) {
  
  if (!is.null(data2)) {
    all.dists=coord2dist(data[, 1:2], is.latlon, lower.tri=FALSE)
    all.dists=all.dists[row(all.dists)!=col(all.dists)]
  }
  else
    all.dists=coord2dist(data[, 1:2], is.latlon)
  
  ## Compute maximum distance
  max.dist=max(all.dists)
  max.extent=max.dist*extent
  
  types=c("semivar", "cov", "pearson", "spearman", "kendall", "moran", "geary")
  type=match.arg(tolower(type), types)
  n.cols=NCOL(data)
  
  if (n.cols > 3)
    is.multivar=TRUE
  else
    is.multivar=FALSE
  
  bins=seq(0, max.extent, length.out=nbins+1)
  grpdata<-cut(all.dists, breaks=bins+1, labels=1:(length(bins)-1))  
  if (is.multivar) {
    glob.mean=NULL
    glob.sd=NULL
    glob.N=NULL
    
    if (is.null(data2)) {
      if (type=="cov")
        vals=cov(t(data[, 3:n.cols]))
      else {
        vals<-cor(t(data[, 3:n.cols]), method=type)
      }
      vals=vals[lower.tri(vals)]
    }
    else {
      if (type=="cov")
        vals=cov(x=t(data[, 3:n.cols]), y=t(data2[, 3:n.cols]))
      else
        vals=cor(x=t(data[, 3:n.cols]), y=t(data2[, 3:n.cols]), method=type)
      vals=vals[row(vals)!=col(vals)]
    }
    
    if (centered)
      vals=vals-mean(vals, na.rm=T)
    
    vario=tapply(vals, grpdata, mean, na.rm=T)
    npoints=tapply(vals, grpdata, FUN=function (x) {length(na.omit(x))})
    bin.dist=tapply(all.dists, grpdata, FUN=mean, na.rm=TRUE)
  }
  else {
    bin.dist=numeric(length(bins)-1)*NA
    vario=numeric(length(bins)-1)*NA
    npoints=numeric(length(bins)-1)*NA
    
    if (!is.null(data2)) {
      all.combs=expand.grid(1:NROW(data), 1:NROW(data2))
      ## Exclude diagonal elements
      all.combs=all.combs[all.combs[,1]!=all.combs[,2], ]
      glob.mean=c(mean(data[,3], na.rm=TRUE), mean(data2[,3], na.rm=TRUE))
      glob.sd=c(sd(data[,3], na.rm=TRUE), sd(data2[,3], na.rm=TRUE))
      glob.N=NROW(data[,3])
    }
    else {
      data2=data
      all.combs=t(combn(NROW(data), 2))
      glob.mean=mean(data[,3], na.rm=TRUE)
      glob.sd=sd(data[,3], na.rm=TRUE)
      glob.N=NROW(data[,3])          
    }
    for (i in 1:(length(bins)-1)) {
      tmp=all.combs[grpdata==i,]
      tmp=tmp[complete.cases(tmp),]
      x=data[tmp[,1], 3:n.cols]
      y=data2[tmp[,2], 3:n.cols]
       
      vario[i]=vario.func(x, y, glob.mean, glob.sd, glob.N, is.multivar, type=type)
      bin.dist[i]=mean(all.dists[grpdata==i], na.rm=T)
      npoints[i]=NROW(x)
    }
    
    if (centered)
      vario=vario-mean(vario)
  }
    
  bins=bins[1:(length(bins))-1]
  col.names=1:length(bins)
  names(bins)=col.names
  names(bin.dist)=col.names
  names(vario)=col.names
  
  if (nrands > 0) {
    if (is.multivar)
      data=vals
    signif <- vario.signif (nrands=nrands, bins=bins, all.combs, grpdata, 
                            data, data2, type, vario, glob.mean, glob.sd, glob.N, 
                            is.multivar, mult.test.corr)
    
    results=list(bins=bins, mean.bin.dist=bin.dist, 
                 vario=vario, npoints=npoints, pvals=signif$pvals, rands=signif$rands, 
                 metric=type, mult.test.corr=signif$mult.test.corr, is.multivar=is.multivar)
  }
  else {
    results=list(bins=bins, mean.bin.dist=bin.dist,
                 vario=vario, npoints=npoints, pvals=NA, rands=NA, metric=type, 
                 mult.test.corr=NA, is.multivar=is.multivar)
  } 
  class(results)="vario"
  return (results)
}
