vario <- function (nbins=20, extent=0.5, data, data2=NULL, is.latlon=TRUE, centered=FALSE,
                   nrands=0, type=c("semivar", "cov", 
                                    "pearson", "spearman", "kendall", "moran", "geary"),
                   alternative=c("one.tailed", "two.tailed"),
                   mult.test.corr=c("none", "holm", "hochberg", "sidak", 
                                    "bonferroni")) {
  
  tails=c("one.tailed", "two.tailed")
  alternative=match.arg(tolower(alternative), tails)
  
  types=c("semivar", "cov", "pearson", "spearman", "kendall", "moran", "geary")
  type=match.arg(tolower(type), types)
  
  mults=c("none", "holm", "hochberg", "sidak", "bonferroni")
  mult.test.corr=match.arg(tolower(mult.test.corr), mults)
  
  n.cols=NCOL(data)
  if (n.cols > 3)
    is.multivar=TRUE
  else
    is.multivar=FALSE
  
  results=vario.aux (nbins=nbins, extent=extent, data=data, data2=data2, 
                     is.latlon=is.latlon, centered=centered, 
                     is.multivar=is.multivar, type=type)
  
  if (nrands > 0) {
    rands=matrix(NA, nrow=nrands+1, ncol=length(results$bins))
    prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    
    if (is.null(data2))
      data2=data
    for (i in 1:nrands) {
      s=sample(results$grpdata)
      
      if (is.multivar) {
        vals=results$vals
        rands[i,]=tapply(vals, s, FUN=mean, na.rm=TRUE)
      }
      else {
        for (j in 1:(length(results$bins))) {        
          tmp=results$all.combs[s==j,]
          tmp=tmp[complete.cases(tmp),]
          rands[i,j]=vario.func(data[tmp[,1], 3], data2[tmp[,2], 3], results$glob.mean, 
                                results$glob.sd, results$glob.N, is.multivar, type)
        }
      }
      setTxtProgressBar(prog.bar, i)
    }
    rands=rands-rowMeans(rands, na.rm=TRUE)
    crit.val=0
    if (!centered) {
      rands=rands+results$regional.mean
      crit.val=results$regional.mean
    }
    rands[nrands+1,]=results$vario
    
    if (alternative == "two.tailed") {
      pvals=apply(rands, MARGIN=2, 
                  function (x) {
                    sum(abs(x - crit.val) >= abs(x[nrands+1] - crit.val))/(nrands+1)
                  })
    }
    else {
      pvals=apply(rands, MARGIN=2, function (x) {
        ifelse (x[nrands+1] > crit.val, 
                sum(x >= x[nrands+1])/(nrands+1), 
                sum(x <= x[nrands+1])/(nrands+1))})
    }

    if (mult.test.corr != "none") {
      pvals=p.adjust(pvals, method=mult.test.corr[1])
    }
    
    colnames(rands)=names(results$bins)
    names(pvals)=names(results$bins)
    
    results$pvals=pvals
    results$rands=rands
    results$alternative=alternative
    results$mult.test.corr=mult.test.corr[1]
  }
  
  ## Remove extraneous elements
  results[8:13]=NULL
  results$is.multivar=is.multivar
  
  class(results)="vario"
  return(results)
}

vario.aux <- function (nbins=20, extent=0.5, data, data2=NULL, is.latlon=TRUE, 
                       centered=FALSE, is.multivar=FALSE,
                       type=c("semivar", "cov", "pearson", "spearman", "kendall", "moran", "geary")) {
  
  n.cols=NCOL(data)
  
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
  
  bins=seq(0, max.extent, length.out=nbins+1)
  grpdata<-cut(all.dists, breaks=bins+1, labels=1:(length(bins)-1))  
  if (is.multivar) {
    glob.mean=NA
    glob.sd=NA
    glob.N=NA
    all.combs=NA
    
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
    
    regional.mean=mean(vals, na.rm=TRUE)
    if (centered) {
      vals=vals-regional.mean
    }
    vario=tapply(vals, grpdata, mean, na.rm=T)
    npoints=tapply(vals, grpdata, FUN=function (x) {length(na.omit(x))})
    bin.dist=tapply(all.dists, grpdata, FUN=mean, na.rm=TRUE)
  }
  else {
    vals=NA
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
    regional.mean=mean(vario, na.rm=TRUE)
    if (centered)
      vario=vario-regional.mean
  }
  
  bins=bins[1:(length(bins))-1]
  col.names=1:length(bins)
  names(bins)=col.names
  names(bin.dist)=col.names
  names(vario)=col.names
  
  return (list(bins=bins, mean.bin.dist=bin.dist,
               vario=vario, npoints=npoints, metric=type, centered=centered,
               regional.mean=regional.mean, all.combs=all.combs, grpdata=grpdata, 
               glob.mean=glob.mean, glob.sd=glob.sd, glob.N=glob.N, vals=vals))
}
