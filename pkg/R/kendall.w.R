kendall.w <- function (data, nrands = 0) {
  data=as.matrix(data)
  kendall.observed=kendall.w.aux (data)
  n=NROW(data)
  M=NCOL(data)
  
  if (nrands == 0) {
    results=kendall.observed
  }
  else {
    rands=numeric(length=nrands+1)
    prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    for (i in 1:nrands) {
      tmp=matrix(sample(data), nrow=n, ncol=M)
      ## Using chi^2 or W is equivalent for permutation purposes 
      rands[i]=kendall.w.aux(tmp)$w.corrected
      setTxtProgressBar(prog.bar, i)  
    }
    rands[i+1]=kendall.observed$w.corrected
    p.val.rand=sum(rands >= kendall.observed$w.corrected)/(nrands+1)
    results=list(w.uncorrected=kendall.observed$w.uncorrected, 
                 w.corrected=kendall.observed$w.corrected, 
                 pval=kendall.observed$p.val, pval.rand=p.val.rand, rands=rands,
                 spearman.corr=kendall.observed$spearman.ranked.corr)
  }
  return (results)
}

kendall.w.aux <- function (data) {
  n=NROW(data)
  M=NCOL(data)
  ## Rank data
  ranks=apply(data, MARGIN=2, rank)
  ties=unlist(apply(ranks, MARGIN=2, FUN=function (x) 
  {tab=table(x); tab[tab >1]}))
  ties=sum(ties^3-ties)
  
  sum.ranks.row=rowSums(ranks)
  squared.sum.ranks=(sum(sum.ranks.row))^2
  sum.squared.ranks=sum(sum.ranks.row^2)
  ## Uncorrected concordance
  w.uncorrected=(sum.squared.ranks - (squared.sum.ranks)/n)/(M^2*(n^3-n)/12)
  ## Corrected concordance
  w.corrected=(sum.squared.ranks - (squared.sum.ranks)/n)/((M^2*(n^3-n)-M*ties)/12)
  spearman.ranked.corr=(M*w.corrected-1)/(M-1)
  chisq.val=M*(n-1)*w.corrected
  ## When n is large enough (>10), chi^2_r approximated as chi^2 with df=n-1 
  p.val=1-pchisq(chisq.val, n-1)
  
  return (list(w.uncorrected=w.uncorrected, 
               w.corrected=w.corrected, 
               pval=p.val, spearman.corr=spearman.ranked.corr))
  } 
