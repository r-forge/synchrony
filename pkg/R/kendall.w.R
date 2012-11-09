kendall.w <-
function (data) {
  n=NROW(data)
  M=NCOL(data)
  ## Rank data
  ranks=apply(as.matrix(data), MARGIN=2, rank)
  ties=unlist(apply(ranks, MARGIN=2, FUN=function (x) {table(x)[table(x) >1]}))
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
               pval=p.val, 
               spearman.corr=spearman.ranked.corr))
}

