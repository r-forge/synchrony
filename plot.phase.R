plot.phase <- function (x, plot.phase=TRUE, at=seq(-pi, by=pi/2, to=pi), len=14,
                        padj=0.7, mod=2, main="", xlab="Phase difference", 
                        ylab="Frequency", col="red", lty=2, lwd=1, ...) {
  if (plot.phase & class(x)=="phase") {
    if (mod==1) {
      p=x$deltaphase$mod_phase_diff_pi
      b=seq(-pi, pi, length=len)
    }
    else {
      p=x$deltaphase$mod_phase_diff_2pi
      at=seq(0, by=pi/2, to=2*pi)
      f=fractions(at/pi)
      labels=paste(expression(f, pi, sep="*"))
      b=seq(0, 2*pi, length=len)
    }
    
    z=hist(p, plot=FALSE, breaks=b)
    hist(p, breaks=z$breaks, sub="", xaxt="n", xlab="", main=main, ylab=ylab)
    box()
    axis(side=1, at=at, labels=labels, padj=padj)
    mtext(side=1, text=xlab, line=3, cex=0.8)
  }
  
  else {
    stop("No permutation data to plot")
  }
}
