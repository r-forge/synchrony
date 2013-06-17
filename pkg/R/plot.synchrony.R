plot.synchrony <- function (x, main="", xlab="Observations from randomizations", ylab="Frequency", col="red", lty=2, lwd=1, ...) {
  if (!is.null(x$rands) & class(x)=="synchrony") {
    if (!is.null(x$w.corrected))
      x$obs=x$w.corrected
    hist(x$rands, main=main, xlab=xlab, ylab=ylab, ...)
    abline(v=x$obs, col=col, lty=lty, lwd=lwd, ...)
    box()
  }
  else {
    stop("No permutation data to plot")
  }
}