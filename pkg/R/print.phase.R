print.phase <- function (x, digits=max(3L, getOption("digits") - 3L), ...) {
  cat("Phases and values of time series 1:\n")
  print(format(x$phases1, digits=digits))
  
  cat("\nPhases and values of time series 2:\n")
  print(format(x$phases2, digits=digits))
  
  cat("\nPhase difference between time series:\n")
  print(format(x$deltaphase, digits=digits))  
}
