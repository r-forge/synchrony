phase.partnered <-
function (n=2000, rho=1, gamma=1, ts.sd=0.10, ts.mean=0) {
  
  ts.part=matrix(nrow=n, ncol=2, NA)
  minVal=0
  maxVal=2*pi
  fs=seq(from=1, to=n/2, by=1)
  epsilon=runif(n=n/2)
  epsilon [epsilon > 0.5]=1
  epsilon [epsilon < 0.5]=-1

  delta=acos(rho)*epsilon
  phi1=minVal+(maxVal-minVal)*runif(n=n/2)
  phi2=phi1+delta
  const.num=1/(fs^(gamma/2))

  for (t in 1:n) {
    ts.part[t,1]=sum(const.num*(sin(2*pi*fs*t/n + phi1)))
    ts.part[t,2]=sum(const.num*(sin(2*pi*fs*t/n + phi2)))
  }

  ts.part=ts.mean+(ts.part/apply(ts.part, 2, sd))*ts.sd
  ts.part=as.data.frame(ts.part)
  colnames(ts.part) <- c("ts1", "ts2")

  return (ts.part)
}
