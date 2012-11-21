data(pisco.data)
## Environmental variogram
d=subset(pisco.data, subset=year==2000, select=c("latitude", "longitude", "upwelling"))
semiv=vario(data=d)
png(file="example_semivar.png", width=800, height=600, pointsize=14)

plot(semiv, xlab="Lag distance (km)")
mod.sph=vario.fit(semiv$vario, semiv$mean.bin.dist)
mod.exp=vario.fit(semiv$vario, semiv$mean.bin.dist, type="expo")
mod.gau=vario.fit(semiv$vario, semiv$mean.bin.dist, type="gauss")
mod.lin=vario.fit(semiv$vario, semiv$mean.bin.dist, type="lin")
lines(semiv$mean.bin.dist, mod.sph$fit, col="red")
lines(semiv$mean.bin.dist, mod.exp$fit, col="black")
lines(semiv$mean.bin.dist, mod.gau$fit, col="blue")
lines(semiv$mean.bin.dist, mod.lin$fit, col="green")
legend(x="topleft", legend=paste(c("Spherical AIC:", "Exponential AIC:", 
                                   "Gaussian AIC:", "Linear AIC:"), 
                                 c(format(mod.sph$AIC, dig=2), 
                                   format(mod.exp$AIC, dig=2), 
                                   format(mod.gau$AIC, dig=2), 
                                   format(mod.lin$AIC, dig=2))), lty=1, col=c("red", "black", "blue", "green"), 
       bty="n")
dev.off()

## Example 2
png(file="example_moran.png", width=800, height=600, pointsize=14)
cover=subset(pisco.data, subset=year==2000, select=c("latitude", "longitude", "mussel_abund"))
moran=vario(data=cover, type="moran")
mod.hol=vario.fit(moran$vario, moran$mean.bin.dist, type="hole", start.vals=list(c0=0.6, a=25, c1=0.01))
mod.per=vario.fit(moran$vario, moran$mean.bin.dist, type="period", start.vals=list(a=1, b=3, c=0))
mod.lin=vario.fit(moran$vario, moran$mean.bin.dist, type="linear")
plot(moran, xlab="Lag distance (km)", ylim=c(-0.6, 0.8))
lines(moran$mean.bin.dist, mod.per$fit, col="red")
lines(moran$mean.bin.dist, mod.hol$fit, col="black")
lines(moran$mean.bin.dist, mod.lin$fit, col="blue")
legend(x="topleft", legend=paste(c("Periodic AIC:", "Hole AIC:", 
                                   "Linear AIC:"), 
                                 c(format(mod.per$AIC, dig=2), 
                                   format(mod.hol$AIC, dig=2), 
                                   format(mod.lin$AIC, dig=2))), 
       lty=1, col=c("red", "black", "blue"), bty="n")
dev.off()

## Example 3
data(pisco.data)
## Compute spatial synchrony
d.upw=subset(pisco.data, select=c("latitude", "longitude", "year", "upwelling"))
d.cov=subset(pisco.data, select=c("latitude", "longitude", "year", "mussel_abund"))
## Reshape the data
d.upw.wide=reshape(data=d.upw, idvar=c("latitude", "longitude"), timevar=c("year"), 
                   direction="wide")
d.cov.wide=reshape(data=d.cov, idvar=c("latitude", "longitude"), timevar=c("year"), 
                   direction="wide")
## Generate variograms
v.upw=vario(nbins=12, data=d.upw.wide, type="pearson", extent=1, nrands=999, alt="two", center=TRUE)
v.cov=vario(nbins=12, data=d.cov.wide, type="pearson", extent=1, nrands=999, alt="two", center=TRUE)
## Fit variograms
v.cov.per=vario.fit(v.cov$vario, v.cov$mean.bin.dist, type="period", 
                    start.vals=list(a=1, b=3, c=0))
v.upw.lin=vario.fit(v.upw$vario, v.upw$mean.bin.dist, type="linear")

png(file="example_synchrony.png", width=600, height=800, pointsize=14)
par(oma=c(0, 0, 0, 1), mar=c(5, 4, 2, 4) + 0.1, mfrow=c(2,1))
plot(v.cov, xlab="Lag distance (km)", bg.sig="red", col.nonsig="red", col.sig="red",
     main="Mussel synchrony", ci=TRUE,
     rug=FALSE, ylim=c(-0.3, 0.3))
lines(v.cov$mean.bin.dist, v.cov.per$fit, col="red")
plot(v.upw, xlab="Lag distance (km)", bg.sig="blue", col.nonsig="blue", col.sig="blue",
     main="Upwelling synchrony", ci=TRUE)
lines(v.upw$mean.bin.dist, v.upw.lin$fit, col="blue")
dev.off()
