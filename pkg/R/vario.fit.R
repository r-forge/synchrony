vario.fit <- function (vario, bins,
                       type=c("spherical", "gaussian", "nugget", "linear",
                              "exponential", "sill", "periodic", "hole"),
                       start.vals=list(c0=min(vario), 
                                       c1=max(vario), 
                                       a=max(bins)/4,
                                       b=0.1,
                                       c=0.1)) {
  
  data=list(vario=vario, bins=bins)
  types=c("spherical", "gaussian", "nugget", "linear", "exponential", "sill", 
          "periodic", "hole")
  type=match.arg(tolower(type), types)
  
  ## Linear model used to determine initial parameters for finicky NLS
  vario.lin=lm(vario ~ bins)
  
  if (type=="nugget") {
    vario.mod=lm(vario ~ 1)
    names=c("nugget")    
    success=TRUE
  }
  else if (type=="linear") {
    vario.mod=vario.lin
    names=c("intercept", "slope")
    success=TRUE
  }
  else if (type=="sill") {
    names=c("c0", "c1", "a", "b")
    vario.mod=try(nls(vario ~ (bins < a)*(c0+b*bins) + (bins >=a)*(c0+c1), 
                      start=list(c0=start.vals$c0,
                                 c1=start.vals$c1,
                                 a=start.vals$a,
                                 b=coef(vario.lin)[2]), 
                      data=data), silent=TRUE)
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        start.vals$c1,
                        start.vals$a,
                        coef(vario.lin)[2]), rmse.sill, control=list(maxit=200000),
                      data=data)
    }
    else {
      success=TRUE
    }
  }
  else if (type=="exponential") {
    names=c("c0", "c1", "a")
    vario.mod=try(nls(vario ~ c0+c1*(1-exp(bins/a)), 
                      start=list(c0=start.vals$c0,
                                 c1=start.vals$c1,
                                 a=start.vals$a), 
                      data=data), silent=TRUE)
  
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        start.vals$c1,
                        start.vals$a), rmse.expo, control=list(maxit=200000),
                      data=data)
    }
    else {
      success=TRUE  
    }
  }
  else if (type=="spherical") {
    names=c("c0", "c1", "a")
    vario.mod=try(nls(vario ~ (bins < a)*(c0+c1*((3*bins)/(2*a)-(1/2)*(bins/a)^3))+(bins >=a)*(c0+c1), 
                      start=list(c0=start.vals$c0,
                                 c1=start.vals$c1,
                                 a=start.vals$a), 
                      data=data), silent=TRUE)
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        start.vals$c1,
                        start.vals$a), rmse.sphere, control=list(maxit=200000),
                      data=data)
    }
    else {
      success=TRUE
    }
  }
  else if (type=="gaussian") {
    names=c("c0", "c1", "a")
    vario.mod=try(nls(vario ~ (bins < a)*(c0+c1*(1-exp(-3*(bins^2)/(a^2))))+(bins >=a)*(c0+c1), 
                      start=list(c0=start.vals$c0,
                                 c1=start.vals$c1,
                                 a=start.vals$a), 
                      data=data), silent=TRUE)
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        start.vals$c1,
                        start.vals$a), rmse.gauss, control=list(maxit=200000),
                      data=data)
    }
    else {
      success=TRUE
    }
  }
  else if (type=="periodic") {
    names=c("a", "b", "c")
    vario.mod=try(nls(vario ~ a*cos(b*pi/max(bins)+c)*bins, 
                      start=list(a=start.vals$a,
                                 b=start.vals$b,
                                 c=start.vals$c), 
                      data=data), silent=TRUE)
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$a, 
                        start.vals$b,
                        start.vals$c), rmse.period, control=list(maxit=200000),
                      data=data)
    }
    else {
      success=TRUE
    }  
  }
  else if (type=="hole") {
    names=c("c0", "c1", "a")
    vario.mod=try(nls(vario ~ c0+c1*(1-(a*sin(bins/a))/bins), 
                      start=list(c0=start.vals$c0,
                                 c1=start.vals$c1,
                                 a=start.vals$a), 
                      data=data), silent=TRUE)
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        start.vals$c1,
                        start.vals$a), rmse.hole, control=list(maxit=200000),
                      data=data)
    }
    else {
      success=TRUE
    }  
  }
  
  opt=vario.stats(data, vario.mod, type, names, success)
  results=list(vario=vario, bins=bins, AIC=opt$AIC, rmse=opt$rmse, 
               params=opt$params, fit=opt$fit, model=type, nls.success=opt$nls.success)
  class(results)="variofit"
  return (results)
}

vario.stats <- function (data, opt, type, names, success) {
  vario=data$vario
  bins=data$bins
  if (success) {
    mod.aic=AIC(opt)
    fit=predict(opt)
    params=coef(opt)
    rmse=sqrt(mean((predict(opt)-vario)^2))
  }
  else {
    N=length(vario)
    mod.logLik=0.5*(-N*(log(2*pi)+1-log(N)+log(opt$value^2*N)))
    mod.aic=-2*mod.logLik+2*(length(opt$par)+1)
    names(opt$par)=names
    params=opt$par
    rmse=opt$value
    
    if (type=="sill") {
      fit=ifelse (bins < opt$par["a"], 
                  opt$par["c0"]+opt$par["b"]*bins,
                  opt$par["c0"]+opt$par["c1"])
    }
    else if (type=="exponential") {
      fit=opt$par["c0"]+opt$par["c1"]*(1-exp(bins/opt$par["a"]))
    }
    else if (type=="spherical") {
      fit=ifelse (bins < opt$par["a"],
                  opt$par["c0"]+
                    opt$par["c1"]*((3*bins)/(2*opt$par["a"])-
                                        0.5*(bins/opt$par["a"])^3),
                  opt$par["c0"]+opt$par["c1"])
    }
    else if (type=="gaussian") {
      fit=ifelse (bins < opt$par["a"],
                  opt$par["c0"]+
                    opt$par["c1"]*(1-exp(-3*(bins)^2/(opt$par["a"]^2))), 
                  opt$par["c0"]+opt$par["c1"])
    }
    else if (type=="periodic") {
      fit=opt$par["a"]*cos(opt$par["b"]*pi*(bins/max(bins))+opt$par["c"])
    }
    else if (type=="hole") {
      fit=opt$par["c0"]+opt$par["c1"]*(1-(opt$par["a"]*sin(bins/opt$par["a"]))/bins)      
    }    
  }
  names(params)=names
  return (list(AIC=mod.aic, rmse=rmse, params=params, fit=fit, nls.success=success))
} 

rmse.period <- function (x, data) {
  a=x[1]; b=x[2]; c=x[3]
  vario=data$vario; bins=data$bins
  
  variohat=a*cos(b*pi*(bins/max(bins))+c)
  rmse=sqrt(sum((vario-variohat)^2)/NROW(bins))  
}

rmse.hole <- function (x, data) {
  c0=x[1]; c1=x[2]; a=x[3]
  vario=data$vario; bins=data$bins
  
  if (a >= 0 & a < max(bins) & c0 >= 0 & c1 >= 0) {
    variohat=c0+c1*(1-(a*sin(bins/a))/bins)
    rmse=sqrt(sum((vario-variohat)^2)/NROW(bins))
  }
  else
    rmse=Inf
}

rmse.sill <- function (x, data) {
  c0=x[1]; c1=x[2]; a=x[3]; b=x[4]
  vario=data$vario; bins=data$bins
  
  if (a < max(bins) & a > 0 & b >= 0 & c1 >= 0 & c0 >= 0) {
    variohat=ifelse (bins < a, 
                     c0+b*bins, 
                     c0+c1)
    rmse=sqrt(sum((vario-variohat)^2)/NROW(bins))
  }
  else
    rmse=Inf
}

rmse.expo <- function (x, data) {
  c0=x[1]; c1=x[2]; a=x[3]
  vario=data$vario; bins=data$bins
  
  if (a < max(bins) & a > 0 & c1 >= 0 & c0 >= 0) {
    variohat=c0+c1*(1-exp(bins/a))
    rmse=sqrt(sum((vario-variohat)^2)/NROW(bins))
  }
  else
    rmse=Inf
}

rmse.sphere <- function (x, data) {
  c0=x[1]; c1=x[2]; a=x[3]
  vario=data$vario; bins=data$bins
  
  if (a < max(bins) & a > 0 & c1 >= 0 & c0 >= 0) {
    variohat=ifelse (bins < a, 
                     c0+c1*(3*bins/(2*a)-0.5*(bins/a)^3), 
                     c0+c1)
    rmse=sqrt(sum((vario-variohat)^2)/NROW(bins))
  }
  else
    rmse=Inf
}

rmse.gauss <- function (x, data) {
  c0=x[1]; c1=x[2]; a=x[3]
  vario=data$vario; bins=data$bins
  
  if (a < max(bins) & a > 0 & c1 >= 0 & c0 >= 0) {
    variohat=ifelse (bins < a, 
                     c0+c1*(1-exp(-3*bins^2/(a^2))), 
                     c0+c1)
    rmse=sqrt(sum((vario-variohat)^2)/NROW(bins))
  }
  else
    rmse=Inf
}
