phase.sync <-
  function (t1, t2, mins = FALSE) {
    ## Find the min/max in both timeseries
    min.max1=find.minmax(t1)
    min.max2=find.minmax(t2)
    ## timesteps, densities, phase
    phases1=matrix(nrow=NROW(t1), ncol=3, NA)
    phases2=matrix(nrow=NROW(t2), ncol=3, NA)
    
    phases1[,1:2]=t1
    phases2[,1:2]=t2
    if (mins) {
      v1=min.max1$mins
      v2=min.max2$mins
    }
    else {
      v1=min.max1$maxs
      v2=min.max2$maxs
    }
    ## Locations of mins/maxs
    locs1=v1$location
    locs2=v2$location
    ## Range of values over which to interpolate
    range1=locs1[1]:locs1[length(locs1)]
    range2=locs2[1]:locs2[length(locs2)]
    
    ## Assign phase values to mins/maxs
    phases1[locs1, 2:3]=cbind(v1[, 2], 
                              seq(from=0, by=2*pi, to=(NROW(v1)-1)*2*pi))
    phases2[locs2, 2:3]=cbind(v2[, 2], 
                              seq(from=0, by=2*pi, to=(NROW(v2)-1)*2*pi))
    ## Interpolate phase values between successive mins/maxs
    phases1[range1, 3]=
      approx(x=phases1[, 1], y=phases1[, 3], 
             n=length(range1))$y
    phases2[range2, 3]=
      approx(x=phases2[,1], y=phases2[,3],
             n=length(range2))$y
    
    deltaphase=matrix(nrow=NROW(t1), ncol=4, NA)
    phase_diff=phases1[,3]-phases2[,3]
    mod_phase_diff1=phase_diff %% (2*pi)
    mod_phase_diff2=((phase_diff+pi) %% (2*pi))-pi
    deltaphase=cbind(1:NROW(t1), phase_diff, mod_phase_diff1, mod_phase_diff2)
    colnames(deltaphase)=c("timestep", "phasediff", "mod_phase_diff_2pi", 
                           "mod_phase_diff_pi")
    colnames(phases1)=c("timestep", "val", "phase")
    colnames(phases2)=c("timestep", "val", "phase")
    return (list(phases1=as.data.frame(phases1), 
                 phases2=as.data.frame(phases2), 
                 deltaphase=as.data.frame(deltaphase)))
  }
