# test script to demonstrate adaptive smoothing

# read in data
curves = read.csv('SimulatedGrowthCurves.csv')

# the first column is time
# ideally it should be regularly spaced

time = seq_along(curves$Time)
X = paste('X', 2, sep='')

xy.in = list(x=time, y=curves[[X]])

smooth.adaptive = function(x, y=NULL, pct=0.1, niter=NULL) {
  # performs jack-knifed / cross-validated smoothing
  
  if (!is.null(x$y)) {
    xy.in = x
  } else {
    xy.in = list(x=x, y=y)
  }
  
  if (is.null(niter)) {
    niter = ceiling(length(xy.in$x)*pct)*10
  }
  
  XY = list(pct=pct, niter=niter)
  for (i in 1:niter) {
    set.seed(i) # this ensures that the same fit happens with repeated runs
    ix = sort(sample(xy.in$x, floor((1-pct)*length(xy.in$x))))
    
    xy = list(x=xy.in$x[ix], y=xy.in$y[ix])
    
    xy.s = smooth.spline(xy)           # smoothing of sampled data
    xy.ss = spline(xy.s, xout=xy.in$x) # spline interpolate to full range
    
    if (i == 1) {
      # first iteration, set trajectory
      XY[['x']]=xy.ss$x
      XY[['y']]=xy.ss$y
    } else {
      # gather the mean trajectory of all iterations
      XY$y = (XY$y + xy.ss$y) / 2
    }
  }
    
  return(XY)
  
}

#xy.out = smooth.adaptive(xy.in)
#plot(xy.out, type='l', log='y', main=sprintf('iter: %d', xy.out$niter))
#points(xy.in, pch=16, cex=0.2, col='red')

curves.s = curves[,-1]
curves.s = lapply(curves.s, function(y){
  y.s = smooth.adaptive(list(x=time, y=y), niter=10)$y
  #y.s.n = y.s / min(y.s[1:floor(length(y.s)*0.05)])
  return(y.s)
})
curves.s = data.frame(curves.s)

#matplot(curves.s, type='l', log='y')

windows(12,8)
par(mar=c(0,0,0,0))
layout(matrix(1:96, nrow=8, ncol=12, byrow=T))
lapply(1:96, function(i){
  plot(time, curves[[i+1]], type='p', pch=16, cex=0.2, col='red', log='y', xaxt='n', yaxt='n', ylim=c(0.005, max(curves[,-1])))
  lines(time, curves.s[[i]], lwd=1.5)
  legend('bottomright', legend=colnames(curves)[i+1], bty='n')
})