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

plot(smooth.adaptive(xy.in), type='l', log='y', main=sprintf('iter: %d', i))
points(xy.in, pch=16, cex=0.2, col='red')
