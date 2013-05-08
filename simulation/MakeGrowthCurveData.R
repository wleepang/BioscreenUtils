# make some fake growth curve data
n.strains = 20 # number of strains being tested - e.g. uniqe rates
n.reps = 10    # number of curve replicates (with added noise)
n.curves = n.wells = n.strains * n.reps

# A Bioscreen plate can hold 100 curves
# A Bioscreen machine can run 2 plates, total of 200 curves
# Of course, one shouldn't be limited to just Bioscreen ...

n.wells = n.curves = 96
n.reps = 4
n.strains = n.wells / n.reps

set.seed(0)

# theses are the values subsequent algorithms are trying to extract
rates.truth = runif(n.strains, min=0.01, max=0.5) # doublings per hour
sats.truth  = runif(n.strains, min=0.5, max=2)    # max od
lags.truth  = runif(n.strains, min=0, max=10)     # lag times

# expand primary growth parameters with noise for replicates
addRepNoise = function(x){rnorm(n.reps, x, runif(1,0.001*x,0.02*x))}
rates = do.call('c', lapply(rates.truth, addRepNoise))
sats  = do.call('c', lapply(sats.truth, addRepNoise))
lags  = do.call('c', lapply(lags.truth, addRepNoise))

# this value should be controlled by the experimenter
inits = rnorm(n.wells, 0.01, 0.001)         # initial ods

# use the generalized logistic DE to generate a growth curve
simulateGrowthCurve = function(tpts, mu, A, tl, X0)
{
  X = rep(NA, length(tpts))
  X[1] = X0
  for (i in seq_along(tpts)[-1]) {
    dt = tpts[i]-tpts[i-1]
    dx = mu*(1 - X[i-1]/A)*X[i-1]*dt
    
    if (tpts[i] < tl) dx = 0
    # the above isn't the best way to produce a delayed logistic
    # but it is effective
    
    X[i] = X[i-1] + dx
  }
  return(X)
}

# timepoints in hours
time = seq(0,48, by=0.25)

curves = sapply(seq_along(rates), function(i){
  x = simulateGrowthCurve(time, rates[i], sats[i], lags[i], inits[i])
  
  # for added realism, add noise to the signal
  x = x + runif(length(x), 0, 0.01)
})
# this returns a matrix with rows for time and columns for wells

# make time the obnoxious h:mm:ss.ss format that Bioscreen uses
hours = floor(time)
mins = floor((time - hours)*60)
secs = time - hours - mins/60
stime = sprintf('%d:%02d:%04.1f', hours, mins, secs)

out = data.frame(Time = stime, curves)
truths = data.frame(rate=rates, sat=sats, lag=lags)
write.csv(out, 'SimulatedGrowthCurves.csv', quote=F, row.names=F)
write.csv(truths, 'SimulatedGrowthCurves_Truth.csv', quote=F)