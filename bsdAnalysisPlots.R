# define layout
lyt = getLayout('./data/20110908/20110908_layout.txt')

# collect data
bsd = bsdProcess('./data/20110908/20110908_data.csv')
# bsd = list(at=bsd$at, # absolute time
# 					 xn=bsd$xn, # normalized smoothed profile
# 					 mu=bsd$mu) # instantaneous specific growth rate

# generate replicate averages
bsd.m = list(
		at = bsd$at,
		xn = bsdMergeReplicates(bsd$xn, lyt)$xbar,
		mu = bsdMergeReplicates(bsd$mu, lyt)$xbar)

# determine global plot ranges
xn.ylim = range(bsd$xn)
mu.ylim = colMeans(t(apply(bsd$mu, 2, range)))*1.1 # robust way to find global plot range for mu

# plot all unique cultures and their replicates
pdf('./data/20110908/20110908_Summary.pdf', width=11, height=8.5)
for (Name in lyt$Names) {
	at = bsd$at
	reps = paste('Well', lyt$Reps[[Name]], sep='.')
	
	# mar = c(b, l, t, r)
	par(mfrow=c(2,2), mar=c(3,3,1,1), oma=c(0,0,2,0))
	
	
	# plot the mean X profile
	plot(at, bsd.m$xn[, Name], type='l', lwd=2, log='y', ann=F, ylim=xn.ylim)
	title(ylab='Normalized Density', line=2)
	
	# plot the replicate X profiles
	matplot(at, bsd$xn[, reps], type='l', lwd=1, lty=1, col=1:length(reps), log='y', ann=F, ylim=xn.ylim)
	matpoints(at, t(t(bsd$x)/bsd$xs[1,])[,reps], col=1:length(reps), pch='.', ann=F)
	legend('bottomright', reps, col=1:length(reps), lwd=2, lty=1, bty='n', cex=0.8)
	
	# plot the mean mu profile
	plot(at[-1], bsd.m$mu[, Name], type='l', lwd=2, ann=F, ylim=mu.ylim)
	title(xlab='Time', ylab='Instantaneous Rate', line=2)
	
	# plot the replicate mu profiles
	matplot(at[-1], bsd$mu[, reps], type='l', lwd=2, lty=1, col=1:length(reps), ann=F, ylim=mu.ylim)
	title(xlab='Time', line=2)
	
	mtext(Name, outer=T)
}
dev.off()