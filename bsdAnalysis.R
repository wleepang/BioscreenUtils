# bsdAnalysis.R
# A simple script for analyzing bioscreen growth data.
# By Lee Pang - Institute for Systems Biology

# History:
# 20090818	Lee Pang	Created
# 20090819	Lee Pang	Massively simplified via modular functions
#						Incorporates replicate handling

# First a few function definitions
source('bsdAnalysisFunc.R')

# Process the data.
bsd = bsdProcess('BioscreenData.csv')
# returns:
#	bsd$at	absolute time in hours
#	bsd$x	raw OD600 values
#	bsd$xs	smoothed raw OD600 values
#	bsd$xn	smoothed normalized OD600 values
#	bsd$mu	instantaneous rate profiles

lyt = getLayout('BioscreenLayout_matrix.txt')
# alternatively:
# lyt = getLayout('BioscreenLayout_list.txt', type='list')
# returns:
#	lyt$Names	the unique culture names defined in the layout
#	lyt$Reps	replicate wells for each culture name

xn.merged = bsdMergeReplicates(bsd$xn, lyt)
mu.merged = bsdMergeReplicates(bsd$mu, lyt)
# each 'merged' variable has the following structure
#	$xbar	matrix of average profiles for all replicates, colnames
#			are unique culture names
#
#	$rmsd	list of root mean squared deviation values for replicates
#			in each culture, a measure of replicate consistency

# calculate the area under the normalized log transform curve to get the
# growth potential factor
xn.auc = apply(bsd$xn, 2, function(y){simp(y, x=bsd$at)})
# alternatively:
# xn.auc.merged = apply(xn.merged$xbar, 2, function(y){simp(y, x=bsd$at)})
# -or-
# xn.auc.merged = bsdMergeReplicates(t(as.matrix(xn.auc)), lyt)
# 
# the same alternates apply for the parameters below


# other useful parameters
# mumax - maximum instantaneous growth rate
mumax = getMuMax(bsd$mu)

# mumaxtime - time when maximum instantaneous growth was reached
mumaxtime = getMuMaxTime(bsd$at, bsd$mu)

# for more complicated growth profiles - e.g. diauxic shifts - it may be
# necessary to get a full list of peaks in the instantaneous growth data
# using:
# mupks = getMuPeaks(bsd$mu, at=bsd$at)

# stitch the results together listed by WellID
result = cbind(xn.auc, mumax, mumaxtime)

# write the results to file
write.table(result, file="result.tsv", sep="\t")