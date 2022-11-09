#############################
##### Read in the Files #####
#############################

###Specify the file path to folder containing the results from Calculate M_i(r y) for Observed Data.R and Simulate the Null Distribution of the Test Statistic.R
path = 'C:/Users/scwatson/Desktop/'
###Read in the proportion of area observed in the data
opr = read.table(paste0(path,'TotalAreaProp.txt'),fill = TRUE)
###Read in the proportion of each circle which falls in observed units in the data
###Each row corresponds to a radius and each column to a centriod
ots.mat = read.table(paste0(path,'ObservedTestStat.txt'),fill = TRUE)
ots.mat = as.matrix(ots.mat)
###Read in the proportion of each circle which falls in observed units from the Monte Carlo Simulations
###Each row corresponds to a Monte Carlo simulation, columns correspond to radii-centriod pairs (radius 1 and centriod 1, then radius 2 and centriod 2, then radius 1 and centriod 3, etc.)
sts.mat  = read.table(paste0(path,'SimulatedProportion.txt'),fill = TRUE)
sts.mat = as.matrix(sts.mat)
###Read in the proportion of area observed in each Monte Carlo Simulation
spr = read.table(paste0(path,'SimulatedTotalProportion.txt'),fill = TRUE)

################################
#####  Specify Parameters  #####
################################

###Number of Monte Carlo simulations performed
mc.reps = dim(sts.mat)[1]
###Number of radii evaluated
r.rad.vec = 10
###Number of observed units
n.points = 40
###Specify quantiles for the rejection of the null hypothesis
lower.quantile = 0.025
upper.quantile = 0.975

###Compute the test statistic for each Monte Carlo Simulation
test.stat.array = array(NA, c(mc.reps,r.rad.vec))
mean.spr.mat = matrix(NA,r.rad.vec,1)
for(r in 1:r.rad.vec){
  mean.spr.mat[r,1] = mean(sts.mat[,((r-1)*n.points +1):(r*n.points)])
}
for(m in 1:mc.reps){
  for(r in 1:r.rad.vec){
    test.stat.array[m,r] = mean(sts.mat[m,((r-1)*n.points +1):(r*n.points)]/spr[m,1]- mean.spr.mat[r] )
  }
}
###Compute the critical values for the hypothesis test
lower.sts.mat = apply(test.stat.array,2,quantile,prob = lower.quantile)
upper.sts.mat = apply(test.stat.array,2,quantile,prob = upper.quantile)

###Compute the Observed Test Statistic
observed.test.stat.array = matrix(NA,1,r.rad.vec)
for(r in 1:r.rad.vec){
  observed.test.stat.array[1,r] = mean(ots.mat[r,]/as.numeric(opr) - mean.spr.mat[r])
}

###Compare the observed test statistic to the critical values at each radius
obs.lower = matrix(NA,1,r.rad.vec)
obs.upper = matrix(NA,1,r.rad.vec)
for(i in 1:r.rad.vec){
  obs.lower[,i] = observed.test.stat.array[,i] <= lower.sts.mat[i]
  obs.upper[,i] = observed.test.stat.array[,i] >= upper.sts.mat[i]
}
###Was the observed test statistic less than the lower critical value?
obs.lower
###Was the observed test statistic greater than the upper critical value?
obs.upper









