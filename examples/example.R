library("MPG")

# set a random seed to generate consistent data sets for testing
set.seed(123)

# number of events to create for each cluster
n_events_clusterA = 500
n_events_clusterB = 500

# number of columns in our data sets
n_parameters = 2

# generate 1st data set
clusterA_events = matrix(rnorm( n_events_clusterA * n_parameters), ncol = n_parameters)
clusterB_events = matrix(rnorm( n_events_clusterA * n_parameters) + 10, ncol = n_parameters)
data_set_1 = rbind(clusterA_events, clusterB_events)

# generate 2nd data set
clusterA_events = matrix(rnorm( n_events_clusterA * n_parameters), ncol = n_parameters)
clusterB_events = matrix(rnorm( n_events_clusterA * n_parameters) + 8, ncol = n_parameters)
data_set_2 = rbind(clusterA_events, clusterB_events)

# shuffle the data sets
data_set_1 = data_set_1[sample(nrow(data_set_1)),]
data_set_2 = data_set_2[sample(nrow(data_set_2)),]

# combine data sets into 1 matrix
Y = rbind(data_set_1, data_set_2)

# map for identifying the data set for each event
C = c(rep(1, nrow(data_set_1)), rep(2, nrow(data_set_2)))

# define our clustering parameters
mcmc = list(nburn = 2000, nsave = 20, nskip = 1, ndisplay = 100, seed = 456)

# run the MPG clustering algorithm
ans = mpg(Y, C, prior=NULL, mcmc)

print(ans$chain$epsilon)

plotDiff(ans, type='weight')
plotDiff(ans, type = "shift")
cal = calibrate(ans)
par(mfrow=c(1,2))
plot(Y, col = C, pch=20, cex=0.5)
plot(cal$Y_cal, col = C, pch=20, cex=0.5)
