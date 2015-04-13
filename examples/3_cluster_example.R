library("MPG")

palette(rainbow(20, start=0.0, end=0.85, alpha=0.5))

# set a random seed to generate consistent data sets for testing
set.seed(123)

# number of events to create for each cluster
n_events_clusterA = 4000
n_events_clusterB = 4000
n_events_clusterC = 2000

# number of columns in our data sets
n_parameters = 2

# generate 1st data set
clusterA_eventsX = rnorm(n_events_clusterA)
clusterA_eventsY = rnorm(n_events_clusterA)
clusterA_events = matrix(c(clusterA_eventsX, clusterA_eventsY), ncol = n_parameters)

clusterB_eventsX = rnorm(n_events_clusterB) + 15
clusterB_eventsY = rnorm(n_events_clusterB) + 7
clusterB_events = matrix(c(clusterB_eventsX, clusterB_eventsY), ncol = n_parameters)

clusterC_eventsX = rnorm(n_events_clusterC) + 7
clusterC_eventsY = rnorm(n_events_clusterC) + 15
clusterC_events = matrix(c(clusterC_eventsX, clusterC_eventsY), ncol = n_parameters)

data_set_1 = rbind(clusterA_events, clusterB_events, clusterC_events)

# generate 2nd data set
clusterA_eventsX = rnorm(n_events_clusterA)
clusterA_eventsY = rnorm(n_events_clusterA)
clusterA_events = matrix(c(clusterA_eventsX, clusterA_eventsY), ncol = n_parameters)

clusterB_eventsX = rnorm(n_events_clusterB) + 13
clusterB_eventsY = rnorm(n_events_clusterB) + 6
clusterB_events = matrix(c(clusterB_eventsX, clusterB_eventsY), ncol = n_parameters)

clusterC_eventsX = rnorm(n_events_clusterC) + 6
clusterC_eventsY = rnorm(n_events_clusterC) + 13
clusterC_events = matrix(c(clusterC_eventsX, clusterC_eventsY), ncol = n_parameters)

data_set_2 = rbind(clusterA_events, clusterB_events, clusterC_events)

# shuffle the data sets
data_set_1 = data_set_1[sample(nrow(data_set_1)),]
data_set_2 = data_set_2[sample(nrow(data_set_2)),]

# combine data sets into 1 matrix
Y = rbind(data_set_1, data_set_2)

# map for identifying the data set for each event
C = c(rep(1, nrow(data_set_1)), rep(2, nrow(data_set_2)))

burn_in_list = c(1000)
#     c(
#     10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
#     200, 300, 400, 500, 600, 700, 800, 900, 1000,
#     2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
#     15000, 20000, 25000, 30000, 35000, 40000  
# )

for (burn_count in burn_in_list) {
    
    # define our clustering parameters
    mcmc = list(nburn = burn_count, nsave = 1, nskip = 1, ndisplay = 100, seed = 456)
    
    # run the MPG clustering algorithm
    ans = mpg(Y, C, prior=NULL, mcmc)
    
    cat("Epsilon: ", ans$chain$epsilon, "\n")
    
    # get the unique clusters
    unique_clusters = unique(ans$chain$Z[ans$mcmc$nsave,1:(nrow(Y))])
    
    # print unique cluster labels from last saved run
    cat("Unique clusters: ", unique_clusters, "\n")

    for (clust in unique_clusters) {
        cat(
            "\t", clust, ": ",
            round(ans$chain$mu_0[1, clust, 1], digits=3),
            ", ",
            round(ans$chain$mu_0[2, clust, 1], digits=3),
            "\n"
        )
    }
    
    par(mfrow=c(2,2))
    
    # plot classifications
    plot(data_set_1, pch=20, cex=0.75, col=ans$chain$Z[1, 1:(nrow(Y)/2)], xlim=c(-5,20), ylim=c(-5,20))
    plot(data_set_2, pch=20, cex=0.75, col=ans$chain$Z[1, (nrow(Y)/2 + 1):(nrow(Y))], xlim=c(-5,20), ylim=c(-5,20))
    
    # align clusters
    cal = calibrate(ans)
    
    # plot both data sets before and after alignment
    plot(Y, col=ans$chain$Z[1, 1:(nrow(Y))], pch=20, cex=0.75, xlim=c(-5,20), ylim=c(-5,20))
    plot(cal$Y_cal, col=ans$chain$Z[1, 1:(nrow(Y))], pch=20, cex=0.75, xlim=c(-5,20), ylim=c(-5,20))
    
    title = paste("burn: ", burn_count)
    mtext(title, side=3, line=-2, outer=TRUE)
    clusters_text = paste(unique_clusters, collapse=", ")
    clusters_text = paste("clusters: ", clusters_text)
    mtext(clusters_text, side=3, line=-3, outer=TRUE)
}
