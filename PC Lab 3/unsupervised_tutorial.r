########################  Load Data  ########################

### Load data
load("rollcall-votes.Rdata")
load("rollcall-members.Rdata")

print('Data loaded.')

##############################################################

print('# Counts of Democrats, Republicans and one special politician')
table(members$party)

print('# Shares of Democrats, Republicans and one special politician')
round(table(members$party)/nrow(members),3)

# Count missing votings for each politician and plot the counts
missings <- rowSums(votes[,(1:ncol(votes))]==0)

# No. politicians who always voted
sum(missings == 0)

# Shares of missing votings
s_missings <- missings/(ncol(votes)-1)

# Histogram with 100 bins
hist(???, breaks = 100)

# Counts - yes and nos
yeas <- rowSums(votes[,(1:ncol(votes))]== ???)
nays <- rowSums(votes[,(1:ncol(votes))]== ???)

# Plots - Party
plot(yeas, nays, col = members$party)
legend('topleft', legend = levels(members$party), col = 1:3,  pch = 1)

# PCA
pr.out = prcomp(??? , center = TRUE, scale = TRUE)

# No of principal components
dim(pr.out$rotation)[2]

# variance explained by each component
pr.var = pr.out$sdev^2

# Proportion of variance explained
pve=pr.var/sum(pr.var)

# Print first 10 PC
pve[1:10]

# Plot the first 10 PC
barplot(pve[1:10], xlab=" Principal Component ", ylab=" Proportion of Variance Explained ", ylim=c(0,1))
barplot(cumsum(pve[1:10]), xlab=" Principal Component ", ylab ="Cumulative Proportion of Variance Explained ", ylim=c(0,1))

# Plot the first two principal components, color the party membership
plot(pr.out$x[,1], pr.out$x[,2], xlab = "PC1", ylab = "PC2", col = members$party, main = "Top two PC directions")
legend('bottomright', legend = levels(members$party), col = 1:3,  pch = 1)

## Far right (very conservative)
head(sort(???))

## Far left (very liberal)
head(sort(???, decreasing=???))

# PC 2
head(sort(???))
# No clear pattern based on party and state information

# Look at the largest loadings in PC2 to discern an interpretation.
loadings <- pr.out$rotation
loadings[order(abs(loadings[,2]), decreasing=TRUE)[1:5],2]

# Analyze voting behavior
table(votes[,1146])
table(votes[,658])
table(votes[,1090])

# Either everyone voted "yea" or missed the voting.
# These votes all correspond to near-unanimous symbolic action.

# Mystery Solved: the second PC is just attendance!
head(sort(rowSums(votes==0), decreasing=TRUE))

set.seed(11122019)

# K-means clustering with 2 clusters
km.out = kmeans(???, 2, nstart = 20)
km.out$cluster

# Tabulate party vs cluster
table(members$party, km.out$cluster)

# How to analyze the optimal number of clusters

sse <- c()
sse[1] <- Inf

for (ind_cl in c(2:20)) {
  set.seed(3)
  km.out = kmeans (votes, ind_cl, nstart = 20)
  sse[ind_cl] = km.out$tot.withinss
}

plot(sse)
# Optimum 4-5 clusters

# Plot the 5 clusters on the PC components graph
set.seed(3)
km.out = kmeans (???, ???, nstart = 20)

# Plot the first two principal components color the party membership
plot(pr.out$x[,1], pr.out$x[,2], xlab = "PC1", ylab = "PC2", col = km.out$cluster, main = "Top two PC directions with 5 clusters")
legend('bottomright', legend = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col = 1:5,  pch = 1)

# Analyzing how the number of starts work
set.seed (3)
print('With nstart = 1')
km.out = kmeans (votes,6, nstart = ???)
km.out$tot.withinss

print('With nstart = 20')
km.out =kmeans (votes,6, nstart = ???)
km.out$tot.withinss
