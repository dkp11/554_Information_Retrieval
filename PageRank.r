# Load the igraph library for processing 
library(igraph)

# Load the input file and form the required adjacency matrix.
# Expected format: adjacency matrix
# Expected file : AM.txt
inputdata <- as.matrix(read.csv("AM.txt", sep =" ", header=FALSE))

# Get the graph from the edge list.
inputGraph <- graph.edgelist(inputdata[,-3],directed=TRUE) 

# Set the dampening factor d = 0.85
d <- 0.85

# Gat the adjacency matrix from the input graph. 
G <- get.adjacency(inputGraph, edges=FALSE, names=FALSE,sparse = TRUE) 

# Due to memory limitation, process only 1000 rows and columns.
# Uncomment the following code only if you are running in a machine with lesser memory.
# Uncommenting this will give a different list altogether.
# G <- G[1:10000,1:10000]

# Get the node count
n <- nrow(G)

###################################################
######## Method 1 begins here #####################
###################################################

# Preprocess/Normalize the data in the graph so that it can be used
# for matrix and Eigen processing.
cvec <- apply(G,2,sum)

# nodes with indegree 0 will cause problems if we divide by 0.
cvec[cvec==0] <- 1 

# Calculate the delta and create the unity matrix same as G
delta <- (1-d)/n
A <- matrix(delta, n, n)

# Calculate the eigen value and vector and the page rank
for (i in 1:n) 
  A[i,] <- A[i,] + d*G[i,]/cvec

# Get the Eigen vectors
x <- Re(eigen(A)$vector[,1])

# Normalize, Sort and display the top 100 results.
PRVect <- order(x/sum(x), decreasing = TRUE)[1:100]
PRVect
x[PRVect]
###################################################
######## Method 1 ends here #######################
###################################################

###################################################
######## Method 2 begins here #####################
###################################################

# Get G as a 2D matrix. 
M <- as.matrix(G)

# Get the node count
n <- nrow(G)

# Normalize matrix.
M <- t(M / rowSums(M))

# create a unity matrix of the same size.
U <- matrix(data=rep(1/n, n^2), nrow=n, ncol=n)

# Calculate the page rank.
A <- d*M+(1-d)*U
rm(M)
rm(U)
A[!is.finite(A)] <- 0
v <- eigen(A)$vector[,1]

# Normalize, Sort and display the top 100 results.
PRVect <- order(as.numeric(v)/sum(as.numeric(v)), decreasing = TRUE)[1:100]
PRVect
v[PRVect]

###################################################
######## Method 2 ends here #######################
###################################################

###################################################
######## Method 3 begins here #######################
###################################################

# This is the precision where convergence will take place.
t = 20

# normalize adjacency matrix by row sum and compute dangling node vector
a = rep(0, n)

# compute row sums (matrix multiplication)
rs = G %*% rep(1,n)
for (i in 1:n) {
  if (rs[i] == 0) {a[i] = 1} 
  else { G[i,] = G[i,] / rs[i] }  
}

# create a unity matrix of the same size.
b <- as.matrix(rep(1/n, n))

# Initialize values.
e = rep(1, n)
pi0 = rep(0, n)
pi1 = rep(1/n, n)
eps = 1/10^t
iter = 0

# Calculate page rank till convergence.
while (sum(abs(pi0 - pi1)) > eps) {
  pi0 = pi1
  pi1 = d * pi1 %*% G + (d * pi1 %*% a + (1 - d) * pi1 %*% e) * b
  iter = iter + 1
} 

# Return the normalized ist, sorted by top 100 page ranks
pi2 = pi1 / sum(pi1)
list(pi = pi2, iter = iter)
PRVect <- (order(pi2, decreasing = TRUE))[1:100]
PRVect

###################################################
######## Method 3  ends here ######################
###################################################

