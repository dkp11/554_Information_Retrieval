###################################################
######## Homework 2 - PageRank Algorithm  #########
#########   Information Retrieval   ###############
###################################################
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
# We can either use the 'power' method or the 'eigen' method.
# We will get negative values in eigen but ranking wont be affected.
method <- 'power' 

# we can either choose fixed iteration or use while loop for 
# convergence. We are going for fixed iterations here as eigen
# has convergence built in.
iter <- 100

# Preprocess/Normalize the data in the graph so that it can be used
# for matrix and Eigen processing.
cvec <- apply(G,2,sum)

# nodes with indegree 0 will cause problems if we divide by 0.
cvec[cvec==0] <- 1 
gvec <- apply(G,1,sum)

# Calculate the delta and create the unity matrix same as G
delta <- (1-d)/n
A <- matrix(delta, n, n)

# Calculate the eigen value and vector and the page rank
for (i in 1:n) 
  A[i,] <- A[i,] + d*G[i,]/cvec

# Get the Eigen vectors
if (method == 'power'){
  x <- rep(1,n)
  for (i in 1: iter) x <- A%*%x
} else { # method == eigen
  x <- Re(eigen(A)$vector[,1])
}

n = 1000
#select the number of results to be displayed
if( n > 100) n <- 100

# Normalize, Sort and display the top 100 results.
PRVect <- order(x/sum(x), decreasing = TRUE)[1:n]
PRVect
abs(x[PRVect])
###################################################
######## Method 1 ends here #######################
###################################################

###################################################
######## Method 2 begins here #######################
# This method should be used for verification purposes only
###################################################

page.rank(inputGraph,directed = TRUE, damping = 0.85)

###################################################
######## Method 2 ends here ######################
###################################################
