

##########################################################
##############  Plot the network  ########################
##########################################################
rm(list=ls()) # Clear memory
load("/Users/sushil/Documents/NetworkA/hypertext.RData") # Load the data
ls() # list the names 
library(igraph)


###### Question 1 ###### 
gra <- graph.adjacency(adj, mode = "undirected") # transform adjacency matrix into igraph object
summary(gra)
dim(gra)

set.seed(12345)
layout <- layout.fruchterman.reingold(gra) # set the positions for the nodes which will be used in all plots


plot.igraph(gra, vertex.label = NA, vertex.size = 1.8*sqrt(rowSums(adj)), vertex.color = c("green","red","yellow","blue"), vertex.frame.color="black", layout = layout)
#plot(gra, vertex.label = NA, vertex.size = 2.5*sqrt(rowSums(adj)),  , layout = layout)
legend("bottomright", c("Group one","Group second","Group three","Group four"), pch = 20, cex = 0.5, pt.cex = 1.2,col = c("green","red","yellow","blue"))
E(gra)
V(gra)

###### Question 2 ######

##########################################################
##############  Statistics of the network  ###############
##########################################################


n_nodes <- nrow(adj)
n_edges <- sum(adj)/2

n_possible_edges <- n_nodes * (n_nodes-1) / 2

n_triangles <- sum(diag(adj %*% adj %*% adj)) / 6
p <- n_edges / n_possible_edges

Degrees <- rowSums(adj)
Degree_freqs <- table(Degrees)
plot(Degree_freqs, type = "h", lwd = 4,col = "blue",panel.first = grid())


###### Question 3 ######


###########################################################
################  Spectral Clustering  ####################
###########################################################

eigen_dec <- eigen(adj) # eigenvalues and eigenvectors
plot(eigen_dec$values, pch = 20,col = "black",panel.first = grid()) # make a decision on K
abline(v=3,col="red",lwd = 3) # add line 
K <- 4 # In this case the choice K = 4 is the most reasonable
embedding <- eigen_dec$vectors[,1:K] # project the nodes into a low-dimensional latent space
memberships <- kmeans(embedding, K, nstart = 100)$cluster # cluster the nodes
plot(gra, vertex.label = NA, vertex.size = 2.5*sqrt(rowSums(adj)), vertex.color = c("green","red","yellow","blue"), layout = layout)
res <- make_clusters(gra, memberships)
set.seed(12)

plot(x = res, y = gra, layout = layout, vertex.label = NA, vertex.size = 2.5*sqrt(rowSums(adj)))
legend("bottomright", c("Group one","Group second","Group three","Group four"), pch = 20, cex = 0.8, pt.cex = 1,col = categorical_pal(8)[1:5])


##### Question 4 ######

###########################################################
################  Stochastic Blockmodels  #################
###########################################################


library(blockmodels) # This package uses S4 classes
sbm <- BM_bernoulli(membership_type = "SBM_sym", # data type
                    adj = adj, # adjacency matrix
                    verbosity = 1, # how much should be printed out while running the algorithm
                    plotting="", # could be used to show how the values of the ICL evolve during the estimation
                    explore_max = 4) # maximum number of clusters to consider

sbm$estimate() # this line runs the VEM on the dataset
K_star <- which.max(sbm$ICL) # extract the best model according to ICL values
soft_clustering <- sbm$memberships[[K_star]]$Z # posterior probabilities for group memberships ("taus")
hard_clustering <- apply(soft_clustering, 1, which.max) # maximum a posteriori clustering
sbm$memberships[[K_star]]$plot() # plot group sizes
sbm$plot_parameters(K_star) # plot connection probabilities
sbm$model_parameters[[K_star]]$pi[1,2]
soft_clustering
hard_clustering


memberships_sbm <- make_clusters(gra, hard_clustering) # create the new communities object according to partition found
plot(memberships_sbm, mark.groups = NULL, edge.color = NULL, gra, vertex.label = NA, vertex.size = 2.9*sqrt(rowSums(adj)), layout = layout)
legend("bottomright", c("Group one","Group second","Group three","Group four"), pch = 20, cex = 0.8, pt.cex = 1,col = categorical_pal(8)[1:5])

#the relative sizes of the group 
colSums(soft_clustering) / sum(colSums(soft_clustering))
table(hard_clustering) / n_nodes

#the probability that an edge appears between a node in group 1 and a node in group 2
sbm$model_parameters[[K_star]]$pi[1,2]

#The groups exhibit a strong community structure and groups do not
sbm$plot_parameters(K_star)# group 1 and 2 exhibit a clear community structure

# The group has the highest expected degree
lambda <- colSums(soft_clustering) / sum(colSums(soft_clustering))
lambda %*% sbm$model_parameters[[K_star]]$pi # Group two has more connections as compare to others.

# The expected number of edges for a node in group 1
n_nodes * lambda %*% sbm$model_parameters[[K_star]]$pi # around 12

table(memberships, hard_clustering)
# the groups with community structure are similar in both partitions, the group of hubs is where the partitions are fundamentally different

