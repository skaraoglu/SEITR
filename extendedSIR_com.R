library(igraph)
library(deSolve)
library(Matrix)

# Model Parameters
Lambda <- 4           # Number of births
beta1 <- 0.8           # S to E rate of transfer
beta2 <- 0.3           # E to I rate of transfer
beta3 <- 0.02          # I to R rate of transfer
alpha1 <- 0.1          # I to Tt rate of transfer (alpha1 - beta3)
alpha2 <- 0.055        # Tt to R rate of transfer
delta_I <- 0.03        # Rate of death due to Infection
delta_T <- 0.03
mu <- 0.02             # Rate of natural death
# Lambda is the number of nodes added at each time step. (1 for this example)
# mu * N + delta_I * I is the number of nodes removed at each time step. (1 + 0.3 for this example)
# With these parameters, total population is expected to decrease slowly.

# Reproduction number
R0 <- (beta1 * beta2) / ((beta3 + mu + delta_I + alpha1) * (beta2 + mu) * (mu + delta_T + alpha2))

# Initial conditions
n <- 200
S <- 170
I <- 10
E <- 20
R <- 0
Tt <- 0
N <- S + E + I + R + Tt

# SEITR model
SEITR <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- Lambda - (beta1 * S * I) / N - mu * S
    dE <- (beta1 * S * I) / N - (beta2 + mu) * E
    dI <- beta2 * E - (beta3 + mu + delta_I + alpha1) * I
    dT <- alpha1 * I - (mu + delta_T + alpha2) * Tt
    dR <- beta3 * I + alpha2 * Tt - mu * R
    dN <- Lambda - N * mu - delta_I * I - delta_T * Tt
    list(c(dS, dE, dI, dT, dR, dN))
  })
}

# Initial state
state <- c(S = S, E = E, I = I, Tt = Tt, R = R, N = n)

# Parameters
parameters <- c(Lambda = Lambda, beta1 = beta1, beta2 = beta2, beta3 = beta3, alpha1 = alpha1, alpha2 = alpha2, delta_I = delta_I, delta_T = delta_T, mu = mu)

# Time steps
times <- seq(0, 100, by = 1)

# Solve ODE
out <- ode(y = state, times = times, func = SEITR, parms = parameters)

# Plot each status separately
par(mfrow = c(3, 2))
plot(out[, "time"], out[, "S"], type = "l", col = 2, xlab = "Time (Days)", ylab = "Susceptibles", main = "Susceptibles")
plot(out[, "time"], out[, "E"], type = "l", col = 3, xlab = "Time (Days)", ylab = "Exposed", main = "Exposed")
plot(out[, "time"], out[, "I"], type = "l", col = 4, xlab = "Time (Days)", ylab = "Infected", main = "Infected")
plot(out[, "time"], out[, "Tt"], type = "l", col = 5, xlab = "Time (Days)", ylab = "Treatment", main = "Treatment")
plot(out[, "time"], out[, "R"], type = "l", col = 6, xlab = "Time (Days)", ylab = "Recovered", main = "Recovered")
plot(out[, "time"], out[, "N"], type = "l", col = 7, xlab = "Time (Days)", ylab = "Total Population", main = "Total Population")

# Create the graph with 4 clusters and 50 nodes in each
npc <-50 # nodes per cluster
n_clust <- 4 # 4 clusters with 25 nodes each
matlist = list()
# Initialize the cluster vector
cluster_vector <- c()

for (i in 1:n_clust){ 
  matlist[[i]] = get.adjacency(erdos.renyi.game(npc, 0.7))
  
  # Assign each node in the current cluster to the cluster
  cluster_vector <- c(cluster_vector, rep(i, npc))
}

# merge clusters into one matrix
mat_clust <- bdiag(matlist)
k <- rowSums(mat_clust) 
print(any(k == 0))

node_vector <- seq(1,npc*n_clust)
for (i in node_vector){
  if (k[i]==0){ # if k=0, connect to something random
    j <- sample(node_vector[-i],1)
    mat_clust[i,j] <- 1
    mat_clust[j,i] <- 1
  }
}

g <- graph_from_adjacency_matrix(mat_clust, mode="undirected", diag=F)
comps <- components(g)
member_vec <- comps$membership
# daisy chain the components together
comp_ids <- seq(1,comps$no-1) # stop short of last one
for (curr_comp in comp_ids){
  # grab random nodes from consecutive components
  i <- sample(which(curr_comp==member_vec),1)
  j <- sample(which((curr_comp+1)==member_vec),1)
  
  # Determine the number of edges to add
  num_edges <- sample(1:5, 1)
  
  for (edge in 1:num_edges) {
    # Check if an edge already exists between the nodes
    if (!are_adjacent(g, i, j)) {
      # join them together
      mat_clust[i,j] <- 1
      mat_clust[j,i] <- 1
    }
    
    # Choose new nodes for the next edge
    i <- sample(which(curr_comp==member_vec),1)
    j <- sample(which((curr_comp+1)==member_vec),1)
  }
}

g <- graph_from_adjacency_matrix(mat_clust, mode="undirected", diag=F)
# Assign the cluster attribute to each vertex
V(g)$cluster <- cluster_vector

# Assign initial conditions to nodes
V(g)$label <- 1:n
V(g)$status <- "S"

# Select a cluster to initialize the infected and exposed nodes
init_cluster <- sample(1:n_clust, 1)

# Get the nodes in the selected cluster
init_cluster_nodes <- V(g)[member_vec == init_cluster]

# Randomly select nodes to be infected and exposed
init_I_nodes <- sample(init_cluster_nodes, I)
init_E_nodes <- sample(setdiff(init_cluster_nodes, init_I_nodes), E)

# Set the status of the selected nodes
V(g)[init_I_nodes]$status <- "I"
V(g)[init_E_nodes]$status <- "E"

# Initialize counters
S_count <- S
E_count <- E
I_count <- I
Tt_count <- Tt
R_count <- R
N_count <- N

# Initialize vectors to store metrics
degree_dist <- vector("list", length(times))
clustering_coeff <- numeric(length(times))
avg_path_length <- numeric(length(times))
largest_comp_size <- numeric(length(times))

# Color mapping for statuses
status_colors <- c("S" = "gray", "E" = "yellow", "I" = "red", "Tt" = "green", "R" = "purple")
par(mfrow = c(1, 1))

# Initial Network plot
plot(g, vertex.color = status_colors[V(g)$status], vertex.size=4, vertex.label="")
title(paste("Initial Network"))
status_counts <- table(V(g)$status)
status_labels <- paste(names(status_colors), " (", status_counts[names(status_colors)], ")", sep = "")
legend("bottomright", legend = status_labels, fill = status_colors, title = "Status")

# Initialize a counter for the node labels
node_counter <- vcount(g)

# Iterate over time steps
for (t in times) {
  # Save old statuses for tracking changes
  old_status <- V(g)$status
  
  # Remove nodes deceased due to Infection (count = delta_I * I)
  nodes_to_remove_count <- delta_I * sum(V(g)$status == "I")
  nodes_with_status_I <- which(V(g)$status == "I")
  
  # Check if the number of nodes to remove is less than 1
  if (nodes_to_remove_count < 1) {
    # Generate a random number between 0 and 1
    random_number <- runif(1)
    
    # Check if the random number is greater than nodes_to_remove_count
    if (random_number < nodes_to_remove_count) {
      node_to_remove <- sample(nodes_with_status_I, 1)
      
      # Print node to be removed
      cat("Time", t, ": Node", V(g)[node_to_remove]$label, "with status", V(g)[node_to_remove]$status, "removed due to Infection\n")
      
      # Remove the selected node
      g <- delete_vertices(g, node_to_remove)
    }
  } else {
    for (i in 1:nodes_to_remove_count){
      node_to_remove <- sample(nodes_with_status_I, 1)
      
      # Print node to be removed
      cat("Time", t, ": Node", V(g)[node_to_remove]$label, "with status", V(g)[node_to_remove]$status, "removed due to Infection\n")
      
      # Remove the selected node
      g <- delete_vertices(g, node_to_remove)
    }
  }
  
  # Remove nodes for deceased (count = mu * N)
  for (i in 1:(mu * N)) {
    node_to_remove <- sample(vcount(g), 1)
    
    # Print node to be removed
    cat("Time", t, ": Node", V(g)[node_to_remove]$label, "with status", V(g)[node_to_remove]$status, "removed\n")
    
    # Remove node from graph
    g <- delete_vertices(g, node_to_remove)
  }
  
  # Add new nodes to the network (count = Lambda)
  for (i in 1:(Lambda)) {
    node_counter <- node_counter + 1
    new_status <- "S"
    new_label <- node_counter
    
    g <- add_vertices(g, 1)
    V(g)[vcount(g)]$status <- new_status
    V(g)[vcount(g)]$label <- new_label
    
    # Assign new node to a cluster
    new_cluster <- sample(1:n_clust, 1)
    V(g)[vcount(g)]$cluster <- new_cluster  # Assign the new node to a cluster
    
    # Get the nodes in the same cluster
    same_cluster_nodes <- V(g)[V(g)$cluster == new_cluster]
    
    # Calculate the probability of attaching to each node in the same cluster
    prob_same_cluster <- degree(g, v = same_cluster_nodes, mode = "all") / sum(degree(g, v = same_cluster_nodes, mode = "all"))
    
    # Determine the number of nodes to attach
    num_nodes_to_attach <- min(sum(prob_same_cluster > 0), round(avg_degree))
    
    # Choose nodes to attach to within the same cluster
    nodes_to_attach_same_cluster <- sample(same_cluster_nodes, size = num_nodes_to_attach, prob = prob_same_cluster)
    
    # Add edges to these nodes
    for (node_to_attach in nodes_to_attach_same_cluster) {
      if (V(g)[node_to_attach]$label != new_label) { # No self edges
        g <- add_edges(g, c(vcount(g), node_to_attach))
      }
    }
    
    # With a small probability, add an edge to a node in a different cluster
    if (runif(1) < 0.1) { # Adjust this probability as needed
      different_cluster_nodes <- V(g)[V(g)$cluster != new_cluster]
      if(length(different_cluster_nodes) > 0) {
        node_to_attach_different_cluster <- sample(different_cluster_nodes, 1)
        g <- add_edges(g, c(vcount(g), node_to_attach_different_cluster))
      }
    }
    
    # Print new node
    cat("Time", t, ": New node", new_label, "added with status", new_status, "\n")
  }
  
  # Iterate over nodes for status changes
  for (i in V(g)) {
    # Get current status
    status <- V(g)[i]$status
    
    # Generate a random number
    rand <- runif(1)
    
    # Update status based on model parameters and current status
    if (status == "S") {
      # Get the statuses of the node's neighbors
      neighbor_statuses <- V(g)[neighbors(g, i)]$status
      # Check if any neighbor is E or I and rnd < beta1 * I / N
      if (any(neighbor_statuses %in% c("E", "I")) && rand < beta1 * I / N) {
        V(g)[i]$status <- "E"
        cat("Time", t, ": Node", V(g)[i]$label, "changed status from S to E\n")
      }
    } # Else if status = E, check if rnd < beta2 
    else if (status == "E" && rand < beta2) {
      V(g)[i]$status <- "I"
      cat("Time", t, ": Node", V(g)[i]$label, "changed status from E to I\n")
    } # Else if status = I, check if rnd < beta3, else if rnd < alpha1  
    else if (status == "I") {
      if (rand < beta3) {
        V(g)[i]$status <- "R"
        cat("Time", t, ": Node", V(g)[i]$label, "changed status from I to R\n")
      } else if (rand < alpha1) {
        V(g)[i]$status <- "Tt"
        cat("Time", t, ": Node", V(g)[i]$label, "changed status from I to Tt\n")
      }
    } # Else if status = Tt, check if rnd < alpha2 
    else if (status == "Tt" && rand < alpha2) {
      V(g)[i]$status <- "R"
      cat("Time", t, ": Node", V(g)[i]$label, "changed status from Tt to R\n")
    }
  }
  
  # Store counts
  S_count[t+1] <- sum(V(g)$status == "S")
  E_count[t+1] <- sum(V(g)$status == "E")
  I_count[t+1] <- sum(V(g)$status == "I")
  Tt_count[t+1] <- sum(V(g)$status == "Tt")
  R_count[t+1] <- sum(V(g)$status == "R")
  N_count[t+1] <- vcount(g)
  
  # Calculate and store metrics
  degree_dist[[t+1]] <- degree_distribution(g)
  clustering_coeff[t+1] <- transitivity(g, type = "global")
  avg_path_length[t+1] <- mean_distance(g, directed = FALSE)
  comps <- components(g)
  largest_comp_size[t+1] <- max(comps$csize)
  
  # Plot network in every 10th time step
  if (t %% 10 == 0) {
    plot(g, vertex.color = status_colors[V(g)$status], vertex.size=4, vertex.label="")
    title(paste("Time Step:", t))
    status_counts <- table(V(g)$status)
    status_labels <- paste(names(status_colors), " (", status_counts[names(status_colors)], ")", sep = "")
    legend("bottomright", legend = status_labels, fill = status_colors, title = "Status")
  }
}

# Set up a single plot for comparison
par(mfrow = c(3, 2))
plot(out[, "time"], out[, "S"], type = "l", col = status_colors["S"], xlab = "Time (Days)", ylab = "Susceptibles", main = "Comparison of Susceptibles")
lines(times, S_count, col = status_colors["S"])
plot(out[, "time"], out[, "E"], type = "l", col = "orange", xlab = "Time (Days)", ylab = "Exposed", main = "Comparison of Exposed")
lines(times, E_count, col = "orange")
plot(out[, "time"], out[, "I"], type = "l", col = status_colors["I"], xlab = "Time (Days)", ylab = "Infected", main = "Comparison of Infected")
lines(times, I_count, col = status_colors["I"])
plot(out[, "time"], out[, "Tt"], type = "l", col = status_colors["Tt"], xlab = "Time (Days)", ylab = "Treatment", main = "Comparison of Treatment")
lines(times, Tt_count, col = status_colors["Tt"])
plot(out[, "time"], out[, "R"], type = "l", col = status_colors["R"], xlab = "Time (Days)", ylab = "Recovered", main = "Comparison of Recovered")
lines(times, R_count, col = status_colors["R"])
plot(out[, "time"], out[, "N"], type = "l", col = "blue", xlab = "Time (Days)", ylab = "Total Population", main = "Comparison of Total Population")
lines(times, N_count, col = "blue")

# Plot metrics
par(mfrow = c(2, 2))
plot(degree_dist[[t+1]], main = paste("Degree Distribution at Time", t), xlab = "Degree", ylab = "Frequency")
plot(times, clustering_coeff, type = "l", xlab = "Time (Days)", ylab = "Clustering Coefficient", main = "Clustering Coefficient")
plot(times, avg_path_length, type = "l", xlab = "Time (Days)", ylab = "Average Path Length", main = "Average Path Length")
plot(times, largest_comp_size, type = "l", xlab = "Time (Days)", ylab = "Size of Largest Component", main = "Size of Largest Component")
