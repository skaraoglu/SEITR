library(igraph)
library(deSolve)

# Model Parameters
Lambda <- 1            # Number of births
beta1 <- 0.8           # S to E rate of transfer
beta2 <- 0.3           # E to I rate of transfer
beta3 <- 0.02          # I to R rate of transfer
alpha1 <- 0.1          # I to Tt rate of transfer (alpha1 - beta3)
alpha2 <- 0.055        # Tt to R rate of transfer
delta_I <- 0.03        # Rate of death due to Infection
delta_T <- 0.03
mu <- 0.01             # Rate of natural death
# Lambda is the number of nodes added at each time step. (1 for this example)
# mu * N + delta_I * I is the number of nodes removed at each time step. (1 + 0.3 for this example)
# With these parameters, total population is expected to decrease slowly.

# Reproduction number
R0 <- (beta1 * beta2) / ((beta3 + mu + delta_I + alpha1) * (beta2 + mu) * (mu + delta_T + alpha2))

# Initial conditions
n <- 100
S <- 70
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

# Create the graph
g <- erdos.renyi.game(n, 0.85)
V(g)$label <- 1:n
V(g)$status <- "S"
V(g)$status[sample(1:n, size = I)] <- "I"
V(g)$status[sample(setdiff(1:n, which(V(g)$status == "I")), size = E)] <- "E"
V(g)$status[sample(setdiff(1:n, which(V(g)$status %in% c("I", "E"))), size = Tt)] <- "Tt"
V(g)$status[sample(setdiff(1:n, which(V(g)$status %in% c("I", "E", "Tt"))), size = R)] <- "R"

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
plot(g, vertex.color = status_colors[V(g)$status], vertex.label = V(g)$label, vertex.size=10)
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
    new_label <- node_counter  # Use the counter for the label
    
    g <- add_vertices(g, 1)
    V(g)[vcount(g)]$status <- new_status
    V(g)[vcount(g)]$label <- new_label  # Assign label to new node
    
    # Get the degree of each node
    degree <- degree(g, mode = "all")
    # Calculate the average degree
    avg_degree <- mean(degree)
    # Calculate the probability of attaching to each node
    prob <- degree / sum(degree)
  
    # Choose nodes to attach to
    nodes_to_attach <- sample(V(g), size = min(round(avg_degree), vcount(g)), prob = prob)
    
    # Add edges to these nodes
    for (node_to_attach in nodes_to_attach) {
      if (V(g)[node_to_attach]$label != new_label) { # No self edges
        g <- add_edges(g, c(vcount(g), node_to_attach))
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
    plot(g, vertex.color = status_colors[V(g)$status], vertex.size=10, vertex.label=V(g)$label)
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

