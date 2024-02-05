# Load the necessary libraries
library(igraph)  # For network analysis and visualization
library(deSolve)  # For solving differential equations
library(ggplot2)
library(dplyr)


# Define the model parameters
Lambda <- 1.1  # Lambda represents the number of births in the population
beta1 <- 0.8  # Beta1 is the rate of transfer from Susceptible (S) to Exposed (E)
beta2 <- 0.18  # Beta2 is the rate of transfer from Exposed (E) to Infected (I)
beta3 <- 0.02  # Beta3 is the rate of transfer from Infected (I) to Recovered (R)
alpha1 <- 0.1  # Alpha1 is the rate of transfer from Infected (I) to Treatment (Tt)
alpha2 <- 0.055  # Alpha2 is the rate of transfer from Treatment (Tt) to Recovered (R)
delta_I <- 0.03  # Delta_I is the rate of death due to Infection
delta_T <- 0.03  # Delta_T is the rate of death during Treatment
mu <- 0.01  # Mu is the rate of natural death

# Calculate the reproduction number (R0)
# R0 is a key parameter in epidemiology, representing the average number of secondary infections produced by a typical case of an infection in a population where everyone is susceptible.
R0 <- (beta1 * beta2) / ((beta3 + mu + delta_I + alpha1) * (beta2 + mu) * (mu + delta_T + alpha2))

# Define the initial conditions
n <- 100  # Total population
S <- 85  # Number of susceptible individuals
I <- 10  # Number of infected individuals
E <- 5  # Number of exposed individuals
R <- 0  # Number of recovered individuals
Tt <- 0  # Number of individuals under treatment
N <- S + E + I + R + Tt  # Total population (should be equal to n)

# Define the SEITR model
# This function represents the system of differential equations for the SEITR model.
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

# Define the initial state
state <- c(S = S, E = E, I = I, Tt = Tt, R = R, N = n)

# Define the parameters
parameters <- c(Lambda = Lambda, beta1 = beta1, beta2 = beta2, beta3 = beta3, alpha1 = alpha1, alpha2 = alpha2, delta_I = delta_I, delta_T = delta_T, mu = mu)

# Define the time steps
times <- seq(0, 100, by = 1)

# Solve the system of ODEs using the 'ode' function from the 'deSolve' package
out <- ode(y = state, times = times, func = SEITR, parms = parameters)

# Plot each status separately
# This section of the code creates six separate plots for each status (S, E, I, Tt, R, N) over time.
par(mfrow = c(3, 2))  # Set up a 3x2 grid for the plots
# Each plot shows the number of individuals in a particular status at each time step.
plot(out[, "time"], out[, "S"], type = "l", col = 2, xlab = "Time (Days)", ylab = "Susceptibles", main = "Susceptibles")
plot(out[, "time"], out[, "E"], type = "l", col = 3, xlab = "Time (Days)", ylab = "Exposed", main = "Exposed")
plot(out[, "time"], out[, "I"], type = "l", col = 4, xlab = "Time (Days)", ylab = "Infected", main = "Infected")
plot(out[, "time"], out[, "Tt"], type = "l", col = 5, xlab = "Time (Days)", ylab = "Treatment", main = "Treatment")
plot(out[, "time"], out[, "R"], type = "l", col = 6, xlab = "Time (Days)", ylab = "Recovered", main = "Recovered")
plot(out[, "time"], out[, "N"], type = "l", col = 7, xlab = "Time (Days)", ylab = "Total Population", main = "Total Population")

# Create the graph
# This section of the code creates an initial network (graph) using the Erdos-Renyi model.
init_p <- 0.9  # The probability of creating an edge between any two nodes
g <- erdos.renyi.game(n, p = init_p)  # Create the graph
# Assign labels and statuses to the nodes
V(g)$label <- 1:n  # Label the nodes with numbers from 1 to n
V(g)$status <- "S"  # Initially, all nodes are Susceptible (S)
# Randomly assign the other statuses (E, I, Tt, R) to the nodes
V(g)$status[sample(1:n, size = I)] <- "I"
V(g)$status[sample(setdiff(1:n, which(V(g)$status == "I")), size = E)] <- "E"
V(g)$status[sample(setdiff(1:n, which(V(g)$status %in% c("I", "E"))), size = Tt)] <- "Tt"
V(g)$status[sample(setdiff(1:n, which(V(g)$status %in% c("I", "E", "Tt"))), size = R)] <- "R"

# Initialize counters for each status
S_count <- S
E_count <- E
I_count <- I
Tt_count <- Tt
R_count <- R
N_count <- N

# Initialize vectors to store metrics
# These vectors will store various network metrics at each time step.
degree_dist <- vector("list", length(times))  # Degree distribution
clustering_coeff <- numeric(length(times))  # Clustering coefficient
avg_path_length <- numeric(length(times))  # Average path length
largest_comp_size <- numeric(length(times))  # Size of the largest connected component

# Define a color mapping for the statuses
status_colors <- c("S" = "gray", "E" = "yellow", "I" = "red", "Tt" = "green", "R" = "purple")

# Plot the initial network
# The nodes are colored according to their status.
par(mfrow = c(1, 1))  # Set up a 1x1 grid for the plot
plot(g, vertex.color = status_colors[V(g)$status], vertex.label = V(g)$label, vertex.size=10)
title(paste("Initial Network"))
# Create a legend showing the number of nodes in each status
status_counts <- table(V(g)$status)
status_labels <- paste(names(status_colors), " (", status_counts[names(status_colors)], ")", sep = "")
legend("bottomright", legend = status_labels, fill = status_colors, title = "Status")

# Initialize a counter for the node labels
node_counter <- vcount(g)  # The current number of nodes in the graph


# Iterate over time steps
# This loop simulates the spread of the disease over time.
for (t in times) {
  # Save old statuses for tracking changes
  # This line saves the current status of each node so that we can track how the status changes over time.
  old_status <- V(g)$status
  
  # Calculate the number of nodes to remove due to Infection
  # This line calculates the number of nodes (individuals) that should be removed from the graph due to Infection.
  nodes_to_remove_count <- delta_I * sum(V(g)$status == "I")
  
  # Separate the integer and decimal parts
  # This line separates the integer and decimal parts of the number of nodes to remove.
  # The integer part represents the number of nodes that will definitely be removed, while the decimal part represents a probability that an additional node will be removed.
  floor_value <- floor(nodes_to_remove_count)
  fractional_part <- nodes_to_remove_count - floor_value
  
  # Get the nodes with status "I"
  # This line gets the IDs of the nodes that are currently Infected.
  nodes_with_status_I <- which(V(g)$status == "I")
  
  # If there are nodes to remove due to Infection
  if (floor_value > 0 && length(nodes_with_status_I) > 0){
    # Remove nodes based on the integer part
    # This loop removes a number of nodes equal to the integer part of the number of nodes to remove.
    for (i in 1:floor_value) {
      node_to_remove <- sample(nodes_with_status_I, 1)  # Select a node to remove at random
      cat("Time", t, ": Node", V(g)[node_to_remove]$label, "with status", V(g)[node_to_remove]$status, "removed due to Infection\n")  # Print a message indicating which node was removed
      g <- delete_vertices(g, node_to_remove)  # Remove the node from the graph
    }
  }
  
  # If there is a probability of removing an additional node due to Infection
  if (fractional_part > 0){
    # Use the decimal part as a probability to remove an additional node
    # This line removes an additional node with a probability equal to the decimal part of the number of nodes to remove.
    if (runif(1) < fractional_part) {
      node_to_remove <- sample(nodes_with_status_I, 1)  # Select a node to remove at random
      cat("Time", t, ": Node", V(g)[node_to_remove]$label, "with status", V(g)[node_to_remove]$status, "removed due to Infection\n")  # Print a message indicating which node was removed
      g <- delete_vertices(g, node_to_remove)  # Remove the node from the graph
    }
  }
  
  # Calculate the number of nodes to remove during Treatment
  # This line calculates the number of nodes (individuals) that should be removed from the graph during Treatment.
  nodes_to_remove_count <- delta_T * sum(V(g)$status == "Tt")
  
  # Separate the integer and decimal parts
  # This line separates the integer and decimal parts of the number of nodes to remove.
  # The integer part represents the number of nodes that will definitely be removed, while the decimal part represents a probability that an additional node will be removed.
  floor_value <- floor(nodes_to_remove_count)
  fractional_part <- nodes_to_remove_count - floor_value
  
  # Get the nodes with status "Tt"
  # This line gets the IDs of the nodes that are currently under Treatment.
  nodes_with_status_T <- which(V(g)$status == "Tt")
  
  # If there are nodes to remove during Treatment
  if (floor_value > 0){
    # Remove nodes based on the integer part
    # This loop removes a number of nodes equal to the integer part of the number of nodes to remove.
    for (i in 1:floor_value) {
      node_to_remove <- sample(nodes_with_status_T, 1)  # Select a node to remove at random
      cat("Time", t, ": Node", V(g)[node_to_remove]$label, "with status", V(g)[node_to_remove]$status, "removed during Treatment\n")  # Print a message indicating which node was removed
      g <- delete_vertices(g, node_to_remove)  # Remove the node from the graph
    }
  }
  
  # If there is a probability of removing an additional node during Treatment
  if (fractional_part > 0){
    # Use the decimal part as a probability to remove an additional node
    # This line removes an additional node with a probability equal to the decimal part of the number of nodes to remove.
    if (runif(1) < fractional_part) {
      node_to_remove <- sample(nodes_with_status_T, 1)  # Select a node to remove at random
      cat("Time", t, ": Node", V(g)[node_to_remove]$label, "with status", V(g)[node_to_remove]$status, "removed during Treatment\n")  # Print a message indicating which node was removed
      g <- delete_vertices(g, node_to_remove)  # Remove the node from the graph
    }
  }
  
  # Calculate the number of nodes to remove due to natural death
  # This line calculates the number of nodes (individuals) that should be removed from the graph due to natural death.
  nodes_to_remove_count <- mu * vcount(g)
  
  # Separate the integer and decimal parts
  # This line separates the integer and decimal parts of the number of nodes to remove.
  # The integer part represents the number of nodes that will definitely be removed, while the decimal part represents a probability that an additional node will be removed.
  floor_value <- floor(nodes_to_remove_count)
  fractional_part <- nodes_to_remove_count - floor_value
  
  # If there are nodes to remove due to natural death
  if (floor_value > 0){
    # Remove nodes based on the integer part
    # This loop removes a number of nodes equal to the integer part of the number of nodes to remove.
    for (i in 1:floor_value) {
      node_to_remove <- sample(vcount(g), 1)  # Select a node to remove at random
      cat("Time", t, ": Node", V(g)[node_to_remove]$label, "with status", V(g)[node_to_remove]$status, "removed\n")  # Print a message indicating which node was removed
      g <- delete_vertices(g, node_to_remove)  # Remove the node from the graph
    }
  }
  
  # If there is a probability of removing an additional node due to natural death
  if (fractional_part > 0){
    # Use the decimal part as a probability to remove an additional node
    # This line removes an additional node with a probability equal to the decimal part of the number of nodes to remove.
    if (runif(1) < fractional_part) {
      node_to_remove <- sample(vcount(g), 1)  # Select a node to remove at random
      cat("Time", t, ": Node", V(g)[node_to_remove]$label, "with status", V(g)[node_to_remove]$status, "removed\n")  # Print a message indicating which node was removed
      g <- delete_vertices(g, node_to_remove)  # Remove the node from the graph
    }
  }
  
  # Calculate the number of nodes to add
  # This line calculates the number of nodes (individuals) that should be added to the graph.
  nodes_to_add_count <- Lambda
  
  # Separate the integer and decimal parts
  # This line separates the integer and decimal parts of the number of nodes to add.
  # The integer part represents the number of nodes that will definitely be added, while the decimal part represents a probability that an additional node will be added.
  floor_value <- floor(nodes_to_add_count)
  fractional_part <- nodes_to_add_count - floor_value
  
  # If there are nodes to add
  if (floor_value > 0){
    # Add nodes based on the integer part
    # This loop adds a number of nodes equal to the integer part of the number of nodes to add.
    for (i in 1:floor_value) {
      node_counter <- node_counter + 1  # Increment the node counter
      new_status <- "S"  # The status of the new nodes is Susceptible (S)
      new_label <- node_counter  # The label of the new nodes is the current value of the node counter
      
      g <- add_vertices(g, 1)  # Add a new node to the graph
      V(g)[vcount(g)]$status <- new_status  # Assign the status to the new node
      V(g)[vcount(g)]$label <- new_label  # Assign the label to the new node
      
      # Choose nodes to attach to randomly
      # This line selects a number of nodes at random to which the new node will be connected.
      # The number of nodes to attach is the minimum of the initial probability times the current number of nodes and the current number of nodes minus one.
      nodes_to_attach <- sample(V(g), size = min(round(init_p * vcount(g)), vcount(g) - 1))
      
      # Add edges to these nodes
      # This loop adds an edge from the new node to each of the selected nodes.
      for (node_to_attach in nodes_to_attach) {
        if (V(g)[node_to_attach]$label != new_label) { # No self edges
          g <- add_edges(g, c(vcount(g), node_to_attach))  # Add an edge from the new node to the selected node
        }
      }
      
      # Print new node
      # This line prints a message indicating that a new node was added and its status.
      cat("Time", t, ": New node", new_label, "added with status", new_status, "\n")
    }
  }
  
  # If there is a probability of adding an additional node
  if (fractional_part > 0){
    # Use the decimal part as a probability to add an additional node
    # This line adds an additional node with a probability equal to the decimal part of the number of nodes to add.
    if (runif(1) < fractional_part) {
      node_counter <- node_counter + 1  # Increment the node counter
      new_status <- "S"  # The status of the new node is Susceptible (S)
      new_label <- node_counter  # The label of the new node is the current value of the node counter
      
      g <- add_vertices(g, 1)  # Add a new node to the graph
      V(g)[vcount(g)]$status <- new_status  # Assign the status to the new node
      V(g)[vcount(g)]$label <- new_label  # Assign the label to the new node
      
      # Choose nodes to attach to randomly
      nodes_to_attach <- sample(V(g), size = min(round(init_p * vcount(g)), vcount(g) - 1))
      
      # Add edges to these nodes
      for (node_to_attach in nodes_to_attach) {
        if (V(g)[node_to_attach]$label != new_label) { # No self edges
          g <- add_edges(g, c(vcount(g), node_to_attach))  # Add an edge from the new node to the selected node
        }
      }
      
      # Print new node
      cat("Time", t, ": New node", new_label, "added with status", new_status, "\n")
    }
  }
  
  # Iterate over nodes for status changes
  # This loop iterates over each node in the graph to update its status based on the model parameters and its current status.
  for (i in V(g)) {
    # Get current status
    status <- V(g)[i]$status
    
    # Generate a random number
    rand <- runif(1)
    
    # Update status based on model parameters and current status
    if (status == "S") {
      # If the node is Susceptible (S) and the random number is less than beta1 * I / N, change the status to Exposed (E).
      if (rand < beta1 * I / N) {
        V(g)[i]$status <- "E"
        cat("Time", t, ": Node", V(g)[i]$label, "changed status from S to E\n")
      }
    } else if (status == "E" && rand < beta2) {
      # If the node is Exposed (E) and the random number is less than beta2, change the status to Infected (I).
      V(g)[i]$status <- "I"
      cat("Time", t, ": Node", V(g)[i]$label, "changed status from E to I\n")
    } else if (status == "I") {
      # If the node is Infected (I), generate another random number.
      rand2 <- runif(1)
      if (rand < beta3) {
        # If the first random number is less than beta3, change the status to Recovered (R).
        V(g)[i]$status <- "R"
        cat("Time", t, ": Node", V(g)[i]$label, "changed status from I to R\n")
      } else if (rand2 < alpha1) {
        # If the second random number is less than alpha1, change the status to Treatment (Tt).
        V(g)[i]$status <- "Tt"
        cat("Time", t, ": Node", V(g)[i]$label, "changed status from I to Tt\n")
      }
    } else if (status == "Tt" && rand < alpha2) {
      # If the node is under Treatment (Tt) and the random number is less than alpha2, change the status to Recovered (R).
      V(g)[i]$status <- "R"
      cat("Time", t, ": Node", V(g)[i]$label, "changed status from Tt to R\n")
    }
  }
  
  # Store counts
  # This section of the code counts the number of nodes in each status at the current time step and stores these counts in the corresponding vectors.
  S_count[t+1] <- sum(V(g)$status == "S")
  E_count[t+1] <- sum(V(g)$status == "E")
  I_count[t+1] <- sum(V(g)$status == "I")
  Tt_count[t+1] <- sum(V(g)$status == "Tt")
  R_count[t+1] <- sum(V(g)$status == "R")
  N_count[t+1] <- vcount(g)
  
  # Calculate and store metrics
  # This section of the code calculates various network metrics at the current time step and stores these metrics in the corresponding vectors.
  degree_dist[[t+1]] <- degree_distribution(g)  # Degree distribution
  clustering_coeff[t+1] <- transitivity(g, type = "global")  # Clustering coefficient
  avg_path_length[t+1] <- mean_distance(g, directed = FALSE)  # Average path length
  comps <- components(g)  # Connected components
  largest_comp_size[t+1] <- max(comps$csize)  # Size of the largest connected component
  
  # Plot network in every 10th time step
  # This loop plots the network at every 10th time step.
  # Each node is colored according to its status, and the plot includes a legend showing the number of nodes in each status.
  if (t %% 10 == 0) {
    plot(g, vertex.color = status_colors[V(g)$status], vertex.size=10, vertex.label=V(g)$label)
    title(paste("Time Step:", t))
    status_counts <- table(V(g)$status)
    status_labels <- paste(names(status_colors), " (", status_counts[names(status_colors)], ")", sep = "")
    legend("bottomright", legend = status_labels, fill = status_colors, title = "Status")
  }
}
# Set up a single plot for comparison
# This section of the code sets up a 3x2 grid of plots for comparing the number of individuals in each status (S, E, I, Tt, R, N) over time.
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
# This section of the code plots various network metrics over time.
# These metrics include the degree distribution, clustering coefficient, average path length, and size of the largest connected component.
par(mfrow = c(2, 2))
plot(degree_dist[[t+1]], main = paste("Degree Distribution at Time", t), xlab = "Degree", ylab = "Frequency")
plot(times, clustering_coeff, type = "l", xlab = "Time (Days)", ylab = "Clustering Coefficient", main = "Clustering Coefficient")
plot(times, avg_path_length, type = "l", xlab = "Time (Days)", ylab = "Average Path Length", main = "Average Path Length")
plot(times, largest_comp_size, type = "l", xlab = "Time (Days)", ylab = "Size of Largest Component", main = "Size of Largest Component")





# Create a data frame from your data
df_ode <- data.frame(time = out[, "time"], 
                     S = out[, "S"], 
                     E = out[, "E"], 
                     I = out[, "I"], 
                     Tt = out[, "Tt"], 
                     R = out[, "R"], 
                     N = out[, "N"],
                     type = "ODE")

df_net <- data.frame(time = times, 
                     S = S_count, 
                     E = E_count, 
                     I = I_count, 
                     Tt = Tt_count, 
                     R = R_count, 
                     N = N_count,
                     type = "Network")

df <- rbind(df_ode, df_net)

# Reshape the data to long format for ggplot
df_long <- tidyr::pivot_longer(df, -c(time, type), names_to = "status", values_to = "count")

# Get peak points for each status and type
peaks <- df_long %>%
  group_by(status, type) %>%
  slice(which.max(count)) %>%
  ungroup() %>%
  mutate(label = paste0("(", round(time, 1), ", ", round(count, 1), ")"))

# Define colors for network solution
network_colors <- c("S" = "gray", "E" = "orange", "I" = "red", "Tt" = "green", "R" = "purple", "N" = "blue")

# Plot
ggplot(df_long, aes(x = time, y = count)) +
  geom_line(data = df_long[df_long$type == "ODE", ], aes(color = type)) +
  geom_line(data = df_long[df_long$type == "Network", ], aes(color = status), linewidth = 1.1) +
  geom_point(data = peaks, aes(fill = type, color = NULL), size = 3, shape = 21) +
  geom_text(data = peaks, aes(label = label), vjust = -1) +
  scale_color_manual(values = c("ODE" = "black", network_colors)) +
  scale_fill_manual(values = c("ODE" = "white", "Network" = network_colors), guide = FALSE) +
  labs(x = "Time", y = "Count") +
  facet_wrap(~ factor(status, levels = c("S", "E", "I", "Tt", "R", "N")), scales = "free_y", ncol = 2) +
  theme_minimal() +
  expand_limits(y = max(df_long$count) * .4) # Expand y-axis limits

