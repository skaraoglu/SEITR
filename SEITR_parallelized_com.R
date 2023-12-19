library(igraph)
library(deSolve)
library(foreach)
library(doParallel)
library(Matrix)

Lambda <- 2
beta1 <- 0.8
beta2 <- 0.3
beta3 <- 0.02
alpha1 <- 0.1
alpha2 <- 0.055
delta_I <- 0.03
delta_T <- 0.03
mu <- 0.01

R0 <- (beta1 * beta2) / ((beta3 + mu + delta_I + alpha1) * (beta2 + mu) * (mu + delta_T + alpha2))

n <- 400
S <- 350
I <- 20
E <- 30
R <- 0
Tt <- 0
N <- S + E + I + R + Tt

state <- c(S = S, E = E, I = I, Tt = Tt, R = R, N = n)
parameters <- c(Lambda = Lambda, beta1 = beta1, beta2 = beta2, beta3 = beta3, alpha1 = alpha1, alpha2 = alpha2, delta_I = delta_I, delta_T = delta_T, mu = mu)

SEITRN <- function(t, state, parameters) {
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

times <- seq(0, 62, by = 1)

out <- ode(y = state, times = times, func = SEITRN, parms = parameters)

par(mfrow = c(3, 2))
plot(out[, "time"], out[, "S"], type = "l", col = 2, xlab = "Time (Days)", ylab = "Susceptibles", main = "Susceptibles")
plot(out[, "time"], out[, "E"], type = "l", col = 3, xlab = "Time (Days)", ylab = "Exposed", main = "Exposed")
plot(out[, "time"], out[, "I"], type = "l", col = 4, xlab = "Time (Days)", ylab = "Infected", main = "Infected")
plot(out[, "time"], out[, "Tt"], type = "l", col = 5, xlab = "Time (Days)", ylab = "Treatment", main = "Treatment")
plot(out[, "time"], out[, "R"], type = "l", col = 6, xlab = "Time (Days)", ylab = "Recovered", main = "Recovered")
plot(out[, "time"], out[, "N"], type = "l", col = 7, xlab = "Time (Days)", ylab = "Total Population", main = "Total Population")

status_colors <- c("S" = "gray", "E" = "yellow", "I" = "red", "Tt" = "green", "R" = "purple")
node_counter <- vcount(g)
status_counts_list <- list()
num_cores <- 6  # Adjust the number of cores based on your machine's capabilities
cl <- makeCluster(num_cores)
registerDoParallel(cl)

status_counts_list <- foreach(run = 1:6, .packages = c("igraph", "deSolve", "Matrix")) %dopar% {
  print(run)
  npc <-80 # nodes per cluster
  n_clust <- 5 # 4 clusters with 25 nodes each
  matlist = list()
  cluster_vector <- c()
  
  for (i in 1:n_clust){ 
    matlist[[i]] = get.adjacency(erdos.renyi.game(npc, 0.8))
    cluster_vector <- c(cluster_vector, rep(i, npc))
  }
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
  comp_ids <- seq(1,comps$no-1) # stop short of last one
  for (curr_comp in comp_ids){
    i <- sample(which(curr_comp==member_vec),1)
    j <- sample(which((curr_comp+1)==member_vec),1)
    num_edges <- sample(1:4, 1)
    
    for (edge in 1:num_edges) {
      if (!are_adjacent(g, i, j)) {
        mat_clust[i,j] <- 1
        mat_clust[j,i] <- 1
      }
      i <- sample(which(curr_comp==member_vec),1)
      j <- sample(which((curr_comp+1)==member_vec),1)
    }
  }
  g <- graph_from_adjacency_matrix(mat_clust, mode="undirected", diag=F)
  V(g)$cluster <- cluster_vector
  V(g)$label <- 1:n
  V(g)$status <- "S"
  init_cluster <- sample(1:n_clust, 1)
  init_cluster_nodes <- V(g)[member_vec == init_cluster]
  init_I_nodes <- sample(init_cluster_nodes, I)
  init_E_nodes <- sample(setdiff(init_cluster_nodes, init_I_nodes), E)
  V(g)[init_I_nodes]$status <- "I"
  V(g)[init_E_nodes]$status <- "E"
  S_count <- S
  E_count <- E
  I_count <- I
  Tt_count <- Tt
  R_count <- R
  N_count <- N
  status_counts <- list()
  plot(g, vertex.color = status_colors[V(g)$status], vertex.size=4, vertex.label="")
  
  for (t in times) {
    old_status <- V(g)$status
    
    nodes_to_remove_count <- delta_I * sum(V(g)$status == "I")
    nodes_with_status_I <- which(V(g)$status == "I")
    if (nodes_to_remove_count < 1) {
      random_number <- runif(1)
      if (random_number < nodes_to_remove_count) {
        node_to_remove <- sample(nodes_with_status_I, 1)
        g <- delete_vertices(g, node_to_remove)
      }
    } else {
      for (i in 1:nodes_to_remove_count){
        node_to_remove <- sample(nodes_with_status_I, 1)
        g <- delete_vertices(g, node_to_remove)
      }
    }
    for (i in 1:(mu * N)) {
      node_to_remove <- sample(vcount(g), 1)
      g <- delete_vertices(g, node_to_remove)
    }
    for (i in 1:(Lambda)) {
      node_counter <- node_counter + 1
      new_status <- "S"
      new_label <- node_counter
      
      g <- add_vertices(g, 1)
      V(g)[vcount(g)]$status <- new_status
      V(g)[vcount(g)]$label <- new_label
      degree <- degree(g, mode = "all")
      avg_degree <- mean(degree)
      prob <- degree / sum(degree)
      nodes_to_attach <- sample(V(g), size = min(round(avg_degree), vcount(g)), prob = prob)
      for (node_to_attach in nodes_to_attach) {
        if (V(g)[node_to_attach]$label != new_label) {
          g <- add_edges(g, c(vcount(g), node_to_attach))
        }
      }
    }
    for (i in V(g)) {
      status <- V(g)[i]$status
      rand <- runif(1)
      if (status == "Tt" && rand < alpha2) {
        V(g)[i]$status <- "R"
      } else if (status == "I") {
        if (rand < beta3) {
          V(g)[i]$status <- "R"
        } else if (rand < alpha1) {
          V(g)[i]$status <- "Tt"
        }
      } else if (status == "E" && rand < beta2) {
        V(g)[i]$status <- "I"
      } else if (status == "S") {
        neighbor_statuses <- V(g)[neighbors(g, i)]$status
        if (any(neighbor_statuses %in% c("E", "I")) && rand < beta1 * I / N) {
          V(g)[i]$status <- "E"
        }
      }
    }
    S_count <- sum(V(g)$status == "S")
    I_count <- sum(V(g)$status == "I")
    E_count <- sum(V(g)$status == "E")
    Tt_count <- sum(V(g)$status == "Tt")
    R_count <- sum(V(g)$status == "R")
    N_count <- vcount(g)
    status_counts[[as.character(t)]] <- c(S_count, I_count, E_count, Tt_count, R_count, N_count)
  }
  status_counts
}

stopCluster(cl)

par(mfrow = c(3, 2))
statuses <- c("S", "E", "I", "Tt", "R", "N")
colors <- c(2, 3, 4, 5, 6, 7)
colors_alpha <- sapply(colors, function(col) rgb(red = col2rgb(col)[1,]/255, green = col2rgb(col)[2,]/255, blue = col2rgb(col)[3,]/255, alpha = 0.5))
for (j in 1:length(statuses)) {
  plot(out[, "time"], out[, statuses[j]], type = "l", col = colors[j], xlab = "Time (Days)", ylab = statuses[j], main = paste("Comparison of", statuses[j]))
  for (i in 1:length(status_counts_list)) {
    status_counts <- status_counts_list[[i]]
    counts <- numeric(length = length(times))
    for (t in times) {
      counts[as.numeric(t)] <- status_counts[[as.character(t)]][j]
    }
    lines(times, counts, col = colors_alpha[j])
  }
}
