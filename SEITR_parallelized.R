library(igraph)
library(deSolve)
library(foreach)
library(doParallel)

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

n <- 200
S <- 150
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

times <- seq(0, 100, by = 1)

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

status_counts_list <- foreach(run = 1:6, .packages = c("igraph", "deSolve")) %dopar% {
  g <- erdos.renyi.game(n, 0.9)
  V(g)$label <- 1:n
  V(g)$status <- "S"
  V(g)$status[sample(1:n, size = I)] <- "I"
  V(g)$status[sample(setdiff(1:n, which(V(g)$status == "I")), size = E)] <- "E"
  V(g)$status[sample(setdiff(1:n, which(V(g)$status %in% c("I", "E"))), size = Tt)] <- "Tt"
  V(g)$status[sample(setdiff(1:n, which(V(g)$status %in% c("I", "E", "Tt"))), size = R)] <- "R"
  S_count <- sum(V(g)$status == "S")
  I_count <- sum(V(g)$status == "I")
  E_count <- sum(V(g)$status == "E")
  Tt_count <- sum(V(g)$status == "Tt")
  R_count <- sum(V(g)$status == "R")
  N_count <- vcount(g)
  status_counts <- list()
  
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
