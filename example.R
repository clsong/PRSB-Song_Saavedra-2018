# load necessary R packages -----------------------------------------------
library(tidyverse)
library(mvtnorm)

# Functions  --------------------------------------------------------------

# compute the interaction strength matrices according to the parameterization given in the methods section
# inputs: web = binary matrix of plant-pollinator interaction, gamma_avg = average mutualistic strength,
# rho = mean field interspecific competition, delta = mutualistic trade-off
# outputs: alpha = interaction strength matrix (full), alphaA = competition among plants, betaP = competition among animals,
# gammaA = mutualistic effect of plants on animals, gammaP = mutualistic effect of animals on plants
construct_interaction_matrix <- function(web, gamma_avg, rho, delta) {
  SA <- nrow(web)
  SP <- ncol(web)
  alphaA <- matrix(rho, SA, SA) + (1 - rho) * diag(rep(1, SA))
  alphaP <- matrix(rho, SP, SP) + (1 - rho) * diag(rep(1, SP))
  gammaA <- diag(rowSums(web) ^ -delta) %*% web
  gammaP <- diag(colSums(web) ^ -delta) %*% t(web)
  f <- sum(gammaA[web == 1] + gammaP[t(web) == 1]) / (2 * sum(web == 1))
  gammaA <- gamma_avg / f * diag(rowSums(web) ^ -delta) %*% web
  gammaP <- gamma_avg / f * diag(colSums(web) ^ -delta) %*% t(web)
  alpha <- rbind(cbind(alphaA, -gammaA), cbind(-gammaP, alphaP))
  
  alpha
}

# stability condition gamma_hat (average mutualistic strength at the stability threshold)
# inputs: web = mutualistic network (binary matrix),
# rho = mean field interspecific competition, delta = mutualistic trade-off
# output: gamma_hat = stability condition
gamma_hat <- function(web, rho, delta) {
  f_eig <- function(gamma_avg, web, rho, delta) {
    alpha <- interaction_matrix(web, gamma_avg, rho, delta)$alpha
    out <- (min(Re(eigen(alpha)$values))) ^ 2
    out
  }
  
  optimize(f_eig, c(0, 1000), web = web, rho = rho, delta = delta)$minimum
}


# compute structural stability as the relative size of the feasibility domain
# input: alpha = interaction strength matrix
# output: Omega = structural stability
Omega <- function(alpha) {
  S <- nrow(alpha)
  Sigma <- solve(t(alpha) %*% alpha)
  m <- matrix(0, S, 1)
  a <- matrix(0, S, 1)
  b <- matrix(Inf, S, 1)
  d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
  
  d[1]^(1/S)
}


# Example -----------------------------------------------------------------
# Load plant-pollinator binary interaction data
web <- read.delim("1996_12.txt") %>%
  as.tibble() %>%
  select(-Day_23) %>% 
  as.matrix()

#generate the interaction matrix
alpha <- construct_interaction_matrix(web, gamma_avg = 0.1, rho = 0, delta = 0.5)

#compute the structural stability
omega <- Omega(alpha) 
