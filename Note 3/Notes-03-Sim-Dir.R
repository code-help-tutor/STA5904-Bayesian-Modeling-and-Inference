# Load packages -----------------------------------------------------------

library(tidyverse)
library(viridis)
library(MCMCpack)
library(MASS)


# Set seed ----------------------------------------------------------------

set.seed(20190131)
ALPHA <- 10

# Make partition ----------------------------------------------------------

lower <- c(-Inf, seq(from = -4, to = 4, length = 1000))
upper <- c(seq(from = -4, to = 4, length = 1000), Inf)

alpha_dir <- ALPHA * (pnorm(upper) - pnorm(lower))

g_dir <- as.numeric(rdirichlet(n = 1, alpha = alpha_dir))

plot(lower, cumsum(g_dir), type = 'l')

plot(pnorm, add = TRUE, col = viridis(10)[2], lwd = 2, xlim = range(lower[-1]))



# Stick breaking ----------------------------------------------------------

set.seed(20190131)

K      <- 100
theta  <- rnorm(K)
V      <- rbeta(K, 1, ALPHA)
V[K]   <- 1
w      <- numeric(K)
w[1]   <- V[1]
w[2:K] <- V[2:K] * cumprod(1 - V[1:(K-1)])

sum(w) ## Check this sums to 1

par(mfrow = c(1,2))
plot(theta, w, type = 'h', col = viridis(10)[2], lwd = 2, ylab = "Probability")
points(theta, w, pch = 20, col = viridis(10)[8], cex = 1)

plot(sort(w, decreasing = TRUE), xlab = "Index", ylab = "Probability")

