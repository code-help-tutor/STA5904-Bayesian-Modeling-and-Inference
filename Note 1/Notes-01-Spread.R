## Load data ---------------------------------------------------------------

spread <- read.table("spread.tab", header = TRUE)
delta  <- spread$favorite - spread$underdog - spread$spread

## Prior -------------------------------------------------------------------

set.seed(1234)

alpha <- 0.5
beta  <- 0.5 * 14^2
m     <- 0
kappa <- 2

sigma_prior <- 1/sqrt(rgamma(1000, alpha, beta))
mu_prior    <- rnorm(1000, m, sigma_prior/sqrt(kappa))

## Update prior ------------------------------------------------------------

N <- length(delta)

alpha_hat <- alpha + N / 2
beta_hat  <- beta + 0.5 * (N - 1) * var(delta) +
  N * kappa * (mean(delta) - m)^2 / (2 * (N + kappa))
kappa_hat <- kappa + N
m_hat <- (sum(delta) + kappa * m) / kappa_hat

## Plot posterior of sigma -------------------------------------------------

set.seed(1235)
library(tidyverse)
library(viridis)

sigma_post <- 1/sqrt(rgamma(1E5, alpha_hat, beta_hat))
mu_post    <- rnorm(1E5, m_hat, sqrt(sigma_post^2 / kappa_hat))

par(mfrow = c(2,2))
hist(sigma_post, freq = FALSE, col = viridis(10)[8], xlab = 'sigma')
hist(mu_post, freq = FALSE, col = viridis(10)[2], xlab = 'mu')

## What happens with flat priors? ------------------------------------------

sigma_flat <- 1/sqrt(rgamma(1E5, (N - 1) / 2, (N - 1) * var(delta) / 2))
mu_flat    <- rnorm(1E5, mean(delta), sigma_flat / sqrt(N))

hist(sigma_flat, freq = FALSE, col = viridis(10)[8], xlab = 'sigma')
hist(mu_flat, freq = FALSE, col = viridis(10)[2], xlab = 'mu')
