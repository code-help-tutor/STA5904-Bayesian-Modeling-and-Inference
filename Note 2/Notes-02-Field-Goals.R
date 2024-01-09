# Load packages -----------------------------------------------------------

library(tidyverse)
library(viridis)
library(MASS)

# Load data ---------------------------------------------------------------

field_goals <- read.csv("Field-Goals.csv")

X <- cbind(1, field_goals$Distance)
y <- field_goals$Success

# Make likelihood function to optimize ------------------------------------

expit <- function(x) 1 / (1 + exp(-x))

neg_log_pi <- function(beta) {
  eta <- X %*% beta
  p <- expit(eta)
  logliks <- y * log(p) + (1 - y) * log(1 - p)
  logliks <- ifelse(is.nan(logliks), 0, logliks)
  loglik <-  sum(logliks)
  loglik <- loglik + sum(dcauchy(beta, scale = 5, log = TRUE))
  return(-loglik)
}

opt_results <- optim(par = c(0,0), fn = neg_log_pi, hessian = TRUE)
theta_hat <- opt_results$par
theta_cov <- solve(opt_results$hessian)

par(mfrow = c(1,2))
plot(function(x) dnorm(x, theta_hat[1], sqrt(theta_cov[1,1])), 
     xlim = c(3.5,8), col = viridis(10)[8], lwd = 4)
plot(function(x) dnorm(x, theta_hat[2], sqrt(theta_cov[2,2])), 
     xlim = c(-.18,-0.05), col = viridis(10)[2], lwd = 4)


# Sampling curves ---------------------------------------------------------

theta_out <- mvrnorm(n = 10000, mu = theta_hat, Sigma = theta_cov)
x_grid <- cbind(1, seq(from = 15, to = 65, length = 1000))

sampled_eta <- theta_out %*% t(x_grid)
sampled_p <- expit(sampled_eta)
p_hat <- colMeans(sampled_p)
quantile_p <- apply(sampled_p, 2, function(x) quantile(x, c(.025, .975)))

par(mfrow = c(1,1))
plot(x_grid[,2], p_hat, type = 'l', ylim = c(0, 1))
lines(x_grid[,2], quantile_p[1,], col = viridis(10)[8])
lines(x_grid[,2], quantile_p[2,], col = viridis(10)[8])

# Using Monte Carlo --------------------------------------------------------

set.seed(12310)
N_samp <- 1E5
beta_samp <- cbind(rcauchy(n = N_samp, scale = 5), rcauchy(n = N_samp, scale = 5))
logliks <- -apply(beta_samp, 1, neg_log_pi)
logsumexp <- function(x) max(x) + log(sum(exp(x - max(x))))
weights <- exp(logliks - logsumexp(logliks))

colSums(weights * beta_samp)
plot(weights)

# Using Importance sampling ------------------------------------------------

library(mvtnorm)

log_weights <- -apply(theta_out, 1, neg_log_pi) - dmvnorm(theta_out, theta_hat, theta_cov, log = TRUE)
weights <- exp(log_weights - logsumexp(log_weights))

plot(weights)
colSums(weights * theta_out)

# Using RWMH ---------------------------------------------------------------

set.seed(12311)
N_rwmh      <- 1E4
beta_0_rwmh <- c(0,0) 
beta_rwmh_samps <- matrix(NA, nrow= N_rwmh, ncol = 2)
for(k in 1:2) {
  for(i in 1:N_rwmh) {
    beta_prop_rwmh <- beta_0_rwmh + runif(2, -1, 1) * c(.4, .01)
    loglik_cur <- -neg_log_pi(beta_0_rwmh)
    loglik_prop <- -neg_log_pi(beta_prop_rwmh)
    if(log(runif(1)) < loglik_prop - loglik_cur) {
      beta_0_rwmh <- beta_prop_rwmh
    }
    beta_rwmh_samps[i,] <- beta_0_rwmh
  }
}

par(mfrow = c(2,2))
plot(beta_rwmh_samps[,1], type = 'l')
plot(beta_rwmh_samps[,2], type = 'l')
plot(density(beta_rwmh_samps[,1]))
plot(density(beta_rwmh_samps[,2]))

mean(beta_rwmh_samps[-1,1] == beta_rwmh_samps[-N_rwmh,1])
colMeans(beta_rwmh_samps)

library(coda)
effectiveSize(beta_rwmh_samps)

# Using STAN ---------------------------------------------------------------

library(rstan)

model_str <-
'
data {
  int<lower = 0> N;
  int<lower = 0> P;
  matrix[N,P] X;
  int<lower = 0, upper = 1> Y[N];
}

parameters {
  vector[P] beta;
}

model {
  beta ~ cauchy(0,5);
  Y ~ bernoulli_logit(X * beta);
}
'

logistic_model <- stan_model(model_code = model_str)
my_data        <- list(Y = y, X = X, N = nrow(X), P = ncol(X))
hmc_samples    <- sampling(logistic_model, data = my_data)

rstan::stan_dens(hmc_samples)
rstan::stan_ac(hmc_samples)
traceplot(hmc_samples)
summary(hmc_samples)[["summary"]]

# Using rjags --------------------------------------------------------------

library(rjags)
library(runjags)
load.module("glm")

jags_code <-
'
model {
  for(p in 1:P) {
    beta[p] ~ dt(0.0, 0.04, 1)
  }
  eta <- X %*% beta
  for(i in 1:N) {
    p[i] <- ilogit(eta[i])
    Y[i] ~ dbern(p[i])
  }
}
'

logistic_jags <- run.jags(model = jags_code, data = my_data, monitor = c("beta"))
plot(logistic_jags)


## Bad mixing for RWMH -----------------------------------------------------


set.seed(12311)
N_rwmh      <- 1E4
beta_0_rwmh <- c(0,0) 
beta_rwmh_samps_2 <- matrix(NA, nrow= N_rwmh, ncol = 2)
beta_rwmh_burn_2 <- matrix(NA, nrow= N_rwmh, ncol = 2)
badwidth <- 0.1
for(i in 1:N_rwmh) {
  beta_prop_rwmh <- beta_0_rwmh + runif(2, -1, 1) * c(badwidth, .01)
  loglik_cur <- -neg_log_pi(beta_0_rwmh)
  loglik_prop <- -neg_log_pi(beta_prop_rwmh)
  if(log(runif(1)) < loglik_prop - loglik_cur) {
    beta_0_rwmh <- beta_prop_rwmh
  }
  beta_rwmh_burn_2[i,] <- beta_0_rwmh
}

for(i in 1:N_rwmh) {
  beta_prop_rwmh <- beta_0_rwmh + runif(2, -1, 1) * c(badwidth, .01)
  loglik_cur <- -neg_log_pi(beta_0_rwmh)
  loglik_prop <- -neg_log_pi(beta_prop_rwmh)
  if(log(runif(1)) < loglik_prop - loglik_cur) {
    beta_0_rwmh <- beta_prop_rwmh
  }
  beta_rwmh_samps_2[i,] <- beta_0_rwmh
}

par(mfrow = c(2,2))
plot(beta_rwmh_burn_2[,1], type = 'l')
plot(beta_rwmh_samps_2[,1], type = 'l')
plot(beta_rwmh_burn_2[,2], type = 'l')
plot(beta_rwmh_samps_2[,2], type = 'l')

mean(beta_rwmh_samps[-1,1] == beta_rwmh_samps[-N_rwmh,1])
colMeans(beta_rwmh_samps)

library(coda)
effectiveSize(beta_rwmh_samps)
