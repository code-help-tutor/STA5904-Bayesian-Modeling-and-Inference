# Load packages ------------------------------------------------------------

library(tidyverse)
library(viridis)

# Simulate a bunch of variables --------------------------------------------

N <- 1E6
set.seed(20181219)
x <- runif(N, -1, 1)
y <- runif(N, -1, 1)
chi <- ifelse(x^2 + y^2 <= 1, 1, 0)

mean_chi <- mean(chi)
var_chi <- var(chi)

4 * mean_chi + 4 * c(-1, 0, 1) * 1.96 * sqrt(var_chi / N)

