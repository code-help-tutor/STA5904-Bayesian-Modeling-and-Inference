## Laplace's Example -------------------------------------------------------

num_female <- 241945
num_total  <- 493472
num_male   <- num_total - num_female

## Make posterior
alpha <- num_female + 1
beta  <- num_male + 1

## Plot posterior ----------------------------------------------------------

xlim <- qbeta(c(0,1) + c(1,-1) * 1E-4, alpha, beta)

plot(function(x) dbeta(x, alpha, beta), xlim = xlim,
     xlab = "theta", ylab = "Posterior Density")

## How likely is more-girls-than-boys under this prior? --------------------

pbeta(0.5, alpha, beta, lower.tail = FALSE)
## [1] 1.146058e-42

## Credible interval -------------------------------------------------------

qbeta(c(.005, .995), alpha, beta)
## [1] 0.4884583 0.4921244

