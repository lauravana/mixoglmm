set.seed(03062019)

## Number of units
N <- 500
## Number of covariates
P <- 3
## Covariates matrix
X <- cbind(Int=1, X1 = rnorm(N), X2 = rnorm(N))
## Coefficients for the binary variable
alpha <- c(0.2, -1, 0.5)
## Variance of the random effects
tau2 <- 0.2
##################################
## Simulation of random effects ##
##################################
u <- rnorm(N, 0, sqrt(tau2))

#############################################
## Simulation of Bernoulli random variable ##
#############################################
## Linear predictor of Bernoulli variable
eta_Be <- X %*% alpha + u

link <- "probit"
pfun <- switch(link,
               "logit"  = plogis,
               "probit" = pnorm,
               "cauchit" = pcauchy)

Be1 <- rbinom(N, 1, pfun(eta_Be))
##########################################
## Simulation of Poison random variable ##
##########################################
eta_Po <- eta_Be + 1 # only the intercept is different
Po1 <- rpois(N, lambda = exp(eta_Po))
###############################################
## Simulation of the normal random variables ##
###############################################
K2 <- 3 # Number of normal random variable
## mean of the normal random variables
m <- cbind(eta_Be - 1, eta_Be, eta_Be + 1) # only the intercept is different
## Correlation matrix
R <- matrix(c(1, 0.9, 0.6, 0.9, 1, 0.5, 0.6, 0.5, 1), ncol = K2)
## Standard deviations of normal variables
omega <- c(0.5, 1, 2)

Omega <- tcrossprod(omega) * R

e <-  MASS::mvrnorm(N, rep.int(0, K2), Omega)
## K2 dimensional normal random variables
y2 <- m + e


data_toy <- cbind.data.frame(y1 = y2[, 1], y2 = y2[, 2], y3 = y2[, 3],
                             Po1 = Po1,
                             Be1 = Be1,
                             X1  = X[, 2],
                             X2  = X[, 3],
                             X3 = factor(sample(c("A", "B", "C"), N, replace = TRUE)))

data_toy$Be1 <- factor(data_toy$Be1, labels = c("ND", "D"))
save(data_toy, file = "data/data_toy.RData")

