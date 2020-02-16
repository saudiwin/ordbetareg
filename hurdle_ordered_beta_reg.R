# new simulation of hurdle ordered beta regression

require(dplyr)
require(tidyr)
require(rstan)
require(loo)

hurdle_ord <- stan_model("beta_logit_infl_v4.stan")

N <- 1000

X <- rnorm(N,runif(1,-2,2),1)
X_miss <- runif(N,-1,1)

X_beta <- 1.5
X_beta_miss <- -0.5
alpha_miss <- -1

eta <- X*X_beta
eta_miss <- alpha_miss + X_miss * X_beta_miss

# ancillary parameter of beta distribution
phi <- 2

# predictor for ordered model
mu1 <- eta
# predictor for beta regression
mu2 <- eta
# probability of 0.5
mu3 <- eta_miss

cutpoints <- c(-2,2)

# probabilities for three possible categories (0, proportion, 1)
low <- 1-plogis(mu2 - cutpoints[1])
middle <- plogis(mu2-cutpoints[1]) - plogis(mu2-cutpoints[2])
high <- plogis(mu2 - cutpoints[2])
point5 <- plogis(mu3)

# we'll assume the same eta was used to generate outcomes

out_beta <- rbeta(N,plogis(mu1) * phi, (1 - plogis(mu1)) * phi) 

# now determine which one we get for each observation
outcomes <- sapply(1:N, function(i) {
  sample(1:3,size=1,prob=c(low[i],middle[i],high[i]))
})

point5_out <- as.numeric(point5>runif(N))

# now combine binary (0/1) with proportion (beta)

final_out <- sapply(1:length(outcomes),function(i) {
  if(point5_out[i]==1) {
    return(0.5)
  } else {
    if(outcomes[i]==1) {
      return(0)
    } else if(outcomes[i]==2) {
      return(out_beta[i])
    } else {
      return(1)
    }
  }
})


