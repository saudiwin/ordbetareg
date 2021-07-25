# helper functions for calculating marginal effects

rbeta_mean <- function(N,mu,phi) {
  rbeta(N, mu * phi, (1 - mu) * phi)
}

# get marginal effects

eps <- 1e-7
setstep <- function(x) {
  x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
}

# functions for doing marginal effects

predict_ordbeta <- function(cutpoints=NULL,X=NULL,
                            X_miss=NULL,X_beta=NULL,X_beta_miss=NULL,
                            combined_out=T) {
  
  # we'll assume the same eta was used to generate outcomes
  eta <- X %*% as.matrix(X_beta)
  
  # probabilities for three possible categories (0, proportion, 1)
  
  
  low <- 1-plogis(eta - cutpoints[1])
  middle <- plogis(eta-cutpoints[1]) - plogis(eta-cutpoints[2])
  high <- plogis(eta - cutpoints[2])
  
  
  # check for whether combined outcome or single outcome
  
  if(combined_out) {
    
    if(!is.null(X_beta_miss)) {
      pr_miss <- X_miss %*% as.matrix(X_beta_miss)
      low*0 + middle*(0.5*pr_miss + (1-pr_miss)*plogis(eta)) + high*1
    } else {
      low*0 + middle*plogis(eta) + high*1
    }
    
  } else {
    
    out_list <- list(pr_zero=low,
         pr_proportion=middle,
         proportion_value=plogis(eta),
         pr_one=high)
    
    if(!is.null(X_beta_miss)) {
      pr_miss <- X_miss %*% as.matrix(X_beta_miss)
      out_list$pr_miss <- pr_miss
    }
    
    out_list
    
  }
  
}

predict_beta <- function(X=NULL,X_beta=NULL) {
  eta <- X %*% as.matrix(X_beta)
  
  prob <- plogis(eta)
  
  prob
}

predict_zoib <- function(coef_g=NULL,coef_a=NULL,coef_m=NULL,
                         alpha1=NULL,
                         alpha2=NULL,
                         alpha3=NULL,
                         X=NULL,
                         combined_out=T) {
  
  # we'll assume the same eta was used to generate outcomes
  psi <- plogis(as.numeric(alpha1) + X %*% coef_a)
  gamma <- plogis(as.numeric(alpha2) + X %*% coef_g)
  eta <- as.numeric(alpha3) + X %*% coef_m
  
  # probabilities for three possible categories (0, proportion, 1)
  low <- psi * (1-gamma)
  middle <- (1-psi)
  high <- psi * gamma
  
  # check for whether combined outcome or single outcome
  
  if(combined_out) {
    low*0 + middle*plogis(eta) + high*1
  } else {
    list(pr_zero=low,
         pr_proportion=middle,
         proportion_value=eta,
         pr_one=high)
  }
  
}
