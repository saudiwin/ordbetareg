# Robert Kubinec
# New York University Abu Dhabi
# January 23, 2020
# let's try to simulate a new distribution for ordered beta regression

require(rstan)
require(bayesplot)
require(dplyr)

N <- 1000

X <- runif(N,-2,2)

# on the logit scale, determines threshold above which we get degenerate (0/1) outcomes
cutpoints <- c(-2,2) 

X_beta <- 2.5


sample_ordbeta <- function(cutpoints=NULL,phi=NULL,X=NULL,X_beta=NULL) {
  
  # we'll assume the same eta was used to generate outcomes
  eta <- X*X_beta
  
  # non-inflated outcome
  out_beta <- rbeta(N,plogis(eta) * phi, (1 - plogis(eta)) * phi) 
  
  # probabilities for three possible categories (0, proportion, 1)
  low <- 1-plogis(eta - cutpoints[1])
  middle <- plogis(eta-cutpoints[1]) - plogis(eta-cutpoints[2])
  high <- plogis(eta - cutpoints[2])
  
  # now determine which one we get for each observation
  
  outcomes <- sapply(1:N, function(i) {
    sample(1:3,size=1,prob=c(low[i],middle[i],high[i]))
  })
  
  # combine and sample from full distribution
  
  sapply(1:length(outcomes),function(i) {
    if(outcomes[i]==1) {
      return(0)
    } else if(outcomes[i]==2) {
      return(out_beta[i])
    } else {
      return(1)
    }
  })
  
}

predict_ordbeta <- function(cutpoints=NULL,phi=NULL,X=NULL,X_beta=NULL,
                         combined_out=T) {
  
  # we'll assume the same eta was used to generate outcomes
  eta <- X*X_beta
  
  # probabilities for three possible categories (0, proportion, 1)
  low <- 1-plogis(eta - cutpoints[1])
  middle <- plogis(eta-cutpoints[1]) - plogis(eta-cutpoints[2])
  high <- plogis(eta - cutpoints[2])
  
  # check for whether combined outcome or single outcome
  
  if(combined_out) {
    low*0 + middle*plogis(eta) + high*1
  } else {
    list(pr_zero=low,
         pr_proportion=middle,
         proportion_value=out_beta,
         pr_one=high)
  }
  
}



# use this function to sample data

final_out <- sample_ordbeta(cutpoints=cutpoints,
                             phi=phi,
                             X=X,
                             X_beta=X_beta)

# we can also use it to calculate "true" marginal effect of X on Y using code from margins package
# i.e., numerical differentiation
# set value of `h` based on `eps` to deal with machine precision

eps <- 1e-7
setstep <- function(x) {
  x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
}

y0 <- predict_ordbeta(cutpoints=cutpoints,
                      phi=phi,
                      X=X - setstep(X),
                      X_beta=X_beta)

y1 <- predict_ordbeta(cutpoints=cutpoints,
                      phi=phi,
                      X=X + setstep(X),
                      X_beta=X_beta)

marg_eff <- (y1-y0)/((X + setstep(X))-(X - setstep(X)))

mean(marg_eff)

# check for pr0 and pr1

y0 <- predict_ordbeta(cutpoints=cutpoints,
                      phi=phi,
                      X=X - setstep(X),
                      X_beta=X_beta,combined_out = F)

y1 <- predict_ordbeta(cutpoints=cutpoints,
                      phi=phi,
                      X=X + setstep(X),
                      X_beta=X_beta,combined_out = F)

marg_eff_0 <- (y1$pr_zero-y0$pr_zero)/((X + setstep(X))-(X - setstep(X)))

marg_eff_1 <- (y1$pr_one-y0$pr_one)/((X + setstep(X))-(X - setstep(X)))

mean(marg_eff_0)
mean(marg_eff_1)

# check what OLS does

ols_fit <- lm(final_out~X)

summary(ols_fit)

# now we need a Stan file

beta_logit <- stan_model("beta_logit.stan")

to_sample <- 400

indices_degen <- sample(1:length(X[final_out %in% c(0,1)]),size=to_sample/2)
indices_prop <- sample(1:length(X[final_out>0 & final_out<1]),size=to_sample/2)

to_bl <- list(N_degen=sum(final_out %in% c(0,1)),
              N_prop=sum(final_out>0 & final_out<1),
              X=1,
              outcome_prop=final_out[final_out>0 & final_out<1],
              outcome_degen=final_out[final_out %in% c(0,1)],
              covar_prop=as.matrix(X[final_out>0 & final_out<1]),
              covar_degen=as.matrix(X[final_out %in% c(0,1)]),
              N_pred_degen=sum(final_out %in% c(0,1)),
              N_pred_prop=sum(final_out>0 & final_out<1),
              indices_degen=1:(sum(final_out %in% c(0,1))),
              indices_prop=1:(sum(final_out>0 & final_out<1)),
              run_gen=1)

fit_model <- sampling(beta_logit,data=to_bl,chains=2,cores=2,iter=1000)

# regenerate data and see how we do capturing it

yrep_ordbeta <- extract(fit_model,"regen_all")[[1]]

# use only a sample of draws

# final_out_prop <- final_out[final_out>0 & final_out<1][indices_prop]
# final_out_degen <- final_out[final_out %in% c(0,1)][indices_degen]

final_out_prop <- final_out[final_out>0 & final_out<1]
final_out_degen <- final_out[final_out %in% c(0,1)]

ppc_ecdf_overlay(y=c(final_out_degen,final_out_prop),yrep=yrep_ordbeta) + 
  ggtitle("Empirical Posterior Predictive Distribution for Ordinal Beta Regression",subtitle="N=1000")

# ppc_dens_overlay(y=c(final_out_degen,final_out_prop),yrep=yrep) + 
#   ggtitle("Posterior Predictive Distribution for Ordinal Beta Regression",subtitle="N=1000")

# looks pretty good

# try with rstanarm betareg where we pre-transform outcome

require(rstanarm)

final_out_scale <- (final_out * (length(final_out) - 1) + .5) / length(final_out)

betareg_fit <- stan_betareg(formula = outcome~X,data=tibble(outcome=final_out_scale,
                                                            X=X),chains=2,cores=2,iter=1000)

summary(betareg_fit)

# X_beta has lost 1/3 its value from original (i.e., effect has been compressed)

# see the distortion re: the original values

yrep_betareg <- posterior_predict(betareg_fit,draws=100)

ppc_dens_overlay(y=final_out,yrep=yrep_betareg) + ggtitle("Posterior Predictive Distribution for Rstanarm Beta Regression",subtitle="N=1000, Transformed to (0,1) Scale")

ppc_dens_overlay(y=final_out_scale,yrep=yrep_betareg) + ggtitle("Posterior Predictive Distribution for Rstanarm Beta Regression",subtitle="N=1000, Transformed to (0,1) Scale")

# Tends to over/under predict middle and outliers (i.e. 0/1)

# now try the ZOIB

zoib_model <- stan_model("zoib.stan")
x <- as.matrix(X)
zoib_fit <- sampling(zoib_model,data=list(n=length(final_out),
                                          y=final_out,
                                          k=ncol(x),
                                          x=x,
                                          run_gen=1),chains=2,cores=2,iter=1000)

# this is X_beta
print(zoib_fit,pars="coef_m")
yrep_zoib <- extract(zoib_fit,"zoib_regen")[[1]]
ppc_dens_overlay(y=final_out,yrep=yrep_zoib) + 
  ggtitle("Empirical Posterior Predictive Distribution for Ordinal Beta Regression",subtitle="N=1000")

# Compare Distributions ---------------------------------------------------

require(loo)

loo_ordbeta <- loo(fit_model,"ord_log")

loo_zoib <- loo(zoib_fit,"zoib_log")

loo_compare(loo_ordbeta,loo_zoib)


# Marginal Effects --------------------------------------------------------

# we can calculate predicted values and marginal effects for a single X

X_beta <- as.matrix(fit_model,"X_beta")
cuts_est <- as.matrix(fit_model,"cutpoints")

# fascinating. we can do numerical differentiation -- cool

pr_y0 <- sapply(1:nrow(X_beta),function(i) {
  b <- X_beta[i,]
  c <- cuts_est[i,1]
  mean(1-plogis(b*X - c)) 
           })

pr_y1 <- sapply(1:nrow(X_beta),function(i) {
  b <- X_beta[i,]
  c <- cuts_est[i,2]
  mean(plogis(b*X - c)) 
})

mu_est <- sapply(1:nrow(X_beta),function(i) {
  b <- X_beta[i,]
  plogis(mean(b*X)) 
})

mean(plogis(mu1))

