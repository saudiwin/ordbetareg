# File with code needed to fit ordered beta regression
# As custom family in R package brms
# See vignette at 

# custom family

ord_beta_reg <- custom_family("ord_beta_reg",
                              dpars=c("mu","phi","cutzero","cutone"),
                              links=c("logit","log",NA,NA),
                              lb=c(NA,0,NA,NA),
                              type="real")

# stan code for density of the model

stan_funs <- "real ord_beta_reg_lpdf(real y, real mu, real phi, real cutzero, real cutone) {

    //auxiliary variables
    real mu_logit = logit(mu);
    vector[2] thresh;
    thresh[1] = cutzero;
    thresh[2] = cutzero + exp(cutone);
    
  if(y==0) {
      return log1m_inv_logit(mu_logit - thresh[1]);
    } else if(y==1) {
      return log_inv_logit(mu_logit  - thresh[2]);
    } else {
      return log(inv_logit(mu_logit  - thresh[1]) - inv_logit(mu_logit - thresh[2])) +
                beta_proportion_lpdf(y|mu,phi);
    }
  }"

stanvars <- stanvar(scode=stan_funs,block="functions")

# For pulling posterior predictions

posterior_predict_ord_beta_reg <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi
  cutzero <- draws$dpars$cutzero
  cutone <- draws$dpars$cutone
  N <- length(phi)
  
  thresh1 <- cutzero
  thresh2 <- cutzero + exp(cutone)
  
  pr_y0 <- 1 - plogis(qlogis(mu) - thresh1)
  pr_y1 <- plogis(qlogis(mu) - thresh2)
  pr_cont <- plogis(qlogis(mu)-thresh1) - plogis(qlogis(mu) - thresh2)
  out_beta <- rbeta(n=N,mu*phi,(1-mu)*phi)
  
  # now determine which one we get for each observation
  outcomes <- sapply(1:N, function(i) {
    sample(1:3,size=1,prob=c(pr_y0[i],pr_cont[i],pr_y1[i]))
  })
  
  final_out <- sapply(1:length(outcomes),function(i) {
    if(outcomes[i]==1) {
      return(0)
    } else if(outcomes[i]==2) {
      return(out_beta[i])
    } else {
      return(1)
    }
  })
  
  final_out
  
}

# for calculating marginal effects/conditional expectations

pp_expect_ord_beta_reg <- function(draws) {
  
  cutzero <- draws$dpars$cutzero
  cutone <- draws$dpars$cutone
  
  mu <- draws$dpars$mu
  
  thresh1 <- cutzero
  thresh2 <- cutzero + exp(cutone)
  
  low <- 1 - plogis(qlogis(mu) - thresh1)
  middle <- plogis(qlogis(mu)-thresh1) - plogis(qlogis(mu) - thresh2)
  high <- plogis(qlogis(mu) - thresh2)
  
  low*0 + middle*mu + high
}

# for calcuating LOO and Bayes Factors

log_lik_ord_beta_reg <- function(i, draws) {
  
  mu <- draws$dpars$mu[,i]
  phi <- draws$dpars$phi
  y <- draws$data$Y[i]
  cutzero <- draws$dpars$cutzero
  cutone <- draws$dpars$cutone
  
  thresh1 <- cutzero
  thresh2 <- cutzero + exp(cutone)
  
  if(y==0) {
    out <- log(1 - plogis(qlogis(mu) - thresh1))
  } else if(y==1) {
    out <- log(plogis(qlogis(mu) - thresh2))
  } else {
    out <- log(plogis(qlogis(mu)-thresh1) - plogis(qlogis(mu) - thresh2)) + dbeta(y,mu*phi,(1-mu)*phi,log=T)
  }
  
  out
  
}

# Feel free to add any other priors / change the priors on b, 
# which represent regression coefficients on the logit
# scale

priors <- set_prior("normal(0,5)",class="b") + 
  prior(constant(0),class="b",coef="Intercept") +
  prior_string("target += normal_lpdf((cutzero + exp(cutone)) - cutzero|0,3) + cutone",check=F) +
  set_prior("exponential(.1)",class="phi")

priors_phireg <- set_prior("normal(0,5)",class="b") + 
  prior(constant(0),class="b",coef="Intercept") +
  prior_string("target += normal_lpdf((cutzero + exp(cutone)) - cutzero|0,3) + cutone",check=F)