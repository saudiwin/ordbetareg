# Robert Kubinec
# New York University Abu Dhabi
# January 23, 2020
# Simulation of 0 - 1 bounded dependent variables
# Note simulation will take some time, approx ~2 days with 3 cores

require(rstan)
require(bayesplot)
require(dplyr)
require(rstanarm)
require(loo)

set.seed(772235)

rstan_options("auto_write" = TRUE)

beta_logit <- stan_model("beta_logit.stan")
zoib_model <- stan_model("zoib_nophireg.stan")

# let's do some simulations

N_rep <- 10000

simul_data <- tibble(N=round(runif(N_rep,100,2000),0),
                     X_beta=runif(N_rep,-2,2),
                     phi=runif(N_rep,0.5,4),
                     cutpoints1=runif(N_rep,-5,-1)) %>% 
              mutate(cutpoints2=cutpoints1+runif(N_rep,0.5,5))

# we can also use it to calculate "true" marginal effect of X on Y using code from margins package
# i.e., numerical differentiation
# set value of `h` based on `eps` to deal with machine precision

eps <- 1e-7
setstep <- function(x) {
  x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
}

# functions for doing marginal effects

predict_ordbeta <- function(cutpoints=NULL,X=NULL,X_beta=NULL,
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
         proportion_value=plogis(eta),
         pr_one=high)
  }
  
}

predict_zoib <- function(coef_g=NULL,coef_a=NULL,coef_m=NULL,
                         alpha1=NULL,
                         alpha2=NULL,
                         alpha3=NULL,
                         X=NULL,
                         combined_out=T) {
  
  # we'll assume the same eta was used to generate outcomes
  psi <- plogis(alpha1 + X %*% as.matrix(coef_a))
  gamma <- plogis(alpha2 + X %*% as.matrix(coef_g))
  eta <- alpha3 + X %*% as.matrix(coef_m)
  
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

r_seeds <- c(6635,2216,8845,9936,3321)

all_simul_data <- parallel::mclapply(1:nrow(simul_data), function(i,simul_data=NULL,r_seeds=NULL) {
  this_data <- slice(simul_data,i)
  cat(file = "simul_status.txt",paste0("Now on row ",i),append = T)
  
  N <- this_data$N
  
  X <- rnorm(N,runif(1,-2,2),1)
  
  X_beta <- this_data$X_beta
  eta <- X*X_beta
  
  # ancillary parameter of beta distribution
  phi <- this_data$phi
  
  # predictor for ordered model
  mu1 <- eta
  # predictor for beta regression
  mu2 <- eta
  
  cutpoints <- c(this_data$cutpoints1,this_data$cutpoints2)
  
  # probabilities for three possible categories (0, proportion, 1)
  low <- 1-plogis(mu2 - cutpoints[1])
  middle <- plogis(mu2-cutpoints[1]) - plogis(mu2-cutpoints[2])
  high <- plogis(mu2 - cutpoints[2])
  
  # we'll assume the same eta was used to generate outcomes
  
  out_beta <- rbeta(N,plogis(mu1) * phi, (1 - plogis(mu1)) * phi) 
  
  print(i)
  
  # now determine which one we get for each observation
  outcomes <- sapply(1:N, function(i) {
    sample(1:3,size=1,prob=c(low[i],middle[i],high[i]))
  })
  
  # now combine binary (0/1) with proportion (beta)
  
  final_out <- sapply(1:length(outcomes),function(i) {
    if(outcomes[i]==1) {
      return(0)
    } else if(outcomes[i]==2) {
      return(out_beta[i])
    } else {
      return(1)
    }
  })
  
  counter <- environment()
  counter$i <- 0
  
  while(((sum(final_out>0 & final_out<1)<5) || (sum(final_out==1)<5) || (sum(final_out==0)<5)) && (counter$i<20)) {
    
    final_out <- sapply(1:length(outcomes),function(i) {
      if(outcomes[i]==1) {
        return(0)
      } else if(outcomes[i]==2) {
        return(out_beta[i])
      } else {
        return(1)
      }
    })
    cat(file = "simul_status.txt","Resampling",append=T)
    counter$i <- counter$i + 1
  }
  
  if((sum(final_out>0 & final_out<1)<5) || (sum(final_out==1)<5) || (sum(final_out==0)<5)) {
    return(this_data)
  }
  
  # calculate "true" marginal effects
  
  X_low <- X - setstep(X)
  X_high <- X + setstep(X)
  
  y0 <- predict_ordbeta(cutpoints=cutpoints,
                        X=X_low,
                        X_beta=X_beta)
  
  y1 <- predict_ordbeta(cutpoints=cutpoints,
                        X=X_high,
                        X_beta=X_beta)
  
  marg_eff <- mean((y1-y0)/((X + setstep(X))-(X - setstep(X))))
  
  # now fit ordinal beta
  
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
  
  fit_model <- sampling(beta_logit,data=to_bl,seed=r_seeds[1],
                        chains=1,cores=1,iter=1000,pars=c("regen_all","X_beta","ord_log","cutpoints"))
  
  
  x <- as.matrix(X)
  zoib_fit <- sampling(zoib_model,data=list(n=length(final_out),
                                            y=final_out,
                                            k=ncol(x),
                                            seed=r_seeds[2],
                                            x=x,
                                            run_gen=1),chains=1,cores=1,iter=1000,pars=c("coef_m","zoib_regen","zoib_log","coef_a","coef_m","coef_g","alpha"))
  
  final_out_scale <- (final_out * (length(final_out) - 1) + .5) / length(final_out)
  
  betareg_fit <- stan_betareg(formula = outcome~X,data=tibble(outcome=final_out_scale,
                                                              X=X),chains=1,cores=1,iter=1000,
                              seed=r_seeds[3])
  yrep_betareg <- posterior_predict(betareg_fit,draws=100)
  # do a second one with just non0/non1 data
  
  betareg_fit2 <- stan_betareg(formula = outcome~X,data=tibble(outcome=final_out[final_out>0 & final_out<1],
                                                               X=X[final_out>0 & final_out<1]),
                               chains=1,seed=r_seeds[4],
                               cores=1,iter=1000)
  
  yrep_betareg2 <- posterior_predict(betareg_fit2,draws=100)
  
  # fit OLS
  
  lm_fit <- stan_lm(formula = outcome~X,data=tibble(outcome=final_out,
                                                    r_seeds[5],
                                                         X=X),chains=1,cores=1,iter=1000,prior=NULL)
  
  yrep_lm <- posterior_predict(lm_fit,draws=100)
  
  # now return the full data frame
  
  X_beta_ord <- as.matrix(fit_model,"X_beta")
  X_beta_zoib <- as.matrix(zoib_fit,"coef_m")
  X_beta_reg <- as.matrix(betareg_fit,pars="X")
  X_beta_reg2 <- as.matrix(betareg_fit2,pars="X")
  X_beta_lm <- as.matrix(lm_fit,pars="X")
  
  # calculate rmse
  
  rmse_ord <- sqrt(mean(apply(as.matrix(fit_model,"regen_all"),1,function(c) { (c - c(final_out[final_out %in% c(0,1)],final_out[final_out>0 & final_out<1]))^2 })))
  rmse_zoib <-sqrt( mean(apply(as.matrix(zoib_fit,"zoib_regen"),1,function(c) { (c - final_out)^2 })))
  rmse_betareg <- sqrt(mean(apply(yrep_betareg,1,function(c) { (c - final_out)^2 })))
  rmse_betareg2 <- sqrt(mean(apply(yrep_betareg2,1,function(c) { (c - final_out[final_out>0 & final_out<1])^2 })))
  rmse_lm <- sqrt(mean(apply(yrep_lm,1,function(c) { (c - final_out)^2 })))
  
  # calculate loo
  
  loo_ordbeta <-try(loo(fit_model,"ord_log"))
  loo_betareg <- loo(betareg_fit)
  loo_betareg2 <- loo(betareg_fit2)
  loo_zoib <- try(loo(zoib_fit,"zoib_log"))
  loo_lm <- loo(lm_fit)
  
  if(any('try-error' %in% c(class(loo_ordbeta),
                            class(loo_zoib)))) {
    comp_loo <- matrix(rep(NA,6),ncol=2)
    row.names(comp_loo) <- c("model3","model2","model1")
    loo_ordbeta <- list(estimates=matrix(c(NA,NA),ncol=2))
    loo_zoib <- list(estimates=matrix(c(NA,NA),ncol=2))
  } else {
    comp_loo <- loo_compare(loo_ordbeta,loo_zoib,loo_lm)
  }
  
  
  
  # calculate marginal effects
  
  cutpoints_est <- as.matrix(fit_model,"cutpoints")
  
  margin_ord <- sapply(1:nrow(X_beta_ord), function(i) {
    y0 <- predict_ordbeta(cutpoints=cutpoints_est[i,],
                        X=X_low,
                        X_beta=X_beta_ord[i,])
  
    y1 <- predict_ordbeta(cutpoints=cutpoints_est[i,],
                            X=X_high,
                            X_beta=X_beta_ord[i,])
    
    marg_eff <- (y1-y0)/(X_high-X_low)
    
    mean(marg_eff)
  })
  
  # now for the ZOIB
  
  coef_a <- as.matrix(zoib_fit,"coef_a")
  coef_g <- as.matrix(zoib_fit,"coef_g")
  alpha <- as.matrix(zoib_fit,"alpha")
  
  margin_zoib <- sapply(1:nrow(X_beta_zoib), function(i) {
    y0 <- predict_zoib(coef_g=coef_g[i],
                       coef_a=coef_a[i],
                       alpha1=alpha[i,1],
                       alpha2=alpha[i,2],
                       alpha3=alpha[i,3],
                       X=X_low,
                       coef_m=X_beta_zoib[i,])
    
    y1 <-  predict_zoib(coef_g=coef_g[i],
                        coef_a=coef_a[i],
                        alpha1=alpha[i,1],
                        alpha2=alpha[i,2],
                        alpha3=alpha[i,3],
                        X=X_high,
                        coef_m=X_beta_zoib[i,])
    
    marg_eff <- (y1-y0)/(X_high-X_low)
    
    mean(marg_eff)
  })
  
  # now betareg
  
  betareg_int <- as.matrix(betareg_fit,pars="(Intercept)")
  
  margin_betareg <- sapply(1:nrow(X_beta_reg), function(i) {
    y0 <- plogis(betareg_int[i,] + X_low*X_beta_reg[i,])
    y1 <- plogis(betareg_int[i,] + X_high*X_beta_reg[i,])
    
    marg_eff <- (y1-y0)/(X_high-X_low)
    
    mean(marg_eff)
  })
  
  betareg2_int <- as.matrix(betareg_fit2,pars="(Intercept)")
  
  margin_betareg2 <- sapply(1:nrow(X_beta_reg2), function(i) {
    y0 <- plogis(betareg2_int[i,] + X_low*X_beta_reg2[i,])
    y1 <- plogis(betareg2_int[i,] + X_high*X_beta_reg2[i,])
    
    marg_eff <- (y1-y0)/(X_high-X_low)
    
    mean(marg_eff)
  })
  
  this_data$marg_eff <- marg_eff
  this_data$true_kurt <- moments::kurtosis(final_out)
  
  bind_cols(purrr::map_dfr(seq_len(5), ~this_data),bind_rows(tibble(model="Ordinal Beta Regression",
                   med_est=mean(X_beta_ord),
                   high=quantile(X_beta_ord,.95),
                   low=quantile(X_beta_ord,.05),
                   sd=sd(X_beta_ord),
                   loo_val=loo_ordbeta$estimates[1,1],
                   win_loo=which(row.names(comp_loo)=="model1"),
                   win_loo_se=comp_loo[which(row.names(comp_loo)=="model1"),2],
                   rmse=rmse_ord,
                   kurt_est=mean(apply(as.matrix(fit_model,"regen_all"),1,moments::kurtosis)),
                   marg_eff_est=mean(margin_ord),
                   high_marg=quantile(margin_ord,.95),
                   low_marg=quantile(margin_ord,.05),
                   sd_marg=sd(margin_ord)),
            tibble(model="ZOIB",
                   med_est=mean(X_beta_zoib),
                   high=quantile(X_beta_zoib,.95),
                   low=quantile(X_beta_zoib,.05),
                   sd=sd(X_beta_zoib),
                   loo_val=loo_zoib$estimates[1,1],
                   kurt_est=mean(apply(as.matrix(zoib_fit,"zoib_regen"),1,moments::kurtosis)),
                   win_loo=which(row.names(comp_loo)=="model2"),
                   win_loo_se=comp_loo[which(row.names(comp_loo)=="model2"),2],
                   rmse=rmse_zoib,
                   marg_eff_est=mean(margin_zoib),
                   high_marg=quantile(margin_zoib,.95),
                   low_marg=quantile(margin_zoib,.05),
                   sd_marg=sd(margin_zoib)),
            tibble(model="Beta Regression - Transformed",
                   med_est=mean(X_beta_reg),
                   high=quantile(X_beta_reg,.95),
                   low=quantile(X_beta_reg,.05),
                   sd=sd(X_beta_reg),
                   loo_val=loo_betareg$estimates[1,1],
                   rmse=rmse_betareg,
                   kurt_est=mean(apply(yrep_betareg,1,moments::kurtosis)),
                   marg_eff_est=mean(margin_betareg),
                   high_marg=quantile(margin_betareg,.95),
                   low_marg=quantile(margin_betareg,.05),
                   sd_marg=sd(margin_betareg)),
            tibble(model="Beta Regression - (0,1)",
                   med_est=mean(X_beta_reg2),
                   high=quantile(X_beta_reg2,.95),
                   low=quantile(X_beta_reg2,.05),
                   sd=sd(X_beta_reg2),
                   loo_val=loo_betareg2$estimates[1,1],
                   kurt_est=mean(apply(yrep_betareg2,1,moments::kurtosis)),
                   rmse=rmse_betareg2,
                    marg_eff_est=mean(margin_betareg2),
                    high_marg=quantile(margin_betareg2,.95),
                    low_marg=quantile(margin_betareg2,.05),
                    sd_marg=sd(margin_betareg2)),
            tibble(model="OLS",
                   med_est=mean(X_beta_lm),
                   high=quantile(X_beta_lm,.95),
                   low=quantile(X_beta_lm,.05),
                   sd=sd(X_beta_lm),
                   kurt_est=mean(apply(yrep_lm,1,moments::kurtosis)),
                   loo_val=loo_lm$estimates[1,1],
                   win_loo=which(row.names(comp_loo)=="lm_fit"),
                   win_loo_se=comp_loo[which(row.names(comp_loo)=="lm_fit"),2],
                   rmse=rmse_lm,
                   marg_eff_est=mean(X_beta_lm),
                   high_marg=quantile(X_beta_lm,.95),
                   low_marg=quantile(X_beta_lm,.05),
                   sd_marg=sd(X_beta_lm))))
  
  
  
},simul_data=simul_data,r_seeds=r_seeds,mc.cores=8)

simul_data_final <- bind_rows(all_simul_data)

saveRDS(simul_data_final,"data/sim_cont_X.rds")
