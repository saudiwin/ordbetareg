# Robert Kubinec
# New York University Abu Dhabi
# October 12, 2021
# Simulation of 0 - 1 bounded dependent variables
# Note simulation will take some time, approx ~2 days with 3 cores

require(cmdstanr)
require(bayesplot)
require(dplyr)
require(brms)
require(loo)
require(posterior)
require(future.apply)
require(faux)

RNGkind("L'Ecuyer-CMRG")

beta_logit <- cmdstanr::cmdstan_model("beta_logit.stan",compile = T)
zoib_model <- cmdstanr::cmdstan_model("zoib_nophireg.stan")
frac_mod <- cmdstanr::cmdstan_model("frac_logit.stan")

set.seed(772235)

# let's do some simulations
  
N_rep <- 10000

simul_data <- tibble(N=round(runif(N_rep,100,4000),0),
                     k=floor(runif(N_rep,1,15)),
                     rho=runif(N_rep,0,.9),
                     phi=runif(N_rep,0.5,30),
                     cutpoints1=runif(N_rep,-10,-1)) %>% 
  mutate(cutpoints2=cutpoints1+runif(N_rep,0.5,10),
         X_beta=sapply(k, function(i) runif(i,-5,5)))

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
  eta <- X%*%matrix(X_beta)[,1]
  
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
  psi <- plogis(alpha1 + X %*% t(coef_a))
  gamma <- plogis(alpha2 + X %*% t(coef_g))
  eta <- alpha3 + X %*% t(coef_m)
  
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

gen_x <- function(k,rho,N_rep) {
  
  # generate true coefs
  
  true_coef <- runif(k, -5, 5)
  
  # generate matrix of correlated Xs
  
  out_x <- sapply(1:length(k), function(i) {
    
    this_x <- 1
      
  })
  
}

# one for each model type

r_seeds <- c(6635,2216,8845,9936,3321,63914)

plan(multicore,workers=40)

all_simul_data <- future_lapply(1:nrow(simul_data), function(i,simul_data=NULL,r_seeds=NULL) {
#all_simul_data <- lapply(1:nrow(simul_data), function(i,simul_data=NULL,r_seeds=NULL) {
  
  this_data <- slice(simul_data,i)
  cat(file = "simul_status.txt",paste0("Now on row ",i),append = T)

# Draw from ordered beta regression ---------------------------------------

  
  
  N <- this_data$N
  
  #X <- rnorm(N,runif(1,-2,2),1)
  
  X_beta <- this_data$X_beta[[1]]
  
  # need to create X
  
  X <- rnorm_multi(n=N,vars=this_data$k,r=this_data$rho,as.matrix=T)
  
  eta <- -2 + X%*%matrix(X_beta)
  
  # ancillary parameter of beta distribution
  phi <- this_data$phi
  
  # predictor for ordered model
  mu1 <- eta[,1]
  # predictor for beta regression
  mu2 <- eta[,1]
  
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
  
  # check for floating point errors
  
  final_out <- ifelse(final_out>(1 - 1e-10) & final_out<1,final_out - 1e-10,
                               final_out)
  
  final_out <- ifelse(final_out<(0 + 1e-10) & final_out>0,final_out + 1e-10,
                               final_out)
  
  counter <- environment()
  counter$i <- 0
  
  # get rid of very rare outcomes, can be difficult for ZOIB
  
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
    
    print(paste0("DGP failed for row ",i,"\n"))
    
    this_data$status <- "dgp_failure"
    
    return(this_data)
  }
  
  # calculate "true" marginal effects
  
  # loop over k
  
  X_low <- X
  X_high <- X
  
  marg_eff <- sapply(1:this_data$k, function(tk)  {
    
    X_low[,tk] <- X[,tk] - setstep(X[,tk])
    X_high[,tk] <- X[,tk] + setstep(X[,tk])
    
    y0 <- predict_ordbeta(cutpoints=cutpoints,
                          X=X_low,
                          X_beta=X_beta)
    
    y1 <- predict_ordbeta(cutpoints=cutpoints,
                          X=X_high,
                          X_beta=X_beta)
    
    mean((y1-y0)/(X_high[,tk]-X_low[,tk]))
    
    
  })
  
  
  

# Fit models --------------------------------------------------------------
  
  # now fit ordinal beta
  
  to_bl <- list(N_degen=sum(final_out %in% c(0,1)),
                N_prop=sum(final_out>0 & final_out<1),
                X=this_data$k,
                outcome_prop=final_out[final_out>0 & final_out<1],
                outcome_degen=final_out[final_out %in% c(0,1)],
                covar_prop=X[final_out>0 & final_out<1,,drop=F],
                covar_degen=X[final_out %in% c(0,1),,drop=F],
                N_pred_degen=sum(final_out %in% c(0,1)),
                N_pred_prop=sum(final_out>0 & final_out<1),
                indices_degen=1:(sum(final_out %in% c(0,1))),
                indices_prop=1:(sum(final_out>0 & final_out<1)),
                run_gen=1)
  
  x <- X
  
  final_out_scale <- (final_out * (length(final_out) - 1) + .5) / length(final_out)
  
  frac_data <- list(y=final_out,
                    x=x,
                    k=ncol(x),
                    run_gen=1,
                    n=nrow(x))
  # fit all models
  
  
  fit_model <-try(beta_logit$sample(data=to_bl,seed=r_seeds[1],
                        chains=1,parallel_chains=1,iter_warmup=500,
                        refresh=0,
                        iter_sampling=500))
  
  
 
  zoib_fit <- try(zoib_model$sample(data=list(n=length(final_out),
                                            y=final_out,
                                            k=ncol(x),
                                            x=x,
                                            run_gen=1), 
                                    seed=r_seeds[2],
                                    refresh=0,chains=1,parallel_chains=1,iter_warmup=500,
                                iter_sampling=500))
  
  X_brms <- X
  colnames(X_brms) <- paste0(rep("Var",ncol(X)),1:ncol(X))
  X_brms <- as_tibble(X_brms)
  
  betareg_fit <- try(brm(formula = outcome~.,data=mutate(X_brms,
                                                         outcome=final_out_scale),
                     chains=1,cores=1,iter=1000,
                              seed=r_seeds[3],
                     silent=0,refresh=0,
                     family="beta",
                     prior=set_prior("normal(0,5)", class = "b"),
                     backend='cmdstanr'))
  
  X_brms_small <- filter(X_brms, final_out>0 & final_out<1)
  
  X_brms_small$outcome <- final_out[final_out>0 & final_out<1]
  
  betareg_fit2 <- try(update(betareg_fit,newdata=X_brms_small,
                             chains=1,cores=1,iter=1000,
                             seed=r_seeds[4],
                             silent=0,refresh=0,
                             family="beta",
                             backend="cmdstanr"))
  
  frac_fit <- try(frac_mod$sample(data=frac_data,seed=r_seeds[6],
                                  chains=1,parallel_chains = 1,refresh=0,
                                  iter_warmup=500,iter_sampling = 500))
  
  lm_fit <- try(brm(formula = outcome~X,data=tibble(outcome=final_out,
                                                    r_seeds[5],
                                                    X=X),chains=1,cores=1,iter=1000,
                    silent=0,refresh=0,
                    backend="cmdstanr",
                    prior=set_prior("normal(0,5)", class = "b")))
  
  if('try-error' %in% c(class(fit_model),
                        class(zoib_model),
                        class(betareg_fit),
                        class(betareg_fit2),
                        class(frac_fit),
                        class(lm_fit))) {
    
    print(paste0("Estimation failed for row ",i,"\n"))
    
    this_data$status <- "estimation_failure"
    
    return(this_data)
    
  }
  
  yrep_ord <- try(as_draws_matrix(fit_model$draws("regen_epred")))
  
  yrep_zoib <- try(as_draws_matrix(zoib_fit$draws("zoib_epred")))
  
  yrep_betareg <- try(posterior_epred(betareg_fit,draws=100))

  yrep_betareg2 <- try(posterior_epred(betareg_fit2,draws=100))

  yrep_lm <- try(posterior_epred(lm_fit,draws=100))
  
  yrep_frac <- try(as_draws_matrix(frac_fit$draws("frac_rep")))
  
  if('try-error' %in% c(class(yrep_ord),
                        class(yrep_zoib),
                        class(yrep_betareg),
                        class(yrep_betareg2),
                        class(yrep_frac),
                        class(yrep_lm))) {
    
    print(paste0("Estimation failed for row ",i,"\n"))
    
    this_data$status <- "estimation_failure"
    
    return(this_data)
    
  }
  
  this_data$status <- "success"
  

# Calculate estimands -----------------------------------------------------

  
  # now return the full data frame
  
  X_beta_ord <- as_draws_matrix(fit_model$draws("X_beta"))
  X_beta_zoib <- as_draws_matrix(zoib_fit$draws("coef_m"))
  X_beta_reg <- as.matrix(betareg_fit,pars=paste0("Var",1:this_data$k))
  X_beta_reg2 <- as.matrix(betareg_fit2,pars=paste0("Var",1:this_data$k))
  X_beta_lm <- as.matrix(lm_fit,pars="X")
  X_beta_frac <- as_draws_matrix(frac_fit$draws(variables="X_beta"))
  
  # calculate rmse
  
  rmse_ord <- sqrt(mean(apply(yrep_ord,1,function(c) { (c - c(final_out[final_out %in% c(0,1)],final_out[final_out>0 & final_out<1]))^2 })))
  rmse_zoib <-sqrt( mean(apply(yrep_zoib,1,function(c) { (c - final_out)^2 })))
  rmse_betareg <- sqrt(mean(apply(yrep_betareg,1,function(c) { (c - final_out)^2 })))
  rmse_betareg2 <- sqrt(mean(apply(yrep_betareg2,1,function(c) { (c - final_out[final_out>0 & final_out<1])^2 })))
  rmse_lm <- sqrt(mean(apply(yrep_lm,1,function(c) { (c - final_out)^2 })))
  rmse_frac <- sqrt(mean(apply(yrep_frac,1,function(c) { (c - final_out)^2 })))
  
  # calculate loo
  
  loo_ordbeta <-try(fit_model$loo("ord_log"))
  loo_betareg <- loo(betareg_fit)
  loo_betareg2 <- loo(betareg_fit2)
  loo_zoib <- try(zoib_fit$loo("zoib_log"))
  loo_lm <- loo(lm_fit)
  loo_frac <- try(frac_fit$loo("frac_log"))
  
  if(any('try-error' %in% c(class(loo_ordbeta),
                            class(loo_zoib),
                            class(loo_frac)))) {
    comp_loo <- matrix(rep(NA,6),ncol=2)
    row.names(comp_loo) <- c("ord","zoib","frac")
    loo_ordbeta <- list(estimates=matrix(c(NA,NA),ncol=2))
    loo_zoib <- list(estimates=matrix(c(NA,NA),ncol=2))
    loo_frac <- list(estimates=matrix(c(NA,NA),ncol=2))
  } else {
    comp_loo <- loo_compare(list(ord=loo_ordbeta,zoib=loo_zoib,
                                 lm=loo_lm,frac=loo_frac))
  }
  
  
  
  # calculate marginal effects
  
  cutpoints_est <- as_draws_matrix(fit_model$draws("cutpoints"))
  
  margin_ord <- lapply(1:ncol(X), function(tk) {
    
    X_low[,tk] <- X[,tk] - setstep(X[,tk])
    X_high[,tk] <- X[,tk] + setstep(X[,tk])
    
    tibble(marg_eff=sapply(1:nrow(X_beta_ord), function(i) {
      y0 <- predict_ordbeta(cutpoints=cutpoints_est[i,],
                            X=X_low,
                            X_beta=c(X_beta_ord[i,]))
      
      y1 <- predict_ordbeta(cutpoints=cutpoints_est[i,],
                            X=X_high,
                            X_beta=c(X_beta_ord[i,]))
      
      marg_eff <- (y1-y0)/(X_high[,tk]-X_low[,tk])
      
      mean(marg_eff)
    }),
    x_col=tk)
    
    }) %>% bind_rows
  
  # now for the ZOIB
  
  coef_a <- as_draws_matrix(zoib_fit$draws("coef_a"))
  coef_g <- as_draws_matrix(zoib_fit$draws("coef_g"))
  alpha <- as_draws_matrix(zoib_fit$draws("alpha"))
  
  margin_zoib <- lapply(1:ncol(X), function(tk) {
    
    X_low[,tk] <- X[,tk] - setstep(X[,tk])
    X_high[,tk] <- X[,tk] + setstep(X[,tk]) 
    
    tibble(marg_eff=    sapply(1:nrow(X_beta_zoib), function(i) {
      y0 <- predict_zoib(coef_g=coef_g[i,],
                         coef_a=coef_a[i,],
                         alpha1=c(alpha[i,1]),
                         alpha2=c(alpha[i,2]),
                         alpha3=c(alpha[i,3]),
                         X=X_low,
                         coef_m=X_beta_zoib[i,])
      
      y1 <-  predict_zoib(coef_g=coef_g[i,],
                          coef_a=coef_a[i,],
                          alpha1=c(alpha[i,1]),
                          alpha2=c(alpha[i,2]),
                          alpha3=c(alpha[i,3]),
                          X=X_high,
                          coef_m=X_beta_zoib[i,])
      
      marg_eff <- (y1-y0)/(X_high[,tk]-X_low[,tk])
      
      mean(marg_eff)
    }),
    x_col=tk)
    
    }) %>% bind_rows
  
  # now betareg
  
  betareg_int <- as.matrix(betareg_fit,pars="(Intercept)")
  
  margin_betareg <- lapply(1:ncol(X), function(tk) {
    
    X_low[,tk] <- X[,tk] - setstep(X[,tk])
    X_high[,tk] <- X[,tk] + setstep(X[,tk]) 
    
    tibble(marg_eff= sapply(1:nrow(X_beta_reg), function(i) {
      y0 <- plogis(betareg_int[i,"b_Intercept"] + X_low%*%X_beta_reg[i,])
      y1 <- plogis(betareg_int[i,"b_Intercept"] + X_high%*%X_beta_reg[i,])
      
      marg_eff <- (y1-y0)/(X_high[,tk]-X_low[,tk])
      
      mean(marg_eff)
    }),
    x_col=tk)
    
    }) %>% bind_rows
  
  betareg2_int <- as.matrix(betareg_fit2,pars="(Intercept)")
  
  margin_betareg2 <- lapply(1:ncol(X), function(tk) {
    
    X_low[,tk] <- X[,tk] - setstep(X[,tk])
    X_high[,tk] <- X[,tk] + setstep(X[,tk]) 
    
    tibble(marg_eff= sapply(1:nrow(X_beta_reg), function(i) {
      y0 <- plogis(betareg2_int[i,"b_Intercept"] + X_low%*%X_beta_reg2[i,])
      y1 <- plogis(betareg2_int[i,"b_Intercept"] + X_high%*%X_beta_reg2[i,])
      
      marg_eff <- (y1-y0)/(X_high[,tk]-X_low[,tk])
      
      mean(marg_eff)
    }),
    x_col=tk)
    
  }) %>% bind_rows
  
  
  # Fractional logit marg effects -------------------------------------------
  
  frac_int <- as_draws_matrix(frac_fit$draws(variables="alpha"))
  
  margin_frac <- lapply(1:ncol(X), function(tk) {
    
    X_low[,tk] <- X[,tk] - setstep(X[,tk])
    X_high[,tk] <- X[,tk] + setstep(X[,tk]) 
    
    tibble(marg_eff= sapply(1:nrow(X_beta_frac), function(i) {
      y0 <- plogis(c(frac_int[i,]) + X_low%*%X_beta_frac[i,,drop=T])
      y1 <- plogis(c(frac_int[i,]) + X_high%*%X_beta_frac[i,,drop=T])
      
      marg_eff <- (y1-y0)/(X_high[,tk]-X_low[,tk])
      
      mean(marg_eff)
    }),
    x_col=tk)
    
    
    }) %>% bind_rows
  
  
  
  this_data$marg_eff <- list(marg_eff)
  this_data$true_kurt <- moments::kurtosis(final_out)
  
  
  
  # Combine estimates -------------------------------------------------------
  
  sum_marg <- function(d,func,...) {
    
   ret_vec <-  arrange(d,x_col) %>% 
      group_by(x_col) %>% 
      summarize(sum_stat=func(marg_eff,...)) %>% 
      pull(sum_stat)
   
   list(ret_vec)
   
  }
  
  try(bind_cols(purrr::map_dfr(seq_len(6), ~this_data),bind_rows(tibble(model="Ordinal Beta Regression",
                                                                    med_est=list(apply(X_beta_ord,2,mean)),
                                                                    high=list(apply(X_beta_ord,2,quantile,.95)),
                                                                    low=list(apply(X_beta_ord,2,quantile,.05)),
                                                                    var_calc=list(apply(X_beta_ord,2,var)),
                                                                    loo_val=loo_ordbeta$estimates[1,1],
                                                                    p_loo=loo_ordbeta$estimates[2,1],
                                                                    bad_k=sum(loo_ordbeta$diagnostics$pareto_k>0.5),
                                                                    win_loo=which(row.names(comp_loo)=="ord"),
                                                                    win_loo_se=comp_loo[which(row.names(comp_loo)=="ord"),2],
                                                                    rmse=rmse_ord,
                                                                    kurt_est=mean(apply(yrep_ord,1,moments::kurtosis)),
                                                                    marg_eff_est=sum_marg(margin_ord,mean),
                                                                    high_marg=sum_marg(margin_ord,quantile,.95),
                                                                    low_marg=sum_marg(margin_ord,quantile,.05),
                                                                    var_marg=sum_marg(margin_ord,var)),
                                                             tibble(model="ZOIB",
                                                                    med_est=list(apply(X_beta_zoib,2,mean)),
                                                                    high=list(apply(X_beta_zoib,2,quantile,.95)),
                                                                    low=list(apply(X_beta_zoib,2,quantile,.05)),
                                                                    var_calc=list(apply(X_beta_zoib,2,var)),
                                                                    loo_val=loo_zoib$estimates[1,1],
                                                                    p_loo=loo_zoib$estimates[2,1],
                                                                    bad_k=sum(loo_zoib$diagnostics$pareto_k>0.5),
                                                                    kurt_est=mean(apply(yrep_zoib,1,moments::kurtosis)),
                                                                    win_loo=which(row.names(comp_loo)=="zoib"),
                                                                    win_loo_se=comp_loo[which(row.names(comp_loo)=="zoib"),2],
                                                                    rmse=rmse_zoib,
                                                                    marg_eff_est=sum_marg(margin_zoib,mean),
                                                                    high_marg=sum_marg(margin_zoib,quantile,.95),
                                                                    low_marg=sum_marg(margin_zoib,quantile,.05),
                                                                    var_marg=sum_marg(margin_zoib,var)),
                                                             tibble(model="Beta Regression - Transformed",
                                                                    med_est=list(apply(X_beta_reg,2,mean)),
                                                                    high=list(apply(X_beta_reg,2,quantile,.95)),
                                                                    low=list(apply(X_beta_reg,2,quantile,.05)),
                                                                    var_calc=list(apply(X_beta_reg,2,var)),
                                                                    loo_val=loo_betareg$estimates[1,1],
                                                                    p_loo=loo_betareg$estimates[2,1],
                                                                    bad_k=sum(loo_betareg$diagnostics$pareto_k>0.5),
                                                                    rmse=rmse_betareg,
                                                                    kurt_est=mean(apply(yrep_betareg,1,moments::kurtosis)),
                                                                    marg_eff_est=sum_marg(margin_betareg,mean),
                                                                    high_marg=sum_marg(margin_betareg,quantile,.95),
                                                                    low_marg=sum_marg(margin_betareg,quantile,.05),
                                                                    var_marg=sum_marg(margin_betareg,var)),
                                                             tibble(model="Beta Regression - (0,1)",
                                                                    med_est=list(apply(X_beta_reg2,2,mean)),
                                                                    high=list(apply(X_beta_reg2,2,quantile,.95)),
                                                                    low=list(apply(X_beta_reg2,2,quantile,.05)),
                                                                    var_calc=list(apply(X_beta_reg2,2,var)),
                                                                    loo_val=loo_betareg2$estimates[1,1],
                                                                    p_loo=loo_betareg2$estimates[2,1],
                                                                    bad_k=sum(loo_betareg2$diagnostics$pareto_k>0.5),
                                                                    kurt_est=mean(apply(yrep_betareg2,1,moments::kurtosis)),
                                                                    rmse=rmse_betareg2,
                                                                    marg_eff_est=sum_marg(margin_betareg2,mean),
                                                                    high_marg=sum_marg(margin_betareg2,quantile,.95),
                                                                    low_marg=sum_marg(margin_betareg2,quantile,.05),
                                                                    var_marg=sum_marg(margin_betareg2,var)),
                                                             tibble(model="OLS",
                                                                    med_est=list(apply(X_beta_lm,2,mean)),
                                                                    high=list(apply(X_beta_lm,2,quantile,.95)),
                                                                    low=list(apply(X_beta_lm,2,quantile,.05)),
                                                                    var_calc=list(apply(X_beta_lm,2,var)),
                                                                    kurt_est=mean(apply(yrep_lm,1,moments::kurtosis)),
                                                                    loo_val=loo_lm$estimates[1,1],
                                                                    p_loo=loo_lm$estimates[2,1],
                                                                    bad_k=sum(loo_lm$diagnostics$pareto_k>0.5),
                                                                    win_loo=which(row.names(comp_loo)=="lm"),
                                                                    win_loo_se=comp_loo[which(row.names(comp_loo)=="lm"),2],
                                                                    rmse=rmse_lm,
                                                                    marg_eff_est=list(apply(X_beta_lm,2,mean)),
                                                                    high_marg=list(apply(X_beta_lm,2,quantile,.95)),
                                                                    low_marg=list(apply(X_beta_lm,2,quantile,.05)),
                                                                    var_marg=list(apply(X_beta_lm,2,var))),
                                                             tibble(model="Fractional",
                                                                    med_est=list(apply(X_beta_frac,2,mean)),
                                                                    high=list(apply(X_beta_frac,2,quantile,.95)),
                                                                    low=list(apply(X_beta_frac,2,quantile,.05)),
                                                                    var_calc=list(apply(X_beta_frac,2,var)),
                                                                    loo_val=loo_frac$estimates[1,1],
                                                                    p_loo=loo_frac$estimates[2,1],
                                                                    bad_k=sum(loo_frac$diagnostics$pareto_k>0.5),
                                                                    win_loo=which(row.names(comp_loo)=="frac"),
                                                                    win_loo_se=comp_loo[which(row.names(comp_loo)=="frac"),2],
                                                                    kurt_est=mean(apply(yrep_frac,1,moments::kurtosis)),
                                                                    rmse=rmse_frac,
                                                                    marg_eff_est=sum_marg(margin_frac,mean),
                                                                    high_marg=sum_marg(margin_frac,quantile,.95),
                                                                    low_marg=sum_marg(margin_frac,quantile,.05),
                                                                    var_marg=sum_marg(margin_frac,var)))))
  
  
#},simul_data=simul_data,r_seeds=r_seeds) 
},simul_data=simul_data,r_seeds=r_seeds,future.seed=TRUE)

#simul_data_final <- bind_rows(all_simul_data)

saveRDS(all_simul_data,"data/sim_cont_X.rds")
all_sim <- all_simul_data
save(all_sim, file="data/sim_cont_X.RData")
