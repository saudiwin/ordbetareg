//
// Ordinal beta regression model for analying experimental outcomes
// with proportion and degenerate responses (i.e. 0 and 1)
// Models 0/1 as ordered categories above/below (0,1) 
// Robert Kubinec
// New York University Abu Dhabi

functions {
  
  // prior from Michael Betancourt for ordered cutpoints
    real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data {
  int<lower=0> N_prop; // number of proportion observations (0,1)
  int<lower=0> N_degen; // number of 0/1 observations
  int X; // number predictors
  vector[N_prop] outcome_prop; // Y in (0,1)
  int outcome_degen[N_degen]; // Y in {0,1}
  matrix[N_prop,X] covar_prop; // covariate X for proportion outcome
  matrix[N_degen,X] covar_degen; // covariate X for degenerate (0,1) outcome
  int N_pred_degen; // number of posterior predictive samples for 0/1
  int N_pred_prop; // number of posterior predictive samples for (0,1)
  int indices_degen[N_pred_degen]; // random row indices to use for posterior predictive calculation of 0/1
  int indices_prop[N_pred_prop]; // random row indices to use for posterior predictive calculation of (0,1)
  int run_gen; // whether to use generated quantities
}
parameters {
  vector[X] X_beta; // common predictor
  real alpha; // common intercept
  ordered[2] cutpoints; // cutpoints on ordered (latent) variable (also stand in as intercepts)
  real<lower=0> kappa; // scale parameter for beta regression
}
transformed parameters {
    // store matrix calculations
  
  vector[N_degen] calc_degen;
  vector[N_prop] calc_prop;
  
  // drop the intercepts so everything is relative to the cutpoints
  if(N_degen>0) {
    calc_degen = alpha + covar_degen*X_beta;
  }
  
  //print(calc_degen[1:10]);
  
  calc_prop = alpha + covar_prop*X_beta;
  
  //print(calc_prop[1:10]);
  
}
model {
  
  // vague priors
  X_beta ~ normal(0,5);
  alpha ~ normal(0,5);
  kappa ~ exponential(.1);
  //cutpoints[2] - cutpoints[1] ~ normal(0,3);
  // induced dirichlet prior on cutpoints:
  
  target += induced_dirichlet_lpdf(cutpoints | rep_vector(1, 3), 0);
  
  // need separate counters for logit (0/1) and beta regression
  if(N_degen>0) {
    for(n in 1:N_degen) {
    if(outcome_degen[n]==0) {
      // Pr(Y==0)
      target += log1m_inv_logit(calc_degen[n] - cutpoints[1]);
    } else {
      //Pr(Y==1)
      target += log_inv_logit(calc_degen[n] - cutpoints[2]);
    }
    }
  }
  
  
  for(n in 1:N_prop) {
    // Pr(Y in (0,1))
    target += log(inv_logit(calc_prop[n] - cutpoints[1]) - inv_logit(calc_prop[n] - cutpoints[2]));
    // Pr(Y==x where x in (0,1))
    outcome_prop[n] ~ beta_proportion(inv_logit(calc_prop[n]),kappa);
  }
  
}

generated quantities {
  
  vector[run_gen==0 ? 0 : N_pred_degen+N_pred_prop] regen_degen; // which model is selected (degenerate or proportional)
  vector[run_gen==0 ? 0 : N_pred_degen+N_pred_prop] regen_all; // final (combined) outcome -- defined as random subset of rows
  vector[run_gen==0 ? 0 : N_pred_degen+N_pred_prop] regen_epred;
  vector[run_gen==0 ? 0 : N_pred_degen+N_pred_prop] ord_log; // store log calculation for loo
  
  if(run_gen==1) {
    
    row_vector[2] all_pr[N_pred_degen+N_pred_prop];
    
    if(N_pred_degen>0) {
     // first do degenerate outcomes 
    // note: these could be *re-generated* as beta/propotions
    for(i in 1:num_elements(indices_degen)) {
      
      // draw an outcome 0 / prop / 1
      regen_degen[i] = ordered_logistic_rng(alpha + covar_degen[indices_degen[i],]*X_beta,cutpoints);
      
      if(outcome_degen[i]==0) {
        ord_log[i] = log1m_inv_logit(calc_degen[indices_degen[i]] - cutpoints[1]);
      } else {
        ord_log[i] = log_inv_logit(calc_degen[indices_degen[i]] - cutpoints[2]);
      }
      
      // don't need zero
      
      all_pr[i] = [log(inv_logit(calc_degen[indices_degen[i]] - cutpoints[1]) - inv_logit(calc_degen[i] - cutpoints[2])) + log_inv_logit(calc_degen[indices_degen[i]]),
                    log_inv_logit(calc_degen[indices_degen[i]] - cutpoints[2])];
      
      if(regen_degen[i]==1) {
        regen_all[i] = 0;
      } else if(regen_degen[i]==3) {
        regen_all[i] = 1;
      } else {
        // did not occur in original data but could re-occur probabilistically
        regen_all[i] = beta_proportion_rng(inv_logit(alpha + covar_degen[indices_degen[i],]*X_beta),kappa);
      }
     }
    }
  if(N_pred_prop>0) {
        // now do originally proportional outcomes
      // can be re-generated as 0s or 1s
      
      int skip = num_elements(indices_degen);
      
      for(i in 1:num_elements(indices_prop)) {
        
        all_pr[i+skip] = [log(inv_logit(calc_prop[indices_prop[i]] - cutpoints[1]) - inv_logit(calc_prop[indices_prop[i]] - cutpoints[2])) + log_inv_logit(calc_prop[indices_prop[i]]),
        log_inv_logit(calc_prop[indices_prop[i]] - cutpoints[2])];
        
        // draw an outcome 0 / prop / 1
        regen_degen[i+skip] = ordered_logistic_rng(alpha + covar_prop[indices_prop[i],]*X_beta,cutpoints);
        
        ord_log[i+skip] = log(inv_logit(calc_prop[indices_prop[i]] - cutpoints[1]) - inv_logit(calc_prop[indices_prop[i]] - cutpoints[2])) +
                        beta_proportion_lpdf(outcome_prop[indices_prop[i]]|inv_logit(calc_prop[indices_prop[i]]),kappa);
        
        if(regen_degen[i+skip]==1) {
          regen_all[i+skip] = 0;
        } else if(regen_degen[i+skip]==3) {
          regen_all[i+skip] = 1;
        } else {
          // did not occur in original data but could re-occur probabilistically
          regen_all[i+skip] = beta_proportion_rng(inv_logit(alpha + covar_prop[indices_prop[i],]*X_beta),kappa);
        }
        
      } 
      }
      
      for(i in 1:(N_pred_degen + N_pred_prop)) {
        
        regen_epred[i] = sum(exp(all_pr[i]));
        
      }
  }
  
}

