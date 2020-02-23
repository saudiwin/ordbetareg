//
// Ordinal beta regression model for analying experimental outcomes
// with proportion and degenerate responses (i.e. 0 and 1)
// Models 0/1 as ordered categories above/below (0,1) 
// Robert Kubinec
// New York University Abu Dhabi
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
  ordered[2] cutpoints; // cutpoints on ordered (latent) variable (also stand in as intercepts)
  real<lower=0> kappa; // scale parameter for beta regression
}
transformed parameters {
    // store matrix calculations
  
  vector[N_degen] calc_degen;
  vector[N_prop] calc_prop;
  
  // drop the intercepts so everything is relative to the cutpoints
  calc_degen = covar_degen*X_beta;
  calc_prop = covar_prop*X_beta;
  
}
model {
  
  // vague priors
  X_beta ~ normal(0,5);
  kappa ~ exponential(1);
  cutpoints[2] - cutpoints[1] ~ normal(0,3);
  
  // need separate counters for logit (0/1) and beta regression
  
  for(n in 1:N_degen) {
    if(outcome_degen[n]==0) {
      // Pr(Y==0)
      target += log1m_inv_logit(calc_degen[n] - cutpoints[1]);
    } else {
      //Pr(Y==1)
      target += log_inv_logit(calc_degen[n] - cutpoints[2]);
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
  vector[run_gen==0 ? 0 : N_pred_degen+N_pred_prop] ord_log; // store log calculation for loo
  
  if(run_gen==1) {
    if(N_pred_degen>0) {
     // first do degenerate outcomes 
    // note: these could be *re-generated* as beta/propotions
    for(i in 1:num_elements(indices_degen)) {
      
      // draw an outcome 0 / prop / 1
      regen_degen[i] = ordered_logistic_rng(covar_degen[indices_degen[i],]*X_beta,cutpoints);
      
      if(outcome_degen[i]==0) {
        ord_log[i] = log1m_inv_logit(calc_degen[i] - cutpoints[1]);
      } else {
        ord_log[i] = log_inv_logit(calc_degen[i] - cutpoints[2]);
      }
      
      
      if(regen_degen[i]==1) {
        regen_all[i] = 0;
      } else if(regen_degen[i]==3) {
        regen_all[i] = 1;
      } else {
        // did not occur in original data but could re-occur probabilistically
        regen_all[i] = beta_proportion_rng(inv_logit(covar_degen[indices_degen[i],]*X_beta),kappa);
      }
      
  }
  
  if(N_pred_prop>0) {
        // now do originally proportional outcomes
      // can be re-generated as 0s or 1s
      
      int skip = num_elements(indices_degen);
      
      for(i in 1:num_elements(indices_prop)) {
        
        // draw an outcome 0 / prop / 1
        regen_degen[i+skip] = ordered_logistic_rng(covar_prop[indices_prop[i],]*X_beta,cutpoints);
        
        ord_log[i+skip] = log(inv_logit(calc_prop[indices_prop[i]] - cutpoints[1]) - inv_logit(calc_prop[indices_prop[i]] - cutpoints[2])) +
                        beta_proportion_lpdf(outcome_prop[indices_prop[i]]|inv_logit(calc_prop[indices_prop[i]]),kappa);
        
        if(regen_degen[i+skip]==1) {
          regen_all[i+skip] = 0;
        } else if(regen_degen[i+skip]==3) {
          regen_all[i+skip] = 1;
        } else {
          // did not occur in original data but could re-occur probabilistically
          regen_all[i+skip] = beta_proportion_rng(inv_logit(covar_prop[indices_prop[i],]*X_beta),kappa);
        }
        
      } 
      }
    }
  }
  
}

