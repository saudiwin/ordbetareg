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
  //int X_miss; // number of predictors for inflated model
  vector[N_prop] outcome_prop; // Y in (0,1)
  real infl_value; // set to value between 0 and 1. If negative, inflation is not used
  int outcome_degen[N_degen]; // Y in {0,1}
  matrix[N_prop,X] covar_prop; // covariate X for proportion outcome
  //matrix[N_prop,X_miss] covar_prop_infl; // covariate X for inflated values
  //matrix[N_degen,X_miss] covar_degen_infl; // covariate X for inflated values
  matrix[N_degen,X] covar_degen; // covariate X for degenerate (0,1) outcome
  int N_pred_degen; // number of posterior predictive samples for 0/1
  int N_pred_prop; // number of posterior predictive samples for (0,1)
  int indices_degen[N_pred_degen]; // random row indices to use for posterior predictive calculation of 0/1
  int indices_prop[N_pred_prop]; // random row indices to use for posterior predictive calculation of (0,1)
  int run_gen; // whether to use generated quantities
}
transformed data {
  int infl_this[N_prop];
  int miss_id[N_prop]; // index for imputed outcome
  int miss_counter = 1;
  
  for(i in 1:N_prop) {
    if(outcome_prop[i]==infl_value) {
      infl_this[i] = 1;
      miss_id[i] = miss_counter;
      miss_counter += 1;
    } else {
      infl_this[i] = 0;
      miss_id[i] = miss_counter;
    }
  }
}
parameters {
  vector[X] X_beta; // common predictor
  vector<lower=0.01,upper=0.99>[miss_counter] prop_miss;
  //vector[!(infl_value<0) ? X_miss : 0] X_beta_miss; // predictor for inflated values
  ordered[2] cutpoints; // cutpoints on ordered (latent) variable (also stand in as intercepts)
  real<lower=0> kappa; // scale parameter for beta regression
  real<lower=0,upper=1> theta; // probability of an observation being missing versus non-missing 
}
transformed parameters {
    // store matrix calculations
  
  vector[N_degen] calc_degen;
  vector[N_prop] calc_prop;
  //vector[!(infl_value<0) ? N_prop : 0] calc_miss;
  //vector[!(infl_value<0) ? N_pred_degen : 0] calc_degen_miss; // must be defined over both distributionss
  
  // drop the intercepts so everything is relative to the cutpoints
  calc_degen = covar_degen*X_beta;
  calc_prop = covar_prop*X_beta;
  
  
  if(!(infl_value<0)) {
   // calc_miss = covar_prop_infl*X_beta_miss;
  //  calc_degen_miss = covar_degen_infl*X_beta_miss;
  }
  
}
model {
  
  // vague priors
  X_beta ~ normal(0,5);
  //X_beta_miss ~ normal(0,5);
  theta ~ normal(0,2); // tight prior around the number of possible missings 
  kappa ~ exponential(1);
  prop_miss ~ normal(0,1); // put the majority of prior miss around 0.5
  cutpoints[2] - cutpoints[1] ~ normal(0,3);
  // for(c in 1:3)
  //   cutpoints[c+1] - cutpoints[c] ~ normal(0,3);
  
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
  
  //target += beta_proportion_lpdf(outcome_prop|inv_logit(calc_prop),kappa);
  
    for(n in 1:N_prop) {
      // - beta_proportion_lcdf(0.49|inv_logit(calc_prop[n]),kappa)
      //- beta_proportion_lccdf(0.51|inv_logit(calc_prop[n]),kappa)
    if(infl_this[n]==infl_value) {
      target += log(inv_logit(calc_prop[n] - cutpoints[1]) - inv_logit(calc_prop[n] - cutpoints[2])) +
                log_mix(theta,
                  beta_proportion_lpdf(prop_miss[miss_id[n]]|inv_logit(calc_prop[n]),kappa),
                    beta_proportion_lpdf(outcome_prop[n]|inv_logit(calc_prop[n]),kappa));
    }  else {
      target += log(inv_logit(calc_prop[n] - cutpoints[1]) - inv_logit(calc_prop[n] - cutpoints[2]));
      outcome_prop[n] ~ beta_proportion(inv_logit(calc_prop[n]),kappa);
    } 
    
    }
  
  
}

generated quantities {
  
  vector[run_gen==0 ? 0 : N_pred_degen+N_pred_prop] regen_degen; // which model is selected (degenerate or proportional)
  vector[run_gen==0 ? 0 : N_pred_degen+N_pred_prop] regen_all; // final (combined) outcome -- defined as random subset of rows
  vector[run_gen==0 ? 0 : N_pred_degen+N_pred_prop] ord_log; // store log calculation for loo
  int infl_gen[(run_gen==0 || infl_value<0) ? 0 : N_pred_prop]; // whether observation belongs to inflated value or not
  
  
  if(run_gen==1) {
    
    if(N_pred_degen>0) {
    
     // first do degenerate outcomes 
    // note: these could be *re-generated* as beta/propotions
    for(i in 1:num_elements(indices_degen)) {
      
      // draw an outcome 0 / prop / 1
      regen_degen[i] = ordered_logistic_rng(calc_degen[i],cutpoints);
      
      if(outcome_degen[i]==1) {
        ord_log[i] = log1m_inv_logit(calc_degen[i] - cutpoints[1]);
      } else {
        ord_log[i] = log_inv_logit(calc_degen[i] - cutpoints[2]);
      }
      
       if(regen_degen[i]==1) {
        regen_all[i] = 0;
      } else if(regen_degen[i]==3) {
        regen_all[i] = 1;
      } else {
          regen_all[i] = beta_proportion_rng(inv_logit(calc_prop[i]),kappa);
      } 
      
      }
    }
  
  if(N_pred_prop>0) {
        // now do originally proportional outcomes
      // can be re-generated as 0s or 1s
      
      int skip = num_elements(indices_degen);
      
      for(i in 1:num_elements(indices_prop)) {
        
        // draw an outcome 0 / prop / 1
        regen_degen[i+skip] = ordered_logistic_rng(calc_prop[i],cutpoints);
        infl_gen[i] = bernoulli_rng(theta);
        
        //ord_log[i+skip] = log(inv_logit(calc_prop[i] - cutpoints[1]) - inv_logit(calc_prop[i] - cutpoints[2]));
        
        if(infl_this[i]==0) {
            ord_log[i+skip] = log(inv_logit(calc_prop[i] - cutpoints[1]) - inv_logit(calc_prop[i] - cutpoints[2]));
            ord_log[i+skip] += beta_proportion_lpdf(outcome_prop[i]|inv_logit(calc_prop[i]),kappa);
        } else {
          ord_log[i+skip] = log(inv_logit(calc_prop[i] - cutpoints[1]) - inv_logit(calc_prop[i] - cutpoints[2]));
          ord_log[i+skip] += log_mix(theta,
                      beta_proportion_lpdf(prop_miss[miss_id[i]]|inv_logit(calc_prop[i]),kappa),
                      beta_proportion_lpdf(outcome_prop[i]|inv_logit(calc_prop[i]),kappa));
        }
        
        if(regen_degen[i+skip]==1) {
          regen_all[i+skip] = 0;
        } else if(regen_degen[i+skip]==3) {
          regen_all[i+skip] = 1;
        } else {
          // did not occur in original data but could re-occur probabilistically
          // check for inflation first
          if(infl_this[i]==infl_value) {
            if(infl_gen[i]==1) {
              regen_all[i+skip] = prop_miss[miss_id[i]];
            } else {
              regen_all[i+skip] = beta_proportion_rng(inv_logit(calc_prop[i]),kappa);
            }
          } else {
            regen_all[i+skip] = beta_proportion_rng(inv_logit(calc_prop[i]),kappa);
          }
          
          
        }
        }

        
    } 
  
  }
  
}
