functions {
  // need a custom cdf
  
  real frac_logit_lcdf(real y, real p) {
    
    return log(1 - 2*p) - (log(2) + log(atan(1 - 2*p)));
    
  }
  
  
}

data {
  int n;
  int k; // number of columns
  matrix[n,k] x; 
  vector<lower=0, upper=1>[n] y;
  int run_gen; // whether to use generated quantities
}
parameters {
  vector[k] X_beta;
  real alpha;
}
model {

  X_beta ~ normal(0,5);
  alpha ~ normal(0,5);
  
  // custom "quasi" likelihood
  
  target += y .* log_inv_logit(alpha + x * X_beta) + (1 - y) .* log1m_inv_logit(alpha + x*X_beta);
  
}
generated quantities {
  vector[run_gen==1 ? n: 0] frac_log;
  vector[run_gen==1 ? n: 0] frac_rep;
  
  if(run_gen==1) {

    frac_log = y .* log_inv_logit(alpha + x * X_beta) + (1 - y) .* log1m_inv_logit(alpha + x*X_beta);
    frac_rep = inv_logit(alpha + x * X_beta);
  }
}
