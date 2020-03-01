data {
  int n;
  int k; // number of columns
  matrix[n,k] x; 
  vector<lower=0, upper=1>[n] y;
  int run_gen; // whether to use generated quantities
}
transformed data {
  int<lower=0, upper=1> is_discrete[n];
  int<lower=-1, upper=1> y_discrete[n];

  // create indicator for whether y is discrete 
  // and an integer value to pass to bernoulli_lpmf for discrete y
  for (i in 1:n) {
    if (y[i] == 0) {
      is_discrete[i] = 1;
      y_discrete[i] = 0;
    } else if (y[i] == 1) {
      is_discrete[i] = 1;
      y_discrete[i] = 1;
    } else {
      is_discrete[i] = 0;
      // hack to ensure that throws error if passed to bernoulli_lpmf
      y_discrete[i] = -1;
    }
  }
}
parameters {
  vector[k] coef_a;
  vector[k] coef_g;
  vector[k] coef_m;
  //vector[k] coef_p;
  real<lower=0> phi;
  vector[3] alpha;
}
transformed parameters {
  vector[n] psi;
  vector[n] gamma;
  vector[n] mu;
  //vector<lower=0>[n] phi;
  
  psi = alpha[1] + x*coef_a;
  gamma = alpha[2] + x*coef_g;
  mu = inv_logit(alpha[3] + x*coef_m);
  //phi = exp(alpha[4] + x*coef_p);
  
}
model {
  coef_a ~ normal(0, 3);
  coef_g ~ normal(0, 3);
  coef_m ~ normal(0, 3);
  //coef_p ~ normal(0, 3);
  phi ~ exponential(.1);
  alpha ~ normal(0,3);
  
  is_discrete ~ bernoulli_logit(psi);
  for (i in 1:n) {
    if (is_discrete[i] == 1) {
      y_discrete[i] ~ bernoulli_logit(gamma[i]);
    } else {
      y[i] ~ beta_proportion(mu[i], phi);
    }
  }
  
}
generated quantities {
  vector[run_gen==1 ? n: 0] zoib_log;
  vector[run_gen==1 ? n: 0] zoib_regen;
  vector[run_gen==1 ? n: 0] is_discrete_regen;
  
  if(run_gen==1) {
    for (i in 1:n) {
      real psit = inv_logit(psi[i]);
      real gammat = inv_logit(gamma[i]);
      if (y[i] == 0) {
        zoib_log[i] = log(psit) + log1m(gammat);
      } else if (y[i] == 1) {
        zoib_log[i] = log(psit) + log(gammat);
      } else {
        zoib_log[i] = log1m(psit) + beta_proportion_lpdf(y[i] | mu[i], phi);
      }
      
      is_discrete_regen[i] = bernoulli_rng(psit);
      
      if(is_discrete_regen[i]==0) {
        zoib_regen[i] = beta_proportion_rng(mu[i], phi);
      } else {
        zoib_regen[i] = bernoulli_rng(gammat);
      }
    }
  
  }
}
