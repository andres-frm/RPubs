
    data{
      int N;
      int N_distrito;
      int N_urbano;
      int N_anticonc;
      array[N] int anticonc;
      array[N] int distrito;
      array[N] int urbano;
      vector[N] edad;
    }
    
    parameters{
      vector[N_urbano] alpha;
      vector[N_distrito] theta;
      real beta;
      real mu;
      real<lower = 0> sigma;
      real<lower = 0> phi;
    }
    
    model{
      vector[N] p;
      alpha ~ normal(mu, sigma);
      mu ~ normal(0, 3);
      sigma ~ exponential(1);
      theta ~ normal(0, phi);
      phi ~ exponential(1);
      beta ~ normal(0.15, 2);
    
      for (i in 1:N) {
        p[i] = alpha[urbano[i]] + theta[distrito[i]] + beta*edad[i];
        p[i] = inv_logit(p[i]);
      }
    
      anticonc ~ binomial(1, p);
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
    
      for (i in 1:N) {
        p[i] = alpha[urbano[i]] + theta[distrito[i]] + beta*edad[i];
        p[i] = inv_logit(p[i]);
      }
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(anticonc[i] | 1, p[i]);
    }
    