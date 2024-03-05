
    data{
      int N;
      int N_parches;
      array[N] int n;
      array[N] int superv;
      array[N] int parche;
      vector[N] prop_veg;
      array[N] int nido;
    }
    
    parameters{
      vector[N] alpha;
      vector[N_parches] psi;
      real mu;
      real beta;
      real<lower = 0> sigma;
      real<lower = 0> tau;
    }
    
    model {
      vector[N] p;
      
      for (i in 1:N) {
        p[i] = alpha[nido[i]] + psi[parche[i]] + beta*prop_veg[i];
        p[i] = inv_logit(p[i]);
      }
    
      alpha ~ normal(mu, sigma);
      psi ~ normal(0, tau);
      mu ~ normal(10, 5);
      beta ~ normal(0.5, 1);
      sigma ~ exponential(1);
      tau ~ exponential(1);
    
      superv ~ binomial(n, p);
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
    
      for (i in 1:N){
        p[i] = alpha[nido[i]] + psi[parche[i]] + beta*prop_veg[i];
        p[i] = inv_logit(p[i]);
      }
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(superv[i] | n[i], p[i]);
    }
    