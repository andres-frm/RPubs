
    data{
      int N;
      int N_parche;
      array[N] int n;
      array[N] int superv;
      array[N] int parche;
      array[N] int nido;
    }
    
    parameters{
      vector[N] alpha;
      vector[N_parche] tau;
      real muA;
      real<lower = 0> sigmaA;
      real<lower = 0> sigmaT;
    }
    
    model{
      vector[N] p;
      alpha ~ normal(muA, sigmaA);
      tau ~ normal(0, sigmaT);
      muA ~ normal(5, 2.5);
      sigmaA ~ exponential(1);
      sigmaT ~ exponential(1);
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]] + tau[parche[i]];
        p[i] = inv_logit(p[i]);
      }
      
      superv ~ binomial(n, p);
    
    }
    
    generated quantities {
      vector[N] log_lik;
      vector[N] p;
    
    for (i in 1:N) {
        p[i] = alpha[nido[i]] + tau[parche[i]];
        p[i] = inv_logit(p[i]);
    }
    
    for (i in 1:N) log_lik[i] = binomial_lpmf(superv[i] | n[i], p[i]);
    }
    