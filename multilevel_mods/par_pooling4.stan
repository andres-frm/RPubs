
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
    
    model{
      vector[N] p;
      mu ~ normal(10, 5);
      alpha ~ normal(0, 1);
      sigma ~ exponential(1);
      psi ~ normal(0, 1);
      tau ~ normal(0, 1);
      beta ~ normal(0.5, 1);
    
      for (i in 1:N) {
        p[i] = mu + alpha[nido[i]]*sigma +
                psi[parche[i]]*tau +
                beta*prop_veg[i];
        p[i] = inv_logit(p[i]);
      }
    
      superv ~ binomial(n, p);
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
      vector[N] nido_P;
      vector[N_parches] parche_P;
      nido_P = mu + alpha*sigma;
      parche_P = psi*tau;
    
      for (i in 1:N) {
        p[i] = mu + alpha[nido[i]]*sigma +
                psi[parche[i]]*tau +
                beta*prop_veg[i];
        p[i] = inv_logit(p[i]);
      }
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(superv[i] | n[i], p[i]);
    }
    