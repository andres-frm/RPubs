
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
      vector[N] z_alpha;
      vector[N_parches] z_psi;
      real mu1;
      real mu2;
      real beta;
      real<lower = 0> sigma;
      real<lower = 0> tau;
    }
    
    model{
      vector[N] p;
      vector[N] alpha;
      vector[N_parches] psi;
      mu1 ~ normal(10, 5);
      mu1 ~ normal(0, 5);
      z_alpha ~ normal(0, 1);
      z_psi ~ normal(0, 1);
      sigma ~ exponential(1);
      tau ~ exponential(1);
      beta ~ normal(0.5, 1);
      alpha = mu1 + z_alpha*sigma;
      psi = mu2 + z_psi*tau;
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]] + psi[parche[i]] + beta*prop_veg[i];
        p[i] = inv_logit(p[i]);
      }
    
      superv ~ binomial(n, p);
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
      vector[N] alpha;
      vector[N_parches] psi;
      alpha = mu1 + z_alpha*sigma;
      psi = mu2 + z_psi*tau;
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]] + psi[parche[i]] + beta*prop_veg[i];
        p[i] = inv_logit(p[i]);
      }
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(superv[i] | n[i], p[i]);
    }
    