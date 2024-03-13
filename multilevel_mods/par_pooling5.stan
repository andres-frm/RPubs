
    data{
      int N;
      int N_distrito;
      int N_urbano;
      int N_anticonc;
      array[N] int anticonc;
      array[N] int distrito;
      array[N] int urbano;
      //vector[N] edad;
    }
    
    parameters{
      vector[N_distrito] z_alpha;
      vector[N_distrito] z_beta;
      real beta1;
      real alpha1;
      real<lower = 0> sigma;
      real<lower = 0> phi;
    }
    
    model{
      vector[N] p;
      vector[N_distrito] alpha;
      vector[N_distrito] beta;
      alpha1 ~ normal(0, 1);
      z_alpha ~ normal(0, 1);
      sigma ~ exponential(1);
      beta1 ~ normal(0, 1);
      z_beta ~ normal(0, 1);
      phi ~ exponential(1);
      alpha = alpha1 + z_alpha*sigma;
      beta = beta1 + z_beta*phi;
    
      for (i in 1:N) {
        p[i] = alpha[distrito[i]] + beta[distrito[i]]*urbano[i]; 
        p[i] = inv_logit(p[i]);
      }
    
      anticonc ~ binomial(1, p);
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
      vector[N_distrito] alpha;
      vector[N_distrito] beta;
      alpha = alpha1 + z_alpha*sigma;
      beta = beta1 + z_beta*phi;
    
      for (i in 1:N) {
        p[i] = alpha[distrito[i]] + beta[distrito[i]]*urbano[i]; 
        p[i] = inv_logit(p[i]);
      }
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(anticonc[i] | 1, p[i]);
    }
    