
    data{
      int N;
      array[N] int surv;
      array[N] int huevos;
      array[N] int nido;
    }
    
    parameters{
      vector[N] alpha;
      real mu;
      real<lower = 0> sigma;
    }
    
    model{
      vector[N] p;
      alpha ~ normal(mu, sigma); // previa para la población de nidos (previa adaptativa)
      mu ~ normal(5, 2.5); // Hiper previa (previa de la previa)
      sigma ~ exponential(1); // Hiper previa (previa de la previa)
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]]; // modelo estratificado por nido
        p[i] = inv_logit(p[i]);
      }
    
    surv ~ binomial(huevos, p);
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
      
      for (i in 1:N) {
        p[i] = alpha[nido[i]];
        p[i] = inv_logit(p[i]);
      }
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(surv[i] | huevos[i], p[i]);
    }
    
    