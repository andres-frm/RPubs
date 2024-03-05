
    data{
      int N;
      array[N] int surv;
      array[N] int huevos;
      array[N] int nido;
    }
    
    parameters{
      vector[N] alpha;
    }
    
    model {
      vector[N] p;
      alpha ~ normal(5, 25); // previa para cada nido
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]]; // estimación estratificada por nido
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
    