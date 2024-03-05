
    data{
      int N;
      array[N] int surv;
      array[N] int huevos;
    }
    
    parameters{
      real alpha;
    }
    
    model{
      real p;
      alpha ~ normal(5, 2.5);
      p = alpha;
      p = inv_logit(p);
      surv ~ binomial(huevos, p);
    }
    
    generated quantities{
      vector[N] log_lik;
      real p;
      p = alpha;
      p = inv_logit(p);
      for (i in 1:N) log_lik[i] = binomial_lpmf(surv[i] | huevos[i], p);
    }
    