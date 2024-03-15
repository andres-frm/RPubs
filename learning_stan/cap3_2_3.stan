
    data {
      int N;
      array[N] int trofeos; 
      array[N] int year;
    }
    
    parameters {
      real alpha;
      real beta;
     
    }
    
    model {
    
      vector[N] lambda;  
    
      for (i in 1:N) {
        
        // M3. versión larga con el muestreo explícito (M 3)
        lambda[i] = alpha + beta*year[i];
        lambda[i] = exp(lambda[i]);
        
      }
    
      trofeos ~ poisson(lambda);  
    
      alpha ~ normal(10, 5); // previa M2
      beta ~ normal(0.5, 1);  // previa M2
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu;
      
      for (i in 1:N) {
      
        mu[i] = alpha + beta*year[i];
        mu[i] = exp(mu[i]);
        
      }
    
      for (i in 1:N) log_lik[i] = poisson_lpmf(trofeos[i] | 
                                               alpha + beta*year[i]);
    }
    