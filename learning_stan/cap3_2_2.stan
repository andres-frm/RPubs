
    data {
      int N;
      array[N] int trofeos; 
      vector[N] year;
    }
    
    parameters {
      real alpha;
      real beta;
     
    }
    
    model {
    
      vector[N] mu;  
    
      for (i in 1:N) {
        
        // M2. versión abreviada con el muestreo explícito (M 2)
        trofeos[1:i] ~ poisson(exp(alpha + beta*year[i]));
        
      }
      
      alpha ~ normal(10, 5); // previa M2
      beta ~ normal(0.5, 1);  // previa M2
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu;
      
      for (i in 1:N) {
      
        mu[1:i] = alpha + beta*year[1:i];
        mu[1:i] = exp(mu[1:i]);
        
      }
    
      for (i in 1:N) log_lik[i] = poisson_lpmf(trofeos[1:i] | 
                                               alpha + beta*year[1:i]);
    }
    