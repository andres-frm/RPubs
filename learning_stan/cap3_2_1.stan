
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
      
      for (i in 1:N) {
        
        // M1. log p de la funcion de distribución de probabilidas
        target += poisson_lpmf(trofeos[i] | 
                               exp(alpha + beta*year[i]));
        
      }
      
      target += normal_lpdf(alpha | 10, 5); // previa M1
      target += normal_lpdf(beta | 0.5, 1);  // previa M1
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu;
      
      for (i in 1:N) {
      
        mu[i] = alpha + beta*year[i];
        mu[i] = exp(mu[i]);
        
      }
    
      for (i in 1:N) log_lik[i] = poisson_lpmf(trofeos[i] | 
                                              exp(alpha + beta*year[i]));
    }
    