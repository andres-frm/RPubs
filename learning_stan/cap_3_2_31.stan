
    data {
      int N;
      int N_x;
      array[N] int trofeos; 
      array[N_x] int year_x;
      array[N] int year;
    }
    
    parameters {
      real alpha;
      real beta;
     
    }
    
    transformed parameters { // se usa para guardar muestras del lambda o para generar variables
                            // que son usadas tanto en el bloque de modelo como en el de  
                            // generated quantities
      
    vector[N] lambda;  
    
      for (i in 1:N) {
        
        // M3. versión larga con el muestreo explícito (M 3)
        lambda[i] = alpha + beta*year[i];
        lambda[i] = exp(lambda[i]); 
        
      }
    
    }
    
    model {

      trofeos ~ poisson(lambda);  
    
      alpha ~ normal(10, 5); // previa M2
      beta ~ normal(0.5, 1);  // previa M2
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu;
      vector[N_x] theta;
      array[N_x] real pred;
      array[N] int ppcheck;
      
      for (i in 1:N) {
      
        mu[i] = alpha + beta*year[i];
        mu[i] = exp(mu[i]);
        
      }
    
      for (i in 1:N) log_lik[i] = poisson_lpmf(trofeos[i] | 
                                               alpha + beta*year[i]);
      
      ppcheck = poisson_rng(mu);
    
      for (i in 1:N_x) {
      
        theta[i] = alpha + beta*year_x[i];
        theta[i] = exp(theta[i]);
        
      }
    
      pred = poisson_rng(theta);
    
    }