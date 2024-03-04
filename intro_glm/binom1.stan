
    data{
      int N;
      array[N] int obs;
      array[N] int edad;
      array[N] int vege;
    }
    
    parameters {
      real alpha;
      real beta;
      real tau;
    }
    
    model{
      vector[N] p;
      
      for (i in 1:N) {
          p[i] = alpha + beta*edad[i] + tau*vege[i];
          p[i] = inv_logit(p[i]);
      }
    
      obs ~ binomial(1, p);  
    
      alpha ~ normal(0, 1.5);
      beta ~ normal(0.5, 1);
      tau ~ normal(0.5, 1);
    }
    
    