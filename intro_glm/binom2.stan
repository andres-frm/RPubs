
    data{
      int N;
      int N_barrio;
      array[N] int barrio;
      array[N] int gatos_norm;
      array[N] int gatos_pois;
    }
    
    parameters{
      vector[N_barrio] alpha_norm;
      real<lower = 0> sigma;
      vector[N_barrio] alpha_pois;
    }
    
    model{
      target += normal_lpdf(gatos_norm | alpha_norm[barrio], sigma);
      target += normal_lpdf(alpha_norm | 7, 5);
      target += exponential_lpdf(sigma | 1);
    
      target += poisson_lpmf(gatos_pois | exp(alpha_pois[barrio]));
      target += normal_lpdf(alpha_pois | 7, 5);
    }
    