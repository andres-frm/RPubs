
    data{
      int N;
      int N_spp;
      int N_sex;
      int N_island;
      int N_year;
      vector[N] flipper_length_mm;
      vector[N] body_mass_g;
      array[N] int species;
      array[N] int year;
      array[N] int island;
      array[N] int sex;
    }
    
    parameters{
      
      vector[N_island] z_I;
      real mu_I;
      real<lower = 0> sigma1;
      matrix[N_spp, N_sex] z_spp;
      real mu_spp;
      real<lower = 0> sigma2;
      vector[N_year] z_Y;
      real mu_Y;
      real<lower = 0> sigma3;
      matrix[N_spp, N_sex] beta;
      real<lower = 0> sigma;
      
    }
    
    model{
      vector[N] mu;
      vector[N_island] I;
      matrix[N_spp, N_sex] spp;
      vector[N_year] Y;
      mu_I ~ normal(0, 1);
      z_I ~ normal(0, 0.5);
      mu_spp ~ normal(200, 50);
      to_vector(z_spp) ~ normal(0, 1);
      mu_Y ~ normal(0, 1);
      z_Y ~ normal(0, 0.5);
      to_vector(beta) ~ normal(1, 1);
      sigma ~ exponential(1);
      sigma1 ~ exponential(1);
      sigma2 ~ exponential(1);
      sigma3 ~ exponential(1);
      I = mu_I + z_I*sigma1;
      spp = mu_spp + z_spp*sigma2;
      Y = mu_Y + z_Y*sigma3;
    
      for (i in 1:N) {
        mu[i] = spp[species[i], sex[i]] +  
                I[island[i]] + Y[year[i]] + 
                beta[species[i], sex[i]]*body_mass_g[i];
      }
      
      flipper_length_mm ~ normal(mu, sigma);
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu;
      vector[N_island] I;
      matrix[N_spp, N_sex] spp;
      vector[N_year] Y;
      I = mu_I + z_I*sigma1;
      spp = mu_spp + z_spp*sigma2;
      Y = mu_Y + z_Y*sigma3;
      array[N] real ppcheck;
    
      for (i in 1:N) {
        mu[i] = spp[species[i], sex[i]] +  
                I[island[i]] + Y[year[i]] + 
                beta[species[i], sex[i]]*body_mass_g[i];
      }
      
      for (i in 1:N) log_lik[i] = normal_lpdf(flipper_length_mm[i] | mu[i], sigma);
    
      ppcheck = normal_rng(mu, sigma);
    }