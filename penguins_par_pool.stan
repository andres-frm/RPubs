
    data{
      int N;
      int N_spp_sex;
      int N_island;
      int N_year;
      vector[N] flipper_length_mm;
      vector[N] body_mass_g;
      array[N] int species_sex;
      array[N] int year;
      array[N] int island;
    }
    
    parameters{
      
      matrix[N_spp_sex, N_island] z_island;  
      cholesky_factor_corr[N_spp_sex] Rho_island;
      vector<lower = 0>[N_spp_sex] sigma_island;
      
      matrix[N_spp_sex, N_year] z_year;
      cholesky_factor_corr[N_spp_sex] Rho_year;
      vector<lower = 0>[N_spp_sex] sigma_year;
      
      vector[N_spp_sex] theta;
      vector[N_spp_sex] beta;
      real<lower = 0> sigma;
    }
    
    transformed parameters{
    
      matrix[N_island, N_spp_sex] alpha;
      matrix[N_year, N_spp_sex] tau;
      alpha = (diag_pre_multiply(sigma_island, Rho_island) * z_island)';
      tau = (diag_pre_multiply(sigma_year, Rho_year) * z_year)';
    }
    
    model{
      vector[N] mu;
      sigma_island ~ exponential(1);
      sigma_year ~ exponential(1);
      Rho_island ~ lkj_corr_cholesky(2);
      Rho_year ~ lkj_corr_cholesky(2);
      to_vector(z_island) ~ normal(0, 1);
      to_vector(z_year) ~ normal(0, 1);
      theta ~ normal(200, 50);
      beta ~ normal(0.5, 0.5);
      sigma ~ exponential(1);
     
      for (i in 1:N) {
        mu[i] = theta[species_sex[i]] +  
                alpha[island[i], species_sex[i]] + tau[year[i], species_sex[i]] + 
                beta[species_sex[i]]*body_mass_g[i];
      }
      
      flipper_length_mm ~ normal(mu, sigma);
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu_;
      array[N] real ppcheck;
      matrix[N_spp_sex, N_spp_sex] Rho_isla_corr;
      matrix[N_spp_sex, N_spp_sex] Rho_year_corr;
      
      Rho_isla_corr = multiply_lower_tri_self_transpose(Rho_island);
      Rho_year_corr = multiply_lower_tri_self_transpose(Rho_year);
    
      for (i in 1:N) {
        mu_[i] = theta[species_sex[i]] +  
                alpha[island[i], species_sex[i]] + tau[year[i], species_sex[i]] + 
                beta[species_sex[i]]*body_mass_g[i];
      }
      
      for (i in 1:N) log_lik[i] = normal_lpdf(flipper_length_mm[i] | mu_[i], sigma);
    
      ppcheck = normal_rng(mu_, sigma);
    }