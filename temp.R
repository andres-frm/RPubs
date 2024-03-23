paquetes <- c("rethinking", "tidyverse", "magrittr", 'patchwork', 'rstan',
              'cmdstanr', 'loo', 'MASS', 'ellipse', 'palmerpenguins')
sapply(paquetes, library, character.only = T)

options(mc.cores = parallel::detectCores())

source('functions_mod_diagnostics.R')

a <- 3.5 # average morning wait time 
b <- (-1) # average difference afternoon wait time
sigma_a <- 1 # std dev in intercepts
sigma_b <- 0.5 # std dev in slopes
rho <- (-0.7) # correlation between intercepts and slopes

Mu <- c( a , b )

# forma 1 de crear la matriz
cov_ab <- sigma_a*sigma_b*rho
Sigma <- matrix( c(sigma_a^2,cov_ab,cov_ab,sigma_b^2) , ncol=2 )

matrix( c(1,2,3,4) , nrow=2 , ncol=2 , byrow = T)

# forma dos de crear la matriz
sigmas <- c(sigma_a,sigma_b) # standard deviations
Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix

# now matrix multiply to get covariance matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)


# simulacion de cafes 
# 
N_cafes <- 20

set.seed(5)# used to replicate example
vary_effects <- mvrnorm( N_cafes , Mu , Sigma )

a_cafe <- vary_effects[,1]
b_cafe <- vary_effects[,2]

plot(a_cafe , b_cafe , col=rangi2 , 
     xlab="intercepts (a_cafe)" , ylab="slopes (b_cafe)" )
# overlay population distribution
for ( l in c(0.1,0.3,0.5,0.8,0.99) )
  lines(ellipse(Sigma,centre=Mu,level=l),col=col.alpha("black",0.2))

set.seed(22)
N_visits <- 10
afternoon <- rep(0:1,N_visits*N_cafes/2)
cafe_id <- rep( 1:N_cafes , each=N_visits )
mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*afternoon
sigma <- 0.5 # std dev within cafes
wait <- rnorm( N_visits*N_cafes , mu , sigma )
d <- list(cafe=cafe_id , afternoon=afternoon , wait=wait, N = length(cafe_id),
          N_cafes = N_cafes)


par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))

for (i in 1:9) {
  
  R <- rlkjcorr(1e3, 
                K = 2, # n?mero de par?metros (dimension de la matriz)
                eta = i) # forma de la distribuci?n
  
  plot(density(R[, 1, 2]),
       main = '', 
       xlab = 'Correlacion'); title(paste('eta=', i))
}

par(mfrow = c(1, 1))

cat(file = 'multilevel_mods_II/par_corr_eg.stan', 
    "
    data{
      int N;
      int N_cafes;
      array[N] int cafe;
      array[N] int<lower = 0> afternoon;
      array[N] real wait;
    }
    
    parameters{
      matrix[N_cafes, 2] S;
      vector[2] mu_bar;
      corr_matrix[2] R;
      vector<lower = 0>[2] sigma_cafe;
      real<lower = 0> sigma;
    }
    
    transformed parameters{
      vector[N_cafes] alpha;
      vector[N_cafes] beta;
      alpha = S[, 1];
      beta = S[, 2];
    }
    
    model{
      vector[N] mu;
      mu_bar ~ normal(0, 5);
      sigma_cafe ~ exponential(1);
      sigma ~ exponential(1);
      R ~ lkj_corr(2);
      for (i in 1:N_cafes) S[i,:] ~ multi_normal(mu_bar, quad_form_diag(R, sigma_cafe));
    
      for (i in 1:N) {
        mu[i] = alpha[cafe[i]] + beta[cafe[i]]*afternoon[i];
      }
      
      wait ~ normal(mu, sigma);
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu;
      array[N] real ppcheck;
      
      for (i in 1:N) {
        mu[i] = alpha[cafe[i]] + beta[cafe[i]]*afternoon[i];
      }
    
      for (i in 1:N) log_lik[i] = normal_lpdf(wait[i] | mu[i], sigma);
    
      ppcheck = normal_rng(mu, sigma);
      
    }")

file <- paste(getwd(), '/multilevel_mods_II/par_corr_eg.stan', sep = '')
fit_eg <- cmdstan_model(file, compile = T)

m_eg <- 
  fit_eg$sample(
    data = d, 
    iter_sampling = 4e3,
    iter_warmup = 500,
    chains = 3, 
    parallel_chains = 3,
    thin = 3, 
    refresh = 500, 
    seed = 123
  )

out_eg <- m_eg$summary()

mod_diagnostics(m_eg, out_eg)

ppcheck_eg <- m_eg$draws(variables = 'ppcheck', 
                         format = 'matrix')

plot(density(d$wait), ylim = c(0, 0.4), main = '')
for (i in 1:100) lines(density(ppcheck_eg[i, ]), lwd = 0.1)
lines(density(d$wait), col = 'red', lwd = 3)

post_eg <- m_eg$draws(variables = c('alpha', 'beta', 'R', 'sigma', 
                                    'sigma_cafe', 'mu_bar'), 
                      format = 'df')

post_eg <- 
  list(mu_bar = post_eg[, grepl('mu_bar', colnames(post_eg))],
       alpha = post_eg[, grepl('alpha', colnames(post_eg))],
       beta = post_eg[, grepl('beta', colnames(post_eg))],
       Rho = post_eg[, grepl('^R', colnames(post_eg))],
       sigma = post_eg[, grepl('sigma', colnames(post_eg))])

plot(density(post_eg$Rho$`R[2,1]`), main = '', xlim = c(-1, 1))
lines(density(rlkjcorr(1e3, 2, 2)[, 1, 2]), col = 'blue')



a1 <- sapply( 1:N_cafes ,
        function(i) mean(wait[cafe_id==i & afternoon==0]) )
b1 <- sapply( 1:N_cafes ,
              function(i) mean(wait[cafe_id==i & afternoon==1]) ) - a1

a2 <- apply( post_eg$alpha , 2 , mean )
b2 <- apply( post_eg$beta , 2 , mean )

# plot both and connect with lines
plot( a1 , b1 , xlab="intercept" , ylab="slope" ,
      pch=16 , col=rangi2 , ylim=c( min(b1)-0.1 , max(b1)+0.1 ) ,
      xlim=c( min(a1)-0.1 , max(a1)+0.1 ) )
points( a2 , b2 , pch=1 )
for ( i in 1:N_cafes ) lines( c(a1[i],a2[i]) , c(b1[i],b2[i]) )

Mu_est <- c(mean(post_eg$mu_bar$`mu_bar[1]`), mean(post_eg$mu_bar$`mu_bar[2]`))
rho_est <- mean( post_eg$Rho$`R[1,2]` )
sa_est <- mean( post_eg$sigma$`sigma_cafe[1]`)
sb_est <- mean( post_eg$sigma$`sigma_cafe[2]`)
cov_ab <- sa_est*sb_est*rho_est
Sigma_est <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )

for ( l in c(0.1,0.3,0.5,0.8,0.99) )
  lines(ellipse(Sigma_est,centre=Mu_est,level=l),
        col=col.alpha("black",0.2))

wait_morning_1 <- (a1)
wait_afternoon_1 <- (a1 + b1)
wait_morning_2 <- (a2)
wait_afternoon_2 <- (a2 + b2)
# plot both and connect with lines
plot( wait_morning_1 , wait_afternoon_1 , xlab="morning wait" ,
      ylab="afternoon wait" , pch=16 , col=rangi2 ,
      ylim=c( min(wait_afternoon_1)-0.1 , max(wait_afternoon_1)+0.1 ) ,
      xlim=c( min(wait_morning_1)-0.1 , max(wait_morning_1)+0.1 ) )
points( wait_morning_2 , wait_afternoon_2 , pch=1 )
for ( i in 1:N_cafes )
  lines( c(wait_morning_1[i],wait_morning_2[i]) ,
         c(wait_afternoon_1[i],wait_afternoon_2[i]) )
abline( a=0 , b=1 , lty=2 )

v <- mvrnorm( 1e4 , Mu_est , Sigma_est )
v[,2] <- v[,1] + v[,2] # calculate afternoon wait
Sigma_est2 <- cov(v)
Mu_est2 <- Mu_est
Mu_est2[2] <- Mu_est[1]+Mu_est[2]

for ( l in c(0.1,0.3,0.5,0.8,0.99) )
  lines(ellipse(Sigma_est2,centre=Mu_est2,level=l),
        col=col.alpha("black",0.5))

#=====================================#



# ===================================#

data(chimpanzees)
d <- chimpanzees
d$block_id <- d$block
d$treatment <- 1 + d$prosoc_left + 2*d$condition
dat <- list(
  L = d$pulled_left,
  tid = d$treatment,
  actor = d$actor,
  block_id = as.integer(d$block_id),
  N = nrow(d),
  N_actor = length(unique(d$actor)),
  N_block = length(unique(d$block_id)), 
  N_tid = length(unique(d$treatment)))

# version centrada
cat(file = 'multilevel_mods_II/par_corr_eg2.stan', 
    '
    data{
      int N;
      int N_actor;
      int N_block;
      int N_tid;
      array[N] int L;
      array[N] int tid;
      array[N] int actor;
      array[N] int block_id;
    }
    
    parameters{
      array[N_actor] vector[N_tid] alpha;
      array[N_block] vector[N_tid] beta;
      vector[N_tid] tau;
      corr_matrix[N_tid] R_actor;
      corr_matrix[N_tid] R_block;
      vector<lower = 0>[N_tid] sigma_actor;
      vector<lower = 0>[N_tid] sigma_block;
    }
    
    model {
      vector[N] p;
      tau ~ normal(0, 1);
      R_actor ~ lkj_corr(2);
      R_block ~ lkj_corr(2);
      sigma_actor ~ exponential(1);
      sigma_block ~ exponential(1);
      alpha ~ multi_normal(rep_vector(0, 4), quad_form_diag(R_actor, sigma_actor));
      beta ~ multi_normal(rep_vector(0, 4), quad_form_diag(R_block, sigma_block));
    
      for (i in 1:N) {
        p[i] = tau[tid[i]] + alpha[actor[i], tid[i]] + beta[block_id[i], tid[i]];
        p[i] = inv_logit(p[i]);
      }  
    
      L ~ binomial(1, p);
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
      array[N] int ppcheck;
    
      for (i in 1:N) {
        p[i] = tau[tid[i]] + alpha[actor[i], tid[i]] + beta[block_id[i], tid[i]];
        p[i] = inv_logit(p[i]);
      }  
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(L[i] | 1, p[i]);
    
      ppcheck = binomial_rng(1, p);
    }
    
    ')

file <- paste(getwd(), '/multilevel_mods_II/par_corr_eg2.stan', sep = '')

fit_eg2 <- cmdstan_model(file, compile = T)

m_eg <- 
  fit_eg2$sample(
    data = dat, 
    iter_sampling = 4e3,
    iter_warmup = 500,
    chains = 3, 
    parallel_chains = 3,
    thin = 3, 
    refresh = 500, 
    seed = 123
  )

out_eg2 <- m_eg$summary()

mod_diagnostics(m_eg, out_eg2)


# version no centrada 
cat(file = 'multilevel_mods_II/par_corr_eg3.stan', 
    "
    data{
      int N;
      int N_actor;
      int N_block;
      int N_tid;
      array[N] int L;
      array[N] int tid;
      array[N] int actor;
      array[N] int block_id;
    }
    
    parameters{
      matrix[N_tid, N_actor] z_alpha;
      matrix[N_tid, N_block] z_beta;
      vector[N_tid] tau;
      cholesky_factor_corr[N_tid] R_actor;
      cholesky_factor_corr[N_tid] R_block;
      vector<lower = 0>[N_tid] sigma_actor;
      vector<lower = 0>[N_tid] sigma_block;
    }
    
    transformed parameters{
      matrix[N_actor, N_tid] alpha;
      matrix[N_block, N_tid] beta;
      alpha = (diag_pre_multiply(sigma_actor, R_actor) * z_alpha)';
      beta = (diag_pre_multiply(sigma_block, R_block) * z_beta)';
    }
    
    model {
      vector[N] p;
      tau ~ normal(0, 1);
      R_actor ~ lkj_corr_cholesky(2);
      R_block ~ lkj_corr_cholesky(2);
      sigma_actor ~ exponential(1);
      sigma_block ~ exponential(1);
      to_vector(z_alpha) ~ normal(0, 1);
      to_vector(z_beta) ~ normal(0, 1);
    
      for (i in 1:N) {
        p[i] = tau[tid[i]] + alpha[actor[i], tid[i]] + beta[block_id[i], tid[i]];
        p[i] = inv_logit(p[i]);
      }  
    
      L ~ binomial(1, p);
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
      array[N] int ppcheck;
      matrix[N_tid, N_tid] R_actor1;
      matrix[N_tid, N_tid] R_block1;
    
      R_actor1 = multiply_lower_tri_self_transpose(R_actor);
      R_actor1 = multiply_lower_tri_self_transpose(R_block);
    
      for (i in 1:N) {
        p[i] = tau[tid[i]] + alpha[actor[i], tid[i]] + beta[block_id[i], tid[i]];
        p[i] = inv_logit(p[i]);
      }  
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(L[i] | 1, p[i]);
    
      ppcheck = binomial_rng(1, p);
    }")

file <- paste(getwd(), '/multilevel_mods_II/par_corr_eg3.stan', sep = '')

fit_eg3 <- cmdstan_model(file, compile = T)

m_eg3 <- 
  fit_eg3$sample(
    data = dat, 
    iter_sampling = 4e3,
    iter_warmup = 500,
    chains = 3, 
    parallel_chains = 3,
    thin = 3, 
    refresh = 500, 
    seed = 123
  )

out_eg3 <- m_eg3$summary()
mod_diagnostics(m_eg3, out_eg3)




# ============================== #


data(bangladesh)
d <- bangladesh

d$district_id <- as.integer(as.factor(d$district))

colnames(d)

dat <- 
  list(
    anticonc = d$use.contraception,
    urbano = as.integer(ifelse(d$urban == 1, 1, 2)),
    distrito = d$district_id,
    edad = d$age.centered,
    N = nrow(d),
    N_distrito = length(unique(d$district_id)),
    N_urbano = 2
  )


cat(file = 'multilevel_mods_II/corr_par1.stan', 
    '
    data{
      int N;
      int N_distrito;
      int N_urbano;
      array[N] int anticonc;
      array[N] int distrito;
      array[N] int urbano;
      //vector[N] edad;
    }
    
    parameters{
      matrix[N_distrito, N_urbano] v;
      vector[2] alpha_bar;
      corr_matrix[2] Rho;
      vector<lower = 0>[2] sigma;
    }
    
    transformed parameters {
      vector[N_distrito] alpha;
      vector[N_distrito] beta;
      alpha = v[, 1];
      beta = v[, 2];
    }
    
    model{
      vector[N] p;
      alpha_bar ~ normal(0, 1);
      Rho ~ lkj_corr(4);
      sigma ~ exponential(1);
      
      for (i in 1:N_distrito) v[i,:] ~ multi_normal(alpha_bar, quad_form_diag(Rho, sigma));
      
      for (i in 1:N) {
        p[i] = alpha[distrito[i]] + beta[distrito[i]]*urbano[i]; 
        p[i] = inv_logit(p[i]);
      }
    
      anticonc ~ binomial(1, p);
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] prob;
      array[N] int sim;

      for (i in 1:N) {
        prob[i] = alpha[distrito[i]] + beta[distrito[i]]*urbano[i];
        prob[i] = inv_logit(prob[i]);
      }

      for (i in 1:N) log_lik[i] = binomial_lpmf(anticonc[i] | 1, prob[i]);
    
      sim = binomial_rng(1, prob);
}
    ')


file <- paste(getwd(), '/multilevel_mods_II/corr_par1.stan', sep = '')

fit_corr_par1 <- cmdstan_model(file, compile = T)

m1 <- 
  fit_corr_par1$sample(
    data = dat,
    iter_sampling = 4e3,
    iter_warmup = 500,
    parallel_chains = 3, 
    chains = 3,
    thin = 3,
    seed = 123,
    refresh = 500
  )


(out_m1 <- m1$summary())


par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
for (i in 1:9) trace_plot(m1, out_m1$variable[i], 3)
par(mfrow = c(1, 1))

mod_diagnostics(m1, out_m1)


cat(file = "multilevel_mods_II/corr_par2.stan", 
    "
    data{
      int N;
      int N_distrito;
      int N_urbano;
      array[N] int anticonc;
      array[N] int distrito;
      array[N] int urbano;
      //vector[N] edad;
    }
    
    parameters{
      matrix[N_urbano, N_distrito] Z;
      vector[2] alpha_bar;
      cholesky_factor_corr[2] L_Rho;
      vector<lower = 0>[2] sigma;
    }
    
    transformed parameters {
      vector[N_distrito] alpha;
      vector[N_distrito] beta;
      matrix[N_distrito, N_urbano] v;
      v = (diag_pre_multiply(sigma, L_Rho) * Z)';
      alpha = alpha_bar[1] + v[, 1];
      beta = alpha_bar[2] + v[, 2];
    }
    
    model{
      vector[N] p;
      alpha_bar ~ normal(0, 1);
      sigma ~ exponential(1);
      L_Rho ~ lkj_corr_cholesky(4);
      to_vector(Z) ~ normal(0, 1);
      
      for (i in 1:N) {
        p[i] = alpha[distrito[i]] + beta[distrito[i]]*urbano[i];
        p[i] = inv_logit(p[i]);
      }
      
      anticonc ~ binomial(1, p);
    }
    
    generated quantities{
      vector[N] log_lik;
      array[N] int sim;
      vector[N] p;
      matrix[2, 2] Rho;
      Rho = multiply_lower_tri_self_transpose(L_Rho);
      
      for (i in 1:N) {
        p[i] = alpha[distrito[i]] + beta[distrito[i]]*urbano[i];
        p[i] = inv_logit(p[i]);
      }
      
      for (i in 1:N) log_lik[i] = binomial_lpmf(anticonc[i] | 1, p[i]);
    
      sim = binomial_rng(1, p);
    }
    ")

file <- paste(getwd(), '/multilevel_mods_II/corr_par2.stan', sep = '')

fit_m2 <- cmdstan_model(file, compile = T)

m2 <- 
  fit_m2$sample(
    data = dat,
    chains = 3,
    parallel_chains = 3, 
    iter_sampling = 1e4,
    iter_warmup = 500,
    thin = 3,
    refresh = 500, 
    seed = 123
  )

out_m2 <- m2$summary()

par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
for (i in 1:9) trace_plot(m2, out_m2$variable[i], 3)
par(mfrow = c(1, 1))

mod_diagnostics(m2, out_m2)


#===============
#

paquetes <- c("tidyverse", "magrittr", 'patchwork', 'rstan',
              'cmdstanr', 'loo', 'palmerpenguins')
sapply(paquetes, library, character.only = T)

options(mc.cores = parallel::detectCores())

source('functions_mod_diagnostics.R')

d <- na.omit(penguins)

d <- d |> dplyr::select(species, year, island, sex, flipper_length_mm, body_mass_g)

d$year <- as.factor(d$year)

dat <- lapply(d, function(x) if(!is.double(x)) as.integer(x) else(x))

dat$N <- length(dat$species)
dat$N_spp <- max(dat$species)
dat$N_sex <- 2
dat$N_island <- max(dat$island)
dat$N_year <- max(dat$year)

cat(file = 'penguins_par_pool.stan', 
    '
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
    }')

file <- paste(getwd(), '/penguins_par_pool.stan', sep = '')

fit_m1 <- cmdstan_model(file, compile = T)

m1 <- 
  fit_m1$sample(
    data = dat,
    iter_sampling = 4e3,
    iter_warmup = 500,
    chains = 3,
    parallel_chains = 3,
    thin = 3,
    refresh = 500,
    seed = 123
  )

out_m1 <- m1$summary()

mod_diagnostics(m1, out_m1)

ppcheck_m1 <- m1$draws(variables = 'ppcheck', format = 'matrix')

plot(NULL, main = 'Partial pooling (intercepts)', 
     xlim = c(160, 245), ylim = c(0, 0.04), xlab = 'mm', 
     ylab = 'Density')
for (i in 1:100) lines(density(ppcheck_m1[i, ]), lwd = 0.1)
lines(density(dat$flipper_length_mm), col = 'red')

post_m1 <- m1$draws(variables = c('spp', 'I', 'Y', 'beta', 'sigma'), 
                    format = 'df')

post_m1 <- 
  list(spp = post_m1[, grepl('spp', colnames(post_m1))],
       I = post_m1[, grepl('^I', colnames(post_m1))],
       Y = post_m1[, grepl('^Y', colnames(post_m1))], 
       beta = post_m1[, grepl('beta', colnames(post_m1))])

x_seq <- seq(min(dat$body_mass_g), max(dat$body_mass_g), length.out = 100)

sex <- c('H', 'M')
vars <- list(x = 1:3, x1 = 4:6)

est_cor_H <- lapply(1:2, FUN = 
                      function(i) {
                        
                        lapply(1:100, FUN = 
                                 function(x) {
                                   
                                   df <- 
                                     sapply(vars[[i]], FUN = 
                                              function(z) {
                                                
                                                post_m1$spp[, z] + # hembras
                                                  apply(post_m1$I, 1, mean) +
                                                  apply(post_m1$Y, 1, mean) +
                                                  post_m1$beta[, z, drop = T]*x_seq[x]
                                                
                                              })
                                   df <- as_tibble(df)
                                   colnames(df) <- levels(d$species)
                                   df <- gather(df)
                                   df$sex <- sex[i]
                                   df$x <- x_seq[x]
                                   df$dat <- x
                                   df
                                   
                                 })
                        
                      })

for (i in 1:2) est_cor_H[[i]] <- do.call('rbind', est_cor_H[[i]])

length(est_cor_H)

est_cor_H <- 
  lapply(est_cor_H, FUN = 
           function(z) {
             
             z$dat <- as.factor(z$dat)
             z |> 
               group_by(key, dat) |> 
               transmute(mu = mean(value), 
                         li = quantile(value, 0.025),
                         ls = quantile(value, 0.975), 
                         x = x, 
                         sex = sex) |> 
               unique()
             
           })

est_cor_H <- do.call('rbind', est_cor_H)

est_cor_H$sex <- ifelse(est_cor_H$sex == 'H', 'female', 'male')
colnames(est_cor_H)[1] <- 'species'

ggplot() +
  geom_point(data = d, aes(body_mass_g, flipper_length_mm, color = sex)) +
  geom_line(data = est_cor_H, 
            aes(x, mu, color = sex)) +
  geom_ribbon(data = est_cor_H, 
              aes(x, mu, fill = sex, 
                  ymin = li, ymax = ls), alpha = 0.3) +
  facet_wrap(~ species)

d1_ <- unique(d[, c("species", "year", "island", "sex")])

d1 <- apply(d1_, 2, FUN = 
              function(x) {
                
                x <- factor(x, labels = 1:length(unique(x)))
                as.integer(x)
                
              })

d1 <- as_tibble(d1)

sigmas <- m1$draws(format = 'df')
sigmas <- sigmas[, grepl('sigma', colnames(sigmas))]

est_islas_year <- lapply(1:nrow(d1), FUN = 
                           function(x) {
                             
                             se <- d1$sex[[x]]
                             s_b <- ifelse(se == 1, d1$species[[x]] + 0, 
                                           d1$species[[x]] + 3)
                             is <- d1$island[[x]]
                             y <- d1$year[[x]]
                             
                             
                             p <- with(post_m1, 
                                       {
                                         spp[, s_b, drop = T] + 
                                           I[, is, drop = T] +
                                           Y[, y, drop = T] +
                                           beta[, s_b, drop = T]*mean(d$body_mass_g)
                                         
                                       })
                             
                             p <- rnorm(1e3, p, sigmas$sigma)
                             
                             cbind(tibble(val = p), d1_[x, ])
                             
                           })

est_islas_year <- as_tibble(do.call('rbind', est_islas_year))

est_islas_year |> 
  ggplot() +
  geom_boxplot(aes(year, val, fill = species)) +
  facet_wrap( ~ island + sex) +
  labs(x = NULL, y = 'Longitud aleta (mm)') +
  theme_bw() +
  theme(legend.position = 'top')

sigmas <- sigmas[, c(1, 3)]
colnames(sigmas) <- c('sigma\n(island)', 'sigma\n(year)')
gather(sigmas) |> 
  ggplot() +
  geom_density(aes(value, fill = key), alpha = 0.5) +
  theme(legend.position = 'top')


# ==============

library(lterdatasampler)

data(and_vertebrates)

?and_vertebrates

d <- and_vertebrates

d %$% aggregate(length_1_mm ~ species + year + 
                  sitecode + section + unittype, FUN = length)


