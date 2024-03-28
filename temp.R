paquetes <- c("rethinking", "tidyverse", "magrittr", 'patchwork', 'rstan',
              'cmdstanr', 'loo', 'MASS', 'ellipse', 'palmerpenguins', 
              'lterdatasampler')

sapply(paquetes, library, character.only = T)

options(mc.cores = parallel::detectCores())

source('functions_mod_diagnostics.R')

alpha <- 0.05
tau <- 0.1
beta <- 1.2
sigma_alpha <- 0.8
sigma_tau <- 1
sigma_beta <- 0.9
rho_alpha_beta <- 0.005
rho_tau_beta <- 0.8

mu1 <- c(alpha, beta)
mu2 <- c(tau, beta)

cov_alpha_beta <- sigma_alpha*sigma_beta*rho_alpha_beta
cov_tau_beta <- sigma_tau*sigma_beta*rho_tau_beta

sig_alpha <- c(sigma_alpha, sigma_beta)
mat_VC_alpha <- matrix(c(1, cov_alpha_beta, cov_alpha_beta, 1), ncol = 2)
mat_VC_alpha <- diag(sig_alpha) %*% mat_VC_alpha %*% diag(sig_alpha)

sig_tau <- c(sigma_tau, sigma_beta)
mat_VC_tau <- matrix(c(1, cov_tau_beta, cov_tau_beta, 1), ncol = 2)
mat_VC_tau <- diag(sig_tau) %*% mat_VC_tau %*% diag(sig_tau)

par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
set.seed(123)
nidos <- rmvnorm(1000, mu1, sigma = mat_VC_alpha)
plot(nidos[, 1], nidos[, 2], xlab = expression(alpha), 
     ylab = expression(beta))
for (i in seq(0.1, 0.9, by = 0.1)) {
  l <- ellipse(mat_VC_alpha, centre = mu1, level = i)
  lines(l[, 1], l[, 2], col = 'tomato', lwd = 2)
}

set.seed(123)
parche <- rmvnorm(1000, mu2, sigma = mat_VC_tau)
plot(parche[, 1], parche[, 2], xlab = expression(tau), 
     ylab = expression(beta))
for (i in seq(0.1, 0.9, by = 0.1)) {
  l <- ellipse(mat_VC_tau, centre = mu2, level = i)
  lines(l[, 1], l[, 2], col = 'tomato', lwd = 2)
}
par(mfrow = c(1, 1))

plot(density(rmvnorm(1e3, mu1, sigma = mat_VC_alpha)[, 2]), 
     main = '', xlab = expression(beta))
lines(density(rmvnorm(1e3, mu2, sigma = mat_VC_tau)[, 2]))

N <- 1000
set.seed(123)
nidada <- rpois(N, 10)
nido_id <- 1:N
parches <- rep(1:20, each = 10 , length.out = N)
bosque <- rep(0:1, each = 5, length.out = N)
alpha_nido <- nidos[, 1]
tau_parche <- parche[, 1]
beta_bosque <- parche[, 2]
mu <- alpha_nido[nido_id] + tau_parche[parches] + beta_bosque[parches]*bosque
mu <- inv_logit(mu)
set.seed(123)
huevos <- rbinom(N, size = nidada, prob = mu)

dat <- 
  list(
    N = N, 
    N_parche = 20,
    N_bosque = 2,
    N_nidos = 1000,
    nido = nido_id,
    nidada = nidada,
    bosque = bosque, 
    huevos = huevos,
    parche = parches
  )

cat(file = 'multilevel_mods_II/eg_cov_0.stan', 
    '
    data{
      int N;
      int N_parche;
      int N_bosque;
      int N_nidos;
      array[N] int nido;
      array[N] int nidada;
      array[N] int bosque;
      array[N] int huevos;
      array[N] int parche;
    }
    
    parameters {
      vector[N_nidos] z_alpha;
      real mu_alpha;
      real<lower = 0> sigma_alpha;
      vector[N_parche] z_tau;
      real mu_tau;
      real<lower = 0> sigma_tau;
      vector[N_parche] z_beta;
      real mu_beta;
      real<lower = 0> sigma_beta;
    }
    
    transformed parameters {
      vector[N_nidos] alpha;
      vector[N_parche] tau;
      vector[N_parche] beta;
      alpha = mu_alpha + z_alpha*sigma_alpha;
      tau = mu_tau + z_tau*sigma_tau;
      beta = mu_beta + z_beta*sigma_beta;
    }
    
    model{
      vector[N] p;
      z_alpha ~ normal(0, 1);
      mu_alpha ~ normal(0, 1);
      sigma_alpha ~ exponential(1);
      z_tau ~ normal(0, 1);
      mu_tau ~ normal(5, 2);
      sigma_tau ~ exponential(1);
      z_beta ~ normal(0, 1);
      mu_beta ~ normal(0, 1);
      sigma_beta ~ exponential(1);
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]] + tau[parche[i]] + beta[parche[i]] * bosque[i];
        p[i] = inv_logit(p[i]);
      }
    
      huevos ~ binomial(nidada, p);
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
      array[N] int ppcheck;
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]] + tau[parche[i]] + beta[parche[i]] * bosque[i];
        p[i] = inv_logit(p[i]);
      }
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(huevos[i] | nidada[i], p[i]);
    
      ppcheck = binomial_rng(nidada, p);
    }
    ')

file <- paste(getwd(), '/multilevel_mods_II/eg_cov_0.stan', sep = '')
fit_m_noCOR <- cmdstan_model(file, compile = T)

m_noCORR <- 
  fit_m_noCOR$sample(
    data = dat,
    iter_sampling = 4e3, 
    iter_warmup = 500,
    chains = 3,
    parallel_chains = 3, 
    thin = 3,
    refresh = 500
  )

m_noCORR$save_object('multilevel_mods_II/eg_cov_0.rds')

m_noCORR <- readRDS('multilevel_mods_II/eg_cov_0.rds')

out_noCOR <- m_noCORR$summary()

mod_diagnostics(m_noCORR, out_noCOR)

ppcheck <- as.matrix(m_noCORR$draws(variables = 'ppcheck', format = 'df'))

plot(density(ppcheck[1, ]), main = '', lwd = 0.1, ylim = c(0, 0.14))
for (i in 1:500) lines(density(ppcheck[i, ]), lwd = 0.1)
lines(density(dat$huevos), col = 'red', lwd = 2)

post <- m_noCORR$draws(variable = c('alpha', 'tau', 'beta'), format = 'df')
post_noCOR <- 
  list(
    alpha = post[, grep('alpha', colnames(post))],
    tau = post[, grep('tau', colnames(post))],
    beta = post[, grep('beta', colnames(post))]
  )


par(mar = c(4, 4, 1, 1))
plot(apply(post_noCOR$tau, 2, mean), 
     apply(post_noCOR$beta, 2, mean), 
     main = '')

# ========
cat(file = 'multilevel_mods_II/eg_cov_par1.stan', 
    '
    data{
      int N;
      int N_parche;
      int N_bosque;
      int N_nidos;
      array[N] int nido;
      array[N] int nidada;
      array[N] int bosque;
      array[N] int huevos;
      array[N] int parche;
    }
    
    parameters{
      matrix[N_parche, N_bosque] M_tau;
      corr_matrix[N_bosque] rho_parche;
      vector<lower = 0>[N_bosque] sigma_parche;
      vector[N_bosque] tau_mu;
    
      matrix[N_nidos, N_bosque] M_alpha;
      corr_matrix[N_bosque] rho_nido;
      vector<lower = 0>[N_bosque] sigma_nido;
    }
    
    transformed parameters{
      vector[N_parche] tau;
      vector[N_parche] beta;
      vector[N_nidos] alpha;
      
      alpha = M_alpha[, 1];
      tau = M_tau[, 1];
      beta = M_tau[, 2];
    }
    
    model{
      vector[N] p;
      sigma_parche ~ exponential(1);
      tau_mu ~ normal(0, 1);
      rho_parche ~ lkj_corr(1);
      sigma_nido ~ exponential(1);
      rho_nido ~ lkj_corr(1);
      for (i in 1:N_parche) M_tau[i, :] ~ multi_normal(tau_mu, 
                                                       quad_form_diag(rho_parche, 
                                                       sigma_parche));
      for (i in 1:N_nidos) M_alpha[i, :] ~ multi_normal(rep_vector(0, N_bosque),
                                                       quad_form_diag(rho_nido,
                                                       sigma_nido));
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]] + tau[parche[i]] + beta[parche[i]]*bosque[i];
        p[i] = inv_logit(p[i]);
      }
      
      huevos ~ binomial(nidada, p);
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
      array[N] int ppcheck;
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]] + tau[parche[i]] + beta[parche[i]]*bosque[i];
        p[i] = inv_logit(p[i]);
      }
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(huevos[i] | nidada[i], p[i]);
    
      ppcheck = binomial_rng(nidada, p);
    }
    ')

file <- paste(getwd(), '/multilevel_mods_II/eg_cov_par1.stan', sep = '')
fit <- cmdstan_model(file, compile = T)

m <- 
  fit$sample(
    data = dat,
    iter_sampling = 4e3,
    iter_warmup = 500, 
    thin = 3, 
    chains = 3, 
    parallel_chains = 3, 
    refresh = 200, 
    seed = 123
  )

m$save_object(paste(getwd(), '/multilevel_mods_II/eg_cov_par1.rds', sep = ''))

m <- readRDS('multilevel_mods_II/eg_cov_par1.rds')

m_out <- m$summary()

par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
for (i in 1:9) trace_plot(m, m_out[m_out$rhat >1.01,][-c(1:5), ]$variable[i+2], 3)
par(mforw = c(1, 1))

mod_diagnostics(m, m_out)

# non centered version

cat(file = 'multilevel_mods_II/eg_cov_par2.stan', 
    "
    data{
      int N;
      int N_parche;
      int N_bosque;
      int N_nidos;
      array[N] int nido;
      array[N] int nidada;
      array[N] int bosque;
      array[N] int huevos;
      array[N] int parche;
    }
    
    parameters{
      matrix[N_bosque, N_parche] Z_tau;
      cholesky_factor_corr[N_bosque] rho_tau;
      vector<lower = 0>[N_bosque] sigma_tau;
      vector[N_bosque] tau_bar;
    
      matrix[N_bosque, N_nidos] Z_alpha;
      cholesky_factor_corr[N_bosque] rho_alpha;
      vector<lower = 0>[N_bosque] sigma_alpha;
      vector[N_bosque] alpha_bar;
    }
    
    transformed parameters{
      vector[N_nidos] alpha;
      vector[N_nidos] beta2;
      vector[N_parche] tau;
      vector[N_parche] beta;
      matrix[N_nidos, N_bosque] M_alpha;
      matrix[N_parche, N_bosque] M_tau;
      M_alpha = (diag_pre_multiply(sigma_alpha, rho_alpha) * Z_alpha)';
      M_tau = (diag_pre_multiply(sigma_tau, rho_tau) * Z_tau)';
      alpha = alpha_bar[1] + M_alpha[, 1];
      beta2 = alpha_bar[2] + M_alpha[, 2];
      tau = tau_bar[1] + M_tau[, 1];
      beta = tau_bar[2] + M_tau[, 2];
    }
    
    model{
      vector[N] p;
      vector[N] p2;
      sigma_alpha ~ exponential(1);
      sigma_tau ~ exponential(1);
      tau_bar ~ normal(0, 1);
      alpha_bar ~ normal(0, 1);
      to_vector(Z_alpha) ~ normal(0, 1);
      to_vector(Z_tau) ~ normal(0, 1);
      rho_alpha ~ lkj_corr_cholesky(2);
      rho_tau ~ lkj_corr_cholesky(2);
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]] + tau[parche[i]] + beta[parche[i]]*bosque[i];
        p[i] = inv_logit(p[i]);
      }
      
      huevos ~ binomial(nidada, p);
    
    for (i in 1:N) {
        p2[i] = alpha[nido[i]] + tau[parche[i]] + beta2[parche[i]]*bosque[i];
        p2[i] = inv_logit(p2[i]);
      }
      
      huevos ~ binomial(nidada, p2);
    
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
      array[N] int ppcheck;
      matrix[N_bosque, N_bosque] rho_parche;
      matrix[N_bosque, N_bosque] rho_nido;
      rho_parche = multiply_lower_tri_self_transpose(rho_tau);
      rho_nido = multiply_lower_tri_self_transpose(rho_alpha);
    
      for (i in 1:N) {
        p[i] = alpha[nido[i]] + tau[parche[i]] + beta[parche[i]]*bosque[i];
        p[i] = inv_logit(p[i]);
      }
    
      for (i in 1:N) log_lik[i] = binomial_lpmf(huevos[i] | nidada[i], p[i]);
    
      ppcheck = binomial_rng(nidada, p);
    }
    ")

file <- paste(getwd(), '/multilevel_mods_II/eg_cov_par2.stan', sep = '')
fit <- cmdstan_model(file, compile = T)

m2 <- 
  fit$sample(
    data = dat,
    iter_sampling = 4e3,
    iter_warmup = 500, 
    thin = 3, 
    chains = 3, 
    parallel_chains = 3, 
    refresh = 200, 
    seed = 123
  )

m2$save_object('multilevel_mods_II/eg_cov_par2.rds')

m2 <- readRDS('multilevel_mods_II/eg_cov_par2.rds')

m2_out <- m2$summary()

mod_diagnostics(m2, m2_out)

ppcheck <- as.matrix(m2$draws(variables = 'ppcheck', format = 'df'))

plot(density(ppcheck[1, ]), xlab = 'Huevos sobrevivientes', 
     ylab = 'Density', main = '', lwd = 0.01, ylim = c(0, 0.12))
for (i in 1:500) lines(density(ppcheck[i, ]), lwd = 0.1)
lines(density(dat$huevos), col = 'red', lwd = 2)

post <- m2$draws(variables = c('alpha', 'tau', 'beta', 'beta2'), format = 'df')
post <- 
  list(alpha = apply(post[, grep('alpha', colnames(post))], 2, mean),
       tau = apply(post[, grep('tau', colnames(post))], 2, mean), 
       beta = apply(post[, grep('beta\\[', colnames(post))], 2, mean), 
       beta2 = apply(post[, grep('beta2', colnames(post))], 2, mean))

rho_parche <- m2$draws(variables = c('rho_parche'), 
                       format = 'df')

sigma_parche <- m2$draws(variables = 'sigma_tau', format = 'df')

sigmas_parche <- 
  matrix(c(mean(sigma_parche$`sigma_tau[1]`)^2,
             mean(rho_parche$`rho_parche[2,1]`),
             mean(rho_parche$`rho_parche[2,1]`),
             mean(sigma_parche$`sigma_tau[2]`)^2), ncol = 2)

mu_parche <- c(mean(post$tau), mean(post$beta))

mu_nido <- c(mean(post$alpha), mean(post$beta2))

sigmas_nido <- m2$draws(variables = 'sigma_alpha', format = 'df')
sigmas_nido <- apply(sigmas_nido, 2, mean)

rho_nido <- m2$draws(variables = 'rho_nido', format = 'df')

sigmas_nido <- 
  matrix(
    c(sigmas_nido[1]^2, 
      mean(rho_nido$`rho_nido[2,1]`), 
      mean(rho_nido$`rho_nido[2,1]`), 
      sigmas_nido[2]^2), ncol = 2
  )


tau_noCORR <- apply(post_noCOR$tau, 2, mean)
beta_noCORR <- apply(post_noCOR$beta, 2, mean)
obs_tau_beta <- parche[1:20, ]

par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
post %$% plot(tau, beta, xlab = expression(tau), 
              ylab = expression(beta), col = 2, 
              xlim = c(-3, 3), ylim = c(-1, 4.5), 
              main = 'Modelo\n par. pool. pars/clusters')
for (i in seq(0.1, 0.9, by = 0.2))
  lines(ellipse(sigmas_parche, centre = mu_parche, level = i), lwd = 0.5)

post %$% plot(tau, beta, xlab = expression(tau), 
              ylab = expression(beta), col = 2, 
              xlim = c(-3, 3), ylim = c(-1, 4.5), 
              main = 'Pars. estimados y reales')
points(tau_noCORR, beta_noCORR, col = 4)
points(obs_tau_beta[, 1], obs_tau_beta[, 2], pch = 16)
par(mfrow = c(1, 1))

parametros <- 
  tibble(beta_cor = post$beta, 
       tau_cor = post$tau,
       beta_noCORR = beta_noCORR, 
       tau_noCORR = tau_noCORR, 
       real_beta = obs_tau_beta[, 2], 
       real_tau = obs_tau_beta[, 1]) 

parametros <- 
  parametros |> 
  mutate(error_corr_beta = real_beta - beta_cor, 
         error_corr_tau = real_tau - tau_cor, 
         error_noCorr_beta = real_beta - beta_noCORR, 
         error_noCorr_tau = real_tau - tau_noCORR, )

par(mfrow = c(2, 2))
parametros %$% 
  plot(tau_noCORR, beta_noCORR, xlim = c(-1.5, 2.5), ylim = c(-1, 4), 
       main = 'Error de estimación\n(par. pool. clusters)', 
       xlab = expression(tau), ylab = expression(beta), col = 4)
parametros %$%
  points(real_tau, real_beta, pch = 16)
parametros %$%
  segments(x0 = tau_noCORR, x1 = real_tau, 
           y0 = beta_noCORR, y1 = real_beta, lty = 3, lwd = 0.8)

parametros %$% 
  plot(tau_cor, beta_cor, xlim = c(-1.5, 1.8), ylim = c(-1, 4),
       main = 'Error de estimación\n(par. pool. par/clusters)', 
       col = 2, xlab = expression(tau), ylab = expression(beta))
parametros %$%
  points(real_tau, real_beta, pch = 16)
parametros %$%
  segments(x0 = tau_cor, x1 = real_tau, 
           y0 = beta_cor, y1 = real_beta, lty = 3, lwd = 0.8)

parametros %$% plot(density(error_corr_beta), col = 2, lwd = 3, 
                    main = 'Error (observado - estimado)', 
                    xlab = expression(beta))
abline(v = median(parametros$error_corr_beta), lty = 3, col = 2, 
       lwd = 3)
parametros %$%
  lines(density(error_noCorr_beta), col = 4, lwd = 3)
abline(v = median(parametros$error_noCorr_beta), lty = 3, col = 4, 
       lwd = 3)

parametros %$% plot(density(error_corr_tau), col = 2, lwd = 3, 
                    main = 'Error (observado - estimado)', 
                    xlab = expression(tau), 
                    xlim = c(-1.5, 0.45))
abline(v = median(parametros$error_corr_tau), lty = 3, col = 2, 
       lwd = 3)
parametros %$%
  lines(density(error_noCorr_tau), col = 4, lwd = 3)
abline(v = median(parametros$error_noCorr_tau), lty = 3, col = 4, 
       lwd = 3)
par(mfrow = c(1, 1))

plot(density(rho_parche$`rho_parche[2,1]`), col =4, lwd = 3, 
     main = expression(rho[beta~tau]), xlab = expression(italic('r')))
abline(v = median(rho_parche$`rho_parche[2,1]`), lty = 3, col = 4, 
       lwd = 3)
abline(v = rho_tau_beta, lty = 3, col = 1, 
       lwd = 3)





plot(density(apply(post$alpha, 2, mean)))


colnames(post)


#
#
#
#
#
#
#
#
#
#
#
#
#
#

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

d <- na.omit(penguins)

d <- d |> dplyr::select(species, year, island, sex, flipper_length_mm, body_mass_g)

d$year <- as.factor(d$year)

d$species_sex <- as.factor(paste(d$species, d$sex, sep = '_'))

dat <- lapply(d, function(x) if(!is.double(x)) as.integer(x) else(x))

dat$N <- length(dat$species)
dat$N_spp <- max(dat$species)
dat$N_sex <- 2
dat$N_island <- max(dat$island)
dat$N_year <- max(dat$year)
dat$N_spp_sex <- max(dat$species_sex)

prior_corr <- rlkjcorr(1e3, 2, eta = 5)

dim(prior_corr)

plot(density(prior_corr[, 2, 1]))

cat(file = 'multilevel_mods_II/penguins_par_pool.stan', 
    "
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
    }")

file <- paste(getwd(), 'multilevel_mods_II/penguins_par_pool.stan', sep = '')

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

par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
for (i in 1:9) trace_plot(m1, out_m1$variable[i], 3)
par(mfrow = c(1, 1))

mod_diagnostics(m1, out_m1)

ppcheck_m1 <- m1$draws(variables = 'ppcheck', format = 'matrix')

plot(NULL, main = 'Partial pooling (intercepts)', 
     xlim = c(160, 245), ylim = c(0, 0.04), xlab = 'mm', 
     ylab = 'Density')
for (i in 1:100) lines(density(ppcheck_m1[i, ]), lwd = 0.1)
lines(density(dat$flipper_length_mm), col = 'red')

post_m1 <- m1$draws(variables = c('theta', 'alpha', 'tau', 'beta', 'sigma'), 
                    format = 'df')

post_m1 <- 
  list(spp = post_m1[, grepl('theta', colnames(post_m1))],
       I = post_m1[, grepl('^alpha', colnames(post_m1))],
       Y = post_m1[, grepl('^tau', colnames(post_m1))], 
       beta = post_m1[, grepl('beta', colnames(post_m1))], 
       sigma = post_m1[, grepl('sigma', colnames(post_m1))])

x_seq <- seq(min(dat$body_mass_g), max(dat$body_mass_g), length.out = 100)

sex <- c('H', 'M')


colnames(post_m1$spp) <- levels(d$species_sex)

nombres <- colnames(post_m1$spp)

est_cor <- lapply(1:ncol(post_m1$spp), FUN = 
                    function(i) {
                      
                      df <- 
                        sapply(x_seq, FUN = 
                                 function(z) {
                                   
                                   post_m1$spp[, i, drop = T] + 
                                     apply(post_m1$I, 1, mean) +
                                     apply(post_m1$Y, 1, mean) +
                                     post_m1$beta[, i, drop = T]*z
                                   
                                 })
                      
                      df <- do.call('rbind',
                                    apply(df, 2, FUN = 
                                            function(j) {
                                              tibble(mu_est = mean(j), 
                                                     li_est = quantile(j, 0.025), 
                                                     ls_est = quantile(j, 0.975))
                                            }, simplify = 'list'))
                      
                      df_pred <- 
                        sapply(x_seq, FUN = 
                                 function(z) {
                                   
                                   mu <- post_m1$spp[, i, drop = T] + 
                                     apply(post_m1$I, 1, mean) +
                                     apply(post_m1$Y, 1, mean) +
                                     post_m1$beta[, i, drop = T]*z
                                   
                                   rnorm(1e3, mu, post_m1$sigma$sigma)
                                   
                                 })
                      
                      df_pred <- do.call('rbind',
                                    apply(df_pred, 2, FUN = 
                                            function(j) {
                                              tibble(li_pred = quantile(j, 0.025), 
                                                     ls_pred = quantile(j, 0.975))
                                            }, simplify = 'list'))
                      
                      df <- cbind(df, df_pred)
                      df$x <- x_seq
                      df$spp <- gsub('^(.*)(_)(.*)$', '\\1', nombres[i])
                      df$sex <- gsub('^(.*)(_)(.*)$', '\\3', nombres[i])
                      df

                    })

est_cor <- do.call('rbind', est_cor)

est_cor <- as_tibble(est_cor)
colnames(est_cor)[7] <- 'species'

ggplot() +
  geom_point(data = d, aes(body_mass_g, flipper_length_mm, color = sex)) +
  geom_line(data = est_cor, 
            aes(x, mu_est, color = sex)) +
  geom_ribbon(data = est_cor, 
              aes(x, mu_est, fill = sex, 
                  ymin = li_pred, ymax = ls_pred), alpha = 0.3) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = 'Masa corporal (g)', y = 'Longitud de la aleta') +
  theme(
    panel.grid = element_blank(),
    legend.position = 'top',
    legend.title = element_blank()
  )

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


rho_isla <- m1$draws(variables = 'Rho_isla_corr', format = 'df')
rho_year <- m1$draws(variables = 'Rho_year_corr', format = 'df')

rho_isla <- lapply(2:6, FUN = 
                     function(i) {
                       
                       df <- sapply(i:6, FUN = 
                                      function(j) {
                                        
                                        rho_isla[, grep(paste('^(.*)', i-1, ',', j,'.$', sep = ''), 
                                                        colnames(rho_isla)), drop = T]
                                        
                                      })
                       name <- i:6
                       df <- as_tibble(df)
                       
                       for (z in seq_along(name)) colnames(df)[z] <- 
                         paste('rho', i-1, name[z], sep = '_')
                       df
                     })
rho_isla <- as_tibble(do.call('cbind', rho_isla))

rho_year <- lapply(2:6, FUN = 
                     function(i) {
                       
                       df <- sapply(i:6, FUN = 
                                      function(j) {
                                        
                                        rho_year[, grep(paste('^(.*)', i-1, ',', j,'.$', sep = ''), 
                                                        colnames(rho_year)), drop = T]
                                        
                                      })
                       name <- i:6
                       df <- as_tibble(df)
                       
                       for (z in seq_along(name)) colnames(df)[z] <- 
                         paste('rho', i-1, name[z], sep = '_')
                       df
                     })
rho_year <- as_tibble(do.call('cbind', rho_year))
  
rho_isla
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(NULL, xlim = c(-1, 1), ylim = c(0, 1.5), xlab = expression('Rho'['isla']), ylab = 'Density')
for (i in 1:ncol(rho_isla)) {
  lines(density(rho_isla[, i, drop = T]), col = i)
  abline(v = mean(rho_isla[, i, drop = T]), lty = 3, col = i)
}

plot(NULL, xlim = c(-1, 1), ylim = c(0, 1.5), xlab = expression('Rho'['year']), ylab = 'Density')
for (i in 1:ncol(rho_year)) {
  lines(density(rho_year[, i, drop = T]), col = i)
  abline(v = mean(rho_year[, i, drop = T]), lty = 3, col = i)
}





# ==============



data(and_vertebrates)

?and_vertebrates

d <- and_vertebrates

d %$% aggregate(length_1_mm ~ species + year + 
                  sitecode + section + unittype, FUN = length)


