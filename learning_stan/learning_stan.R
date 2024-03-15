# ======== Learnin stan =======
# 
# 
# ========Chapter 3=======

paquetes <- c("tidyverse", "magrittr", 'patchwork', 'rstan',
              'cmdstanr', 'loo', 'palmerpenguins')
sapply(paquetes, 
       function(x) {
         if(require(x, character.only = T)) library(x, character.only = T)
         else install.packages(x); library(x, character.only = T)
       })

options(mc.cores = parallel::detectCores(), mar = c(4, 4, 1, 1))

# ===== functions ======

trace_plot <- function(fit, par, n_chains) {
  
  d <- as_mcmc.list(fit)
  
  d <- lapply(d, FUN = 
                function(x) {
                  x <- as.data.frame(unclass(as.matrix(x)))
                  x$iter <- 1:nrow(x)
                  x
                })
  
  for (i in 1:n_chains) d[[i]]$chain <- i
  
  plot(d[[1]][, 'iter', drop = T], d[[1]][, par, drop = T], 
       type = 'l', col = 1, main = '', ylab = par, 
       xlab = 'Iteration')
  for (i in 2:n_chains) {
    lines(d[[i]][, 'iter', drop = T], d[[i]][, par, drop = T], 
          col = i, main = '', ylab = par, 
          xlab = 'Iteration')
  }
  
}

mod_diagnostics <- function(model, output) {
  
  l <- model$loo()
  
  pk <- l$diagnostics$pareto_k
  
  pareto <- tibble(x = 1:length(pk), 
                   pareto_k = pk)
  
  diags <- na.omit(tibble(x = 1:nrow(output), 
                          rhat = output$rhat, 
                          ess_bulk = output$ess_bulk, 
                          ess_tail = output$ess_tail))
  
  par(mfrow = c(1, 3), mar = c(4, 4, 1, 1))
  diags %$% plot(rhat ~ x, pch = 16, col = 'tomato3', xlab = 'Parameters', 
                 ylab = 'Rhat')
  abline(h = 1.1, lty = 3, col = 'red')
  diags %$% plot(sort(ess_bulk) ~ x, pch = 16, col = 'cyan4', xlab = 'Parameters', 
                 ylab = 'ESS')
  diags %$% points(sort(ess_tail) ~ x, pch = 16, col = 'tan1')
  abline(h = 1000, lty = 3, col = 'red')
  pareto %$% plot(pareto_k ~ x, pch = 16, col = 'purple', xlab = 'Observations', 
                  ylab = 'Pareto-k', ylim = c(-0.5, 1.2))
  abline(h = c(0.5, 0.7, 1), lty = 3, col = 'red')
  text(x = rep(nrow(pareto)/2, 4), y = c(0.25, 0.6, 0.8, 1.2), 
       labels = c('Perfect', 'ok', 'high', 'too hight'))
  par(mfrow = c(1, 1))
  
}

# ======================================== #


cat(file = 'learning_stan/chapter3_1.stan', 
    '
    data {
        int N;
        vector[N] est;
        array[N] int est1;
        vector[N] est2;
    }
    
    parameters {
        real mu;
        real<lower = 0> sigma;
        real lambda;
        real mu2;
        real<lower = 0> sigma1;
    }
    
    model {
        
        for (i in 1:N) {
          est[i] ~ normal(mu, sigma);
        }
    
        for (i in 1:N) {
          est1[i] ~ poisson(exp(lambda));
        }
    
        for (i in 1:N) {
          target += normal_lpdf(est2[i] | mu2, sigma1);
        }
    
    
        mu ~ normal(120, 30);
        sigma ~ exponential(1);
        lambda ~ normal(5, 1);
        target += normal_lpdf(mu2 | 120, 30);
        target += exponential_lpdf(sigma1 | 1);
    
    // este es un comentario
    
    /* este también es un comentario,
    pero puede distribuirse en diferentes 
    lines */
    
    }
  
    ')

writeLines(readLines('learning_stan/chapter3_1.stan'))

m3_1_file <- paste(getwd(), '/learning_stan/chapter3_1.stan', sep = '')

m3_1 <- cmdstan_model(m3_1_file, compile = T)

d <- list(est = rnorm(100, 140, 20), 
          est1 = rpois(100, 120),
          est2 = rnorm(100, 140, 20),
          N = 100)

fit_me_1 <- 
  m3_1$sample(data = d,
              seed = 123, 
              refresh = 500, 
              parallel_chains = 3, 
              chains = 3, 
              iter_sampling = 2e3, 
              iter_warmup = 500)

fit_me_1_out <- fit_me_1$summary()

par(mfrow = c(3, 2), mar = c(4, 4, 1, 1))
for (i in 1:6) trace_plot(fit_me_1, fit_me_1_out$variable[i], 3)
par(mfrow = c(1, 1))



# ====== chapter 4 =======

# first LM

# simulating years of experience and anual income 


N <- 30

year <- rpois(N, lambda = 15)

income <- rnorm(N, (15 + 2.5*year), 2) + rnorm(N, 3, 2)

trofeos <- rpois(N, (2 + 0.5*year))

plot(income ~ year)
plot(trofeos ~ year)

# version 1

cat(file = 'learning_stan/cap3_2_1.stan', 
    '
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
    ')

path.file <- paste(getwd(), '/learning_stan/cap3_2_1.stan', sep = '')

fit_cap3_2_1 <- cmdstan_model(path.file, compile = T)

m_cap3_2_1 <- 
  fit_cap3_2_1$sample(data = 
                      list(N = N, year = year, 
                           trofeos = trofeos), 
                    chains = 3, 
                    parallel_chains = 3, 
                    iter_sampling = 4e3, 
                    iter_warmup = 500, 
                    thin = 3,
                    seed = 123)

m_cap3_21_out <- m_cap3_2_1$summary()

par(mfrow = c(4, 3))
for (i in 1:10) trace_plot(m_cap3_2_1, m_cap3_21_out$variable[i], 3)
par(mfrow = c(1, 1))

mod_diagnostics(m_cap3_2_1, m_cap3_21_out)

# version 2

cat(file = 'learning_stan/cap3_2_3.stan', 
    '
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
    ')

path.file <- paste(getwd(), '/learning_stan/cap3_2_3.stan', sep = '')

fit_cap3_2_3 <- cmdstan_model(path.file, compile = T)

m_cap3_2_3 <- 
  fit_cap3_2_3$sample(data = 
                        list(N = N, year = year,
                             trofeos = trofeos), 
                      chains = 3, 
                      parallel_chains = 3, 
                      iter_sampling = 4e3, 
                      iter_warmup = 500, 
                      thin = 3,
                      seed = 123)

m_cap3_23_out <- m_cap3_2_3$summary()

par(mfrow = c(4, 3))
for (i in 1:10) trace_plot(m_cap3_2_3, m_cap3_23_out$variable[i], 3)
par(mfrow = c(1, 1))

mod_diagnostics(m_cap3_2_3, m_cap3_23_out)

loo_compare(m_cap3_2_1$loo(), m_cap3_2_3$loo())

m_cap3_21_out[2:3,]
m_cap3_23_out[2:3,]



#### transformed parameters and generated quantities 


cat(file = 'learning_stan/cap_3_2_31.stan', 
    '
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
    
    }')

path.file <- paste(getwd(), '/learning_stan/cap_3_2_31.stan', sep = '')

fit_cap3_2_31 <- cmdstan_model(path.file, compile = T)

m_cap3_2_31 <- 
  fit_cap3_2_31$sample(data = 
                        list(N = N, year = year,
                             trofeos = trofeos, 
                             year_x = seq(min(year), max(year), by = 1),
                             N_x = length(seq(min(year), max(year), by = 1))), 
                      chains = 3, 
                      parallel_chains = 3, 
                      iter_sampling = 4e3, 
                      iter_warmup = 500, 
                      thin = 3,
                      seed = 123)

m_cap3_231_out <- m_cap3_2_31$summary()

par(mfrow = c(4, 3), mar = c(4, 4, 1, 1))
for (i in 1:10) trace_plot(m_cap3_2_31, m_cap3_231_out$variable[i], 3)
par(mfrow = c(1, 1))

post <- m_cap3_2_31$draws(format = 'df')

pred <- post[, grepl('pred', colnames(post))]

pred <- 
  apply(pred, 2, 
      function(x) {
        tibble(mu = mean(x), 
               li = quantile(x, 0.025), 
               ls = quantile(x, 0.975))
      })

pred <- do.call('rbind', pred)

pred$x <- seq(min(year), max(year), by = 1)

est <- post[, grepl('theta', colnames(post))]

est <- 
  apply(est, 2, 
        function(x) {
          tibble(mu = mean(x), 
                 li = quantile(x, 0.025), 
                 ls = quantile(x, 0.975))
        })

est <- do.call('rbind', est)

est$x <- seq(min(year), max(year), by = 1)

ggplot() +
  geom_line(data = pred, aes(x = x, y = mu)) +
  geom_ribbon(data = pred, aes(x = x, ymin = li, ymax = ls), 
              alpha = 0.3) +
  geom_ribbon(data = est, aes(x = x, ymin = li, ymax = ls), 
              alpha = 0.3) +
  geom_point(data = tibble(x = year, y = trofeos), aes(x, y))


ppcheck <- as.matrix(post[, grepl('ppcheck', colnames(post))])

plot(density(trofeos), main = '', xlab = 'Trofeos', ylim = c(0, 0.20),
     xlim = c(-5, 20))
for (i in 1:200) lines(density(ppcheck[i, ]), lwd = 0.1)
lines(density(trofeos), col = 'red', lwd = 2)










