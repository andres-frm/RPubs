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

options(mc.cores = parallel::detectCores())

cat(file = 'stan/chapter3_1.stan', 
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

writeLines(readLines('stan/chapter3_1.stan'))

m3_1_file <- paste(getwd(), '/stan/', 'chapter3_1.stan', sep = '')

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

print(fit_me_1)

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


cat(file = 'stan/cap3_2.stan', 
    '
    data {
      int N;
      vector[N] income;
      array[N] int trofeos; 
      array[N] int trofeos2;
      array[N] int trofeos3; 
      vector[N] year;
    }
    
    parameters {
      real alpha;
      real beta;
      real<lower = 0> sigma;
      real alpha1;
      real beta1;
      real alpha2;
      real beta2;
      real alpha3;
      real beta3;
    }
    
    model {
    
      vector[N] lambda; 
      
      for (i in 1:N) {
        
        // M1. log de la funcion de distribución de probabilidas
        target += normal_lpdf(income[1:i] | 
                              alpha + beta*year[1:i], sigma);
        
        // M2. versión abreviada con el muestreo explícito (M 2)
        trofeos[1:i] ~ poisson(exp(alpha1 + beta1*year[1:i]));
    
        // M3. versión larga con el muestreo explícito (M 3)
        lambda[i] = alpha2 + beta2*year[i];
        lambda[i] = exp(lambda[i]);
        
      }
      
      // M4. versión vectorizada del log de la función de distribución de probabilidad
      target += poisson_lpmf(trofeos3 | exp(alpha3 + beta3*year));
      target += normal_lpdf(alpha3 | 10, 5); // previa M4
      target += normal_lpdf(beta3 | 0.5, 1); // previa M4
    
      target += normal_lpdf(alpha | 50, 15); // previa M1
      target += normal_lpdf(beta | 0.5, 1);  // previa M1
      target += exponential_lpdf(sigma | 1); // previa M1
      
      alpha1 ~ normal(10, 5); // previa M2
      beta1 ~ normal(0.5, 1); // previa M2
    
      alpha2 ~ normal(10, 5); // previa M3
      beta2 ~ normal(0.5, 1); // previa M3
      trofeos2 ~ poisson(lambda);
    
    }
    ')

path.file <- paste(getwd(), '/stan/', 'cap3_2.stan', sep = '')

fit_cap3_2 <- cmdstan_model(path.file, compile = T)

m_cap3_2 <- 
  fit_cap3_2$sample(data = 
                      list(N = N, year = year, 
                           income = income, 
                           trofeos = trofeos,
                           trofeos2 = trofeos,
                           trofeos3 = trofeos), 
                    chains = 3, 
                    parallel_chains = 3, 
                    iter_sampling = 2e3, 
                    iter_warmup = 500, 
                    thin = 3,
                    seed = 123)

m_cap3_2_out <- m_cap3_2$summary()



par(mfrow = c(4, 3))
for (i in 1:10) trace_plot(m_cap3_2, m_cap3_2_out$variable[i], 3)
par(mfrow = c(1, 1))

post_m_cap3_2 <- m_cap3_2$draws(format = 'df')

post_m_cap3_2 <- post_m_cap3_2[, -c(1:4)]

alpha <- grep('alpha', colnames(post_m_cap3_2))
beta <- grep('beta', colnames(post_m_cap3_2))

par(mfrow = c(1, 2))
plot(density(post_m_cap3_2[, alpha[1], drop = T]), 
     xlab = 'Alpha', main = '', xlim = c(0.8, 3))
for (i in 1:3) lines(density(post_m_cap3_2[, alpha[i], drop = T]), 
                    xlab = 'Alpha', main = '', col = i)

plot(density(post_m_cap3_2[, beta[1], drop = T]), 
     xlab = 'Beta', main = '', xlim = c(-0.02, 0.07))
for (i in 1:3) lines(density(post_m_cap3_2[, beta[i], drop = T]), 
                     xlab = 'Beta', main = '', col = i)
par(mfrow = c(1, 1))









