
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
    