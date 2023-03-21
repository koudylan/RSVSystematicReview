## estimation of mean of antibody level using mixture model assuming SIW
rm(list=ls(all=TRUE))
library(tidyverse)
library(rstan)
library(ggmcmc)
library(patchwork)

# model (SIW catalytic model)
# sigma follow student_t
modelString <- "data {
  int N;
  real Y[N];
  real Age[N];
  int N_old;
  real Y_old[N_old];
  real Age_old[N_old];
  int N_new;
  real X_new[N_new];
}

parameters {
  positive_ordered[2] mu;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real<lower=0> foi1;
  real<lower=0> foi2;
  real<lower=0> foi3;
  real<lower=0> decay;
  real<lower=0> boost;
  real<lower=0> sigma3;
}

model {
  foi1 ~ gamma(64, 80);
  foi2 ~ gamma(56, 38);
  foi3 ~ gamma(56, 38);
  decay ~ gamma(25,250);
  //boost ~ exponential(0.2);
  boost ~ gamma(4,0.8);
  sigma1 ~ student_t(4, 0, 10);
  sigma2 ~ student_t(4, 0, 10);
  sigma3 ~ student_t(4, 0, 10);
  for (n in 1:N)
    if(Age[n] < 1){
    target += log_sum_exp(
      log((((1-boost)*foi1/((1-boost)*foi1-decay))*(exp(-(boost*foi1+decay)*(Age[n]-0.5))-exp(-foi1*(Age[n]-0.5)))+
(boost*foi1/(boost*foi1+decay))*(1-exp(-(boost*foi1+decay)*(Age[n]-0.5))))) + normal_lpdf(Y[n] | mu[2], sigma1),
      log1m((((1-boost)*foi1/((1-boost)*foi1-decay))*(exp(-(boost*foi1+decay)*(Age[n]-0.5))-exp(-foi1*(Age[n]-0.5)))+
(boost*foi1/(boost*foi1+decay))*(1-exp(-(boost*foi1+decay)*(Age[n]-0.5))))) + normal_lpdf(Y[n] | mu[1], sigma2)
      );
    } else {
    target += log_sum_exp(
      log((((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(Age[n]-1))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5))*(exp(-(boost*foi2+decay)*(Age[n]-1))-exp(-foi2*(Age[n]-1)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(Age[n]-1))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*(1-exp(-(boost*foi2+decay)*(Age[n]-1))))) + normal_lpdf(Y[n] | mu[2], sigma1),
      log1m((((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(Age[n]-1))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5))*(exp(-(boost*foi2+decay)*(Age[n]-1))-exp(-foi2*(Age[n]-1)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(Age[n]-1))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*(1-exp(-(boost*foi2+decay)*(Age[n]-1))))) + normal_lpdf(Y[n] | mu[1], sigma2)
      );   
    }
   for (n in 1:N_old){
     Y_old[n] ~ normal(mu[2]*(((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(Age_old[n]-5))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5)-(boost*foi3+decay)*(Age_old[n]-5))*(exp(-(boost*foi2+decay)*(5-1))-exp(-foi2*(5-1)))+
((1-boost)*foi3/((1-boost)*foi3-decay))*exp(-foi1*(1-0.5)-foi2*(5-1))*(exp(-(boost*foi3+decay)*(Age_old[n]-5))-exp(-foi3*(Age_old[n]-5)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(Age_old[n]-5))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*exp(-(boost*foi3+decay)*(Age_old[n]-5))*(1-exp(-(boost*foi2+decay)*(5-1)))+
(boost*foi3/(boost*foi3+decay))*(1-exp(-(boost*foi3+decay)*(Age_old[n]-5)))) + 
    mu[1]*(1 - (((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(Age_old[n]-5))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5)-(boost*foi3+decay)*(Age_old[n]-5))*(exp(-(boost*foi2+decay)*(5-1))-exp(-foi2*(5-1)))+
((1-boost)*foi3/((1-boost)*foi3-decay))*exp(-foi1*(1-0.5)-foi2*(5-1))*(exp(-(boost*foi3+decay)*(Age_old[n]-5))-exp(-foi3*(Age_old[n]-5)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(Age_old[n]-5))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*exp(-(boost*foi3+decay)*(Age_old[n]-5))*(1-exp(-(boost*foi2+decay)*(5-1)))+
(boost*foi3/(boost*foi3+decay))*(1-exp(-(boost*foi3+decay)*(Age_old[n]-5))))), sigma3);
   }
}

generated quantities{
  real Y_pre_new[N_new];
  real infected_new[N_new];
  real log_lik[N+N_old];
  // predicted interval of antibody levels
  for (n in 1:N_new)
   if(X_new[n] < 1){
   Y_pre_new[n] = (((1-boost)*foi1/((1-boost)*foi1-decay))*(exp(-(boost*foi1+decay)*(X_new[n]-0.5))-exp(-foi1*(X_new[n]-0.5)))+
(boost*foi1/(boost*foi1+decay))*(1-exp(-(boost*foi1+decay)*(X_new[n]-0.5))))*normal_rng(mu[2], sigma1) + 
   (1-(((1-boost)*foi1/((1-boost)*foi1-decay))*(exp(-(boost*foi1+decay)*(X_new[n]-0.5))-exp(-foi1*(X_new[n]-0.5)))+
(boost*foi1/(boost*foi1+decay))*(1-exp(-(boost*foi1+decay)*(X_new[n]-0.5)))))*normal_rng(mu[1], sigma2);
   } else if(X_new[n] >= 1 && X_new[n] < 5){
   Y_pre_new[n] = (((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(X_new[n]-1))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5))*(exp(-(boost*foi2+decay)*(X_new[n]-1))-exp(-foi2*(X_new[n]-1)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(X_new[n]-1))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*(1-exp(-(boost*foi2+decay)*(X_new[n]-1))))*normal_rng(mu[2], sigma1) + 
  (1-(((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(X_new[n]-1))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5))*(exp(-(boost*foi2+decay)*(X_new[n]-1))-exp(-foi2*(X_new[n]-1)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(X_new[n]-1))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*(1-exp(-(boost*foi2+decay)*(X_new[n]-1)))))*normal_rng(mu[1], sigma2);
   } else {
   Y_pre_new[n]=normal_rng(mu[2]*(((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(X_new[n]-5))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5)-(boost*foi3+decay)*(X_new[n]-5))*(exp(-(boost*foi2+decay)*(5-1))-exp(-foi2*(5-1)))+
((1-boost)*foi3/((1-boost)*foi3-decay))*exp(-foi1*(1-0.5)-foi2*(5-1))*(exp(-(boost*foi3+decay)*(X_new[n]-5))-exp(-foi3*(X_new[n]-5)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(X_new[n]-5))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*exp(-(boost*foi3+decay)*(X_new[n]-5))*(1-exp(-(boost*foi2+decay)*(5-1)))+
(boost*foi3/(boost*foi3+decay))*(1-exp(-(boost*foi3+decay)*(X_new[n]-5)))) + 
    mu[1]*(1 - (((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(X_new[n]-5))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5)-(boost*foi3+decay)*(X_new[n]-5))*(exp(-(boost*foi2+decay)*(5-1))-exp(-foi2*(5-1)))+
((1-boost)*foi3/((1-boost)*foi3-decay))*exp(-foi1*(1-0.5)-foi2*(5-1))*(exp(-(boost*foi3+decay)*(X_new[n]-5))-exp(-foi3*(X_new[n]-5)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(X_new[n]-5))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*exp(-(boost*foi3+decay)*(X_new[n]-5))*(1-exp(-(boost*foi2+decay)*(5-1)))+
(boost*foi3/(boost*foi3+decay))*(1-exp(-(boost*foi3+decay)*(X_new[n]-5))))), sigma3);
   }
  // credible interval of infected proportion 
  for (n in 1:N_new)
   if(X_new[n] < 1){
   infected_new[n] = ((1-boost)*foi1/((1-boost)*foi1-decay))*(exp(-(boost*foi1+decay)*(X_new[n]-0.5))-exp(-foi1*(X_new[n]-0.5)))+
(boost*foi1/(boost*foi1+decay))*(1-exp(-(boost*foi1+decay)*(X_new[n]-0.5)));
   } else if(X_new[n] >= 1 && X_new[n] < 5){
  infected_new[n] = ((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(X_new[n]-1))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5))*(exp(-(boost*foi2+decay)*(X_new[n]-1))-exp(-foi2*(X_new[n]-1)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(X_new[n]-1))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*(1-exp(-(boost*foi2+decay)*(X_new[n]-1)));
  } else {
  infected_new[n] = ((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(X_new[n]-5))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5)-(boost*foi3+decay)*(X_new[n]-5))*(exp(-(boost*foi2+decay)*(5-1))-exp(-foi2*(5-1)))+
((1-boost)*foi3/((1-boost)*foi3-decay))*exp(-foi1*(1-0.5)-foi2*(5-1))*(exp(-(boost*foi3+decay)*(X_new[n]-5))-exp(-foi3*(X_new[n]-5)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(X_new[n]-5))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*exp(-(boost*foi3+decay)*(X_new[n]-5))*(1-exp(-(boost*foi2+decay)*(5-1)))+
(boost*foi3/(boost*foi3+decay))*(1-exp(-(boost*foi3+decay)*(X_new[n]-5)));  
  }
  
  // waic
  for (n in 1:N)
   if(Age[n] < 1){
   log_lik[n] = (((1-boost)*foi1/((1-boost)*foi1-decay))*(exp(-(boost*foi1+decay)*(Age[n]-0.5))-exp(-foi1*(Age[n]-0.5)))+
(boost*foi1/(boost*foi1+decay))*(1-exp(-(boost*foi1+decay)*(Age[n]-0.5))))*normal_lpdf(Y[n] | mu[2], sigma1) + 
   (1-(((1-boost)*foi1/((1-boost)*foi1-decay))*(exp(-(boost*foi1+decay)*(Age[n]-0.5))-exp(-foi1*(Age[n]-0.5)))+
(boost*foi1/(boost*foi1+decay))*(1-exp(-(boost*foi1+decay)*(Age[n]-0.5)))))*normal_lpdf(Y[n] | mu[1], sigma2);
   } else {
   log_lik[n] = (((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(Age[n]-1))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5))*(exp(-(boost*foi2+decay)*(Age[n]-1))-exp(-foi2*(Age[n]-1)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(Age[n]-1))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*(1-exp(-(boost*foi2+decay)*(Age[n]-1))))*normal_lpdf(Y[n] | mu[2], sigma1) + 
  (1-(((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(Age[n]-1))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5))*(exp(-(boost*foi2+decay)*(Age[n]-1))-exp(-foi2*(Age[n]-1)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(Age[n]-1))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*(1-exp(-(boost*foi2+decay)*(Age[n]-1)))))*normal_lpdf(Y[n] | mu[1], sigma2);
   }
   for (t in 1:N_old){
   log_lik[N+t]=normal_lpdf(Y_old[t] | mu[2]*(((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(Age_old[t]-5))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5)-(boost*foi3+decay)*(Age_old[t]-5))*(exp(-(boost*foi2+decay)*(5-1))-exp(-foi2*(5-1)))+
((1-boost)*foi3/((1-boost)*foi3-decay))*exp(-foi1*(1-0.5)-foi2*(5-1))*(exp(-(boost*foi3+decay)*(Age_old[t]-5))-exp(-foi3*(Age_old[t]-5)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(Age_old[t]-5))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*exp(-(boost*foi3+decay)*(Age_old[t]-5))*(1-exp(-(boost*foi2+decay)*(5-1)))+
(boost*foi3/(boost*foi3+decay))*(1-exp(-(boost*foi3+decay)*(Age_old[t]-5)))) + 
    mu[1]*(1 - (((1-boost)*foi1/((1-boost)*foi1-decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(Age_old[t]-5))*(exp(-(boost*foi1+decay)*(1-0.5))-exp(-foi1*(1-0.5)))+
((1-boost)*foi2/((1-boost)*foi2-decay))*exp(-foi1*(1-0.5)-(boost*foi3+decay)*(Age_old[t]-5))*(exp(-(boost*foi2+decay)*(5-1))-exp(-foi2*(5-1)))+
((1-boost)*foi3/((1-boost)*foi3-decay))*exp(-foi1*(1-0.5)-foi2*(5-1))*(exp(-(boost*foi3+decay)*(Age_old[t]-5))-exp(-foi3*(Age_old[t]-5)))+
(boost*foi1/(boost*foi1+decay))*exp(-(boost*foi2+decay)*(5-1)-(boost*foi3+decay)*(Age_old[t]-5))*(1-exp(-(boost*foi1+decay)*(1-0.5)))+
(boost*foi2/(boost*foi2+decay))*exp(-(boost*foi3+decay)*(Age_old[t]-5))*(1-exp(-(boost*foi2+decay)*(5-1)))+
(boost*foi3/(boost*foi3+decay))*(1-exp(-(boost*foi3+decay)*(Age_old[t]-5))))), sigma3);
   }
}
"
# excluding aged < 0.5 year
d <- read.csv(file = 'data_under5.txt')
d <- d %>% filter(age_days > 183)
d <- d %>% mutate(age_days/365)
d <- d %>% mutate(log(IgG_PreF))

d_old <- read.csv(file = 'data_over5.txt')
d_old <- d_old %>% mutate(log(Prefusion_F))

X_new <- c(0.504, 0.6, 0.7, 0.8, 0.9, 1,1.5, 2, 3, 4, 5.08, 7, 12, 20, 40, 60, 80) 

data <- list(N=nrow(d), Y=d$`log(IgG_PreF)`, Age=d$`age_days/365`, 
             N_old=nrow(d_old), Y_old=d_old$`log(Prefusion_F)`, Age_old=d_old$age, N_new=length(X_new), X_new=X_new) # age as years

stanmodel <- stan_model(model_code = modelString) 

fit <- sampling(
  stanmodel,
  data=data,
  init = function(){
    list(sigma1=log(100), sigma2=log(100), sigma3=log(100))
  },
  seed = 1234,
  chains=4, iter=5000, warmup=1000, thin=2
)

ms <- rstan::extract(fit)

quantile(ms$foi1, probs = c(0.025, 0.5, 0.975))
quantile(ms$foi2, probs = c(0.025, 0.5, 0.975))
quantile(ms$foi3, probs = c(0.025, 0.5, 0.975))
quantile(ms$decay, probs = c(0.025, 0.5, 0.975))
quantile(ms$decay2, probs = c(0.025, 0.5, 0.975))
quantile(ms$mu[,1], probs = c(0.025, 0.5, 0.975))
quantile(ms$mu[,2], probs = c(0.025, 0.5, 0.975))

write.table(data.frame(summary(fit)$summary),
            file = 'fit-summary.txt', sep = '\t', quote = FALSE, col.names = NA)

ggmcmc(ggs(fit, inc_warmup=TRUE, stan_include_auxiliar=TRUE),
       file='fit-traceplot.pdf', plot='traceplot')

ggmcmc(ggs(fit), file='fit-ggmcmc.pdf')


# plot of observed and predicted antibody levels
data.frame.quantile.mcmc <- function(x, y_mcmc, probs=c(2.5, 25, 50, 75, 97.5)/100) {
  qua <- apply(y_mcmc, 2, quantile, probs=probs)
  d <- data.frame(X=x, t(qua))
  colnames(d) <- c('X', paste0('p', probs*100))
  return(d)
}

ggplot.5quantile <- function(data){
  p <- ggplot(data = data, aes(x=X, y=p50))
  p <- p + theme_classic(base_size = 18)
  p <- p + geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill='black', alpha=1/6)
  p <- p + geom_ribbon(aes(ymin=p25, ymax=p75), fill='black', alpha=2/6)
  p <- p + geom_line(size=1)
  return(p)
}

d_est <- data.frame.quantile.mcmc(x=X_new, y_mcmc = exp(ms$Y_pre_new))

# plot for aged < 5years
d_est1 <- d_est %>% filter(X < 7)
p1 <- ggplot.5quantile(data = d_est1)
p1 <- p1 + geom_point(data=d, aes(x=age_days/365, y=IgG_PreF),
                    alpha=0.4,
                    shape=21,
                    size=1.5)
p1 <- p1 + scale_y_log10(breaks = 10^(-1:5))
p1 <- p1 + labs(x='Age', y='Prefusion F IgG(AU/mL)')
print(p1)

# plot for aged >= 5 years
d_est2 <- d_est %>% filter(X >= 7)
p2 <- ggplot.5quantile(data = d_est2)
p2 <- p2 + geom_point(data=d_old, aes(x=age, y=Prefusion_F),
                      alpha=0.4,
                      shape=21,
                      size=1.5)
p2 <- p2 + scale_y_log10(breaks= 10^(-1:5))
p2 <- p2 + expand_limits(y=c(0.1:10000))
p2 <- p2 + labs(x='Age', y='Prefusion F IgG(AU/mL)')
print(p2)

# plot for proportion of infected
d_est3 <- data.frame.quantile.mcmc(x=X_new, y_mcmc = ms$infected_new)
p3 <- ggplot.5quantile(data = d_est3)
p3 <- p3 + labs(x='Age', y='Proportion of infected')
print(p3)

p3 + p1 + p2
