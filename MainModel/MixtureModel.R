## part 2: estimation of mean of antibody level using mixture model

rm(list=ls(all=TRUE))
library(tidyverse)
library(rstan)
library(ggmcmc)
library(patchwork)

d <- read.csv(file = 'data_under5.txt')
d <- d %>% filter(age_days > 183)
d <- d %>% mutate(age_days/365)
d <- d %>% mutate(log(IgG_PreF))

d_old <- read.csv(file = 'data_over5.txt')
d_old <- d_old %>% mutate(log(Prefusion_F))

X_new <- c(0.504, 0.6, 0.7, 0.8, 0.9, 1,1.5, 2, 3, 4, 5.08, 7, 12, 20, 40, 60, 80) 

data <- list(N=nrow(d), Y=d$`log(IgG_PreF)`, Age=d$`age_days/365`, 
             N_old=nrow(d_old), Y_old=d_old$`log(Prefusion_F)`, Age_old=d_old$age, N_new=length(X_new), X_new=X_new) # age as years

stanmodel <- stan_model(file = 'MixtureModel.stan') ## SIS model
stanmodel <- stan_model(file = 'SensitivityAnalysisMixtureModel.stan') ## boosting model

# SIS model
fit <- sampling(
  stanmodel,
  data=data,
  init = function(){
    list(sigma1=log(100), sigma2=log(100), sigma3=log(100))
  },
  seed = 1234,
  chains=4, iter=5000, warmup=500, thin=1
)


# SIW model
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
