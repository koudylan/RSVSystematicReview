# estimation of FoIs using seroepidemiological studies (ELISA/F/N/Neutr/CF)
rm(list=ls(all=TRUE))
library(tidyverse)
library(rstan)
library(ggmcmc)
library(ggplot2)
library(patchwork)

# Model (catalytic MSIS model)
modelString <-"data {
  int N_Nyiro;
  int Y_Nyiro[N_Nyiro];
  real age_Nyiro[N_Nyiro];
  int M_Nyiro[N_Nyiro];
  int N_Aran;
  int Y_Aran[N_Aran];
  real age_Aran[N_Aran];
  int M_Aran[N_Aran];
  int N_Zhang;
  int Y_Zhang[N_Zhang];
  real age_Zhang[N_Zhang];
  int M_Zhang[N_Zhang];
  int N_Sastre;
  int Y_Sastre[N_Sastre];
  real age_Sastre[N_Sastre];
  int M_Sastre[N_Sastre];
  int N_Lu;
  int Y_Lu[N_Lu];
  real age_Lu[N_Lu];
  int M_Lu[N_Lu];
  int N_Leo;
  int Y_Leo[N_Leo];
  real age_Leo[N_Leo];
  int M_Leo[N_Leo];
  int N_Jen;
  int Y_Jen[N_Jen];
  real age_Jen[N_Jen];
  int M_Jen[N_Jen];
  int N_Mos;
  int Y_Mos[N_Mos];
  real age_Mos[N_Mos];
  int M_Mos[N_Mos];
  int N_Gol;
  int Y_Gol[N_Gol];
  real age_Gol[N_Gol];
  int M_Gol[N_Gol];
  int N_Mad;
  int Y_Mad[N_Mad];
  real age_Mad[N_Mad];
  int M_Mad[N_Mad];
}

parameters {
  real<lower=0> foi1_Nyiro;
  real<lower=0> foi2_Nyiro;
  real<lower=0> foi3_Nyiro;
  real<lower=0> foi1_Aran;
  real<lower=0> foi2_Aran;
  real<lower=0> foi3_Aran;
  real<lower=0> foi1_Zhang;
  real<lower=0> foi2_Zhang;
  real<lower=0> foi3_Zhang;
  real<lower=0, upper=5> decay;
  real<lower=0> foi1_Sastre;
  real<lower=0> foi2_Sastre;
  real<lower=0> foi3_Sastre;
  real<lower=0, upper=5> decay_F;
  real<lower=0> foi1_Lu;
  real<lower=0> foi2_Lu;
  real<lower=0> foi3_Lu;
  real<lower=0, upper=5> decay_N;
  real<lower=0> foi1_Leo;
  real<lower=0> foi2_Leo;
  real<lower=0> foi3_Leo;
  real<lower=0, upper=5> decay_Neu;
  real<lower=0> foi1_Jen;
  real<lower=0> foi2_Jen;
  real<lower=0> foi3_Jen;
  real<lower=0> foi1_Mos;
  real<lower=0> foi2_Mos;
  real<lower=0> foi3_Mos;
  real<lower=0> foi1_Gol;
  real<lower=0> foi2_Gol;
  real<lower=0> foi3_Gol;
  real<lower=0> foi1_Mad;
  real<lower=0> foi2_Mad;
  real<lower=0> foi3_Mad;
  real<lower=0, upper=5> decay_CF;
}

transformed parameters {
  real q_Nyiro[N_Nyiro];
  real q_Aran[N_Aran];
  real q_Zhang[N_Zhang];
  real q_Sastre[N_Sastre];
  real q_Lu[N_Lu];
  real q_Leo[N_Leo];
  real q_Jen[N_Jen];
  real q_Mos[N_Mos];
  real q_Gol[N_Gol];
  real q_Mad[N_Mad];
  
  for (n in 1:N_Nyiro)
    if (age_Nyiro[n] < 1){
      q_Nyiro[n] = (foi1_Nyiro / (foi1_Nyiro + decay)) * (1 - exp(-(foi1_Nyiro + decay)*(age_Nyiro[n]-0.5)));
    } else if (age_Nyiro[n] >= 1 && age_Nyiro[n] < 5) {
      q_Nyiro[n] = ((foi1_Nyiro / (foi1_Nyiro + decay)) * (1 - exp(-(foi1_Nyiro + decay)*(1-0.5))) - 
                (foi2_Nyiro / (foi2_Nyiro + decay))) * exp(-(foi2_Nyiro + decay)*(age_Nyiro[n] - 1)) +
        (foi2_Nyiro / (foi2_Nyiro + decay));
    } else {
      q_Nyiro[n] = ((((foi1_Nyiro / (foi1_Nyiro + decay)) * (1 - exp(-(foi1_Nyiro + decay)*(1-0.5))) - 
                  (foi2_Nyiro / (foi2_Nyiro + decay))) * exp(-(foi2_Nyiro + decay)*(5 - 1)) + 
                 decay/(foi3_Nyiro+decay) - decay/(foi2_Nyiro+decay))*exp(-(foi3_Nyiro + decay)*(age_Nyiro[n]-5)) + 
                foi3_Nyiro/(foi3_Nyiro+decay));
    }
    
    for (n in 1:N_Aran)
    if (age_Aran[n] < 1){
      q_Aran[n] = (foi1_Aran / (foi1_Aran + decay)) * (1 - exp(-(foi1_Aran + decay)*(age_Aran[n]-0.5)));
    } else if (age_Aran[n] >= 1 && age_Aran[n] < 5) {
      q_Aran[n] = ((foi1_Aran / (foi1_Aran + decay)) * (1 - exp(-(foi1_Aran + decay)*(1-0.5))) - 
                (foi2_Aran / (foi2_Aran + decay))) * exp(-(foi2_Aran + decay)*(age_Aran[n] - 1)) +
        (foi2_Aran / (foi2_Aran + decay));
    } else {
      q_Aran[n] = ((((foi1_Aran / (foi1_Aran + decay)) * (1 - exp(-(foi1_Aran+ decay)*(1-0.5))) - 
                  (foi2_Aran / (foi2_Aran + decay))) * exp(-(foi2_Aran + decay)*(5 - 1)) + 
                 decay/(foi3_Aran+decay) - decay/(foi2_Aran+decay))*exp(-(foi3_Aran + decay)*(age_Aran[n]-5)) + 
                foi3_Aran/(foi3_Aran+decay));
    }
    
    for (n in 1:N_Zhang)
    if (age_Zhang[n] < 1){
      q_Zhang[n] = (foi1_Zhang / (foi1_Zhang + decay)) * (1 - exp(-(foi1_Zhang + decay)*(age_Zhang[n]-0.5)));
    } else if (age_Zhang[n] >= 1 && age_Zhang[n] < 5) {
      q_Zhang[n] = ((foi1_Zhang / (foi1_Zhang + decay)) * (1 - exp(-(foi1_Zhang + decay)*(1-0.5))) - 
                (foi2_Zhang / (foi2_Zhang + decay))) * exp(-(foi2_Zhang + decay)*(age_Zhang[n] - 1)) +
        (foi2_Zhang / (foi2_Zhang + decay));
    } else {
      q_Zhang[n] = ((((foi1_Zhang / (foi1_Zhang + decay)) * (1 - exp(-(foi1_Zhang+ decay)*(1-0.5))) - 
                  (foi2_Zhang / (foi2_Zhang + decay))) * exp(-(foi2_Zhang + decay)*(5 - 1)) + 
                 decay/(foi3_Zhang+decay) - decay/(foi2_Zhang+decay))*exp(-(foi3_Zhang + decay)*(age_Zhang[n]-5)) + 
                foi3_Zhang/(foi3_Zhang+decay));
    }
    
    for (n in 1:N_Sastre)
    if (age_Sastre[n] < 1){
      q_Sastre[n] = (foi1_Sastre / (foi1_Sastre + decay_F)) * (1 - exp(-(foi1_Sastre + decay_F)*(age_Sastre[n]-0.5)));
    } else if (age_Sastre[n] >= 1 && age_Sastre[n] < 5) {
      q_Sastre[n] = ((foi1_Sastre / (foi1_Sastre + decay_F)) * (1 - exp(-(foi1_Sastre + decay_F)*(1-0.5))) - 
                (foi2_Sastre / (foi2_Sastre + decay_F))) * exp(-(foi2_Sastre + decay_F)*(age_Sastre[n] - 1)) +
        (foi2_Sastre / (foi2_Sastre + decay_F));
    } else {
      q_Sastre[n] = ((((foi1_Sastre / (foi1_Sastre + decay_F)) * (1 - exp(-(foi1_Sastre+ decay_F)*(1-0.5))) - 
                  (foi2_Sastre / (foi2_Sastre + decay_F))) * exp(-(foi2_Sastre + decay_F)*(5 - 1)) + 
                 decay_F/(foi3_Sastre+decay_F) - decay_F/(foi2_Sastre+decay_F))*exp(-(foi3_Sastre + decay_F)*(age_Sastre[n]-5)) + 
                foi3_Sastre/(foi3_Sastre+decay_F));
    }
    
    for (n in 1:N_Lu)
    if (age_Lu[n] < 1){
      q_Lu[n] = (foi1_Lu / (foi1_Lu + decay_N)) * (1 - exp(-(foi1_Lu + decay_N)*(age_Lu[n]-0.5)));
    } else if (age_Lu[n] >= 1 && age_Lu[n] < 5) {
      q_Lu[n] = ((foi1_Lu / (foi1_Lu + decay_N)) * (1 - exp(-(foi1_Lu + decay_N)*(1-0.5))) - 
                (foi2_Lu / (foi2_Lu + decay_N))) * exp(-(foi2_Lu + decay_N)*(age_Lu[n] - 1)) +
        (foi2_Lu / (foi2_Lu + decay_N));
    } else {
      q_Lu[n] = ((((foi1_Lu / (foi1_Lu + decay_N)) * (1 - exp(-(foi1_Lu+ decay_N)*(1-0.5))) - 
                  (foi2_Lu / (foi2_Lu + decay_N))) * exp(-(foi2_Lu + decay_N)*(5 - 1)) + 
                 decay_N/(foi3_Lu+decay_N) - decay_N/(foi2_Lu+decay_N))*exp(-(foi3_Lu + decay_N)*(age_Lu[n]-5)) + 
                foi3_Lu/(foi3_Lu+decay_N));
    }
    
     for (n in 1:N_Leo)
    if (age_Leo[n] < 1){
      q_Leo[n] = (foi1_Leo / (foi1_Leo + decay_Neu)) * (1 - exp(-(foi1_Leo + decay_Neu)*(age_Leo[n]-0.5)));
    } else if (age_Leo[n] >= 1 && age_Leo[n] < 5) {
      q_Leo[n] = ((foi1_Leo / (foi1_Leo + decay_Neu)) * (1 - exp(-(foi1_Leo + decay_Neu)*(1-0.5))) - 
                (foi2_Leo / (foi2_Leo + decay_Neu))) * exp(-(foi2_Leo + decay_Neu)*(age_Leo[n] - 1)) +
        (foi2_Leo / (foi2_Leo + decay_Neu));
    } else {
      q_Leo[n] = ((((foi1_Leo / (foi1_Leo + decay_Neu)) * (1 - exp(-(foi1_Leo+ decay_Neu)*(1-0.5))) - 
                  (foi2_Leo / (foi2_Leo + decay_Neu))) * exp(-(foi2_Leo + decay_Neu)*(5 - 1)) + 
                 decay_Neu/(foi3_Leo+decay_Neu) - decay_Neu/(foi2_Leo+decay_Neu))*exp(-(foi3_Leo + decay_Neu)*(age_Leo[n]-5)) + 
                foi3_Leo/(foi3_Leo+decay_Neu));
    }
    
    for (n in 1:N_Jen)
    if (age_Jen[n] < 1){
      q_Jen[n] = (foi1_Jen / (foi1_Jen + decay_CF)) * (1 - exp(-(foi1_Jen + decay_CF)*(age_Jen[n]-0.5)));
    } else if (age_Jen[n] >= 1 && age_Jen[n] < 5) {
      q_Jen[n] = ((foi1_Jen / (foi1_Jen + decay_CF)) * (1 - exp(-(foi1_Jen + decay_CF)*(1-0.5))) - 
                (foi2_Jen / (foi2_Jen + decay_CF))) * exp(-(foi2_Jen + decay_CF)*(age_Jen[n] - 1)) +
        (foi2_Jen / (foi2_Jen + decay_CF));
    } else {
      q_Jen[n] = ((((foi1_Jen / (foi1_Jen + decay_CF)) * (1 - exp(-(foi1_Jen+ decay_CF)*(1-0.5))) - 
                  (foi2_Jen / (foi2_Jen + decay_CF))) * exp(-(foi2_Jen + decay_CF)*(5 - 1)) + 
                 decay_CF/(foi3_Jen+decay_CF) - decay_CF/(foi2_Jen+decay_CF))*exp(-(foi3_Jen + decay_CF)*(age_Jen[n]-5)) + 
                foi3_Jen/(foi3_Jen+decay_CF));
    }
    
    for (n in 1:N_Mos)
    if (age_Mos[n] < 1){
      q_Mos[n] = (foi1_Mos / (foi1_Mos + decay_CF)) * (1 - exp(-(foi1_Mos + decay_CF)*(age_Mos[n]-0.5)));
    } else if (age_Mos[n] >= 1 && age_Mos[n] < 5) {
      q_Mos[n] = ((foi1_Mos / (foi1_Mos + decay_CF)) * (1 - exp(-(foi1_Mos + decay_CF)*(1-0.5))) - 
                (foi2_Mos / (foi2_Mos + decay_CF))) * exp(-(foi2_Mos + decay_CF)*(age_Mos[n] - 1)) +
        (foi2_Mos / (foi2_Mos + decay_CF));
    } else {
      q_Mos[n] = ((((foi1_Mos / (foi1_Mos + decay_CF)) * (1 - exp(-(foi1_Mos+ decay_CF)*(1-0.5))) - 
                  (foi2_Mos / (foi2_Mos + decay_CF))) * exp(-(foi2_Mos + decay_CF)*(5 - 1)) + 
                 decay_CF/(foi3_Mos+decay_CF) - decay_CF/(foi2_Mos+decay_CF))*exp(-(foi3_Mos + decay_CF)*(age_Mos[n]-5)) + 
                foi3_Mos/(foi3_Mos+decay_CF));
    }
    
    for (n in 1:N_Gol)
    if (age_Gol[n] < 1){
      q_Gol[n] = (foi1_Gol / (foi1_Gol+ decay_CF)) * (1 - exp(-(foi1_Gol + decay_CF)*(age_Gol[n]-0.5)));
    } else if (age_Gol[n] >= 1 && age_Gol[n] < 5) {
      q_Gol[n] = ((foi1_Gol / (foi1_Gol + decay_CF)) * (1 - exp(-(foi1_Gol + decay_CF)*(1-0.5))) - 
                (foi2_Gol / (foi2_Gol + decay_CF))) * exp(-(foi2_Gol + decay_CF)*(age_Gol[n] - 1)) +
        (foi2_Gol / (foi2_Gol + decay_CF));
    } else {
      q_Gol[n] = ((((foi1_Gol / (foi1_Gol + decay_CF)) * (1 - exp(-(foi1_Gol+ decay_CF)*(1-0.5))) - 
                  (foi2_Gol / (foi2_Gol + decay_CF))) * exp(-(foi2_Gol + decay_CF)*(5 - 1)) + 
                 decay_CF/(foi3_Gol+decay_CF) - decay_CF/(foi2_Gol+decay_CF))*exp(-(foi3_Gol + decay_CF)*(age_Gol[n]-5)) + 
                foi3_Gol/(foi3_Gol+decay_CF));
    }
    
    for (n in 1:N_Mad)
    if (age_Mad[n] < 1){
      q_Mad[n] = (foi1_Mad / (foi1_Mad+ decay_CF)) * (1 - exp(-(foi1_Mad + decay_CF)*(age_Mad[n]-0.5)));
    } else if (age_Mad[n] >= 1 && age_Mad[n] < 5) {
      q_Mad[n] = ((foi1_Mad / (foi1_Mad + decay_CF)) * (1 - exp(-(foi1_Mad + decay_CF)*(1-0.5))) - 
                (foi2_Mad / (foi2_Mad + decay_CF))) * exp(-(foi2_Mad + decay_CF)*(age_Mad[n] - 1)) +
        (foi2_Mad / (foi2_Mad + decay_CF));
    } else {
      q_Mad[n] = ((((foi1_Mad / (foi1_Mad + decay_CF)) * (1 - exp(-(foi1_Mad+ decay_CF)*(1-0.5))) - 
                  (foi2_Mad / (foi2_Mad + decay_CF))) * exp(-(foi2_Mad + decay_CF)*(5 - 1)) + 
                 decay_CF/(foi3_Mad+decay_CF) - decay_CF/(foi2_Mad+decay_CF))*exp(-(foi3_Mad + decay_CF)*(age_Mad[n]-5)) + 
                foi3_Mad/(foi3_Mad+decay_CF));
    }
    
}

model {
  foi1_Nyiro ~ gamma(64, 80);
  foi2_Nyiro ~ gamma(64, 40);
  foi3_Nyiro ~ gamma(64, 40);
  foi1_Aran ~ gamma(64, 80);
  foi2_Aran ~ gamma(64, 40);
  foi3_Aran ~ gamma(64, 40);
  foi1_Zhang ~ gamma(64, 80);
  foi2_Zhang ~ gamma(64, 40);
  foi3_Zhang ~ gamma(64, 40);
  
  foi1_Sastre ~ gamma(64, 80);
  foi2_Sastre ~ gamma(64, 40);
  foi3_Sastre ~ gamma(64, 40);
  
  foi1_Lu ~ gamma(64, 80);
  foi2_Lu ~ gamma(64, 40);
  foi3_Lu ~ gamma(64, 40);
  
  foi1_Leo ~ gamma(64, 80);
  foi2_Leo ~ gamma(64, 40);
  foi3_Leo ~ gamma(64, 40);
  
  foi1_Jen ~ gamma(64, 80);
  foi2_Jen ~ gamma(64, 40);
  foi3_Jen ~ gamma(64, 40);
  foi1_Mos ~ gamma(64, 80);
  foi2_Mos ~ gamma(64, 40);
  foi3_Mos ~ gamma(64, 40);
  foi1_Gol ~ gamma(64, 80);
  foi2_Gol ~ gamma(64, 40);
  foi3_Gol ~ gamma(64, 40);
  foi1_Mad ~ gamma(64, 80);
  foi2_Mad ~ gamma(64, 40);
  foi3_Mad ~ gamma(64, 40);
  
  for (n in 1:N_Nyiro)
    Y_Nyiro[n] ~ binomial(M_Nyiro[n], q_Nyiro[n]);
  for (n in 1:N_Aran)
    Y_Aran[n] ~ binomial(M_Aran[n], q_Aran[n]);
  for (n in 1:N_Zhang)
    Y_Zhang[n] ~ binomial(M_Zhang[n], q_Zhang[n]);
  
  for (n in 1:N_Sastre)
    Y_Sastre[n] ~ binomial(M_Sastre[n], q_Sastre[n]);
    
  for (n in 1:N_Lu)
    Y_Lu[n] ~ binomial(M_Lu[n], q_Lu[n]);  
    
  for (n in 1:N_Leo)
    Y_Leo[n] ~ binomial(M_Leo[n], q_Leo[n]);
  
  for (n in 1:N_Jen)
    Y_Jen[n] ~ binomial(M_Jen[n], q_Jen[n]);
  for (n in 1:N_Mos)
    Y_Mos[n] ~ binomial(M_Mos[n], q_Mos[n]);
  for (n in 1:N_Gol)
    Y_Gol[n] ~ binomial(M_Gol[n], q_Gol[n]);
  for (n in 1:N_Mad)
    Y_Mad[n] ~ binomial(M_Mad[n], q_Mad[n]);
}
"
# Read in processed data
datComb <- readRDS("cleanedData.RDS")

# Subset data by study
datNyiro <- subset(datComb, author == "datNyiro")
datAran <- subset(datComb, author == "datAran")
datZhang <- subset(datComb, author == "datZhang")
datSastre <- subset(datComb, author == "datSastre")
datLu <- subset(datComb, author == "datLu")
datLeo <- subset(datComb, author == "datLeo")
datJen <- subset(datComb, author == "datJen")
datMos <- subset(datComb, author == "datMos")
datGol <- subset(datComb, author == "datGol")
datMad <- subset(datComb, author == "datMad")

data <- list(N_Nyiro=nrow(datNyiro), Y_Nyiro=datNyiro$N_positive, age_Nyiro=datNyiro$ageMid, M_Nyiro=datNyiro$N,
             N_Aran=nrow(datAran), Y_Aran=datAran$N_positive, age_Aran=datAran$ageMid, M_Aran=datAran$N,
             N_Zhang=nrow(datZhang), Y_Zhang=datZhang$N_positive, age_Zhang=datZhang$ageMid, M_Zhang=datZhang$N,
             N_Sastre=nrow(datSastre), Y_Sastre=datSastre$N_positive, age_Sastre=datSastre$ageMid, M_Sastre=datSastre$N,
             N_Lu=nrow(datLu), Y_Lu=datLu$N_positive, age_Lu=datLu$ageMid, M_Lu=datLu$N,
             N_Leo=nrow(datLeo), Y_Leo=datLeo$N_positive, age_Leo=datLeo$ageMid, M_Leo=datLeo$N,
             N_Jen=nrow(datJen), Y_Jen=datJen$N_positive, age_Jen=datJen$ageMid, M_Jen=datJen$N,
             N_Mos=nrow(datMos), Y_Mos=datMos$N_positive, age_Mos=datMos$ageMid, M_Mos=datMos$N,
             N_Gol=nrow(datGol), Y_Gol=datGol$N_positive, age_Gol=datGol$ageMid, M_Gol=datGol$N,
             N_Mad=nrow(datMad), Y_Mad=datMad$N_positive, age_Mad=datMad$ageMid, M_Mad=datMad$N)

fit <- stan(model_code =modelString, data=data, seed=1234)

ggmcmc(ggs(fit, inc_warmup=TRUE, stan_include_auxiliar=TRUE),
       file='fit-traceplot.pdf', plot='traceplot')
ggmcmc(ggs(fit), file='fit-ggmcmc.pdf')

ms <- rstan::extract(fit)

# scatter plot of FoIs and decay (Nyiro)
d <- data.frame(foi1=ms$foi1_Nyiro, foi2=ms$foi2_Nyiro, foi3=ms$foi3_Nyiro, decay=ms$decay)
N_col <- ncol(d)
ggp <- ggpairs(d, upper='blank', diag='blank', lower='blank')

label_list <- list(foi1='foi1', foi2='foi2', foi3='foi3', decay='decay')
for(i in 1:N_col) {
  x <- d[,i]
  bw <- (max(x)-min(x))/10
  p <- ggplot(data.frame(x), aes(x))
  p <- p + theme_bw(base_size=14)
  p <- p + theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
  p <- p + geom_histogram(binwidth=bw, fill='white', color='grey5')
  p <- p + geom_line(eval(bquote(aes(y=..count..*.(bw)))), stat='density')
  p <- p + geom_label(data=data.frame(x=-Inf, y=Inf, label=label_list[[colnames(d)[i]]]), aes(x=x, y=y, label=label), hjust=0, vjust=1)
  ggp <- putPlot(ggp, p, i, i)
}

zcolat <- seq(-1, 1, length=81)
zcolre <- c(zcolat[1:40]+1, rev(zcolat[41:81]))

for(i in 1:(N_col-1)) {
  for(j in (i+1):N_col) {
    x <- as.numeric(d[,i])
    y <- as.numeric(d[,j])
    r <- cor(x, y, method='spearman', use='pairwise.complete.obs')
    zcol <- lattice::level.colors(r, at=zcolat, col.regions=grey(zcolre))
    textcol <- ifelse(abs(r) < 0.4, 'grey20', 'white')
    ell <- ellipse::ellipse(r, level=0.95, type='l', npoints=50, scale=c(.2, .2), centre=c(.5, .5))
    p <- ggplot(data.frame(ell), aes(x=x, y=y))
    p <- p + theme_bw() + theme(
      plot.background=element_blank(),
      panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
      panel.border=element_blank(), axis.ticks=element_blank()
    )
    p <- p + geom_polygon(fill=zcol, color=zcol)
    p <- p + geom_text(data=NULL, x=.5, y=.5, label=100*round(r, 2), size=6, col=textcol)
    ggp <- putPlot(ggp, p, i, j)
  }
}

for(j in 1:(N_col-1)) {
  for(i in (j+1):N_col) {
    x <- d[,j]
    y <- d[,i]
    p <- ggplot(data.frame(x, y), aes(x=x, y=y))
    p <- p + theme_bw(base_size=14)
    p <- p + theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
    p <- p + geom_hex()
    p <- p + scale_fill_gradientn(colours=gray.colors(7, start=0.1, end=0.9))
    ggp <- putPlot(ggp, p, i, j)
  }
}
ggp


# credible interval based on mcmc sample
data.frame.quantile.mcmc <- function(x, y_mcmc, probs=c(2.5, 50, 97.5)/100) {
  qua <- apply(y_mcmc, 2, quantile, probs=probs)
  d <- data.frame(X=x, t(qua))
  colnames(d) <- c('ageMid', paste0('p', probs*100))
  return(d)
}

# plot of observed seropositivity and prediction (ELISA)
N_mcmc <- length(ms$lp__)
N_X <- length(X_new)
X_new <- 1:80

q_Nyiro_mcmc <- as.data.frame(matrix(nrow=N_mcmc, ncol=N_X))
q_Aran_mcmc <- as.data.frame(matrix(nrow=N_mcmc, ncol=N_X))
q_Zhang_mcmc <- as.data.frame(matrix(nrow=N_mcmc, ncol=N_X))
for (i in 1:N_X) 
  if(X_new[i] < 1){
    q_Nyiro_mcmc[,i] <- (ms$foi1_Nyiro / (ms$foi1_Nyiro + ms$decay)) * (1 - exp(-(ms$foi1_Nyiro + ms$decay)*(X_new[i]-0.5)))
  } else if (X_new[i] >= 1 && X_new[i] < 5) {
    q_Nyiro_mcmc[,i] <- ((ms$foi1_Nyiro / (ms$foi1_Nyiro + ms$decay)) * (1 - exp(-(ms$foi1_Nyiro + ms$decay)*(1-0.5))) - 
                    (ms$foi2_Nyiro / (ms$foi2_Nyiro + ms$decay))) * exp(-(ms$foi2_Nyiro + ms$decay)*(X_new[i] - 1)) +
      (ms$foi2_Nyiro / (ms$foi2_Nyiro + ms$decay));
  } else {
    q_Nyiro_mcmc[,i] <- ((((ms$foi1_Nyiro / (ms$foi1_Nyiro + ms$decay)) * (1 - exp(-(ms$foi1_Nyiro + ms$decay)*(1-0.5))) - 
                      (ms$foi2_Nyiro / (ms$foi2_Nyiro + ms$decay))) * exp(-(ms$foi2_Nyiro + ms$decay)*(5 - 1)) + 
                     ms$decay/(ms$foi3_Nyiro+ms$decay) - ms$decay/(ms$foi2_Nyiro+ms$decay))*exp(-(ms$foi3_Nyiro + ms$decay)*(X_new[i]-5)) + 
                    ms$foi3_Nyiro/(ms$foi3_Nyiro+ms$decay));
  }
for (i in 1:N_X) 
  if(X_new[i] < 1){
    q_Aran_mcmc[,i] <- (ms$foi1_Aran / (ms$foi1_Aran + ms$decay)) * (1 - exp(-(ms$foi1_Aran + ms$decay)*(X_new[i]-0.5)))
  } else if (X_new[i] >= 1 && X_new[i] < 5) {
    q_Aran_mcmc[,i] <- ((ms$foi1_Aran / (ms$foi1_Aran + ms$decay)) * (1 - exp(-(ms$foi1_Aran + ms$decay)*(1-0.5))) - 
                           (ms$foi2_Aran / (ms$foi2_Aran + ms$decay))) * exp(-(ms$foi2_Aran + ms$decay)*(X_new[i] - 1)) +
      (ms$foi2_Aran / (ms$foi2_Aran + ms$decay));
  } else {
    q_Aran_mcmc[,i] <- ((((ms$foi1_Aran / (ms$foi1_Aran + ms$decay)) * (1 - exp(-(ms$foi1_Aran + ms$decay)*(1-0.5))) - 
                             (ms$foi2_Aran / (ms$foi2_Aran + ms$decay))) * exp(-(ms$foi2_Aran + ms$decay)*(5 - 1)) + 
                            ms$decay/(ms$foi3_Aran+ms$decay) - ms$decay/(ms$foi2_Aran+ms$decay))*exp(-(ms$foi3_Aran + ms$decay)*(X_new[i]-5)) + 
                           ms$foi3_Aran/(ms$foi3_Aran+ms$decay));
  }
for (i in 1:N_X) 
  if(X_new[i] < 1){
    q_Zhang_mcmc[,i] <- (ms$foi1_Zhang / (ms$foi1_Zhang + ms$decay)) * (1 - exp(-(ms$foi1_Zhang + ms$decay)*(X_new[i]-0.5)))
  } else if (X_new[i] >= 1 && X_new[i] < 5) {
    q_Zhang_mcmc[,i] <- ((ms$foi1_Zhang / (ms$foi1_Zhang + ms$decay)) * (1 - exp(-(ms$foi1_Zhang + ms$decay)*(1-0.5))) - 
                          (ms$foi2_Zhang / (ms$foi2_Zhang + ms$decay))) * exp(-(ms$foi2_Zhang + ms$decay)*(X_new[i] - 1)) +
      (ms$foi2_Zhang / (ms$foi2_Zhang + ms$decay));
  } else {
    q_Zhang_mcmc[,i] <- ((((ms$foi1_Zhang / (ms$foi1_Zhang + ms$decay)) * (1 - exp(-(ms$foi1_Zhang + ms$decay)*(1-0.5))) - 
                            (ms$foi2_Zhang / (ms$foi2_Zhang + ms$decay))) * exp(-(ms$foi2_Zhang + ms$decay)*(5 - 1)) + 
                           ms$decay/(ms$foi3_Zhang+ms$decay) - ms$decay/(ms$foi2_Zhang+ms$decay))*exp(-(ms$foi3_Zhang + ms$decay)*(X_new[i]-5)) + 
                          ms$foi3_Zhang/(ms$foi3_Zhang+ms$decay));
  }

d_Nyiro <- data.frame.quantile.mcmc(x=X_new, y_mcmc=q_Nyiro_mcmc )
d_Aran <- data.frame.quantile.mcmc(x=X_new, y_mcmc=q_Aran_mcmc )
d_Zhang <- data.frame.quantile.mcmc(x=X_new, y_mcmc=q_Zhang_mcmc )
p <- ggplot(d_Nyiro, aes(x = ageMid, y = p50)) +
  geom_ribbon(alpha = 0.2,  aes(ymin = p2.5, ymax = p97.5), fill = "#558C8C")+
  geom_line(color = "#558C8C")+
  geom_point(data = datNyiro, aes(x = ageMid, y = mean), color = "#558C8C")+
  geom_linerange(data = datNyiro, aes(y = mean, ymin = lower, ymax = upper), color = "#558C8C") +
  geom_ribbon(data= d_Aran, alpha = 0.2,  aes(x=ageMid, y= p50, ymin = p2.5, ymax = p97.5), fill = "#C05746")+
  geom_line(data=d_Aran, color = "#C05746")+
  geom_point(data = datAran, aes(x = ageMid, y = mean), color = "#C05746")+
  geom_linerange(data = datAran, aes(y = mean, ymin = lower, ymax = upper), color = "#C05746") +
  geom_ribbon(data= d_Zhang, alpha = 0.2,  aes(x=ageMid, y= p50, ymin = p2.5, ymax = p97.5), fill = "#075E9D")+
  geom_line(data=d_Zhang, color = "#075E9D")+
  geom_point(data = datZhang, aes(x = ageMid, y = mean), color = "#075E9D")+
  geom_linerange(data = datZhang, aes(y = mean, ymin = lower, ymax = upper), color = "#075E9D") +
  xlab("Age (years)") + ylab("Proportion Seropositive") +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  ylim(0, 1) +
  ggtitle("ELISA") +
  theme_classic()
p

# plot of observed seropositivity and prediction (F)
X_new_f <- 1:90
N_X_f <- length(X_new_f)

q_Sastre_mcmc <- as.data.frame(matrix(nrow=N_mcmc, ncol=N_X_f))
for (i in 1:N_X_f) 
  if(X_new_f[i] < 1){
    q_Sastre_mcmc[,i] <- (ms$foi1_Sastre / (ms$foi1_Sastre + ms$decay_F)) * (1 - exp(-(ms$foi1_Sastre + ms$decay_F)*(X_new_f[i]-0.5)))
  } else if (X_new_f[i] >= 1 && X_new_f[i] < 5) {
    q_Sastre_mcmc[,i] <- ((ms$foi1_Sastre / (ms$foi1_Sastre + ms$decay_F)) * (1 - exp(-(ms$foi1_Sastre + ms$decay_F)*(1-0.5))) - 
                           (ms$foi2_Sastre / (ms$foi2_Sastre + ms$decay_F))) * exp(-(ms$foi2_Sastre + ms$decay_F)*(X_new_f[i] - 1)) +
      (ms$foi2_Sastre / (ms$foi2_Sastre + ms$decay_F));
  } else {
    q_Sastre_mcmc[,i] <- ((((ms$foi1_Sastre / (ms$foi1_Sastre + ms$decay_F)) * (1 - exp(-(ms$foi1_Sastre + ms$decay_F)*(1-0.5))) - 
                             (ms$foi2_Sastre / (ms$foi2_Sastre + ms$decay_F))) * exp(-(ms$foi2_Sastre + ms$decay_F)*(5 - 1)) + 
                            ms$decay_F/(ms$foi3_Sastre+ms$decay_F) - ms$decay_F/(ms$foi2_Sastre+ms$decay_F))*exp(-(ms$foi3_Sastre + ms$decay_F)*(X_new_f[i]-5)) + 
                           ms$foi3_Sastre/(ms$foi3_Sastre+ms$decay_F));
  }


d_Sastre <- data.frame.quantile.mcmc(x=X_new_f, y_mcmc=q_Sastre_mcmc )
p_F <- ggplot(d_Sastre, aes(x = ageMid, y = p50)) +
  geom_ribbon(alpha = 0.2,  aes(ymin = p2.5, ymax = p97.5), fill = "#558C8C")+
  geom_line(color = "#558C8C")+
  geom_point(data = datSastre, aes(x = ageMid, y = mean), color = "#558C8C")+
  geom_linerange(data = datSastre, aes(y = mean, ymin = lower, ymax = upper), color = "#558C8C") +
  xlab("Age (years)") + ylab("Proportion Seropositive") +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  ylim(0, 1) +
  ggtitle("ELISA_F") +
  theme_classic()
p_F

# plot of observed seropositivity and prediction (N)
q_Lu_mcmc <- as.data.frame(matrix(nrow=N_mcmc, ncol=N_X))
for (i in 1:N_X) 
  if(X_new[i] < 1){
    q_Lu_mcmc[,i] <- (ms$foi1_Lu / (ms$foi1_Lu + ms$decay_N)) * (1 - exp(-(ms$foi1_Lu + ms$decay_N)*(X_new[i]-0.5)))
  } else if (X_new[i] >= 1 && X_new[i] < 5) {
    q_Lu_mcmc[,i] <- ((ms$foi1_Lu / (ms$foi1_Lu + ms$decay_N)) * (1 - exp(-(ms$foi1_Lu + ms$decay_N)*(1-0.5))) - 
                           (ms$foi2_Lu / (ms$foi2_Lu + ms$decay_N))) * exp(-(ms$foi2_Lu + ms$decay_N)*(X_new[i] - 1)) +
      (ms$foi2_Lu / (ms$foi2_Lu + ms$decay_N));
  } else {
    q_Lu_mcmc[,i] <- ((((ms$foi1_Lu / (ms$foi1_Lu + ms$decay_N)) * (1 - exp(-(ms$foi1_Lu + ms$decay_N)*(1-0.5))) - 
                             (ms$foi2_Lu / (ms$foi2_Lu + ms$decay_N))) * exp(-(ms$foi2_Lu + ms$decay_N)*(5 - 1)) + 
                            ms$decay_N/(ms$foi3_Lu +ms$decay_N) - ms$decay_N/(ms$foi2_Lu+ms$decay_N))*exp(-(ms$foi3_Lu + ms$decay_N)*(X_new[i]-5)) + 
                           ms$foi3_Lu/(ms$foi3_Lu+ms$decay_N));
  }

d_Lu <- data.frame.quantile.mcmc(x=X_new, y_mcmc=q_Lu_mcmc )
p_N <- ggplot(d_Lu, aes(x = ageMid, y = p50)) +
  geom_ribbon(alpha = 0.2,  aes(ymin = p2.5, ymax = p97.5), fill = "#558C8C")+
  geom_line(color = "#558C8C")+
  geom_point(data = datLu, aes(x = ageMid, y = mean), color = "#558C8C")+
  geom_linerange(data = datLu, aes(y = mean, ymin = lower, ymax = upper), color = "#558C8C") +
  xlab("Age (years)") + ylab("Proportion Seropositive") +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  ylim(0, 1) +
  ggtitle("ELISA_N") +
  theme_classic()
p_N

p + p_F + p_N


# plot of observed seropositivity and prediction (Neutralization)
N_mcmc <- length(ms$lp__)
X_new <- 1:80
N_X <- length(X_new)

q_Leo_mcmc <- as.data.frame(matrix(nrow=N_mcmc, ncol=N_X))
for (i in 1:N_X) 
  if(X_new[i] < 1){
    q_Leo_mcmc[,i] <- (ms$foi1_Leo / (ms$foi1_Leo + ms$decay_Neu)) * (1 - exp(-(ms$foi1_Leo + ms$decay_Neu)*(X_new[i]-0.5)))
  } else if (X_new[i] >= 1 && X_new[i] < 5) {
    q_Leo_mcmc[,i] <- ((ms$foi1_Leo / (ms$foi1_Leo + ms$decay_Neu)) * (1 - exp(-(ms$foi1_Leo + ms$decay_Neu)*(1-0.5))) - 
                        (ms$foi2_Leo / (ms$foi2_Leo + ms$decay_Neu))) * exp(-(ms$foi2_Leo + ms$decay_Neu)*(X_new[i] - 1)) +
      (ms$foi2_Leo / (ms$foi2_Leo + ms$decay_Neu));
  } else {
    q_Leo_mcmc[,i] <- ((((ms$foi1_Leo / (ms$foi1_Leo + ms$decay_Neu)) * (1 - exp(-(ms$foi1_Leo + ms$decay_Neu)*(1-0.5))) - 
                          (ms$foi2_Leo / (ms$foi2_Leo + ms$decay_Neu))) * exp(-(ms$foi2_Leo + ms$decay_Neu)*(5 - 1)) + 
                         ms$decay_Neu/(ms$foi3_Leo +ms$decay_Neu) - ms$decay_Neu/(ms$foi2_Leo+ms$decay_Neu))*exp(-(ms$foi3_Leo + ms$decay_Neu)*(X_new[i]-5)) + 
                        ms$foi3_Leo/(ms$foi3_Leo+ms$decay_Neu));
  }

d_Leo <- data.frame.quantile.mcmc(x=X_new, y_mcmc=q_Leo_mcmc )
p_Neu <- ggplot(d_Leo, aes(x = ageMid, y = p50)) +
  geom_ribbon(alpha = 0.2,  aes(ymin = p2.5, ymax = p97.5), fill = "#558C8C")+
  geom_line(color = "#558C8C")+
  geom_point(data = datLeo, aes(x = ageMid, y = mean), color = "#558C8C")+
  geom_linerange(data = datLeo, aes(y = mean, ymin = lower, ymax = upper), color = "#558C8C") +
  xlab("Age (years)") + ylab("Proportion Seropositive") +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  ylim(0, 1) +
  ggtitle("Neutralization") +
  theme_classic()
p_Neu

# plot of observed seropositivity and prediction (CF)
q_Jen_mcmc <- as.data.frame(matrix(nrow=N_mcmc, ncol=N_X))
q_Mos_mcmc <- as.data.frame(matrix(nrow=N_mcmc, ncol=N_X))
q_Gol_mcmc <- as.data.frame(matrix(nrow=N_mcmc, ncol=N_X))
q_Mad_mcmc <- as.data.frame(matrix(nrow=N_mcmc, ncol=N_X))
for (i in 1:N_X) 
  if(X_new[i] < 1){
    q_Jen_mcmc[,i] <- (ms$foi1_Jen / (ms$foi1_Jen + ms$decay_CF)) * (1 - exp(-(ms$foi1_Jen + ms$decay_CF)*(X_new[i]-0.5)))
  } else if (X_new[i] >= 1 && X_new[i] < 5) {
    q_Jen_mcmc[,i] <- ((ms$foi1_Jen / (ms$foi1_Jen + ms$decay_CF)) * (1 - exp(-(ms$foi1_Jen + ms$decay_CF)*(1-0.5))) - 
                           (ms$foi2_Jen / (ms$foi2_Jen + ms$decay_CF))) * exp(-(ms$foi2_Jen + ms$decay_CF)*(X_new[i] - 1)) +
      (ms$foi2_Jen / (ms$foi2_Jen + ms$decay_CF));
  } else {
    q_Jen_mcmc[,i] <- ((((ms$foi1_Jen / (ms$foi1_Jen + ms$decay_CF)) * (1 - exp(-(ms$foi1_Jen + ms$decay_CF)*(1-0.5))) - 
                             (ms$foi2_Jen / (ms$foi2_Jen + ms$decay_CF))) * exp(-(ms$foi2_Jen + ms$decay_CF)*(5 - 1)) + 
                            ms$decay_CF/(ms$foi3_Jen+ms$decay_CF) - ms$decay_CF/(ms$foi2_Jen+ms$decay_CF))*exp(-(ms$foi3_Jen + ms$decay_CF)*(X_new[i]-5)) + 
                           ms$foi3_Jen/(ms$foi3_Jen+ms$decay_CF));
  }
for (i in 1:N_X) 
  if(X_new[i] < 1){
    q_Mos_mcmc[,i] <- (ms$foi1_Mos / (ms$foi1_Mos + ms$decay_CF)) * (1 - exp(-(ms$foi1_Mos + ms$decay_CF)*(X_new[i]-0.5)))
  } else if (X_new[i] >= 1 && X_new[i] < 5) {
    q_Mos_mcmc[,i] <- ((ms$foi1_Mos / (ms$foi1_Mos + ms$decay_CF)) * (1 - exp(-(ms$foi1_Mos + ms$decay_CF)*(1-0.5))) - 
                         (ms$foi2_Mos / (ms$foi2_Mos + ms$decay_CF))) * exp(-(ms$foi2_Mos + ms$decay_CF)*(X_new[i] - 1)) +
      (ms$foi2_Mos / (ms$foi2_Mos + ms$decay_CF));
  } else {
    q_Mos_mcmc[,i] <- ((((ms$foi1_Mos / (ms$foi1_Mos + ms$decay_CF)) * (1 - exp(-(ms$foi1_Mos + ms$decay_CF)*(1-0.5))) - 
                           (ms$foi2_Mos / (ms$foi2_Mos + ms$decay_CF))) * exp(-(ms$foi2_Mos + ms$decay_CF)*(5 - 1)) + 
                          ms$decay_CF/(ms$foi3_Mos+ms$decay_CF) - ms$decay_CF/(ms$foi2_Mos+ms$decay_CF))*exp(-(ms$foi3_Mos + ms$decay_CF)*(X_new[i]-5)) + 
                         ms$foi3_Mos/(ms$foi3_Mos+ms$decay_CF));
  }
for (i in 1:N_X) 
  if(X_new[i] < 1){
    q_Gol_mcmc[,i] <- (ms$foi1_Gol / (ms$foi1_Gol + ms$decay_CF)) * (1 - exp(-(ms$foi1_Gol + ms$decay_CF)*(X_new[i]-0.5)))
  } else if (X_new[i] >= 1 && X_new[i] < 5) {
    q_Gol_mcmc[,i] <- ((ms$foi1_Gol / (ms$foi1_Gol + ms$decay_CF)) * (1 - exp(-(ms$foi1_Gol + ms$decay_CF)*(1-0.5))) - 
                         (ms$foi2_Gol / (ms$foi2_Gol + ms$decay_CF))) * exp(-(ms$foi2_Gol + ms$decay_CF)*(X_new[i] - 1)) +
      (ms$foi2_Gol / (ms$foi2_Gol + ms$decay_CF));
  } else {
    q_Gol_mcmc[,i] <- ((((ms$foi1_Gol / (ms$foi1_Gol + ms$decay_CF)) * (1 - exp(-(ms$foi1_Gol + ms$decay_CF)*(1-0.5))) - 
                           (ms$foi2_Gol / (ms$foi2_Gol + ms$decay_CF))) * exp(-(ms$foi2_Gol + ms$decay_CF)*(5 - 1)) + 
                          ms$decay_CF/(ms$foi3_Gol+ms$decay_CF) - ms$decay_CF/(ms$foi2_Gol+ms$decay_CF))*exp(-(ms$foi3_Gol + ms$decay_CF)*(X_new[i]-5)) + 
                         ms$foi3_Gol/(ms$foi3_Gol+ms$decay_CF));
  }
for (i in 1:N_X) 
  if(X_new[i] < 1){
    q_Mad_mcmc[,i] <- (ms$foi1_Mad / (ms$foi1_Mad + ms$decay_CF)) * (1 - exp(-(ms$foi1_Mad + ms$decay_CF)*(X_new[i]-0.5)))
  } else if (X_new[i] >= 1 && X_new[i] < 5) {
    q_Mad_mcmc[,i] <- ((ms$foi1_Mad / (ms$foi1_Mad + ms$decay_CF)) * (1 - exp(-(ms$foi1_Mad + ms$decay_CF)*(1-0.5))) - 
                         (ms$foi2_Mad / (ms$foi2_Mad + ms$decay_CF))) * exp(-(ms$foi2_Mad + ms$decay_CF)*(X_new[i] - 1)) +
      (ms$foi2_Mad / (ms$foi2_Mad + ms$decay_CF));
  } else {
    q_Mad_mcmc[,i] <- ((((ms$foi1_Mad / (ms$foi1_Mad + ms$decay_CF)) * (1 - exp(-(ms$foi1_Mad + ms$decay_CF)*(1-0.5))) - 
                           (ms$foi2_Mad / (ms$foi2_Mad + ms$decay_CF))) * exp(-(ms$foi2_Mad + ms$decay_CF)*(5 - 1)) + 
                          ms$decay_CF/(ms$foi3_Mad+ms$decay_CF) - ms$decay_CF/(ms$foi2_Mad+ms$decay_CF))*exp(-(ms$foi3_Mad + ms$decay_CF)*(X_new[i]-5)) + 
                         ms$foi3_Mad/(ms$foi3_Mad+ms$decay_CF));
  }

d_Jen <- data.frame.quantile.mcmc(x=X_new, y_mcmc=q_Jen_mcmc )
d_Mos <- data.frame.quantile.mcmc(x=X_new, y_mcmc=q_Mos_mcmc )
d_Gol <- data.frame.quantile.mcmc(x=X_new, y_mcmc=q_Gol_mcmc )
d_Mad <- data.frame.quantile.mcmc(x=X_new, y_mcmc=q_Mad_mcmc )
p_CF <- ggplot(d_Jen, aes(x = ageMid, y = p50)) +
  geom_ribbon(alpha = 0.2,  aes(ymin = p2.5, ymax = p97.5), fill = "#558C8C")+
  geom_line(color = "#558C8C")+
  geom_point(data = datJen, aes(x = ageMid, y = mean), color = "#558C8C")+
  geom_linerange(data = datJen, aes(y = mean, ymin = lower, ymax = upper), color = "#558C8C") +
  geom_ribbon(data= d_Mos, alpha = 0.2,  aes(x=ageMid, y= p50, ymin = p2.5, ymax = p97.5), fill = "#C05746")+
  geom_line(data=d_Mos, color = "#C05746")+
  geom_point(data = datMos, aes(x = ageMid, y = mean), color = "#C05746")+
  geom_linerange(data = datMos, aes(y = mean, ymin = lower, ymax = upper), color = "#C05746") +
  geom_ribbon(data= d_Gol, alpha = 0.2,  aes(x=ageMid, y= p50, ymin = p2.5, ymax = p97.5), fill = "#075E9D")+
  geom_line(data=d_Gol, color = "#075E9D")+
  geom_point(data = datGol, aes(x = ageMid, y = mean), color = "#075E9D")+
  geom_linerange(data = datGol, aes(y = mean, ymin = lower, ymax = upper), color = "#075E9D") +
  geom_ribbon(data= d_Mad, alpha = 0.2,  aes(x=ageMid, y= p50, ymin = p2.5, ymax = p97.5), fill = "#A03E99")+
  geom_line(data=d_Mad, color = "#A03E99")+
  geom_point(data = datMad, aes(x = ageMid, y = mean), color = "#A03E99")+
  geom_linerange(data = datMad, aes(y = mean, ymin = lower, ymax = upper), color = "#A03E99") +
  xlab("Age (years)") + ylab("Proportion Seropositive") +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  ylim(0, 1) +
  ggtitle("CF") +
  theme_classic()
p_CF

p_Neu + p_CF

##plot of mcmc samples (foi1 of ELISA)
N_mcmc <- length(ms$lp__)
d_est <- data.frame(1:N_mcmc, ms$foi1_Nyiro, ms$foi1_Aran, ms$foi1_Zhang)
colnames(d_est) <- c('mcmc', paste0('foi1', 1:3))

d_mode <- data.frame(t(apply(d_est[-1], 2, function(x) {
  dens <- density(x)
  mode_i <- which.max(dens$y)
  mode_x <- dens$x[mode_i]
  mode_y <- dens$y[mode_i]
  c(mode_x, mode_y)
})))
colnames(d_mode) <- c('X', 'Y')
d_melt <- tidyr::pivot_longer(d_est, cols=!mcmc, names_to = "foi1", values_to = "value")

p_1 <- ggplot() +
  theme_bw(base_size=18)+ 
  geom_density(data=d_melt, aes(x=value, fill=foi1), alpha=0.15)+
  geom_segment(data=d_mode, aes(x=X, xend=X, y=Y, yend=0), color='black', linetype='dashed', alpha=0.6)+
  geom_rug(data=d_mode, aes(x=X), sides='b')+ 
  labs(x='value', y='density')+ 
  scale_x_continuous(breaks=seq(from=0, to=3, by=1))+ 
  scale_fill_discrete(labels=c("Nyiro", "Arankalle", "Zhang")) +
  guides(fill = FALSE)+
  theme_classic()
p_1 

##plot of mcmc samples (foi2 of ELISA)
d_est <- data.frame(1:N_mcmc, ms$foi2_Nyiro, ms$foi2_Aran, ms$foi2_Zhang)
colnames(d_est) <- c('mcmc', paste0('foi2', 1:3))

d_mode <- data.frame(t(apply(d_est[-1], 2, function(x) {
  dens <- density(x)
  mode_i <- which.max(dens$y)
  mode_x <- dens$x[mode_i]
  mode_y <- dens$y[mode_i]
  c(mode_x, mode_y)
})))
colnames(d_mode) <- c('X', 'Y')
d_melt <- tidyr::pivot_longer(d_est, cols=!mcmc, names_to = "foi2", values_to = "value")

p_2 <- ggplot() +
  theme_bw(base_size=18)+ 
  geom_density(data=d_melt, aes(x=value, fill=foi2), alpha=0.15)+
  geom_segment(data=d_mode, aes(x=X, xend=X, y=Y, yend=0), color='black', linetype='dashed', alpha=0.6)+
  geom_rug(data=d_mode, aes(x=X), sides='b')+ 
  labs(x='value', y='density')+ 
  scale_x_continuous(breaks=seq(from=0, to=3, by=1))+ 
  scale_fill_discrete(labels=c("Nyiro", "Arankalle", "Zhang")) +
  guides(fill = guide_legend(title = NULL))+
  theme_classic()+
  theme(legend.position = "bottom")
p_2 


##plot of mcmc samples (foi3 of ELISA)
d_est <- data.frame(1:N_mcmc, ms$foi3_Nyiro, ms$foi3_Aran, ms$foi3_Zhang)
colnames(d_est) <- c('mcmc', paste0('foi3', 1:3))

d_mode <- data.frame(t(apply(d_est[-1], 2, function(x) {
  dens <- density(x)
  mode_i <- which.max(dens$y)
  mode_x <- dens$x[mode_i]
  mode_y <- dens$y[mode_i]
  c(mode_x, mode_y)
})))
colnames(d_mode) <- c('X', 'Y')
d_melt <- tidyr::pivot_longer(d_est, cols=!mcmc, names_to = "foi3", values_to = "value")

p_3 <- ggplot() +
  theme_bw(base_size=18)+ 
  geom_density(data=d_melt, aes(x=value, fill=foi3), alpha=0.15)+
  geom_segment(data=d_mode, aes(x=X, xend=X, y=Y, yend=0), color='black', linetype='dashed', alpha=0.6)+
  geom_rug(data=d_mode, aes(x=X), sides='b')+ 
  labs(x='value', y='density')+ 
  scale_x_continuous(breaks=seq(from=0, to=3, by=1))+ 
  scale_fill_discrete(labels=c("Nyiro", "Arankalle", "Zhang")) +
  guides(fill = FALSE)+
  theme_classic()
p_3 

p_1 + p_2 + p_3

