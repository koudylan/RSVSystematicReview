# estimation of foi (ELISA/F/N/Neutr/CF) using non-informative prior
rm(list=ls(all=TRUE))
library(tidyverse)
library(rstan)
library(ggmcmc)
library(ggplot2)
library(patchwork)

# model
modelString <- "data {
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
  real<lower=0, upper=5> foi1_Nyiro;
  real<lower=0, upper=5> foi2_Nyiro;
  real<lower=0, upper=5> foi3_Nyiro;
  real<lower=0, upper=5> foi1_Aran;
  real<lower=0, upper=5> foi2_Aran;
  real<lower=0, upper=5> foi3_Aran;
  real<lower=0, upper=5> foi1_Zhang;
  real<lower=0, upper=5> foi2_Zhang;
  real<lower=0, upper=5> foi3_Zhang;
  real<lower=0, upper=5> decay;
  real<lower=0, upper=5> foi1_Sastre;
  real<lower=0, upper=5> foi2_Sastre;
  real<lower=0, upper=5> foi3_Sastre;
  real<lower=0, upper=5> decay_F;
  real<lower=0, upper=5> foi1_Lu;
  real<lower=0, upper=5> foi2_Lu;
  real<lower=0, upper=5> foi3_Lu;
  real<lower=0, upper=5> decay_N;
  real<lower=0, upper=5> foi1_Leo;
  real<lower=0, upper=5> foi2_Leo;
  real<lower=0, upper=5> foi3_Leo;
  real<lower=0, upper=5> decay_Neu;
  real<lower=0, upper=5> foi1_Jen;
  real<lower=0, upper=5> foi2_Jen;
  real<lower=0, upper=5> foi3_Jen;
  real<lower=0, upper=5> foi1_Mos;
  real<lower=0, upper=5> foi2_Mos;
  real<lower=0, upper=5> foi3_Mos;
  real<lower=0, upper=5> foi1_Gol;
  real<lower=0, upper=5> foi2_Gol;
  real<lower=0, upper=5> foi3_Gol;
  real<lower=0, upper=5> foi1_Mad;
  real<lower=0, upper=5> foi2_Mad;
  real<lower=0, upper=5> foi3_Mad;
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
}"

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

fit <- stan(model_code = modelString, data=data, seed=1234)

