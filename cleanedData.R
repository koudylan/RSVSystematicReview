## Cleaning seroprevalence data
## Removing under aged 0.5 year

library(tidyverse)
library(binom)

dat <- read_csv("Data_seroprevalence.csv")

######################################
### ELISA IgG
######################################

###################
### Nyiro
###################

datNyiro <- dat %>%
  filter(PMID == 28531224) %>%
  select(ageStart,ageEnd,ageMid, N,N_positive) %>% 
  mutate(author = "datNyiro")

datNyiro <- datNyiro %>%
  mutate(seropostive = N_positive/N) %>% 
  filter(ageStart > 0.5)

# Adding confidence intervals to serosamples
datNyiro[,c("mean","lower","upper")] <- binom.confint(datNyiro$N_positive, datNyiro$N,method="exact")[,c("mean","lower","upper")]

###################
### Arankalle
###################

datAran <- dat %>%
  filter(PMID == 31012488) %>%
  select(ageStart,ageEnd,ageMid,N,N_positive)

datAran <- datAran %>%
  mutate(seropostive = N_positive/N) %>% 
  mutate(author = "datAran") %>% 
  filter(ageStart > 0.5)

# Adding confidence intervals to serosamples
datAran[,c("mean","lower","upper")] <- binom.confint(datAran$N_positive, datAran$N,method="exact")[,c("mean","lower","upper")]

###################
### Arankalle
###################

datZhang <- dat %>%
  filter(PMID == 19080178) %>%
  select(ageStart,ageEnd,ageMid,N,N_positive)

datZhang <- datZhang %>%
  mutate(seropostive = N_positive/N) %>% 
  mutate(author = "datZhang") %>% 
  filter(ageStart > 0.5)

# Adding confidence intervals to serosamples
datZhang[,c("mean","lower","upper")] <- binom.confint(datZhang$N_positive, datZhang$N,method="exact")[,c("mean","lower","upper")]

######################################
### ELISA F
######################################

###################
### Sastre
###################

datSastre<- dat %>%
  filter(PMID == 22748150) %>%
  select(ageStart,ageEnd,ageMid, N,N_positive) %>% 
  mutate(author = "datSastre")

datSastre <- datSastre %>%
  mutate(seropostive = N_positive/N) %>% 
  filter(ageStart > 0.5)

# Adding confidence intervals to serosamples
datSastre[,c("mean","lower","upper")] <- binom.confint(datSastre$N_positive, datSastre$N,method="exact")[,c("mean","lower","upper")]

######################################
### ELISA N
######################################

###################
### Lu
###################

datLu<- dat %>%
  filter(PMID == 21310026) %>%
  select(ageStart,ageEnd,ageMid, N,N_positive) %>% 
  mutate(author = "datLu")

datLu <- datLu %>%
  mutate(seropostive = N_positive/N) %>% 
  filter(ageStart > 0.5)

# Adding confidence intervals to serosamples
datLu[,c("mean","lower","upper")] <- binom.confint(datLu$N_positive, datLu$N,method="exact")[,c("mean","lower","upper")]

######################################
### Neutralizatioin
######################################

###################
### Leogrande
###################

datLeo<- dat %>%
  filter(PMID == 1335108) %>%
  select(ageStart,ageEnd,ageMid, N,N_positive) %>% 
  mutate(author = "datLeo")

datLeo <- datLeo %>%
  mutate(seropostive = N_positive/N) %>% 
  filter(ageStart > 0.5)

# Adding confidence intervals to serosamples
datLeo[,c("mean","lower","upper")] <- binom.confint(datLeo$N_positive, datLeo$N,method="exact")[,c("mean","lower","upper")]

######################################
### CF
######################################

###################
### Jennings
###################

datJen<- dat %>%
  filter(PMID == 4341999) %>%
  select(ageStart,ageEnd,ageMid, N,N_positive) %>% 
  mutate(author = "datJen")

datJen <- datJen %>%
  mutate(seropostive = N_positive/N) %>% 
  filter(ageStart > 0.5)

# Adding confidence intervals to serosamples
datJen[,c("mean","lower","upper")] <- binom.confint(datJen$N_positive, datJen$N,method="exact")[,c("mean","lower","upper")]

###################
### Moss
###################

datMos<- dat %>%
  filter(PMID == 13936211) %>%
  select(ageStart,ageEnd,ageMid, N,N_positive) %>% 
  mutate(author = "datMos")

datMos <- datMos %>%
  mutate(seropostive = N_positive/N) %>% 
  filter(ageStart > 0.5)

# Adding confidence intervals to serosamples
datMos[,c("mean","lower","upper")] <- binom.confint(datMos$N_positive, datMos$N,method="exact")[,c("mean","lower","upper")]

###################
### Golubjatnikov
###################

datGol<- dat %>%
  filter(PMID == 165720) %>%
  select(ageStart,ageEnd,ageMid, N,N_positive) %>% 
  mutate(author = "datGol")

datGol <- datGol %>%
  mutate(seropostive = N_positive/N) %>% 
  filter(ageStart > 0.5)

# Adding confidence intervals to serosamples
datGol[,c("mean","lower","upper")] <- binom.confint(datGol$N_positive, datGol$N,method="exact")[,c("mean","lower","upper")]

###################
### Madhavan
###################

datMad<- dat %>%
  filter(PMID == 4468943) %>%
  select(ageStart,ageEnd,ageMid, N,N_positive) %>% 
  mutate(author = "datMad")

datMad <- datMad %>%
  mutate(seropostive = N_positive/N) %>% 
  filter(ageStart > 0.5)

# Adding confidence intervals to serosamples
datMad[,c("mean","lower","upper")] <- binom.confint(datMad$N_positive, datMad$N,method="exact")[,c("mean","lower","upper")]


###################
###################

datComb <- datNyiro %>% 
  bind_rows(datAran) %>% 
  bind_rows(datZhang) %>% 
  bind_rows(datSastre) %>% 
  bind_rows(datLu) %>% 
  bind_rows(datLeo) %>% 
  bind_rows(datJen) %>% 
  bind_rows(datMos) %>% 
  bind_rows(datGol) %>% 
  bind_rows(datMad) 

saveRDS(datComb, "cleanedData.RDS")
