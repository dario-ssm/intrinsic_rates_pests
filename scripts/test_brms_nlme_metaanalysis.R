
#### Test meta-analysis nonlinear function ####

# - SCRIPT INFO ----
#     Authors: Dario San Segundo Molina, Sara Villen P?rez, Ignacio Morales Castilla
#     Title: brms & nlme for hierarchical meta-analysis
#     Aim: test model fitting to pooled studies to obtain ecologically-informative parameters and to assess heterogeneity
#     Date: January 2022
#     

# - 1. Prepare script and load dataset of simulated data from our database ----
rm(list=ls())
library(tidyverse)
library(nlme)
library(brms)
#load data without controversial papers:

## NOTE THAT THE FINAL SIMULATIONS MUST BE RETAKEN FROM TEST_LOOP__NONLINEAR_REGRESSION
intrapest_test <-read_csv("/Users/dario-ssm/Documents/Dario Investigacion/IntRaPest/intrinsic_rates_pests/data/IR_data_complete_sim.csv" ) %>% 
  filter(id != 19 &
         id != 47 &
         order == "Acari>Prostigmata" |
         order =="Acari>Trombidiformes") %>% 
  mutate(order = "Acari") %>% 
  glimpse()


# - 2. nlme package ----
#  .... a) lme() ----

###random-effects
re_lme_acari <- lme(r ~ 1,
                    random = ~ 1|id,
                    data = intrapest_test)
summary(re_lme_acari)

#  .... b) nlme() ----

nlme_br1_all_varExp_sim_sens <- nlme(model= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                     start = starts_all_sim_sens,
                                     fixed = a+Tmin+Tmax ~1,
                                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                                     #random =a+Tmin+Tmax~1|id,
                                     data = IR_data_sim_sens,
                                     weights = varExp(),
                                     na.action=na.exclude,
                                     control = nlmeControl(pnlsTol = 1,
                                                           msMaxIter = 100,
                                                           msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1



# - 3. lmer package ----
#  .... a) lme() ----

nlme_br1_all_varExp_sim_sens <- nlme(model= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                     start = starts_all_sim_sens,
                                     fixed = a+Tmin+Tmax ~1,
                                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                                     #random =a+Tmin+Tmax~1|id,
                                     data = IR_data_sim_sens,
                                     weights = varExp(),
                                     na.action=na.exclude,
                                     control = nlmeControl(pnlsTol = 1,
                                                           msMaxIter = 100,
                                                           msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1




# - 4. brms paclage ----

