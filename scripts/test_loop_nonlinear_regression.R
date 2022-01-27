#### Meta-Analysis: nonlinear mixed effects meta-regression  ####

#### SCRIPT INFO ####     
#     Authors: Dar?o San Segundo Molina, Sara Vill?n P?rez, Ignacio Morales Castilla
#     Title: Bri?re-1 Modeling
#     Aim: test model fitting to each study to obtain ecologically-informative parameters (traits)
#     Date: December 2021
#     
#     Workflow: 
# if you have "IR_data_complete_sim.csv", jump directly to section 4 after running the function sections
# otherwise:
#       1. Run Section 1
#       2. Run Section 2
#       3. Run Section 3
#        
#  In both cases, continue running 4 c)
#_________________________ ####


#### 1. Load Dataset ####
rm(list=ls())
library(tidyverse)
library(svglite)
library(nlstools)
library(nls2)
library(msm)
library(magrittr)
library(cowplot)
library(car)
library(nlme)
library(brms)
#load data:
IR_data <-read_csv("/Users/dario-ssm/Documents/Dario Investigacion/IntRaPest/intrinsic_rates_pests/data/IR_metaanalysis_suitable.csv") %>%
  filter(Filter_3 == "yes") %>% #select only those which have passed the filter 3 in the source csv
  #mutate(growth_rate = replace(growth_rate, growth_rate <= 0,0))%>% #we assign 0 to negative values of intrinsic rates (biological nonsense?)
  filter(as.numeric(temperature)<50) %>% #exclude possible missleading points
  dplyr::select(Authors,Year,title, DOI,Filter_3,Approach,Subapproach,temperature,growth_rate,
                error,n_1,order,family,genus,species,feeding_guild,h_p_family,
                diet_family,RH,daylength,lat,lon)%>%
  glimpse()

#### 2. Define fitting function ####
# We use a Briere-1 model (Briere et al., 1999).This equation represents a trade-off 
# between model complexity (only 3 parameters) and non-linearity modeling effort, and 
# it has been widely used in the physiological responses literature to obtain 
# thermal performance curves, especially in development curves but also in r.

briere1 <- function(a, temp, Tmin, Tmax){
  a*temp*(temp-Tmin)*(Tmax-temp)^(1/2)
}
# Since Topt is not a direct parameter of the model, it can be derived from Tmin and Tmax
# according to Marchioro & Foerster, 2011:
Topt <- function(Tmin,Tmax,m){
  Topt=((2*m*Tmax+(m+1)*Tmin)+sqrt(4*(m^2)*(Tmax^2)+((m+1)^2)*(Tmin^2)-4*(m^2)*Tmin*Tmax))/(4*m+2)
  return(Topt)
}

#### 3. Data cleaning and manipulation ####
#### _ _ 3.1. Dataframe preparation ####
IR_data_sd <- IR_data %>%
  mutate(stdev=error*sqrt(n_1))
IR_data_title <- IR_data_sd %>%
  distinct(title)%>%
  mutate(id=row_number()) #one id per distinct paper in a new dataframe
IR_data_all <- inner_join(IR_data_sd,IR_data_title,by='title')
## now we will avoid errors by summarising mean for those papers who have repeated temperatures into one temperature : 1 row
# ids: 13, 34, 35, 41, 44, 55
problematic_ids <- c(13,34,35,41,44,55)
num_vars_except_temp <- c("id","growth_rate","error","stdev","n_1","lon","lat")
### 13
IR_data_13_numvars <- IR_data_all %>% #select numeric variables and summarise mean for each one of them
  filter(id==13) %>%
  group_by(temperature,.drop=FALSE) %>%
  summarise_all(mean) %>%
  select(temperature,all_of(num_vars_except_temp))
IR_data_13_extra <-  IR_data_all %>% #recall the categorical variables
  filter(id==13) %>%
  select(-num_vars_except_temp)%>%
  group_by(temperature)%>%
  summarise_all(unique) %>%
  select(-temperature)
IR_data_13 <- IR_data_13_numvars %>% #bind both
  bind_cols(IR_data_13_extra)
### 34
IR_data_34_numvars <- IR_data_all %>%
  filter(id==34) %>%
  group_by(temperature,.drop=FALSE) %>%
  summarise_all(mean) %>%
  select(temperature,all_of(num_vars_except_temp))
IR_data_34_extra <-  IR_data_all %>%
  filter(id==34) %>%
  select(-num_vars_except_temp,-RH)%>%
  group_by(temperature)%>%
  summarise_all(unique) %>%
  select(-temperature)
IR_data_34 <- IR_data_34_numvars %>%
  bind_cols(IR_data_34_extra)
### 35
species_35 <- IR_data_all %>% filter(id==35) %>% distinct(species) %>% select(species) #only requires to separate species
names_species_35 <- species_35$species
coordinates_35 <- IR_data_all %>%
  filter(id == 35) %>%
  group_by(lon, lat) %>% 
  summarise(lon=unique(lon),
            lat=unique(lat))
  
IR_data_35_all <- IR_data_all %>%
  filter(id==35) %>%
  mutate(species=rep(names_species_35,5))

IR_data_35 <- IR_data_35_all %>%
  filter(species == "urticae")%>%
  mutate(lon=coordinates_35$lon[6],
         lat=coordinates_35$lat[6])

IR_data_57 <- IR_data_35_all %>%
  filter(species == "ludeni") %>%
  mutate(id=57) %>% 
  mutate(lon=coordinates_35$lon[4],
         lat=coordinates_35$lat[4])

IR_data_58 <- IR_data_35_all %>%
  filter(species == "phaselus")%>%
  mutate(id=58)%>% 
  mutate(lon=coordinates_35$lon[3],
         lat=coordinates_35$lat[3])

IR_data_59 <- IR_data_35_all %>%
  filter(species == "piercei")%>%
  mutate(id=59)%>% 
  mutate(lon=coordinates_35$lon[2],
         lat=coordinates_35$lat[2])

IR_data_60 <- IR_data_35_all %>%
  filter(species == "truncatus")%>%
  mutate(id=60)%>% 
  mutate(lon=coordinates_35$lon[1],
         lat=coordinates_35$lat[1])
### 41
IR_data_41_numvars <- IR_data_all %>%
  filter(id==41) %>%
  group_by(temperature,.drop=FALSE) %>%
  summarise_all(mean) %>%
  select(temperature,all_of(num_vars_except_temp),RH)
IR_data_41_extra <-  IR_data_all %>%
  filter(id==41) %>%
  select(-num_vars_except_temp,-RH)%>%
  group_by(temperature)%>%
  summarise_all(unique) %>%
  select(-temperature)
IR_data_41 <- IR_data_41_numvars %>%
  bind_cols(IR_data_41_extra)
### 41
IR_data_44_numvars <- IR_data_all %>%
  filter(id==44) %>%
  group_by(temperature) %>%
  summarise_all(mean) %>%
  mutate(lon = 138, lat=36)%>%
  select(temperature,all_of(num_vars_except_temp))
IR_data_44_extra <-  IR_data_all %>%
  filter(id==44) %>%
  select(-num_vars_except_temp)%>%
  group_by(temperature)%>%
  summarise_all(unique) %>%
  select(-temperature)
IR_data_44 <- IR_data_44_numvars %>%
  bind_cols(IR_data_44_extra)
### 55
IR_data_55_numvars <- IR_data_all %>%
  filter(id==55 &
           species == "fragariae") %>%
  group_by(temperature) %>%
  summarise_all(mean) %>%
  select(temperature,all_of(num_vars_except_temp))
IR_data_55_extra <-  IR_data_all %>%
  filter(id==55 & species == "fragariae") %>%
  select(-num_vars_except_temp)%>%
  group_by(temperature)%>%
  summarise_all(unique) %>%
  select(-temperature)
IR_data_55 <- IR_data_55_numvars %>%
  bind_cols(IR_data_55_extra)
IR_data_61 <- IR_data_all %>%
  filter(id==55 &
           species == "miscanthi")%>%
  mutate(id=61)
# now we ensemble all these subsets into the main dataset
IR_data_all_rev <- IR_data_all %>%
  filter(id != 13 &
           id!= 34 &
           id!= 35 &
           id!= 41 &
           id!= 44 &
           id!= 55) %>%
  bind_rows(IR_data_13,
            IR_data_34,
            IR_data_35,
            IR_data_41,
            IR_data_44,
            IR_data_55,
            IR_data_57,
            IR_data_58,
            IR_data_59,
            IR_data_60,
            IR_data_61)%>%
  arrange(id)
IR_data_all <- IR_data_all_rev #recover original name of the tibble

#### _ _ 3.2. Sensitivity analyses ####
# Here we run the models after excluding
# IDs 19 and 47 (see Appendix I). Then we compare the output with the whole dataset (including those)
# from the different approaches explored along this script.
IR_data_all_without_papers <- IR_data_all %>% 
  filter(id != 19 &
         id != 47) %>% 
  glimpse()
## plot it
plot_without_paperes <- ggplot(IR_data_all_without_papers, aes(x = temperature, y = growth_rate))+
  geom_point(aes(color = as.factor(id)),
             alpha = .5)+
  theme_bw()+
  labs(x = "Temperature",
       y = "Intrinsic Rate of Increase (r)",
       title = "Effect sizes",
       subtitle = "Excluding papers 19 and 47; color = Study")+
  geom_smooth(method = "loess",
              color = "lightcoral",
              fill = "lightcoral")+
  theme(legend.position = "none")
plot_without_paperes
#### _ _ _ _ a) nlme raw ####
#### _ _ _ _ _ _ i) without papers + varExp ####
grid_br1_sens <- expand.grid(list(a=seq(3e-05,1.5e-04,by=1e-05),
                                   Tmin=seq(8,10,by=0.5),
                                   Tmax=seq(33,48,by=0.5)))
fitted_br1_sens_brute<- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                              data = IR_data_all_without_papers,
                              start = grid_br1_sens,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_brute_sens <- summary(fitted_br1_sens_brute)
starts_all_sens <- sum_brute_sens$parameters[,1]
nlme_br1_all_varExp_sens <- nlme(model= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     start = starts_all_sens,
                     fixed = a+Tmin+Tmax ~1,
                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                     #random =a+Tmin+Tmax~1|id,
                     data = IR_data_all_without_papers,
                     weights = varExp(),
                     na.action=na.exclude,
                     control = nlmeControl(pnlsTol = 10,
                                           msMaxIter = 100,
                                           msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
sum_nlme_all_sens_varExp_ <- summary(nlme_br1_all_varExp_sens)
AIC_varExp_without <- AIC(nlme_br1_all_varExp_sens)
loglik_varExp_without <- logLik(nlme_br1_all_varExp_sens)[1]
params_varExp_without <- sum_nlme_all_sens_varExp_$coefficients$fixed
varExp_coef_without <- sum_nlme_all_sens_varExp_$modelStruct$varStruct[1]
overview_nlme_all_varExp_without <- tibble(a = params_varExp_without[1],
                                    Tmin = params_varExp_without[2],
                                    Tmax = params_varExp_without[3],
                                    varExp_coef = varExp_coef_without,
                                    AIC = AIC_varExp_without,
                                    log_likelihood = loglik_varExp_without)
#### _ _ _ _ _ _ ii) without papers + varPower ####
nlme_br1_all_sens_varPow <- nlme(model= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                          start = starts_all_sens,
                          fixed = a+Tmin+Tmax ~1,
                          groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                          #random =a+Tmin+Tmax~1|id,
                          data = IR_data_all_without_papers,
                          weights = varPower(),
                          na.action=na.exclude,
                          control = nlmeControl(pnlsTol = 10,
                                                msMaxIter = 100,
                                                msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
sum_nlme_all_sens_varPow <- summary(nlme_br1_all_sens_varPow)
AIC_varPower_without <- AIC(nlme_br1_all_sens_varPow)
loglik_varPower_without <- logLik(nlme_br1_all_sens_varPow)[1]
params_varPower_without <- sum_nlme_all_sens_varPow$coefficients$fixed
varPow_coef_without <- sum_nlme_all_sens_varPow$modelStruct$varStruct[1]
overview_nlme_all_varPow_without <- tibble(a = params_varPower_without[1],
                                    Tmin = params_varPower_without[2],
                                    Tmax = params_varPower_without[3],
                                    varExp_coef = varPow_coef_without,
                                    AIC = AIC_varPower_without,
                                    log_likelihood = loglik_varPower_without)
#### _ _ _ _ _ _ iii) with papers + varExp ####
nlme_br1_all_varExp <- nlme(model= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                            start = starts_all_sens,
                            fixed = a+Tmin+Tmax ~1,
                            groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                            #random =a+Tmin+Tmax~1|id,
                            data = IR_data_all,
                            weights = varExp(),
                            na.action=na.exclude,
                            control = nlmeControl(pnlsTol = 10,
                                                  msMaxIter = 100,
                                                  msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
sum_nlme_all_varExp <- summary(nlme_br1_all_varExp)
AIC_varExp <- AIC(nlme_br1_all_varExp)
loglik_varExp <- logLik(nlme_br1_all_varExp)[1]
params_varExp <- sum_nlme_all_varExp$coefficients$fixed
varExp_coef <- sum_nlme_all_varExp$modelStruct$varStruct[1]
overview_nlme_all_varExp <- tibble(a = params_varExp[1],
                                   Tmin = params_varExp[2],
                                   Tmax = params_varExp[3],
                                   varExp_coef = varExp_coef,
                                   AIC = AIC_varExp,
                                   log_likelihood = loglik_varExp)

#### _ _ _ _ _ _ iv) with papers + varPower ####
nlme_br1_all_varPow <- nlme(model= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                                 start = starts_all_sens,
                                 fixed = a+Tmin+Tmax ~1,
                                 groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                                 #random =a+Tmin+Tmax~1|id,
                                 data = IR_data_all,
                                 weights = varPower(),
                                 na.action=na.exclude,
                                 control = nlmeControl(pnlsTol = 10,
                                                       msMaxIter = 100,
                                                       msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
sum_nlme_all_varPow <- summary(nlme_br1_all_varPow)
AIC_varPower <- AIC(nlme_br1_all_sens_varPow)
loglik_varPower <- logLik(nlme_br1_all_varPow)[1]
params_varPower <- sum_nlme_all_varPow$coefficients$fixed
varPow_coef <- sum_nlme_all_varPow$modelStruct$varStruct[1]
overview_nlme_all_varPow <- tibble(a = params_varPower[1],
                                           Tmin = params_varPower[2],
                                           Tmax = params_varPower[3],
                                           varExp_coef = varPow_coef,
                                           AIC = AIC_varPower,
                                           log_likelihood = loglik_varPower)

#### _ _ _ _ _ _ v) compare them ####
first_sensitivity_analysis <- overview_nlme_all_varExp_without %>% 
  bind_rows(overview_nlme_all_varExp,
            overview_nlme_all_varPow_without,
            overview_nlme_all_varPow) %>% 
  mutate(analysis = c("varExp, excluding papers",
                      "varExp, including papers",
                      "varPower, excluding papers",
                      "varPower, including papers"))
View(first_sensitivity_analysis)
#### _ _ _ _ b) nlme simul ####
#### _ _ _ _ _ _ i) without papers + varExp ####
IR_data_sim_sens <- read_csv("/Users/dario-ssm/Documents/Dario Investigacion/IntRaPest/intrinsic_rates_pests/data/IR_data_complete_sim.csv") %>% 
  filter(id != 19 &
         id != 47)
grid_br1_sim_sens <- expand.grid(list(a=seq(9e-05,2e-04,by=1e-05),
                                  Tmin=seq(-5,10,by=0.5),
                                  Tmax=seq(33,48,by=0.5)))
fitted_br1_sim_sens_brute<- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                             data = IR_data_sim_sens,
                             start = grid_br1_sim_sens,
                             algorithm = "brute-force",
                             trace = FALSE)
sum_brute_sim_sens <- summary(fitted_br1_sim_sens_brute)
starts_all_sim_sens <- sum_brute_sim_sens$parameters[,1]
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
sum_nlme_all_sim_sens_varExp <- summary(nlme_br1_all_varExp_sim_sens)
AIC_varExp_sim_without <- AIC(nlme_br1_all_varExp_sim_sens)
loglik_varExp_sim_without <- logLik(nlme_br1_all_varExp_sim_sens)[1]
params_varExp_sim_without <- sum_nlme_all_sim_sens_varExp$tTable[,1]
varExp_coef_sim_without <- sum_nlme_all_sim_sens_varExp$modelStruct$varStruct[1]
overview_nlme_all_varExp_sim_without <- tibble(a = params_varExp_sim_without[1],
                                               se_a = sum_nlme_all_sim_sens_varExp$tTable[1,2],
                                               Tmin = params_varExp_sim_without[2],
                                               se_Tmin = sum_nlme_all_sim_sens_varExp$tTable[2,2],
                                               Tmax = params_varExp_sim_without[3],
                                               se_Tmax = sum_nlme_all_sim_sens_varExp$tTable[3,2],
                                               Topt =  Topt(Tmin = params_varExp_sim_without[2],
                                                            Tmax = params_varExp_sim_without[3],
                                                            m=2),
                                               Topt_se<- deltaMethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                                                                         coef(nlme_br1_all_varExp_sim_sens), vcov(nlme_br1_all_varExp_sim_sens)),
                                               deltaMethod(nlme_br1_all_varExp_sim_sens, g. = "((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2)",
                                                           vcov. = vcov(nlme_br1_all_varExp_sim_sens, complete=FALSE),
                                                           parameterNames = names(fixef(nlme_br1_all_varExp_sim_sens)), envir=parent.frame())                                               
                                               )
                                               ### complete it with Topt_se deltaMethod
                                               varExp_coef = varExp_coef_sim_without,
                                               AIC = AIC_varExp_sim_without,
                                               log_likelihood = loglik_varExp_sim_without)
#### _ _ _ _ _ _ ii) without papers + varPower ####

nlme_br1_all_varPow_sim_sens <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                                     start = starts_all_sim_sens,
                                     fixed = a+Tmin+Tmax ~1,
                                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                                     #random =a+Tmin+Tmax~1|id,
                                     data = IR_data_sim_sens,
                                     weights = varPower(),
                                     control = nlmeControl(pnlsTol = 1,
                                                           msMaxIter = 100,
                                                           msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
## ERROR convergence
#### _ _ _ _ _ _ iii) with papers + varExp ####
IR_data_sim_all <- read_csv("/Users/Ecologia/Desktop/DAR?O_actualizada septiembre 2021/Intrinsic_metaanalysis/IR_data_complete_sim.csv")

grid_br1_sim_all <- expand.grid(list(a=seq(9e-05,2e-04,by=1e-05),
                                      Tmin=seq(-5,10,by=0.5),
                                      Tmax=seq(33,48,by=0.5)))
fitted_br1_sim_all_brute<- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                                 data = IR_data_sim_all,
                                 start = grid_br1_sim_all,
                                 algorithm = "brute-force",
                                 trace = FALSE)
sum_brute_sim_all <- summary(fitted_br1_sim_all_brute)
starts_all_sim_all<- sum_brute_sim_all$parameters[,1]
nlme_br1_all_varExp_sim_all <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                            start = c(a = 0.00005, # if we use starts_all_sim_all, it does not converge
                                      Tmin = -5,
                                      Tmax = 48),
                            fixed = a+Tmin+Tmax ~1,
                            groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                            #random =a+Tmin+Tmax~1|id,
                            data = IR_data_sim_all,
                            weights = varExp(),
                            control = nlmeControl(pnlsTol = 1e-01,
                                                  msMaxIter = 100,
                                                  msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
sum_nlme_all_varExp_sim_all  <- summary(nlme_br1_all_varExp_sim_all)
AIC_varExp_sim_all <- AIC(nlme_br1_all_varExp_sim_all)
loglik_varExp_sim_all <- logLik(nlme_br1_all_varExp_sim_all)[1]
params_varExp_sim_all <- sum_nlme_all_varExp_sim_all$coefficients$fixed
varExp_coef_sim_all <- sum_nlme_all_varExp_sim_all$modelStruct$varStruct[1]
overview_nlme_all_varExp_sim_all <- tibble(a = params_varExp_sim_all[1],
                                               se_a = sum_nlme_all_varExp_sim_all$tTable[1,2],
                                               Tmin = params_varExp_sim_all[2],
                                               se_Tmin = sum_nlme_all_varExp_sim_all$tTable[2,2],
                                               Tmax = params_varExp_sim_all[3],
                                               se_Tmax = sum_nlme_all_varExp_sim_all$tTable[3,2],
                                               varExp_coef = varExp_coef_sim_all,
                                               AIC = AIC_varExp_sim_all,
                                               log_likelihood = loglik_varExp_sim_all)

#### _ _ _ _ _ _ iv) with papers + varPower ####

nlme_br1_all_varPow_sim_all <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                                    start = starts_all_sim_all,
                                    fixed = a+Tmin+Tmax ~1,
                                    groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                                    #random =a+Tmin+Tmax~1|id,
                                    data = IR_data_sim_all,
                                    weights = varPower(),
                                    control = nlmeControl(pnlsTol = 100,
                                                          msMaxIter = 100,
                                                          msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
## ERROR of CONVERGENCE

#### _ _ _ _ _ _ v) compare them ####
second_sensitivity_analysis <- overview_nlme_all_varExp_sim_without %>% 
  bind_rows(overview_nlme_all_varExp_sim_all) %>% 
  mutate(analysis = c("varExp, excluding papers",
                      "varExp, including papers")) %>% 
  write_csv("comparison_including_vs_exluding_ids47&19.csv")
View(second_sensitivity_analysis)

#### _ _ _ _ c) Conclusion: ####
first_sensitivity_analysis
second_sensitivity_analysis
# Although with rough data inluding studies 19 and 47 slightly improves the model
# (see first_sensitivity_analysis), the simulated data fitting
# dramatically improves when we exclude these papers (see second_sensitivity_analysis)











### 4. Loop for simulation csv generator & export ####
setwd("C:/Users/Ecologia/Desktop/DAR?O_actualizada septiembre 2021/Intrinsic_metaanalysis/simulations") #to save plots here together
distinct_ids <- IR_data_all %>%  
  distinct(id)
#### _ _ _ 4.1. Data simulations ####
## simulated_dataframe
ID <- 1:length(distinct_ids$id)
for (i in ID){ #para cada valor entre los que est?n en nuestro vector ("Study" or "id"; es decir, del 1 al 56)
  
  print(ID[i])  #primero ver por qu? iteraci?n vamos
  IR_data_ID <- IR_data_all %>% #subset con el conjunto que hacemos en cada iteraci?n
    filter(id==i)
  temp_ID <- IR_data_ID %>% select(temperature)  # first we assume normal distribution
  simul_ID <- tibble(id =rep(i,as.numeric(IR_data_ID %>% filter(id==i) %>% select(n_1) %>% summarise(n_1=sum(n_1)))),
                     "temp" = 0,
                     "r" = 0 ,
                     "sd" = 0)
  enes <-IR_data_ID %>% filter(id==i) %>% select(n_1)
  position <- cumsum(enes)-enes[1]+1
  iter_temp <- 1:length(temp_ID$temperature)
  for (num in iter_temp){
    simul_ID$temp[position$n_1[num]:(position$n_1[num]-1+enes$n_1[num])] <- rep(temp_ID$temperature[num], each = enes$n_1[num])
  }
  simul_ID
  
  for(t in 1:length(temp_ID$temperature)){
    
    temper <- as.numeric(temp_ID[t,])
    n <-  as.numeric(IR_data_ID %>% filter(id==i & temperature==temper)%>% select(n_1))
    mu <- as.numeric(IR_data_ID %>% filter(id==i & temperature==temper)%>% select(growth_rate))
    sd <-as.numeric( IR_data_ID %>% filter(id==i & temperature==temper)%>% select(stdev))
    sim_r <- rnorm(n,mu,sd)
    rep_stdev <- rep(sd, n)
    simul_ID[simul_ID$temp == temper,"r"] <- tibble(sim_r) 
    simul_ID[simul_ID$temp == temper,"sd"] <- tibble(rep_stdev)
  }
  simul_ID
  png(filename = paste0("/Users/Ecologia/Desktop/DAR?O_actualizada septiembre 2021/Intrinsic_metaanalysis/simulations/simulated_",i,".png"))
  plot(IR_data_ID$temperature,IR_data_ID$growth_rate) #ploteamos y guardamos un plot para tenerlos todos luego
  dev.off()
  write_csv(simul_ID,file = paste0("simulated_data_study_id_",i,".csv"))
} 

#now read each csv and assign other variables
## What's the number of rows for each simulation?
repetitions <- rep(NA,length(ID))
for(i in ID) {
  enes <- IR_data_all %>% filter(id==i) %>% select(n_1)
  position <- cumsum(enes)
  num_reps <- tail(position$n_1, n=1)
  repetitions[i] <- num_reps
}
repetitions
dimensions <- tail(cumsum(repetitions),n=1)
dimensions_position <- cumsum(repetitions)
## now we generate the empty dataframe
simulations_empty <- tibble(id = rep(0, dimensions),
                            temp = rep(0, dimensions),
                            r = rep(0, dimensions),
                            sd = rep(0, dimensions)) #empty to fill with the loop

for(i in ID) {
  simulated_ID <- read_csv(file = paste0("simulated_data_study_id_",i,".csv"))#read each year csv previously generated
  simulations_empty[((dimensions_position[i]-repetitions[i]+1):dimensions_position[i]),] <- simulated_ID
}

simulations_numeric <- simulations_empty
#### _ _ _ 4.2. Prepare covariates df to bind ####
#and now repeat that along length of the simulations
covariate_IR_data <- IR_data_all %>% 
  dplyr::select(id,
                order, 
                family, 
                genus, 
                species, 
                feeding_guild, 
                lat, 
                lon, 
                title, 
                Authors)
simulations_others_empty<- tibble(id=rep(0,dimensions),
                                  order=rep("letras",dimensions),
                                  family=rep("letras",dimensions),
                                  genus=rep("letras",dimensions),
                                  species=rep("letras",dimensions),
                                  feeding_guild=rep("letras",dimensions),
                                  lat=rep(0,dimensions),
                                  lon=rep(0,dimensions),
                                  title=rep("letras",dimensions),
                                  Authors=rep("letras",dimensions)
) #empty to fill with the loop
##  add
for(i in ID) {
  simulations_ID <- covariate_IR_data %>%
    filter(id==i) %>% 
    mutate_all(unique)
  simulations_ID_rep <- simulations_ID %>%
    slice(rep(row_number(1), repetitions[i]))  #slice is useful to take the values of a row
  simulations_others_empty[(dimensions_position[i]-repetitions[i]+1):dimensions_position[i],] <- simulations_ID_rep
print(paste0(i/61*100, " %"))
}
simulations_covars <- simulations_others_empty
#### _ _ _ 4.3. Ensemble sim dataset ####
IR_data_complete_sim <- simulations_numeric %>% 
  select(-id)%>%
  bind_cols(simulations_covars)%>%
  write_csv(file = "IR_data_complete_sim.csv")
#let's see sd values
boxplot(IR_data_complete_sim$sd) # some values that we have to take appart...
exclusion_limit <- quantile(IR_data_complete_sim$sd,probs = .95)
#and clean
#simulated_IR_data_complete <- IR_data_complete_sim %>% 
#  filter(abs(r) <= 0.8 &
#           temp!=0 &
#           r >= 0 &
#           sd <= exclusion_limit)
#




plot_simul <- ggplot(IR_data_complete_sim, aes(x = temp, y = r))+
  geom_point(aes(color = as.factor(id)),
             alpha = .1,
             #position = position_jitter(width = 4)
             )+
  theme_bw()+
  labs(x = "Temperature",
       y = "Intrinsic Rate of Increase (r)",
       title = "Simulated effect sizes",
       subtitle = "r ~ N(mu,sigma); color = Study",
       color = "Study")+
  geom_smooth(method = "loess",
              color = "lightpink4")+
  theme(legend.position = "none")
plot_simul
ggsave("plot_simulated_IR_data_all.png",
       dpi = 300,
       width = 20,
       height = 25,
       units = "cm")


#### 5. nlme model covariates ####
#### _ _ a) Order, all ####
acari_data <- IR_data_sim_sens %>%
  filter(order == "Acari>Prostigmata" | 
           order == "Acari>Trombidiformes")%>%
  mutate(order = "Acari")
IR_data_all_non_acari <- IR_data_sim_sens %>%
  filter(order != "Acari>Prostigmata" & 
           order != "Acari>Trombidiformes") 
IR_data_sim_sens_rev <- IR_data_all_non_acari %>%
  bind_rows(acari_data) 
  
starts_all_sim_sens
startvals_order_subset <- c(starts_all_sim_sens[1],
                       rep(0, 6),
                       starts_all_sim_sens[2],
                       rep(0, 6),
                       starts_all_sim_sens[3],
                       rep(0, 6))

nlme_br1_order_subset <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = list(fixed = startvals_order_subset),
                     fixed = a+Tmin+Tmax ~ as.factor(order),
                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                     data = IR_data_sim_sens_rev,
                     weights = varExp(),
                     control = nlmeControl(pnlsTol = 1,
                                           msMaxIter = 50,
                                           maxIter = 50,
                                           msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
nlme_order_all <- nlme_br1_order_subset

#### _ _ b) Order, subset ####
## let's try with well represented order subset (lepidoptera, hemiptera, acari)
IR_sim_ord_subset <- IR_data_sim_sens_rev %>% 
  filter(order == "Acari" |
         order == "Lepidoptera" |
         order == "Hemiptera") %>% 
  filter(id != 37)
IR_sim_acari <- IR_sim_ord_subset %>% 
  filter(order == "Acari")
IR_sim_hemiptera <- IR_sim_ord_subset %>% 
  filter(order == "Hemiptera")
IR_sim_lepidoptera <- IR_sim_ord_subset %>% 
  filter(order == "Lepidoptera")
#### _ _ _ _ b.1. starting values ####
## Acari:
grid_br1_sim_order_acari <- expand.grid(list(a=seq(5e-05,2e-04,by=1e-05),
                                     Tmin=seq(0,15,by=0.5),
                                     Tmax=seq(33,48,by=0.5)))
fitted_br1_sim_acari <- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                                data = IR_sim_acari,
                                start = grid_br1_sim_order_acari,
                                algorithm = "brute-force",
                                trace = FALSE)
sum_brute_sim_acari <- summary(fitted_br1_sim_acari)
starts_sim_acari <- sum_brute_sim_acari$parameters[,1]

## Hemiptera:
grid_br1_sim_order_hemiptera <- expand.grid(list(a=seq(1e-05,2e-04,by=1e-05),
                                             Tmin=seq(5,15,by=0.5),
                                             Tmax=seq(33,48,by=0.5)))
fitted_br1_sim_hemiptera <- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                             data = IR_sim_hemiptera,
                             start = grid_br1_sim_order_hemiptera,
                             algorithm = "brute-force",
                             trace = FALSE)
sum_brute_sim_hemiptera <- summary(fitted_br1_sim_hemiptera)
starts_sim_hemiptera <- sum_brute_sim_hemiptera$parameters[,1]

## Lepidoptera:
grid_br1_sim_order_lepidoptera <- expand.grid(list(a=seq(1e-05,2e-04,by=1e-05),
                                                 Tmin=seq(10,15,by=0.5),
                                                 Tmax=seq(33,48,by=0.5)))
fitted_br1_sim_lepidoptera <- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                                 data = IR_sim_lepidoptera,
                                 start = grid_br1_sim_order_lepidoptera,
                                 algorithm = "brute-force",
                                 trace = FALSE)
sum_brute_sim_lepidoptera <- summary(fitted_br1_sim_lepidoptera)
starts_sim_lepidoptera <- sum_brute_sim_lepidoptera$parameters[,1]

## ensemble all startings:
startvals_order_subset <- c(starts_sim_acari[1],
                            starts_sim_hemiptera[1],
                            starts_sim_lepidoptera[1],
                            starts_sim_acari[2],
                            starts_sim_hemiptera[2],
                            starts_sim_lepidoptera[2],
                            starts_sim_acari[3],
                            starts_sim_hemiptera[3],
                            starts_sim_lepidoptera[3])
#### _ _ _ _ b.2. fit model ####

nlme_br1_order_subset <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                              start = c(fixed = startvals_order_subset),
                              fixed = a+Tmin+Tmax ~ as.factor(order),
                              groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                              data = IR_sim_ord_subset,
                              weights = varExp(),
                              control = nlmeControl(pnlsTol = 10,
                                                    msMaxIter = 50,
                                                    maxIter = 50,
                                                    msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
sum_nlme_order_subset <- summary(nlme_br1_order_subset)
##:( error of convergence)
#let's try different start values combinations
startvals_order_subset
startvals_order_invent <- c(0.00010, 0, 0,
                            6.5,0, 0,
                            43, 0, 0)

nlme_br1_order_subset_invent <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                              start = c(fixed = startvals_order_invent),
                              fixed = a+Tmin+Tmax ~ as.factor(order),
                              groups = ~id,
                              data = IR_sim_ord_subset,
                              weights = varExp(),
                              control = nlmeControl(pnlsTol = 10,
                                                    msMaxIter = 50,
                                                    maxIter = 50,
                                                    msVerbose = TRUE))


#### _ _ b) Feeding guild, all ####
miners_data <- IR_data_sim_sens %>%
  filter(feeding_guild == "miner" | 
         feeding_guild == "leafminer")%>%
  mutate(feeding_guild = "miner")
IR_data_all_non_miners <- IR_data_sim_sens %>%
  filter(feeding_guild != "miner" & 
         feeding_guild != "leafminer")
IR_data_sim_sens_rev <- IR_data_all_non_miners %>%
  bind_rows(miners_data) 
IR_data_sim_sens_rev %>% 
  distinct(feeding_guild)

## borers:
borer_rev_data <- IR_data_sim_sens_rev %>% 
  filter(feeding_guild == "borer")
grid_br1_sim_fg_borer <- expand.grid(list(a=seq(1e-05,2e-04,by=1e-05),
                                           Tmin=seq(13,18,by=0.5),
                                           Tmax=seq(33,48,by=0.5)))
fitted_br1_sim_borer <- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                              data = borer_rev_data,
                              start = grid_br1_sim_fg_borer,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_brute_sim_borer  <- summary(fitted_br1_sim_borer)
starts_sim_borer <- sum_brute_sim_borer$parameters[,1]
## miner:
miner_rev_data <- IR_data_sim_sens_rev %>% 
  filter(feeding_guild == "miner")

grid_br1_sim_fg_miner <- expand.grid(list(a=seq(1e-05,2e-04,by=1e-05),
                                           Tmin=seq(8,15,by=0.5),
                                           Tmax=seq(33,48,by=0.5)))
fitted_br1_sim_miner <- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                              data = miner_rev_data,
                              start = grid_br1_sim_fg_miner,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_brute_sim_miner  <- summary(fitted_br1_sim_miner)
starts_sim_miner <- sum_brute_sim_miner$parameters[,1]
## sucker:
sucker_rev_data <- IR_data_sim_sens_rev %>% 
  filter(feeding_guild == "sucker")

grid_br1_sim_fg_sucker <- expand.grid(list(a=seq(1e-05,2e-04,by=1e-05),
                                          Tmin=seq(8,15,by=0.5),
                                          Tmax=seq(33,48,by=0.5)))
fitted_br1_sim_sucker <- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                             data = sucker_rev_data,
                             start = grid_br1_sim_fg_sucker,
                             algorithm = "brute-force",
                             trace = FALSE)
sum_brute_sim_sucker  <- summary(fitted_br1_sim_sucker)
starts_sim_sucker <- sum_brute_sim_sucker$parameters[,1]
## chewer:
chewer_rev_data <- IR_data_sim_sens_rev %>% 
  filter(feeding_guild == "chewer")

grid_br1_sim_fg_chewer <- expand.grid(list(a=seq(1e-05,2e-04,by=1e-05),
                                           Tmin=seq(8,15,by=0.5),
                                           Tmax=seq(33,48,by=0.5)))
fitted_br1_sim_chewer <- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                              data = chewer_rev_data,
                              start = grid_br1_sim_fg_chewer,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_brute_sim_chewer  <- summary(fitted_br1_sim_chewer)
starts_sim_chewer <- sum_brute_sim_chewer$parameters[,1]

startvals_all_fg <- c(starts_sim_borer[1],
                      starts_sim_chewer[1],
                      starts_sim_miner[1],
                      starts_sim_sucker[1],
                      starts_sim_borer[2],
                      starts_sim_chewer[2],
                      starts_sim_miner[2],
                      starts_sim_sucker[2],
                      starts_sim_borer[3],
                      starts_sim_chewer[3],
                      starts_sim_miner[3],
                      starts_sim_sucker[3]) 
startvals_all_fg

nlme_br1_sim_fg <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = c(fixed = startvals_all_fg),
                     fixed = a+Tmin+Tmax ~ as.factor(feeding_guild),
                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                     data = IR_data_sim_sens_rev,
                     weights = varExp(),
                     control = nlmeControl(pnlsTol = 1,
                                           msMaxIter = 50,
                                           maxIter = 100,
                                           msVerbose = TRUE,
                                           ))#to avoid error of singularity in backsolve at level 0; block 1
# Error in nlme.formula(model = r ~ briere1(a, temp, Tmin, Tmax), start = c(fixed = startvals_idea_fg),  : 
#                         maximum number of iterations (maxIter = 100) reached without convergence
#                       In addition: There were 50 or more warnings (use warnings() to see the first 50)

#### 6. Bayesian brms ####
#### _ _  a) example for a small subset ####
acari_data <- IR_data %>% 
  filter(order == "Acari>Prostigmata" |
         order == "Acari>Trombidiformes") %>% 
  mutate(order = "Acari")

nls_acari <- nls(r ~ briere1(a,temperature,Tmin,Tmax),
                 data = acari_data,
                 start = c(a = 2e-04, Tmin = 6.5, Tmax = 38))
summary(nls_acari)
prior_1 <- c(set_prior("normal(1.317e-04, (2.855e-05)^2)", nlpar = "a"),
             set_prior("normal(8.158e+00, (2.546e+00)^2)", nlpar = "Tmin"),
             set_prior("normal(4.097e+01, (1.979e+00)^2)", nlpar = "Tmax"))

briere_formula <- brmsformula(formula =  r ~ a*temp*(temp-Tmin)*(Tmax-temp)^(1/2)+1|id,
                     # Nonlinear variables
                     a+Tmin+Tmax ~ 1,
                     # Nonlinear fit
                     nl = TRUE)
bayes_fit <- brm(
  briere_formula,
  family=gaussian(), 
  data = acari_data,
  prior = prior_1)
summary(bayes_fit)



#### 7. Loop for individual study model fitting ####
# let's compile standard deviations
IR_data_sd <- IR_data %>%
  mutate(stdev=error*sqrt(n_1)) %>% 
  filter(word(Authors,1) != "Vangansbeke," &
         word(Authors,1) != "Xie,") #problems: that study has only two temperature treatments
# and let's ad an id to each unique paper (although some papers with different treatments within
 # e.g. for different species in the same study, will be treated as different studies)
IR_data_title <- IR_data_sd %>%
  distinct(title)%>%
  mutate(id=row_number()) #one id per distinct paper in a new dataframe
IR_data_all <- inner_join(IR_data_sd,IR_data_title,by='title')

## now we will avoid errors by sumarising mean for those papers who have repeated temperatures into one temperature : 1 row
# first let's see problematic studies:
# let's check out if duplicate temperature treatments have been summarised
subset_problems <- IR_data_all %>% 
  group_by(id, temperature) %>%
  count() %>% 
  filter(n >1) %>% 
  print()
troublemakers <- subset_problems %>% 
  ungroup() %>% 
  select(id) %>% 
  distinct(id) %>% 
  glimpse()
# ids: 10, 27, 32, 42, 50, 52
num_vars_except_temp <- c("id","growth_rate","error","stdev","n_1","lon","lat")
### 10
IR_data_10_numvars <- IR_data_all %>% #select numeric variables and summarise mean for each one of them
  filter(id==10) %>%
  group_by(temperature,.drop=FALSE) %>%
  summarise_all(mean) %>%
  select(temperature,all_of(num_vars_except_temp))
IR_data_10_extra <-  IR_data_all %>% #recall the categorical variables
  filter(id==10) %>%
  select(-num_vars_except_temp)%>%
  group_by(temperature)%>%
  summarise_all(unique) %>%
  select(-temperature)
IR_data_10 <- IR_data_10_numvars %>% #bind both
  bind_cols(IR_data_10_extra)

### 27
IR_data_27_numvars <- IR_data_all %>%
  filter(id==27) %>%
  group_by(temperature,.drop=FALSE) %>%
  summarise_all(mean) %>%
  select(temperature,all_of(num_vars_except_temp),RH)
IR_data_27_extra <-  IR_data_all %>%
  filter(id==27) %>%
  select(-num_vars_except_temp,-RH)%>%
  group_by(temperature)%>%
  summarise_all(unique) %>%
  select(-temperature)
IR_data_27 <- IR_data_27_numvars %>%
  bind_cols(IR_data_27_extra)

### 32
IR_data_32_numvars <- IR_data_all %>%
  filter(id==32) %>%
  group_by(temperature,.drop=FALSE) %>%
  summarise_all(mean) %>%
  select(temperature,all_of(num_vars_except_temp))
IR_data_32_extra <-  IR_data_all %>%
  filter(id==32) %>%
  select(-num_vars_except_temp,-RH)%>%
  group_by(temperature)%>%
  summarise_all(unique) %>%
  select(-temperature)
IR_data_32 <- IR_data_32_numvars %>%
  bind_cols(IR_data_32_extra)

### 50
species_50 <- IR_data_all %>% filter(id==50) %>% distinct(species) %>% select(species) #only requires to separate species
names_species_50 <- species_50$species
coordinates_50 <- IR_data_all %>%
  filter(id == 50) %>%
  group_by(lon, lat) %>% 
  summarise(lon=unique(lon),
            lat=unique(lat))

IR_data_50_all <- IR_data_all %>%
  filter(id==50) %>%
  mutate(species=rep(names_species_50,5))

IR_data_50 <- IR_data_50_all %>%
  filter(species == "urticae")%>%
  mutate(lon=coordinates_50$lon[6],
         lat=coordinates_50$lat[6])

IR_data_55 <- IR_data_50_all %>%
  filter(species == "ludeni") %>%
  mutate(id=55) %>% 
  mutate(lon=coordinates_50$lon[4],
         lat=coordinates_50$lat[4])

IR_data_56 <- IR_data_50_all %>%
  filter(species == "phaselus")%>%
  mutate(id=56)%>% 
  mutate(lon=coordinates_50$lon[3],
         lat=coordinates_50$lat[3])

IR_data_57 <- IR_data_50_all %>%
  filter(species == "piercei")%>%
  mutate(id=57)%>% 
  mutate(lon=coordinates_50$lon[2],
         lat=coordinates_50$lat[2])

IR_data_58 <- IR_data_50_all %>%
  filter(species == "truncatus")%>%
  mutate(id=58)%>% 
  mutate(lon=coordinates_50$lon[1],
         lat=coordinates_50$lat[1])

### 52
IR_data_52_numvars <- IR_data_all %>%
  filter(id==52) %>%
  group_by(temperature) %>%
  summarise_all(mean) %>%
  mutate(lon = 138, lat=36)%>%
  select(temperature,all_of(num_vars_except_temp))
IR_data_52_extra <-  IR_data_all %>%
  filter(id==52) %>%
  select(-num_vars_except_temp)%>%
  group_by(temperature)%>%
  summarise_all(unique) %>%
  select(-temperature)
IR_data_52 <- IR_data_52_numvars %>%
  bind_cols(IR_data_52_extra)

### 42
IR_data_42_numvars <- IR_data_all %>%
  filter(id==42 &
           species == "fragariae") %>%
  group_by(temperature) %>%
  summarise_all(mean) %>%
  select(temperature,all_of(num_vars_except_temp))
IR_data_42_extra <-  IR_data_all %>%
  filter(id==42 & species == "fragariae") %>%
  select(-num_vars_except_temp)%>%
  group_by(temperature)%>%
  summarise_all(unique) %>%
  select(-temperature)
IR_data_42 <- IR_data_42_numvars %>%
  bind_cols(IR_data_42_extra)
IR_data_59 <- IR_data_all %>%
  filter(id==42 &
           species == "miscanthi")%>%
  mutate(id=59)

# now we ensemble all these subsets into the main dataset
print(troublemakers)
IR_data_all_rev <- IR_data_all %>%
  filter(id != 10 &
         id != 27 &
         id != 32 &
         id != 42 &
         id != 50 &
         id != 52) %>%
  bind_rows(IR_data_10,
            IR_data_27,
            IR_data_32,
            IR_data_42,
            IR_data_50,
            IR_data_52,
            IR_data_55,
            IR_data_56,
            IR_data_57,
            IR_data_58,
            IR_data_59)%>%
  arrange(id)  #problems with that paper which only has two treatments (we need three to gnls operation)
# let's check out if duplicate temperature treatments have been summarised
remaining_troublemakers <- IR_data_all_rev %>% 
  group_by(id, temperature) %>%
  count() %>% 
  filter(n >1) %>% 
  print()
#none is problematic :_)
# thus we can rename the dataset and correct NAs at year variable
IR_data_year_corr<- IR_data_all_rev %>%  #recover original name of the tibble
  mutate(across(Year, ~replace_na(.,2021))) %>%  #it is 2021
  glimpse()

#now we summarise an unique median standard deviation for each study for later weighting
IR_data_all <- IR_data_year_corr %>% 
  group_by(id) %>% 
  mutate(sd = median(stdev)) %>% 
  rename(sd_treat = stdev,
         sd_median = sd) %>% 
  print()

distinct_ids <- IR_data_all %>%  
  distinct(id)
ID <- rep(1:length(distinct_ids$id), each = 100)
myList <- tibble(a_est = rep(NA,length(ID)), #create a list to use as replacement of NAs in dplyr format
                 a_se = rep(NA,length(ID)),
                 Tmin_est = rep(NA,length(ID)),
                 Tmin_se = rep(NA,length(ID)), 
                 Tmax_est = rep(NA,length(ID)),
                 Tmax_se = rep(NA,length(ID)), 
                 Topt_est = rep(NA,length(ID)),
                 Topt_se = rep(NA,length(ID)),
                 starting_a = rep(NA,length(ID)),
                 starting_Tmin = rep(NA,length(ID)),
                 starting_Tmax = rep(NA,length(ID)))

params_br1_individual <- tibble(id=ID,myList)

for (i in unique(ID)){ 
  IR_data_ID <- IR_data_all %>% 
    filter(id==i)
  png(filename = paste0("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/synchro_github_ir/intrinsic_rates_pests/data_",i,".png"))
  plot(IR_data_ID$temperature,IR_data_ID$growth_rate) 
  dev.off()
  iterative_simul_ID <- params_br1_individual %>% 
    filter(id == i)
    
  for (nrep in 1:100){
    cat(paste(paste("Study",i,"/",length(unique(ID))),
              paste("simulation",nrep,"/",100),
              paste("total progress:",(100*unique(ID)[i-1]+nrep/100)/length(unique(ID)),"%"),
              sep = "\n"))
    temp_ID <- IR_data_ID %>% filter(id==i) %>% select(temperature)  # first we assume normal distribution
    simul_ID <- tibble(id =rep(i,IR_data_ID %>% filter(id==i) %>% select(n_1) %>% summarise(n_1=sum(n_1))),
                       "temp"=0,
                       "r"=0)
    enes <-IR_data_ID %>% filter(id==i) %>% select(n_1)
    position <- cumsum(enes)-enes[1]+1
    iter_temp <- 1:length(temp_ID$temperature)
    for (num in iter_temp){
      simul_ID$temp[position$n_1[num]:(position$n_1[num]-1+enes$n_1[num])] <- rep(temp_ID$temperature[num], each = enes$n_1[num])
    }
    simul_ID
    
    for(t in 1:length(temp_ID$temperature)){
      
      temper <- as.numeric(temp_ID[t,])
      n <-  as.numeric(IR_data_ID %>% filter(id==i & temperature==temper)%>% select(n_1))
      mu <- as.numeric(IR_data_ID %>% filter(id==i & temperature==temper)%>% select(growth_rate))
      sd <-as.numeric(IR_data_ID %>% filter(id==i & temperature==temper)%>% select(stdev))
      sim_r <- rnorm(n,mu,sd)
      simul_ID[simul_ID$temp == temper,"r"] <- tibble(sim_r) 
    }
    simul_ID
    #write_csv(simul_ID,file = paste0("/Users/Ecologia/Desktop/DAR?O_actualizada septiembre 2021/Intrinsic_metaanalysis/simulations/simulated_data_study_id_",i,".csv"))
    grid_br1_ID <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                                    Tmin=seq(-5,21.5,by=1),
                                    Tmax=seq(25.5,48,by=1)))
    capture.output(type="message",
                   fitted_br1_ID_brute <- try(nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                               data = simul_ID,
                                               start = grid_br1_ID,
                                               algorithm = "brute-force",
                                               trace = FALSE),
                                    silent=TRUE)
    )
    sum_grid_ID <- summary(fitted_br1_ID_brute) #save the summary of this first scan
    starVals_ID <- sum_grid_ID$parameters[,1] #these are the starting values for
    print("fitting model begins")
    skip_to_next <- FALSE
    tryCatch({
      fitted_br1_ID_gnls <- gnls(r ~ briere1(a,temp = temp,Tmin,Tmax),
                                 data = simul_ID,
                                 start = starVals_ID,
                                 weights = varExp(),
                                 control = gnlsControl(nlsTol = 1e-2))},
      error = function(e){ skip_to_next <<- TRUE})
    if(skip_to_next){ next }
    sum_br1_ID <- summary(fitted_br1_ID_gnls)
    coefs_ID <- data.frame(coef(sum_br1_ID)[,1:2],row.names = NULL)
    Topt_est_ID <- Topt(Tmin=coefs_ID[2,1],
                        Tmax=coefs_ID[3,1],
                        m=2)
    Topt_se_ID <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                              coef(fitted_br1_ID_gnls), vcov(fitted_br1_ID_gnls))
    id <- i
    a_est <-     coef(sum_br1_ID)[1,1]
    a_se <-      coef(sum_br1_ID)[1,2]
    Tmin_est <-  coef(sum_br1_ID)[2,1]
    Tmin_se <-   coef(sum_br1_ID)[2,2]
    Tmax_est <-  coef(sum_br1_ID)[3,1]
    Tmax_se <-   coef(sum_br1_ID)[3,2]
    Topt_est <-  Topt_est_ID
    Topt_se <-   Topt_se_ID
    starting_a <-  starVals_ID[1]
    starting_Tmin <- starVals_ID[2]
    starting_Tmax <-  starVals_ID[3]
    myList <- tibble(id = i,
                     a_est, #create a list to use as replacement of NAs in dplyr format
                     a_se,
                     Tmin_est,
                     Tmin_se, 
                     Tmax_est,
                     Tmax_se, 
                     Topt_est,
                     Topt_se,
                     starting_a,
                     starting_Tmin,
                     starting_Tmax)
    iterative_simul_ID[nrep,] <- tibble(myList)
    
    
  }
  
  params_br1_individual[((i-1)*100+1):(i*100),] <- iterative_simul_ID
  }

params_br1_individual 
View(params_br1_individual)
write_csv(params_br1_individual,"repeated_simul_parameters.csv")



## do not run from here

#### _ _ 7.2. Dataset traits ensemble ####
# let's check if id number is correct
IR_data_all %>% distinct(id) #yes it is
# let's extract non-numneric vars and group acari into one order
acari_data <- IR_data_all %>% 
  filter(order == "Acari>Prostigmata" |
           order == "Acari>Trombidiformes") %>% 
  mutate(order = "Acari")
non_acari <- IR_data_all %>% 
  filter(order != "Acari>Prostigmata" &
           order != "Acari>Trombidiformes")
IR_data_covs <- acari_data %>% 
  bind_rows(non_acari) 

data4params <- IR_data_covs %>%
  select(Authors,order,family,genus,species,feeding_guild,
         lat,lon,id) %>%
  mutate(spp = paste(genus, species)) %>% 
  group_by_all() %>%
  summarise(id=unique(id)) %>%
  arrange(id) %>% 
  mutate(Authors = word(Authors,1,2)) %>% 
  relocate(id, Authors, order, family, spp, feeding_guild, lat, lon) %>% 
  glimpse()

thermal_traits_indiv <- data4params %>%
  select(-id) %>% 
  bind_cols(params_br1_individual) %>% 
  filter(id != 19 &
         id != 45) %>% 
  glimpse()
write_csv(thermal_traits_indiv,"thermal_traits_individual.csv")
#### _ _ 7.3. exploratory ####
boxplot(thermal_traits_indiv$Tmin_est,thermal_traits_indiv$Tmax_est)
range(thermal_traits_indiv$Tmin_est)
range(thermal_traits_indiv$Tmax_est)
lat_tmin <- ggplot(thermal_traits_indiv, aes(abs(lat), Tmin_est))+
  geom_point(color = "lightblue2", alpha = 0.7)+
  geom_smooth(method = "lm")+
  theme_classic()

# remove the outliers
thermal_traits_indiv_filter <- thermal_traits_indiv %>% 
  filter(Tmin_est >-10 &
         Tmax_est < 50)
boxplot(thermal_traits_indiv_filter$Tmin_est,thermal_traits_indiv_filter$Tmax_est)
range(thermal_traits_indiv_filter$Tmin_est)
range(thermal_traits_indiv_filter$Tmax_est)
lat_tmin_filter <- ggplot(thermal_traits_indiv_filter, aes(abs(lat), Tmin_est))+
  geom_point(color = "lightblue2", alpha = 0.7)+
  geom_smooth(method = "lm", 
              color = "lightblue4",
              fill = "lightblue1")+
  theme_classic()
lat_tmin_filter
tmin_to_lat_filter_lm <- lm(Tmin_est ~ abs(lat),
                            data = thermal_traits_indiv_filter)
summary(tmin_to_lat_filter_lm)
lat_tmax_filter <- ggplot(thermal_traits_indiv_filter, aes(abs(lat), Tmax_est))+
  geom_point(color = "firebrick2", alpha = 0.7)+
  geom_smooth(method = "lm", 
              color = "firebrick4",
              fill = "firebrick1")+
  theme_classic()
lat_tmax_filter
tmax_to_lat_filter_lm <- lm(Tmax_est ~ abs(lat),
                            data = thermal_traits_indiv_filter)
summary(tmax_to_lat_filter_lm)


#### _ _ extra: possibilities for error at individual loop ####

#  for (i in 1:1000){
#    while(TRUE){
#      df <- try(downloadfnc("URL", file = i), silent=TRUE)
#      if(!is(df, 'try-error')) break
#    }
#    table[i,] <- df
#  }

#tryCatch
# https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r/12195574#12195574 



#### Appendix I. Decisions to exclude papers ####
#ID 19 <- no explanation nor discussion of such those high values in rm.
#         methods seem adequate, n's are high, seems like typo error. Almost
#         two orders of magnitude above r values.
# http://www.scielo.org.ar/pdf/rsea/v78n4/v78n4a05.pdf 

# ID 37 <- no apparent problems (there is only one slightly high value of 
# error in the 25?C treatment, but within the same order of magnitude)
# https://journals.flvc.org/flaent/article/view/83895/80785

# ID 47 <- some problems:
#           1) three order-magnitude higher errors for two first treatments 
#           2) zero-rounded error for other variables (...)
#           3) they do not discuss the errors
#  mat & meth seem okay but I would remove this paper...

