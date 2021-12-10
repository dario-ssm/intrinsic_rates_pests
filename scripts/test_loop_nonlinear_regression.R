#### Meta-Analysis: nonlinear mixed effects meta-regression  ####

#### SCRIPT INFO ####     
#     Authors: Darío San Segundo Molina, Sara Villén Pérez, Ignacio Morales Castilla
#     Title: Brière-1 Modeling
#     Aim: test model fitting to each study to obtain ecologically-informative parameters (traits)
#     Date: December 2021
#     
#     Workflow: 
# if you have "IR_data_complete_sim.csv", jump directly to section 4
# otherwise:
#       1. Run Section 1
#       2. Run Section 2
#       3. Run Section 3
#        
#  In both cases, continue running 4 c)
#_________________________ ####


#### 1. Load Dataset ####
rm(list=ls())
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(svglite)
library(nlstools)
library(tidyr)
library(nls2)
library(msm)
library(magrittr)
library(cowplot)
library(car)
library(nlme)
library(brms)
#load data:
IR_data <-read_csv("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/IR_metaanalysis_suitable.csv") %>%
  filter(Filter_3 == "yes") %>% #select only those which have passed the filter 3 in the source csv
  #mutate(growth_rate = replace(growth_rate, growth_rate <= 0,0))%>% #we assign 0 to negative values of intrinsic rates (biological nonsense?)
  filter(as.numeric(temperature)<50) %>% #exclude possible missleading points
  dplyr::select(Authors,title,Filter_3,Approach,Subapproach,temperature,growth_rate,
                error,n_1,order,family,genus,species,feeding_guild,h_p_family,
                diet_family,RH,daylength,lat,lon)%>%
  glimpse()

#### 2. Define fitting function ####
# We use a Brière-1 model (Briere et al., 1999).This equation represents a trade-off 
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

#### 3. Loop for simulation csv generator & export ####
#### _ _ _ 3.1. Data cleaning and manipulation ####
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
#now set wd:
setwd("C:/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/simulations") #to save plots here together
distinct_ids <- IR_data_all %>%  
  distinct(id)
#### _ _ _ 3.2. Data simulations ####
## simulated_dataframe
ID <- 1:length(distinct_ids$id)
for (i in ID){ #para cada valor entre los que están en nuestro vector ("Study" or "id"; es decir, del 1 al 56)
  
  print(ID[i])  #primero ver por qué iteración vamos
  IR_data_ID <- IR_data_all %>% #subset con el conjunto que hacemos en cada iteración
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
  png(filename = paste0("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/simulations/simulated_",i,".png"))
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
#### _ _ _ 3.3. Prepare covariates df to bind ####
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
#### _ _ _ 3.4. Ensemble sim dataset ####
IR_data_complete_sim <- simulations_numeric %>% 
  select(-id)%>%
  bind_cols(simulations_covars)%>%
  write_csv(file = "IR_data_complete_sim.csv")
#let's see sd values
boxplot(simulated_IR_data$sd) # some values that we have to take appart...
exclusion_limit <- quantile(IR_data_complete_sim$sd,probs = .95)
#and clean
simulated_IR_data_complete <- IR_data_complete_sim %>% 
  filter(abs(r) <= 0.8 &
           temp!=0 &
           r >= 0 &
           sd <= exclusion_limit)





plot_simul <- ggplot(simulated_IR_data_complete, aes(x = temp, y = r))+
  geom_point(aes(color = as.factor(id)),
             alpha = .1,
             position = position_jitter(width = 4))+
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


#### 4. nlme model to entire simulated data ####
#### _ _ a) no negative values, manual starting values ####
# IR_data_complete_sim <- read_csv("IR_data_complete_sim.csv")

simulated_IR_data <- IR_data_complete_sim %>% 
  filter(abs(r) <= 0.8 &
           temp!=0 & 
           r >=0)
## look for appropriate starting values
grid_br1_simul <- expand.grid(list(a=seq(3e-05,1.5e-04,by=1e-05),
                                Tmin=seq(0,10,by=0.5),
                                Tmax=seq(33,48,by=0.5)))
fitted_br1_simul_brute<- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                            data = simul_IR_data,
                            start = grid_br1_simul,
                            algorithm = "brute-force",
                            trace = FALSE)
sum_brute <- summary(fitted_br1_simul_brute)
starts_all <- sum_brute$parameters[,1]
nlme_br1_all <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = starts_all,
                     fixed = a+Tmin+Tmax ~1,
                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                     #random =a+Tmin+Tmax~1|id,
                     data = simul_IR_data,
                     weights = varExp(),
                     control = nlmeControl(pnlsTol = 10,
                                           msMaxIter = 50,
                                           msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
summary(nlme_br1_all)
gnls_br1_all <- gnls(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = starts_all,
                     data = simul_IR_data,
                     weights = varExp(),
                     control = gnlsControl())
summary(gnls_br1_all)
#### _ _ b) also negative values, manual starting values ####
simulated_IR_data <- simulations_empty %>% 
  filter(abs(r) <= 0.8 &
           temp!=0)

simul_IR_data <- simulated_IR_data
## look for appropriate starting values
grid_br1_simul <- expand.grid(list(a=seq(1e-05,2e-04,by=1e-05),
                                   Tmin=seq(-5,18,by=1),
                                   Tmax=seq(30,48,by=1)))
fitted_br1_simul_brute<- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                              data = simul_IR_data,
                              start = grid_br1_simul,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_brute <- summary(fitted_br1_simul_brute)
starts_all <- sum_brute$parameters[,1]
nlme_br1_all <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = starts_all,
                     fixed = a+Tmin+Tmax ~1,
                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                     #random =a+Tmin+Tmax~1|id,
                     data = simul_IR_data,
                     weights = varExp(),
                     control = nlmeControl(pnlsTol = 10,
                                           msMaxIter = 1000,
                                           msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
summary(nlme_br1_all)
gnls_br1_all <- gnls(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = starts_all,
                     data = simul_IR_data,
                     weights = varExp(),
                     control = gnlsControl())

#### _ _ c) exclude high stdev. rows, take also negative values, manual starting values ####
simulated_IR_data <- simulations_empty %>% 
  filter(abs(r) <= 0.8 &
           temp!=0 & 
           sd <= exclusion_limit)
simul_IR_data <- simulated_IR_data
## look for appropriate starting values
grid_br1_simul <- expand.grid(list(a=seq(3e-05,1.5e-04,by=1e-05),
                                   Tmin=seq(-8,10,by=0.5),
                                   Tmax=seq(33,48,by=0.5)))
fitted_br1_simul_brute<- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                              data = simul_IR_data,
                              start = grid_br1_simul,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_brute <- summary(fitted_br1_simul_brute)
starts_all <- sum_brute$parameters[,1]
nlme_br1_all <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = starts_all,
                     fixed = a+Tmin+Tmax ~1,
                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                     #random =a+Tmin+Tmax~1|id,
                     data = simul_IR_data,
                     weights = varExp(),
                     control = nlmeControl(pnlsTol = 100,
                                           msMaxIter = 100,
                                           msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
summary(nlme_br1_all)
gnls_br1_all <- gnls(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = starts_all,
                     data = simul_IR_data,
                     weights = varExp(),
                     control = gnlsControl())
summary(gnls_br1_all)

#### _ _ d) exclude high stdev. rows, only non negative values, manual starting values ####
simulated_IR_data <- simulations_empty %>% 
  filter(abs(r) <= 0.8 &
           temp!=0 &
           r >= 0 &
           sd <= exclusion_limit)
simul_IR_data <- simulated_IR_data
## look for appropriate starting values
grid_br1_simul <- expand.grid(list(a=seq(3e-05,1.5e-04,by=1e-05),
                                   Tmin=seq(-8,10,by=0.5),
                                   Tmax=seq(33,48,by=0.5)))
fitted_br1_simul_brute<- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                              data = simul_IR_data,
                              start = grid_br1_simul,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_brute <- summary(fitted_br1_simul_brute)
starts_all <- sum_brute$parameters[,1]
nlme_br1_all <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = starts_all,
                     fixed = a+Tmin+Tmax ~1,
                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                     #random =a+Tmin+Tmax~1|id,
                     data = simul_IR_data,
                     weights = varExp(),
                     control = nlmeControl(pnlsTol = 10,
                                           msMaxIter = 100,
                                           msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
summary(nlme_br1_all)
gnls_br1_all <- gnls(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = starts_all,
                     data = simul_IR_data,
                     weights = varExp(),
                     control = gnlsControl())
summary(gnls_br1_all)


#### e) Other variables
##let's add other variables...
acari_data <- IR_data_all %>%
  filter(order == "Acari>Prostigmata" | 
           order == "Acari>Trombidiformes")%>%
  mutate(order = "Acari")
IR_data_all_non_acari <- IR_data_all %>%
  filter(order != "Acari>Prostigmata" & 
           order != "Acari>Trombidiformes") 
IR_data_all_rev <- IR_data_all_non_acari %>%
  bind_rows(acari_data) %>%
  select(order,family,genus,species,feeding_guild,lat,lon,title,Authors,id) %>% 
  group_by(id) %>% 
  summarise_all(unique)
  

#temp<-length(coef(lm(r ~ temp + order,IR_data_complete_sim)))
nlme_br1_all <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = list(fixed = c(starts_all,rep(0,18)),
                                  random = c(starts_all)),
                     fixed = a+Tmin+Tmax ~ as.factor(order),
                     random =a+Tmin+Tmax~1|id,
                     data = simulated_IR_data_complete,
                     weights = varExp(),
                     control = nlmeControl(pnlsTol = 1e-1, msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
nlme_br1_all <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = c(fixed = startvals_idea),
                     fixed = a+Tmin+Tmax ~ as.factor(order),
                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                     data = simulated_IR_data_complete,
                     weights = varExp(),
                     control = nlmeControl(pnlsTol = 1,
                                           msMaxIter = 100,
                                           msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
plot(simulated_IR_data_complete$temp,simulated_IR_data_complete$sd)
summary(nlme_br1_all)
starts_all

## let's try feeding_guild
unique(IR_data_complete_sim$feeding_guild)
startvals_idea_fg <- c(starts_all[1],
                       rep(0, 4),
                       starts_all[2],
                       rep(0, 4),
                       starts_all[3],
                       rep(0, 4))
length(startvals_idea_fg) #15 = 5 values* 3 parameters
nlme_br1_all <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                     start = c(fixed = startvals_idea_fg),
                     fixed = a+Tmin+Tmax ~ as.factor(feeding_guild),
                     groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                     data = simulated_IR_data_complete,
                     weights = varExp(),
                     control = nlmeControl(pnlsTol = 10,
                                           msMaxIter = 100,
                                           maxIter = 100,
                                           msVerbose = TRUE,
                                           ))#to avoid error of singularity in backsolve at level 0; block 1
# Error in nlme.formula(model = r ~ briere1(a, temp, Tmin, Tmax), start = c(fixed = startvals_idea_fg),  : 
#                         maximum number of iterations (maxIter = 100) reached without convergence
#                       In addition: There were 50 or more warnings (use warnings() to see the first 50)
startvals_idea <- c(0.00004,0,0,0,0,0,-5,0,0,0,0,0,0,46,0,0,0,0,0,0,0)




#### 5. Bayesian brms ####
#### _ _  a) example for a small subset ####
acari_data <- IR_data %>% 
  filter(order == "Acari>Prostigmata" |
         order == "Acari>Trombidiformes") %>% 
  mutate(order = "Acari")

nls_acari <- nls(growth_rate ~ briere1(a,temperature,Tmin,Tmax),
                 data = acari_data,
                 start = c(a = 2e-04, Tmin = 6.5, Tmax = 38))
summary(nls_acari)
prior_1 <- c(set_prior("normal(1.317e-04, (2.855e-05)^2)", nlpar = "a"),
             set_prior("normal(8.158e+00, (2.546e+00)^2)", nlpar = "Tmin"),
             set_prior("normal(4.097e+01, (1.979e+00)^2)", nlpar = "Tmax"))

briere_formula <- bf(growth_rate ~ a*temperature*(temp-Tmin)*(Tmax-temp)^(1/2),
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



#### 5. Loop for individual study model fitting ####
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
#now set wd:
setwd("C:/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/simulations") #to save plots here together
distinct_ids <- IR_data_all %>%  
  distinct(id)
ID <- 1:length(distinct_ids$id)
myList <- tibble(a_est = rep(0,length(ID)), #create a list to use as replacement of NAs in dplyr format
                 a_se = rep(0,length(ID)),
                 Tmin_est = rep(0,length(ID)),
                 Tmin_se = rep(0,length(ID)), 
                 Tmax_est = rep(0,length(ID)),
                 Tmax_se = rep(0,length(ID)), 
                 Topt_est = rep(0,length(ID)),
                 Topt_se = rep(0,length(ID)),
                 starting_a = rep(0,length(ID)),
                 starting_Tmin = rep(0,length(ID)),
                 starting_Tmax = rep(0,length(ID)))
params_br1_individual <- tibble(id=ID,myList)

for (i in ID){ #para cada valor entre los que están en nuestro vector ("Study" or "id"; es decir, del 1 al 56)
  
  print(ID[i])  #primero ver por qué iteración vamos
  IR_data_ID <- IR_data_all %>% #subset con el conjunto que hacemos en cada iteración
    filter(id==i)
  png(filename = paste0("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/Explore/explore_",i,".png"))
  plot(IR_data_ID$temperature,IR_data_ID$growth_rate) #ploteamos y guardamos un plot para tenerlos todos luego
  dev.off()
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
    sd <-as.numeric( IR_data_ID %>% filter(id==i & temperature==temper)%>% select(stdev))
    sim_r <- rnorm(n,mu,sd)
    simul_ID[simul_ID$temp == temper,"r"] <- tibble(sim_r) 
  }
  
  write_csv(simul_ID,file = paste0("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/simulations/simulated_data_study_id_",i,".csv"))
  
  print("brute-force starVals searching along grid begins")
  grid_br1_ID <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                                  Tmin=seq(-5,21.5,by=1),
                                  Tmax=seq(25.5,48,by=1)))
  fitted_br1_ID_brute <- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                              data = simul_ID,
                              start = grid_br1_ID,
                              algorithm = "brute-force",
                              trace = FALSE)
  sum_grid_ID <- summary(fitted_br1_ID_brute) #save the summary of this first scan
  starVals_ID <- sum_grid_ID$parameters[,1] #these are the starting values for
  print("fitting model begins")
  print(i)
  fitted_br1_ID_gnls <- gnls(r ~ briere1(a,temp = temp,Tmin,Tmax),
                              data = simul_ID,
                              start = starVals_ID,
                              weights = varExp(),
                             control = gnlsControl(nlsTol = 1e-2))
  sum_br1_ID <- summary(fitted_br1_ID_gnls)
  coefs_ID <- data.frame(coef(sum_br1_ID)[,1:2],row.names = NULL)
  Topt_est_ID <- Topt(Tmin=coefs_ID[2,1],
                      Tmax=coefs_ID[3,1],
                      m=2)
  Topt_se_ID <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                            coef(fitted_br1_ID_gnls), vcov(fitted_br1_ID_gnls))
  
  #output -->
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
  myList <- tibble(a_est, #create a list to use as replacement of NAs in dplyr format
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
  params_br1_individual[params_br1_individual$id==i,-1] <- tibble(myList)
}
params_br1_individual

                                 
                             