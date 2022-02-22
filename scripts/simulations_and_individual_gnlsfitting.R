
#- Simulations & generalized non-linear fitting ----------------------------

#     Authors: Dario San Segundo Molina, Sara Villen Perez, Ignacio Morales Castilla
#     Title: Briere-1 Modeling
#     Aim: test model fitting to each study to obtain ecologically-informative parameters (traits)
#     Date: February 2022
#     


# 1. Load dataset & config options ----------------------------------------
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
setwd("~/Dario Investigacion/IntRaPest/intrinsic_rates_pests/data")
IR_data <-read_csv("IR_metaanalysis_suitable.csv") %>%
  filter(Filter_3 == "yes") %>% #select only those which have passed the filter 3 in the source csv
  #mutate(growth_rate = replace(growth_rate, growth_rate <= 0,0))%>% #we assign 0 to negative values of intrinsic rates (biological nonsense?)
  filter(as.numeric(temperature)<50) %>% #exclude possible missleading points
  dplyr::select(Authors,Year,title, DOI,Filter_3,Approach,Subapproach,temperature,growth_rate,
                error,n_1,order,family,genus,species,feeding_guild,h_p_family,
                diet_family,RH,daylength,lat,lon)%>%
  glimpse()

# 2. Define fitting function for modelling -----------------------------------

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

# 3. Model fitting ------------------------------

# .... a) Data Cleaning & Manipulation -----------------------------------


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
  mutate(across(Year & where(is.numeric),
                ~ case_when(DOI != "10.25085/rsea.780405" ~ replace_na(.,2021)))) %>% # poner 2021 en el primero
  mutate(across(Year, ~ replace_na(.,2019))) %>% #poner 2019 en el otro
  glimpse()

#now we summarise an unique median standard deviation for each study for later weighting
IR_data_all <- IR_data_year_corr %>% 
  group_by(id) %>% 
  mutate(sd = median(stdev)) %>% 
  rename(sd_treat = stdev,
         sd_median = sd,
         vi = error) %>% 
  filter(id != 14 &
           id != 36) %>%  #two papers excluded a priori for non-sense standard errors
  ungroup() %>% 
  glimpse()
write_csv(IR_data_all,"IR_data_all_clean.csv")
IR_data_all_probs <- IR_data_year_corr %>% 
  filter(DOI == "10.25085/rsea.780405" |
           DOI == "10.1007/s10340-008-0198-9") %>% 
  distinct(id)

# .... b) Individual Fitting Approach -----------------------------------

# ............ i) Loop code -----------------------------------

distinct_ids <- IR_data_all %>%  
  distinct(id) #to justify why those ids, since those DOIs are problematic ones
params_br1_individual <- tibble(a_est = rep(NULL,length(distinct_ids$id)), #create a list to use as replacement of NULLs in dplyr format
                                a_se = rep(NULL,length(distinct_ids$id)),
                                Tmin_est = rep(NULL,length(distinct_ids$id)),
                                Tmin_se = rep(NULL,length(distinct_ids$id)), 
                                Tmax_est = rep(NULL,length(distinct_ids$id)),
                                Tmax_se = rep(NULL,length(distinct_ids$id)), 
                                Topt_est = rep(NULL,length(distinct_ids$id)),
                                Topt_se = rep(NULL,length(distinct_ids$id)),
                                starting_a = rep(NULL,length(distinct_ids$id)),
                                starting_Tmin = rep(NULL,length(distinct_ids$id)),
                                starting_Tmax = rep(NULL,length(distinct_ids$id)),
)
current_nrep <- 100*(which(distinct_ids$id == i)/length(distinct_ids$id)-1/length(distinct_ids$id)) 
to_next_nrep <- 100*(which(distinct_ids$id == i)/length(distinct_ids$id)/length(distinct_ids$id)) 

set.seed(654)
for (i in distinct_ids$id){ 
  set.seed(654)
  IR_data_ID <- IR_data_all %>% 
    filter(id==i) 
  png(filename = paste0("/Users/Ecologia/Desktop/DAR?O_actualizada septiembre 2021/Intrinsic_metaanalysis/synchro_github_ir/intrinsic_rates_pests/data_",i,".png"))
  plot(IR_data_ID$temperature,IR_data_ID$growth_rate) 
  dev.off()
  myList_i <- tibble(a_est = rep(NULL,length(distinct_ids$id)), #create a list to use as replacement of NULLs in dplyr format
                     a_se = rep(NULL,length(distinct_ids$id)),
                     Tmin_est = rep(NULL,length(distinct_ids$id)),
                     Tmin_se = rep(NULL,length(distinct_ids$id)), 
                     Tmax_est = rep(NULL,length(distinct_ids$id)),
                     Tmax_se = rep(NULL,length(distinct_ids$id)), 
                     Topt_est = rep(NULL,length(distinct_ids$id)),
                     Topt_se = rep(NULL,length(distinct_ids$id)),
                     starting_a = rep(NULL,length(distinct_ids$id)),
                     starting_Tmin = rep(NULL,length(distinct_ids$id)),
                     starting_Tmax = rep(NULL,length(distinct_ids$id))
  )
  for (nrep in 1:100){
    cat(paste(paste("Study",i,"/",length(distinct_ids$id)),
              paste("simulation",nrep,"/",100),
              paste("total progress:",(current_nrep + (nrep/100)* to_next_nrep)),"%"),
        paste("=)"),
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
      sd <-as.numeric(IR_data_ID %>% filter(id==i & temperature==temper)%>% select(sd_treat))
      sim_r <- rnorm(n,mu,sd)
      simul_ID[simul_ID$temp == temper,"r"] <- tibble(sim_r) 
    }
    simul_ID
    #write_csv(simul_ID,file = paste0("simulated_data_study_id_",i,"_rep_",nrep,".csv"))
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
    print("fitting model ends")
    skip_to_next <- FALSE
    tryCatch({
      fitted_br1_ID_gnls <- gnls(r ~ briere1(a,temp = temp,Tmin,Tmax),
                                 data = simul_ID,
                                 start = starVals_ID,
                                 weights = varExp(form = ~temp),
                                 control = gnlsControl(nlsTol = 1e-07))},
      error = function(e){skip_to_next <<- TRUE})
    if(skip_to_next | is.null(fitted_br1_ID_gnls)){ next }
    sum_br1_ID <- summary(fitted_br1_ID_gnls)
    coefs_ID <- as_tibble(coef(sum_br1_ID)[,1:2])
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
    myList_nrep <- tibble(id = i,
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
                          starting_Tmax
    )
    myList_i <- myList_i %>% 
      bind_rows(myList_nrep)
  }
  params_br1_individual <- params_br1_individual %>% 
    bind_rows(myList_i)
}
params_br1_individual

# ............ ii) Dataset ensamblage from loop outputs -----------------------------------

## and count NAs
thermal_traits_raw <- read_csv("/Users/dario-ssm/Documents/Dario Investigacion/IntRaPest/intrinsic_rates_pests/data/repeated_simul_parameters.csv") %>% 
  print()
counting_nas <- thermal_traits_raw %>% 
  group_by(id) %>% 
  summarise(counting = sum(!is.na(a_est))) %>% 
  print()
convergence_number <- rep(counting_nas$counting, each = 100)
# add convergence number to the dataset
thermal_traits_simulations <- thermal_traits_raw %>% 
  bind_cols(studies_converging = convergence_number)


# let's check if id number is correct
IR_data_all %>% distinct(id) #yes it is
# let's extract non-numneric vars and group acari into one order
IR_data_covs <- IR_data_all %>%
  select(Authors,Year,DOI,order,family,genus,species,feeding_guild,
         lat,lon,sd_median,id) %>%
  mutate(spp = paste(genus, species)) %>% 
  group_by_all() %>%
  summarise(id=unique(id)) %>%
  arrange(id) %>% 
  mutate(Authors = word(Authors,1,2)) %>% 
  rename(year = Year, authors = Authors) %>% 
  relocate(id, authors, year, feeding_guild, order, family, genus, species, spp, lat, lon, sd_median) %>% 
  glimpse()

# now repeat 100 times each row
data4params <- IR_data_covs %>% 
  slice(rep(1:n(), each =100)) %>% 
  arrange(id) %>% 
  relocate() %>%
  ungroup() %>% 
  print()

#and bind datasets
glimpse(thermal_traits_simulations)
glimpse(data4params)

thermal_traits_indiv <- data4params %>%
  select(-id) %>% 
  bind_cols(thermal_traits_simulations) %>% 
  relocate(id) %>% 
  glimpse() 

write_csv(thermal_traits_indiv,"/Users/dario-ssm/Documents/Dario Investigacion/IntRaPest/intrinsic_rates_pests/data/thermal_traits_individual.csv")

# .... c) Pooled Data Approach -----------------------------------
# ........ i) Simulations loop (100 x rnorm x study) -----------------------------------
# Write again simulations using the same looping process before model fitting; we set the seed to
# ensure simulations return the same values as befores.
distinct_ids <- IR_data_all %>%  
  distinct(id) #to justify why those ids, since those DOIs are problematic ones
set.seed(654)
total_simulations <- tibble(id = rep(NULL, 0),
                            temp = rep(NULL,0),
                            r = rep(NULL, 0))
empty_simulations <- tibble(id = rep(NULL, 0),
                            temp = rep(NULL,0),
                            r = rep(NULL, 0))
for (i in distinct_ids$id){ 
  set.seed(654)
  IR_data_ID <- IR_data_all %>% 
    filter(id==i) 
  current_nrep <- 100*(which(distinct_ids$id == i)/length(distinct_ids$id)-1/length(distinct_ids$id)) 
  to_next_nrep <- 100*(which(distinct_ids$id == i)/length(distinct_ids$id)/length(distinct_ids$id)) 
  
  for (nrep in 1:100){
    cat(paste(paste("Study",i,"/",length(distinct_ids$id)),
              paste("simulation",nrep,"/",100),
              paste("total progress:",(current_nrep + (nrep/100)* to_next_nrep)),"%"),
              paste("=)"),
              sep = "\n")
    
    temp_ID <- IR_data_ID %>% filter(id==i) %>% select(temperature)  # first we assume normal distribution
    simul_ID <- tibble(id =rep(i,IR_data_ID %>% filter(id==i) %>% select(n_1) %>% summarise(n_1=sum(n_1))),
                       "temp"= 0,
                       "r"= 0,
                       "nrep" = 0)
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
      sd <-as.numeric(IR_data_ID %>% filter(id==i & temperature==temper)%>% select(sd_treat))
      sim_r <- rnorm(n,mu,sd)
      simul_ID[simul_ID$temp == temper,"r"] <- tibble(sim_r) 
    }
    empty_simulations <- empty_simulations %>% 
      bind_rows(simul_ID) 
  }
  filled_simulations_i <- empty_simulations
  
  total_simulations <- total_simulations %>% 
    bind_rows(filled_simulations_i)
}

write_csv(total_simulations, "simulations_pooled.csv")
# and add repetition number
#rep_values <- IR_data_all %>% 
#  group_by(id) %>% 
#  summarise(reps = sum(n_1))
#study_rows <- as.numeric(rep_values$reps)

# ........ ii) Simulations loop (rnorm x study) -----------------------------------

distinct_ids <- IR_data_all %>%  
  distinct(id) #to justify why those ids, since those DOIs are problematic ones
set.seed(654)
less_simulations <- tibble(id = rep(NULL, 0),
                            temp = rep(NULL,0),
                            r = rep(NULL, 0))
empty_less_simulations <- tibble(id = rep(NULL, 0),
                            temp = rep(NULL,0),
                            r = rep(NULL, 0))
for (i in distinct_ids$id){ 
  set.seed(654)
  IR_data_ID <- IR_data_all %>% 
    filter(id==i) 
  cat(paste(paste("Study",i,"/",length(distinct_ids$id)),
            paste("total progress:", 100*which(distinct_ids$id == i)/length(distinct_ids$id),"%"),
      sep = "\n"))
  
  temp_ID <- IR_data_ID %>% filter(id==i) %>% select(temperature)  # first we assume normal distribution
  simul_ID <- tibble(id =rep(i,IR_data_ID %>% filter(id==i) %>% select(n_1) %>% summarise(n_1=sum(n_1))),
                     "temp"= 0,
                     "r"= 0,
                     "nrep" = 0)
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
    sd <-as.numeric(IR_data_ID %>% filter(id==i & temperature==temper)%>% select(sd_treat))
    sim_r <- rnorm(n,mu,sd)
    simul_ID[simul_ID$temp == temper,"r"] <- tibble(sim_r) 
  }
  empty_less_simulations <- empty_less_simulations %>% 
    bind_rows(simul_ID) 
}
simulations_less_dataset <- empty_less_simulations
write_csv(simulations_less_dataset,"simulations_less_dataset.csv")


# Appendix I. Decisions to exclude papers ---------------
#ID ?? <- no explanation nor discussion of such those high values in rm.
#         methods seem adequate, n's are high, seems like typo error. Almost
#         two orders of magnitude above r values.
# http://www.scielo.org.ar/pdf/rsea/v78n4/v78n4a05.pdf 

# ID ?? <- no apparent problems (there is only one slightly high value of 
# error in the 25?C treatment, but within the same order of magnitude)
# https://journals.flvc.org/flaent/article/view/83895/80785

# ID ?? <- some problems:
#           1) three order-magnitude higher errors for two first treatments 
#           2) zero-rounded error for other variables (...)
#           3) they do not discuss the errors
#  mat & meth seem okay but I would remove this paper...

