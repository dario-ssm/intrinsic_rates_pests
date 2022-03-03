#- Pooled simulations dataset model fitting ----------------------------

#     Authors: Dario San Segundo Molina, Sara Villen Perez, Ignacio Morales Castilla
#     Title: Briere-1 Modeling
#     Aim: test model fitting to each study to obtain ecologically-informative parameters (traits)
#     Date: February 2022
#     


# 1. Load dataset & config options ----------------------------------------
#rm(list=ls())
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
library(ggdark)
# ...... a) load data --------------------------------------------

setwd("~/Dario Investigacion/IntRaPest/intrinsic_rates_pests/data")
IR_data_all <- read_csv("IR_data_all_clean.csv") %>% 
  rename(r = growth_rate,
         temp = temperature)

# ...... b) simulate (Papadimitropoulou 2019) --------------------------------------------

intrapest_raw <- read_csv("IR_data_all_clean.csv") %>% 
  filter(sd_treat < 0.6) #exclude one with unusual large error treatment at id = 29
intrapest_rep <- data.frame(study = rep(intrapest_raw$id, intrapest_raw$n_1),
                            r = rep(intrapest_raw$growth_rate, intrapest_raw$n_1),
                            stdev = rep(intrapest_raw$sd_treat, intrapest_raw$n_1),
                            temp = rep(intrapest_raw$temperature, intrapest_raw$n_1),
                            order = rep(intrapest_raw$order, intrapest_raw$n_1),
                            fg = rep(intrapest_raw$feeding_guild, intrapest_raw$n_1),
                            lat = rep(abs(intrapest_raw$lat),intrapest_raw$n_1),
                            year = rep(intrapest_raw$Year, intrapest_raw$n_1),
                            vi = rep(intrapest_raw$vi, intrapest_raw$n_1))
set.seed(28022022)
intrapest_rep$r_sim = rnorm(n = nrow(intrapest_rep),
                            0,
                            1)
intrapest_rep$r_sim = rnorm(n = nrow(intrapest_rep),
                            mean = mean(intrapest_rep$r),
                            sd = mean(intrapest_rep$stdev))


sum_intrapest <- intrapest_rep %>% 
  group_by(temp, study) %>% 
  summarise(r_sum = mean(r_sim),
            sd_sum = sd(r_sim)) 

simulated_intrapest <- inner_join(intrapest_rep,sum_intrapest) %>% 
  mutate(int_rate = r + (r_sim - r_sum)*(stdev/sd_sum)) %>% 
  select(study, year, order, fg, lat, temp, int_rate, vi) %>% 
  as_tibble() %>% 
  print()

# ...... c) formulas  --------------------------------------------

briere1 <- function(a, temp, Tmin, Tmax){
  a*temp*(temp-Tmin)*(Tmax-temp)^(1/2)
}
# Since Topt is not a direct parameter of the model, it can be derived from Tmin and Tmax
# according to Marchioro & Foerster, 2011:
Topt <- function(Tmin,Tmax,m){
  Topt=((2*m*Tmax+(m+1)*Tmin)+sqrt(4*(m^2)*(Tmax^2)+((m+1)^2)*(Tmin^2)-4*(m^2)*Tmin*Tmax))/(4*m+2)
  return(Topt)
}

# 2. Exploratory Data Analysis --------------------------------------------
# ...... a) simulated dataset --------------------------------------------
#ranges
explore_r1 <- simulated_intrapest %>% 
  select(int_rate) %>% 
  as_vector()
range(explore_r1)
quantile(explore_r1, na.rm = TRUE)
scatter_sims <- ggplot(simulated_intrapest, aes(temp,int_rate))+
  geom_point(aes(color = as_factor(study)), 
             alpha = 0.32,
             position = position_jitter(width = 5))+
  ggdark::dark_theme_light()+
  labs(x = "Temperature (?C)",
       y = "Intrinsic rate of increase (r)")+
  theme(legend.position = "none")
scatter_sims
# ...... c) without simulations --------------------------------------------
#ranges
explore_r3 <- IR_data_all %>% 
  select(r) %>% 
  as_vector()
range(explore_r3)
quantile(explore_r3)

scatter_no_sims <- ggplot(IR_data_all, aes(temp,r))+
  geom_point(aes(color = as_factor(id)), 
             alpha = 0.5,
             position = position_jitter(width = 2))+
  ggdark::dark_theme_light()+
  labs(x = "Temperature (?C)",
       y = "Intrinsic rate of increase (r)")+
  theme(legend.position = "none")
scatter_no_sims
lm3 <- lm(r ~ temp, data = IR_data_all)
performance::check_model(lm3)


# 3. Model fitting: Summary Effects --------------------------------------------------------

# ...... a) simulated dataset --------------------------------------------
sims_inh <- read_csv("simulations_less_dataset.csv")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(5,21.5,by=1),
                             Tmax=seq(28,48,by=1)))

fitted_br1_ID_brute <- nls2::nls2(formula= int_rate ~ briere1(a,temp,Tmin,Tmax),
                                  data = simulated_intrapest,
                                  start = grid_br1,
                                  algorithm = "brute-force",
                                  trace = FALSE)
sum_nls2_fit3 <- summary(fitted_br1_ID_brute)
starts_nls2 <- c(a = 6e-5, Tmin = 5, Tmax = 42)

## and fit the model:

nlme_fit_sims <- nlme(int_rate ~ briere1(a,temp = temp,Tmin,Tmax),
                  start = starts_nls2,
                  fixed = a+Tmin+Tmax ~ 1, 
                  groups = ~ as.factor(study), #study level
                  weights = varComb(varFixed(value = ~vi),
                                    varExp(form = ~ temp)), # heteroscedasticity accounted for
                  data = simulated_intrapest, 
                  na.action = na.exclude,
                  control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                        maxIter =  100, #recommendation by the console
                                        pnlsTol = 10, #to achieve convergence
                                        sigma = 1 # to avoid within-studies variance modelling (for meta-analysis) 
                  )
)

sum_nlme_fit_sims <-summary(nlme_fit_sims)
traits_nlme_fit3 <- as_tibble(sum_nlme_fit_sims$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(nlme_fit_sims)[1:3,2]),
         AIC = rep(sum_nlme_fit_sims$AIC),
         BIC = rep(sum_nlme_fit_sims$BIC),
         log_lik = rep(sum_nlme_fit_sims$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_nlme_fit3,"traits_nlme_poooled_sim.txt")
# ...... c) without simulations --------------------------------------------
# ............. i) nlme() --------------------------------------------
## let's look for starting values
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                                Tmin=seq(-5,21.5,by=1),
                                Tmax=seq(25.5,48,by=1)))

fitted_br1_ID_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                  data = IR_data_all,
                                  start = grid_br1,
                                  algorithm = "brute-force",
                                  trace = FALSE)
sum_nls2_fit3 <- summary(fitted_br1_ID_brute)
starts_nls2 <- coef(sum_nls2_fit3)[,1]

## and fit the model:
nlme_fit3 <- nlme(r ~ briere1(a,temp = temp,Tmin,Tmax),
                   start = starts_nls2,
                   fixed = a+Tmin+Tmax ~ 1, 
                   groups = ~ as.factor(id), #study level
                   weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                     varExp(form = ~ temp)), # heteroscedasticity accounted for
                   data = IR_data_all, 
                   control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                         maxIter =  100, #recommendation by the console
                                         pnlsTol = 1, #to achieve convergence
                                         sigma = 1 # to avoid within-studies variance modelling (for meta-analysis) 
                                         )
)

sum_nlme_fit3 <-summary(nlme_fit3)
traits_nlme_fit3 <- as_tibble(sum_nlme_fit3$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(nlme_fit3)[1:3,2]),
         AIC = rep(sum_nlme_fit3$AIC),
         BIC = rep(sum_nlme_fit3$BIC),
         log_lik = rep(nlme_fit3$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_nlme_fit3,"traits_nlme_poooled_raw.txt")

## without random-effects
gnls_fit3 <- gnls(r ~ briere1(a, temp = temp, Tmin, Tmax),
                  start = starts_nls2,
                  weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                    varExp(form = ~ temp)), # heteroscedasticity accounted for
                  data = IR_data_all,
                  control = gnlsControl(nlsTol = 100))
gnls_fit3 #fatality

# 4. Model fitting: Covariates --------------------------------------------------------

# ...... a) simulated dataset --------------------------------------------
# ............. i) order --------------------------------------------
## let's look for starting values for each order (only Acari, Lepidoptera, Hemiptera and Diptera)
counts_order <- IR_data_all %>%
  group_by(order) %>% 
  summarise(n= length(unique(id))) %>% 
  print()
#### acari
acari <- IR_data_all %>% 
  filter(order == "Acari")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(-5,21.5,by=1),
                             Tmax=seq(25.5,48,by=1)))

acari_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                              data = acari,
                              start = grid_br1,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_acari_nls2 <- summary(acari_fit_brute)
starts_nls2_acari <- coef(sum_acari_nls2)[,1]

#### hemiptera
hemiptera <- IR_data_all %>% 
  filter(order == "Hemiptera")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(10,15,by=1),
                             Tmax=seq(32,48,by=1)))

hemiptera_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                  data = hemiptera,
                                  start = grid_br1,
                                  algorithm = "brute-force",
                                  trace = FALSE)
sum_hemiptera_nls2 <- summary(hemiptera_fit_brute)
starts_nls2_hemiptera <- coef(sum_hemiptera_nls2)[,1]

#### lepidoptera
lepidoptera <- IR_data_all %>% 
  filter(order == "Lepidoptera")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(8,21.5,by=1),
                             Tmax=seq(25.5,48,by=1)))

lepidoptera_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                    data = lepidoptera,
                                    start = grid_br1,
                                    algorithm = "brute-force",
                                    trace = FALSE)
sum_lepidoptera_nls2 <- summary(lepidoptera_fit_brute)
starts_nls2_lepidoptera <- coef(sum_lepidoptera_nls2)[,1]

## and give them appropriate starting values:
#### since there are 3 orders and 3 parameters:
#### c(a*acari, a*lepi, a*hemi, Tmin_acari, Tmin_lepi, etc)
starts_ord <- c(starts_nls2_acari[1],
                starts_nls2_lepidoptera[1],
                starts_nls2_hemiptera[1],
                starts_nls2_acari[2],
                starts_nls2_lepidoptera[2],
                starts_nls2_hemiptera[2],
                starts_nls2_acari[3],
                starts_nls2_lepidoptera[3],
                starts_nls2_hemiptera[3])

order_subset <- simulated_intrapest %>% 
  filter(order == "Acari" |
           order == "Lepidoptera" |
           order == "Hemiptera") %>% 
  glimpse()
starts_ord 
nlme_fit_sims_order <- nlme(int_rate ~ briere1(a,temp = temp,Tmin,Tmax),
                      start = starts_ord,
                      fixed = a+Tmin+Tmax ~ as.factor(order), 
                      groups = ~ as.factor(study), #study level
                      weights = varComb(varFixed(value = ~vi),
                                        varExp(form = ~ temp)), # heteroscedasticity accounted for
                      data = order_subset, 
                      na.action = na.exclude,
                      control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                            maxIter =  100, #recommendation by the console
                                            pnlsTol = 1, #to achieve convergence
                                            sigma = 1 # to avoid within-studies variance modelling (for meta-analysis) 
                      )
)
# ............. ii) feeding guild --------------------------------------------
## let's look for starting values for each order (only Borers, Suckers, Chewers)
counts_feedguild <- IR_data_all %>%
  group_by(feeding_guild) %>% 
  summarise(n= length(unique(id))) %>% 
  print()
#### borers
borers <- IR_data_all %>% 
  filter(feeding_guild == "borer")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(-5,21.5,by=1),
                             Tmax=seq(25.5,48,by=1)))

borers_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                               data = borers,
                               start = grid_br1,
                               algorithm = "brute-force",
                               trace = FALSE)
sum_borers_nls2 <- summary(borers_fit_brute)
starts_nls2_borers <- coef(sum_borers_nls2)[,1]

#### suckers
suckers <- IR_data_all %>% 
  filter(feeding_guild == "sucker")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(10,15,by=1),
                             Tmax=seq(32,48,by=1)))

suckers_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                data = suckers,
                                start = grid_br1,
                                algorithm = "brute-force",
                                trace = FALSE)
sum_suckers_nls2 <- summary(suckers_fit_brute)
starts_nls2_suckers <- coef(sum_suckers_nls2)[,1]

#### chewers
chewers <- IR_data_all %>% 
  filter(feeding_guild == "chewer")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(8,21.5,by=1),
                             Tmax=seq(25.5,48,by=1)))

chewers_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                data = chewers,
                                start = grid_br1,
                                algorithm = "brute-force",
                                trace = FALSE)
sum_chewers_nls2 <- summary(chewers_fit_brute)
starts_nls2_chewers<- coef(sum_chewers_nls2)[,1]

## now we subset for those orders
feedguild_subset <- simulated_intrapest %>% 
  filter(feeding_guild == "sucker" |
           feeding_guild == "borer" |
           feeding_guild == "chewer") %>% 
  glimpse()

## and give them appropriate starting values:
#### since there are 4 orders and 3 parameters:
#### c(a*acari, a*lepi, a*hemi, int_Tmin, etc)
starts_feedguild <- c(starts_nls2_borers[1],
                      starts_nls2_chewers[1],
                      starts_nls2_suckers[1],
                      starts_nls2_borers[2],
                      starts_nls2_chewers[2],
                      starts_nls2_suckers[2],
                      starts_nls2_borers[3],
                      starts_nls2_chewers[3],
                      starts_nls2_suckers[3])
starts_feedguild
nlme_feedguild_sim <- nlme(int_rate ~ briere1(a,temp = temp,Tmin,Tmax),
                       start = list(fixed = starts_feedguild),
                       fixed = a+Tmin+Tmax ~ as.factor(feeding_guild), 
                       groups = ~ as.factor(study), #study level
                       weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                         varExp(form = ~ temp)), # heteroscedasticity accounted for
                       data = simulated_intrapest, 
                       na.action = na.exclude,
                       control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                             maxIter =  100, #recommendation by the console
                                             pnlsTol = 10, #to achieve convergence
                                             sigma = 1 # to avoid within-studies variance modelling (for meta-analysis) 
                       )
)
summary(nlme_feedguild_sim)
# ............. iii) latitude --------------------------------------------
starts_latitude <- c(starts_nls2[1],
                     0,
                     starts_nls2[2],
                     0,
                     starts_nls2[3],
                     0)

nlme_lat <- nlme(int_rate ~ briere1(a,temp = temp,Tmin,Tmax),
                 start = list(fixed = starts_latitude),
                 fixed = a+Tmin+Tmax ~ abs(lat), 
                 groups = ~ as.factor(study), #study level
                 weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                   varExp(form = ~ temp)), # heteroscedasticity accounted for
                 data = simulated_intrapest, 
                 na.action = na.exclude,
                 control = nlmeControl(msMaxIter = 50, #recommendation by the console
                                       maxIter =  50, #recommendation by the console
                                       pnlsTol = 10, #to achieve convergence
                                       sigma = 1,
                                       msVerbose = FALSE
                                       # to avoid within-studies variance modelling (for meta-analysis) 
                 )
)
summary(nlme_lat)  #nothing is significant



# ...... c) without simulations --------------------------------------------

# ............. i) order --------------------------------------------
## let's look for starting values for each order (only Acari, Lepidoptera, Hemiptera and Diptera)
counts_order <- IR_data_all %>%
  group_by(order) %>% 
  summarise(n= length(unique(id))) %>% 
  print()
#### acari
acari <- IR_data_all %>% 
  filter(order == "Acari")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(-5,21.5,by=1),
                             Tmax=seq(25.5,48,by=1)))

acari_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                  data = acari,
                                  start = grid_br1,
                                  algorithm = "brute-force",
                                  trace = FALSE)
sum_acari_nls2 <- summary(acari_fit_brute)
starts_nls2_acari <- coef(sum_acari_nls2)[,1]

#### hemiptera
hemiptera <- IR_data_all %>% 
  filter(order == "Hemiptera")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(10,15,by=1),
                             Tmax=seq(32,48,by=1)))

hemiptera_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                              data = hemiptera,
                              start = grid_br1,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_hemiptera_nls2 <- summary(hemiptera_fit_brute)
starts_nls2_hemiptera <- coef(sum_hemiptera_nls2)[,1]

#### lepidoptera
lepidoptera <- IR_data_all %>% 
  filter(order == "Lepidoptera")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(8,21.5,by=1),
                             Tmax=seq(25.5,48,by=1)))

lepidoptera_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                  data = lepidoptera,
                                  start = grid_br1,
                                  algorithm = "brute-force",
                                  trace = FALSE)
sum_lepidoptera_nls2 <- summary(lepidoptera_fit_brute)
starts_nls2_lepidoptera <- coef(sum_lepidoptera_nls2)[,1]

## now we subset for thoes orders
order_subset <- IR_data_all %>% 
  filter(order == "Acari" |
         order == "Lepidoptera" |
         order == "Hemiptera") %>% 
  glimpse()

## and give them appropriate starting values:
#### since there are 3 orders and 3 parameters:
#### c(int_a, a*acari, a*lepi, a*hemi, int_Tmin, etc)
starts_ord <- c(starts_nls2_acari[1],
                starts_nls2_lepidoptera[1],
                starts_nls2_hemiptera[1],
                starts_nls2_acari[2],
                starts_nls2_lepidoptera[2],
                starts_nls2_hemiptera[2],
                starts_nls2_acari[3],
                starts_nls2_lepidoptera[3],
                starts_nls2_hemiptera[3])

nlme_fit3 <- nlme(r ~ briere1(a,temp = temp,Tmin,Tmax),
                  start = list(fixed = starts_ord),
                  fixed = a+Tmin+Tmax ~ as.factor(order), 
                  groups =  ~ id, #study level
                  weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                    varExp(form = ~ temp)), # heteroscedasticity accounted for
                  data = order_subset, 
                  control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                        maxIter =  100, #recommendation by the console
                                        pnlsTol = 10,
                                        sigma = 1 # to avoid within-studies variance modelling (for meta-analysis) 
                  )
)
ACF(nlme_fit3)
plot(ACF(nlme_fit3, maxLag = 10),alpha = 0.05 ) #0.41606872
nlme_fit3 <- nlme(r ~ briere1(a,temp = temp,Tmin,Tmax),
                  start = list(fixed = starts_ord),
                  fixed = a+Tmin+Tmax ~ as.factor(order), 
                  groups =  ~ id, #study level
                  weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                    varExp(form = ~ temp)), # heteroscedasticity accounted for
                  data = order_subset,
                  correlation = corAR1(0.41606872),
                  control = nlmeControl(msMaxIter = 1000, #recommendation by the console
                                        maxIter =  1000, #recommendation by the console
                                        pnlsTol = 1,
                                        sigma = 1 # to avoid within-studies variance modelling (for meta-analysis) 
                  )
)

summary(nlme_fit3)
#pdDiag(list(i ~ 1, A ~ 1))
# ............. ii) feeding guild --------------------------------------------
## let's look for starting values for each order (only Borers, Suckers, Chewers)
counts_feedguild <- IR_data_all %>%
  group_by(feeding_guild) %>% 
  summarise(n= length(unique(id))) %>% 
  print()
#### borers
borers <- IR_data_all %>% 
  filter(feeding_guild == "borer")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(-5,21.5,by=1),
                             Tmax=seq(25.5,48,by=1)))

borers_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                              data = borers,
                              start = grid_br1,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_borers_nls2 <- summary(borers_fit_brute)
starts_nls2_borers <- coef(sum_borers_nls2)[,1]

#### suckers
suckers <- IR_data_all %>% 
  filter(feeding_guild == "sucker")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(10,15,by=1),
                             Tmax=seq(32,48,by=1)))

suckers_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                  data = suckers,
                                  start = grid_br1,
                                  algorithm = "brute-force",
                                  trace = FALSE)
sum_suckers_nls2 <- summary(suckers_fit_brute)
starts_nls2_suckers <- coef(sum_suckers_nls2)[,1]

#### chewers
chewers <- IR_data_all %>% 
  filter(feeding_guild == "chewer")
grid_br1 <- expand.grid(list(a=seq(1e-05,6e-04,by=1e-05),
                             Tmin=seq(8,21.5,by=1),
                             Tmax=seq(25.5,48,by=1)))

chewers_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                                    data = chewers,
                                    start = grid_br1,
                                    algorithm = "brute-force",
                                    trace = FALSE)
sum_chewers_nls2 <- summary(chewers_fit_brute)
starts_nls2_chewers<- coef(sum_chewers_nls2)[,1]

## now we subset for those orders
feedguild_subset <- IR_data_all %>% 
  filter(feeding_guild == "sucker" |
           feeding_guild == "borer" |
           feeding_guild == "chewer") %>% 
  glimpse()

## and give them appropriate starting values:
#### since there are 4 orders and 3 parameters:
#### c(a*acari, a*lepi, a*hemi, int_Tmin, etc)
starts_feedguild <- c(starts_nls2_borers[1],
                      starts_nls2_chewers[1],
                      starts_nls2_suckers[1],
                      starts_nls2_borers[2],
                      starts_nls2_chewers[2],
                      starts_nls2_suckers[2],
                      starts_nls2_borers[3],
                      starts_nls2_chewers[3],
                      starts_nls2_suckers[3])
starts_feedguild
nlme_feedguild <- nlme(r ~ briere1(a,temp = temp,Tmin,Tmax),
                  start = list(fixed = starts_feedguild),
                  fixed = a+Tmin+Tmax ~ as.factor(feeding_guild), 
                  groups = ~ as.factor(id), #study level
                  weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                    varExp(form = ~ temp)), # heteroscedasticity accounted for
                  data = feedguild_subset, 
                  control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                        maxIter =  100, #recommendation by the console
                                        pnlsTol = 10, #to achieve convergence
                                        sigma = 1 # to avoid within-studies variance modelling (for meta-analysis) 
                  )
)
summary(nlme_feedguild)

# ............. iii) latitude --------------------------------------------
starts_latitude <- c(starts_nls2[1],
                      0,
                      starts_nls2[2],
                      0,
                      starts_nls2[3],
                      0)

nlme_lat <- nlme(r ~ briere1(a,temp = temp,Tmin,Tmax),
                       start = list(fixed = starts_latitude),
                       fixed = a+Tmin+Tmax ~ abs(lat), 
                       groups = ~ as.factor(id), #study level
                       weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                         varExp(form = ~ temp)), # heteroscedasticity accounted for
                       data = IR_data_all, 
                       control = nlmeControl(msMaxIter = 1000, #recommendation by the console
                                             maxIter =  1000, #recommendation by the console
                                             pnlsTol = 100, #to achieve convergence
                                             sigma = 1 # to avoid within-studies variance modelling (for meta-analysis) 
                       )
)
summary(nlme_lat)

## no hay manera...

# 5. Model fitting: Subsets --------------------------------------------------------
# ...... a) Each Order --------------------------------------------
# ............. i) Acari --------------------------------------------
sim_acari <- simulated_intrapest %>% 
  filter(order == "Acari")
raw_acari <- IR_data_all %>% 
  filter(order == "Acari")

sim_acari_sumef <- nlme(int_rate ~ briere1(a,temp,Tmin,Tmax),
                        fixed = a + Tmin + Tmax ~ 1,
                        groups = ~ as.factor(study),
                        start = starts_nls2_acari,
                        weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                          varExp(form = ~ temp)), # heteroscedasticity accounted for
                        data = sim_acari, 
                        na.action = na.exclude,
                        control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                              maxIter =  100, #recommendation by the console
                                              pnlsTol = 11, #to achieve convergence
                                              sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
                        )
sum_sim_acari_sumef <-summary(sim_acari_sumef)
traits_acari_sim <- as_tibble(sum_sim_acari_sumef$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(sum_sim_acari_sumef)[1:3,2]),
         AIC = rep(sum_sim_acari_sumef$AIC),
         BIC = rep(sum_sim_acari_sumef$BIC),
         log_lik = rep(sum_sim_acari_sumef$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_acari_sim, "traits_acari_sim.txt")

raw_acari_sumef <- nlme(r ~ briere1(a,temp,Tmin,Tmax),
                        fixed = a + Tmin + Tmax ~ 1,
                        groups = ~ as.factor(id),
                        start = starts_nls2_acari,
                        weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                          varExp(form = ~ temp)), # heteroscedasticity accounted for
                        data = raw_acari, 
                        control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                              maxIter =  100, #recommendation by the console
                                              pnlsTol = 1e-01, #to achieve convergence
                                              sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
)
raw_sim_acari_sumef <-summary(raw_acari_sumef)
traits_acari_raw <- as_tibble(raw_sim_acari_sumef$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(raw_sim_acari_sumef)[1:3,2]),
         AIC = rep(raw_sim_acari_sumef$AIC),
         BIC = rep(raw_sim_acari_sumef$BIC),
         log_lik = rep(raw_sim_acari_sumef$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_acari_raw, "traits_acari_raw.txt")


sim_acari_lat <- nlme(int_rate ~ briere1(a,temp,Tmin,Tmax),
                      fixed = a + Tmin + Tmax ~ abs(lat),
                      groups = ~ as.factor(study),
                      start = c(starts_nls2_acari[1],0,
                                 starts_nls2_acari[2],0,
                                 starts_nls2_acari[3],0),
                      weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                        varExp(form = ~ temp)), # heteroscedasticity accounted for
                      data = sim_acari,
                      na.action = na.exclude,
                      control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                                maxIter =  100, #recommendation by the console
                                                pnlsTol = 1e-01, #to achieve convergence
                                                sigma = 1 )
)
sum_sim_acari_lat_sumef <- summary(sim_acari_lat)
traits_acari_lat_sim <- as_tibble(sum_sim_acari_lat_sumef$tTable) %>% 
  mutate(parameter = c("a_intercept", "a_slope",
                       "Tmin_intercept", "Tmin_slope",
                       "Tmax_intercept", "Tmax_slope"),
         tau_2 = as.numeric(VarCorr(sum_sim_acari_lat_sumef)[1:6,2]),
         AIC = rep(sum_sim_acari_lat_sumef$AIC),
         BIC = rep(sum_sim_acari_lat_sumef$BIC),
         log_lik = rep(sum_sim_acari_lat_sumef$logLik)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_acari_lat_sim, "traits_acari_lat_sim.txt")

raw_acari_lat <- nlme(r ~ briere1(a,temp,Tmin,Tmax),
                      fixed = a + Tmin + Tmax ~ abs(lat),
                      groups = ~ as.factor(id),
                      start = c(starts_nls2_acari[1],0,
                                 starts_nls2_acari[2],0,
                                 starts_nls2_acari[3],0),
                      weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                        varExp(form = ~ temp)), # heteroscedasticity accounted for
                      data = raw_acari, 
                      control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                            maxIter =  100, #recommendation by the console
                                            pnlsTol = 10, #to achieve convergence
                                            sigma = 1 )
) # not converging

# ............. ii) Lepidoptera --------------------------------------------
sim_lepidoptera <- simulated_intrapest %>% 
  filter(order == "Lepidoptera")
raw_lepidoptera <- IR_data_all %>% 
  filter(order == "Lepidoptera")

sim_lepidoptera_sumef <- nlme(int_rate ~ briere1(a,temp,Tmin,Tmax),
                        fixed = a + Tmin + Tmax ~ 1,
                        groups = ~ as.factor(study),
                        start = starts_nls2_lepidoptera,
                        weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                          varExp(form = ~ temp)), # heteroscedasticity accounted for
                        data = sim_lepidoptera, 
                        na.action = na.exclude,
                        control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                              maxIter =  100, #recommendation by the console
                                              pnlsTol = 10, #to achieve convergence
                                              sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
)
sum_sim_lepidoptera_sumef <- summary(sim_lepidoptera_sumef)
traits_lepidoptera_sim <- as_tibble(sum_sim_lepidoptera_sumef$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(sum_sim_lepidoptera_sumef)[1:3,2]),
         AIC = rep(sum_sim_lepidoptera_sumef$AIC),
         BIC = rep(sum_sim_lepidoptera_sumef$BIC),
         log_lik = rep(sum_sim_lepidoptera_sumef$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_lepidoptera_sim, "traits_lepidoptera_sim.txt")


raw_lepidoptera_sumef <- nlme(r ~ briere1(a,temp,Tmin,Tmax),
                        fixed = a + Tmin + Tmax ~ 1,
                        groups = ~ as.factor(id),
                        start = starts_nls2_lepidoptera,
                        weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                          varExp(form = ~ temp)), # heteroscedasticity accounted for
                        data = raw_lepidoptera, 
                        control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                              maxIter =  100, #recommendation by the console
                                              pnlsTol = 10, #to achieve convergence
                                              sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
)
#ojo porque muchísima correlacion...
sum_raw_lepidoptera_sumef <- summary(raw_lepidoptera_sumef)
traits_lepidoptera_raw <- as_tibble(sum_raw_lepidoptera_sumef$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(sum_raw_lepidoptera_sumef)[1:3,2]),
         AIC = rep(sum_raw_lepidoptera_sumef$AIC),
         BIC = rep(sum_raw_lepidoptera_sumef$BIC),
         log_lik = rep(sum_raw_lepidoptera_sumef$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_lepidoptera_raw, "traits_lepidoptera_raw.txt")


sim_lepidoptera_lat <- nlme(int_rate ~ briere1(a,temp,Tmin,Tmax),
                          fixed = a + Tmin + Tmax ~ abs(lat),
                          groups = ~ as.factor(study),
                          start = c(starts_nls2_lepidoptera[1],0,
                                     starts_nls2_lepidoptera[2],0,
                                     starts_nls2_lepidoptera[3],0),
                          weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                            varExp(form = ~ temp)), # heteroscedasticity accounted for
                          data = sim_lepidoptera, 
                          na.action = na.exclude,
                          control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                                maxIter =  100, #recommendation by the console
                                                pnlsTol = 1e-02, #to achieve convergence
                                                sigma = 1 )
)
sum_sim_lepidoptera_lat <- summary(sim_lepidoptera_lat)
traits_lepidoptera_lat_sim <- as_tibble(sum_sim_lepidoptera_lat$tTable) %>% 
  mutate(parameter = c("a_intercept", "a_slope",
                       "Tmin_intercept", "Tmin_slope",
                       "Tmax_intercept", "Tmax_slope"),
         tau_2 = as.numeric(VarCorr(sum_sim_lepidoptera_lat)[1:6,2]),
         AIC = rep(sum_sim_lepidoptera_lat$AIC),
         BIC = rep(sum_sim_lepidoptera_lat$BIC),
         log_lik = rep(sum_sim_lepidoptera_lat$logLik)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_sim_lepidoptera_lat, "traits_sim_lepidoptera_lat.txt")


raw_lepidoptera_lat <- nlme(r ~ briere1(a,temp,Tmin,Tmax),
                      fixed = a + Tmin + Tmax ~ abs(lat),
                      groups = ~ as.factor(id),
                      start = c(starts_nls2_lepidoptera[1],0,
                                 starts_nls2_lepidoptera[2],0,
                                 starts_nls2_lepidoptera[3],0),
                      weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                        varExp(form = ~ temp)), # heteroscedasticity accounted for
                      data = raw_lepidoptera, 
                      control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                            maxIter =  100, #recommendation by the console
                                            pnlsTol = 100, #to achieve convergence
                                            sigma = 1 )
                      )
#not converging
# ............. iii) Hemiptera --------------------------------------------
sim_hemiptera <- simulated_intrapest %>% 
  filter(order == "Hemiptera")
raw_hemiptera <- IR_data_all %>% 
  filter(order == "Hemiptera")
sim_hemiptera_sumef <- nlme(int_rate ~ briere1(a,temp,Tmin,Tmax),
                            fixed = a + Tmin + Tmax ~ 1,
                            groups = ~ as.factor(study),
                            start = starts_nls2_hemiptera,
                            weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                              varExp(form = ~ temp)), # heteroscedasticity accounted for
                            data = sim_hemiptera, 
                            na.action = na.exclude,
                            control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                                  maxIter =  100, #recommendation by the console
                                                  pnlsTol = 1e-01, #to achieve convergence
                                                  sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
)   
sum_sim_hemiptera_sumef <- summary(sim_hemiptera_sumef)
traits_hemiptera_sim <- as_tibble(sum_sim_hemiptera_sumef$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(sum_sim_hemiptera_sumef)[1:3,2]),
         AIC = rep(sum_sim_hemiptera_sumef$AIC),
         BIC = rep(sum_sim_hemiptera_sumef$BIC),
         log_lik = rep(sum_sim_hemiptera_sumef$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_hemiptera_sim, "traits_hemiptera_sim.txt")


raw_hemiptera_sumef <- nlme(r ~ briere1(a,temp,Tmin,Tmax),
                            fixed = a + Tmin + Tmax ~ 1,
                            groups = ~ as.factor(id),
                            start = c(a=8e-05,Tmin= 8,Tmax= 40),
                            weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                              varExp(form = ~ temp)), # heteroscedasticity accounted for
                            data = raw_hemiptera, 
                            control = nlmeControl(msMaxIter = 1000, #recommendation by the console
                                                  maxIter =  100, #recommendation by the console
                                                  pnlsTol = 100, #to achieve convergence
                                                  sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
)
#not converging
sim_hemiptera_lat <- nlme(int_rate ~ briere1(a,temp,Tmin,Tmax),
                          fixed = a + Tmin + Tmax ~ abs(lat),
                          groups = ~ as.factor(study),
                          start = c(starts_nls2_hemiptera[1],0,
                                     starts_nls2_hemiptera[2],0,
                                     starts_nls2_hemiptera[3],0),
                          weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                            varExp(form = ~ temp)), # heteroscedasticity accounted for
                          data = sim_hemiptera,
                          na.action = na.exclude,
                          control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                                maxIter =  100, #recommendation by the console
                                                pnlsTol = 100, #to achieve convergence
                                                sigma = 1)
                          )
#not converging

raw_hemiptera_lat <- nlme(r ~ briere1(a,temp,Tmin,Tmax),
                      fixed = a + Tmin + Tmax ~ abs(lat),
                      groups = ~ as.factor(id),
                      start = c(starts_nls2_hemiptera[1],0,
                                 starts_nls2_hemiptera[2],0,
                                 starts_nls2_hemiptera[3],0),
                      weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                        varExp(form = ~ temp)), # heteroscedasticity accounted for
                      data = raw_hemiptera, 
                      control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                            maxIter =  100, #recommendation by the console
                                            pnlsTol = 10, #to achieve convergence
                                            sigma = 1 )
)
# not converging

# ...... a) Each Feeding Guild --------------------------------------------
# ............. i) borers --------------------------------------------
sim_borer <- simulated_intrapest %>% 
  filter(fg == "borer")
raw_borer <- IR_data_all %>% 
  filter(feeding_guild == "borer")

sim_borer_sumef <- nlme(int_rate ~ briere1(a,temp,Tmin,Tmax),
                        fixed = a + Tmin + Tmax ~ 1,
                        groups = ~ as.factor(study),
                        start = starts_nls2_borers,
                        weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                          varExp(form = ~ temp)), # heteroscedasticity accounted for
                        data = sim_borer, 
                        na.action = na.exclude,
                        control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                              maxIter =  100, #recommendation by the console
                                              pnlsTol = 1e-02, #to achieve convergence
                                              sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
)
sum_sim_borer_sumef <- summary(sim_borer_sumef)
#high correlation...
traits_borer_sim <- as_tibble(sum_sim_borer_sumef$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(sum_sim_borer_sumef)[1:3,2]),
         AIC = rep(sum_sim_borer_sumef$AIC),
         BIC = rep(sum_sim_borer_sumef$BIC),
         log_lik = rep(sum_sim_borer_sumef$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_borer_sim, "traits_borer_sim.txt")


raw_borer_sumef <- nlme(r ~ briere1(a,temp,Tmin,Tmax),
                        fixed = a + Tmin + Tmax ~ 1,
                        groups = ~ as.factor(id),
                        start = starts_nls2_borers,
                        weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                          varExp(form = ~ temp)), # heteroscedasticity accounted for
                        data = raw_borer, 
                        control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                              maxIter =  100, #recommendation by the console
                                              pnlsTol = 1e-01, #to achieve convergence
                                              sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
)
#ojo porque completa correlación...

# ............. ii) chewers --------------------------------------------
sim_chewer <- simulated_intrapest %>% 
  filter(fg == "chewer")
raw_chewer <- IR_data_all %>% 
  filter(feeding_guild == "chewer")

sim_chewer_sumef <- nlme(int_rate ~ briere1(a,temp,Tmin,Tmax),
                         fixed = a + Tmin + Tmax ~ 1,
                         groups = ~ as.factor(study),
                         start = starts_nls2_chewers,
                         weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                           varExp(form = ~ temp)), # heteroscedasticity accounted for
                         data = sim_chewer,
                         na.action = na.exclude,
                         control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                               maxIter =  100, #recommendation by the console
                                               pnlsTol = 10, #to achieve convergence
                                              sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
)
sum_sim_chewer_sumef <- summary(sim_chewer_sumef)
#ojo alta correlacion
traits_chewer_sim <- as_tibble(sum_sim_chewer_sumef$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(sum_sim_borer_sumef)[1:3,2]),
         AIC = rep(sum_sim_borer_sumef$AIC),
         BIC = rep(sum_sim_borer_sumef$BIC),
         log_lik = rep(sum_sim_borer_sumef$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_chewer_sim, "traits_chewer_sim.txt")


raw_chewer_sumef <- nlme(r ~ briere1(a,temp,Tmin,Tmax),
                        fixed = a + Tmin + Tmax ~ 1,
                        groups = ~ as.factor(id),
                        start = starts_nls2_chewers,
                        weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                          varExp(form = ~ temp)), # heteroscedasticity accounted for
                        data = raw_chewer, 
                        control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                              maxIter =  100, #recommendation by the console
                                              pnlsTol = 1e-02, #to achieve convergence
                                              sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
)
sum_raw_chewers <- summary(raw_chewer_sumef) #high correlations
traits_chewer_raw <- as_tibble(sum_raw_chewers$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(sum_raw_chewers)[1:3,2]),
         AIC = rep(sum_raw_chewers$AIC),
         BIC = rep(sum_raw_chewers$BIC),
         log_lik = rep(sum_raw_chewers$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_chewer_raw, "traits_chewer_raw.txt")


# ............. iii) suckers --------------------------------------------
sim_sucker <- simulated_intrapest %>% 
  filter(fg == "sucker")
raw_sucker <- IR_data_all %>% 
  filter(feeding_guild == "sucker")

sim_sucker_sumef <- nlme(int_rate ~ briere1(a,temp,Tmin,Tmax),
                        fixed = a + Tmin + Tmax ~ 1,
                        groups = ~ as.factor(study),
                        start = starts_nls2_suckers,
                        weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                          varExp(form = ~ temp)), # heteroscedasticity accounted for
                        data = sim_sucker, 
                        na.action = na.exclude,
                        control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                              maxIter =  100, #recommendation by the console
                                              pnlsTol = 10, #to achieve convergence
                                              sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
)
sum_sim_sucker_sumef <-summary(sim_sucker_sumef)
traits_sucker_sim <- as_tibble(sum_sim_sucker_sumef$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(sum_sim_sucker_sumef)[1:3,2]),
         AIC = rep(sum_sim_sucker_sumef$AIC),
         BIC = rep(sum_sim_sucker_sumef$BIC),
         log_lik = rep(sum_sim_sucker_sumef$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_sucker_sim, "traits_sucker_sim.txt")





raw_sucker_sumef <- nlme(r ~ briere1(a,temp,Tmin,Tmax),
                        fixed = a + Tmin + Tmax ~ 1,
                        groups = ~ as.factor(id),
                        start = starts_nls2_suckers,
                        weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                          varExp(form = ~ temp)), # heteroscedasticity accounted for
                        data = raw_sucker, 
                        control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                              maxIter =  100, #recommendation by the console
                                              pnlsTol = 1, #to achieve convergence
                                              sigma = 1 )# to avoid within-studies variance modelling (for meta-analysis) 
                        
)
sum_raw_sucker_sumef <- summary(raw_sucker_sumef)
traits_sucker_raw <- as_tibble(sum_raw_sucker_sumef$tTable) %>% 
  mutate(parameter = c("a", "Tmin", "Tmax"),
         tau_2 = as.numeric(VarCorr(sum_raw_sucker_sumef)[1:3,2]),
         AIC = rep(sum_raw_sucker_sumef$AIC),
         BIC = rep(sum_raw_sucker_sumef$BIC),
         log_lik = rep(sum_raw_sucker_sumef$logLik,3)) %>% 
  rename(value = Value, 
         se = Std.Error, 
         df = DF, 
         t = `t-value`,
         p = `p-value`) %>% 
  relocate(parameter, value, se, tau_2, AIC, BIC, log_lik, df, t, p) %>% 
  print()
write_csv(traits_sucker_raw, "traits_sucker_raw.txt")






# 6. No temperature: --------------------------------------------------------
# .... a) r ~ order --------------------------------------------------------
r_order_scatterplot <-ggplot(IR_data_all, aes(abs(lat), r))+
  geom_point(aes(color = order))+
  geom_smooth(method = "lm", color = "turquoise3")+
  facet_wrap(.~order)+
  theme_half_open()
r_order_scatterplot_pooled <-ggplot(IR_data_all, aes(abs(lat), r))+
  geom_point(aes(color = order))+
  geom_smooth(method = "lm", color = "turquoise3")+
  theme_half_open()
# higher latitutes, larger body-sizes
lme_r_lat <- lme(r ~ abs(lat),
                 random = ~ 1|id, #study level
                weights = varFixed(value = ~vi), # inverse-variance weighting
                data = IR_data_all,
                control = lmeControl(sigma = 1))
lme_r_lat_varExp <- lme(r ~ abs(lat),
                        random = ~ 1|id, #study level
                        weights = varComb(varFixed(value = ~vi),
                                          varExp(form= ~temp)),# inverse-variance weighting
                        data = IR_data_all,
                        control = lmeControl(sigma = 1))

lme_r_lat_temp_varExp <- lme(r ~ abs(lat) + temp,
                             random = ~ 1|id, #study level
                             weights = varComb(varFixed(value = ~vi),
                                               varExp(form= ~temp)),# inverse-variance weighting
                             data = IR_data_all,
                             control = lmeControl(sigma = 1))
                                
lme_r_lat_temp_varExp_order <- lme(r ~ abs(lat) + temp + order,
                             random = ~ 1|id, #study level
                             weights = varComb(varFixed(value = ~vi),
                                               varExp(form= ~temp)),# inverse-variance weighting
                             data = IR_data_all,
                             control = lmeControl(sigma = 1))

anova(lme_r_lat, lme_r_lat_varExp)

# 7. Traits plotting ----------------
# .... a) Tmin & Tmax by order ----------------
# first we'll need to assemble de dataset with model parameters:
counts_order <- simulated_intrapest %>% 
  filter(order == "Acari" |
         order == "Lepidoptera" |
         order == "Hemiptera") %>% 
  group_by(order) %>% 
  summarise(n = length(unique(study))) %>% 
  select(n) %>% 
  as_vector() %>% 
  print()

thermal_traits_order <- traits_acari_sim %>% 
  bind_rows(traits_hemiptera_sim) %>% 
  bind_rows(traits_lepidoptera_sim) %>% 
  mutate(order = rep(c("Acari","Hemiptera","Lepidoptera"), each =3),
         n = rep(counts_order, each = 3)) %>% 
  filter(parameter != "a") %>% 
  print()

plot_thermal_traits_order <- ggplot(thermal_traits_order, aes(x = order, y = value))+
  geom_point(aes(color = parameter,
                  size = n))+
  geom_pointrange(aes(x = order,
                      y=value,
                      ymin = value-se,
                      ymax = value + se,
                      color = parameter))+
  theme_classic()
plot_thermal_traits_order
ggsave("plot_thermal_traits_order.png",
       dpi = 300,
       width = 20,
       height = 15,
       units = "cm")

# .... b) Tmin ~ lat for acari ----------------
sim_acari_lat$groups


latvals_acari <- sim_acari %>%
  group_by(study) %>% 
  summarise(lat = mean(lat)) 
weights_acari <- sim_acari %>% 
  group_by(study) %>% 
  summarise(weights = mean(1/vi))


sum_eff_acari_lat <- as_tibble(sum_sim_acari_lat_sumef$tTable[,1:2]) %>%
  t()
colnames(sum_eff_acari_lat) <-c("a_int","a_slope",
                                    "tmin_int","tmin_slope",
                                    "tmax_int","tmax_slope")
#random_sumeff_acari_lat <- VarCorr(sim_acari_lat)[1:6,2]
#names(random_sumeff_acari_lat) <- c("random_a_int","random_a_slope",
 #                                      "random_tmin_int","random_tmin_slope",
  #                                     "random_tmax_int","random_tmax_slope")
coefs_lat_acari <- coefs_acari %>%
  mutate(lat = latvals_acari$lat,
         weights = weights_acari$weights,
         study = weights_acari$study) %>%
  bind_cols(ranefs_acari) %>% 
  print()

fixeff_sum <- intervals(sim_acari_lat,which = "fixed")
sum_effs <- tibble(estimate = fixeff_sum$fixed[,2],
                   lower_ci = fixeff_sum$fixed[,1],
                   upper_ci = fixeff_sum$fixed[,3],
                   parameter = c("a_int","a_slope",
                                 "tmin_int","tmin_slope",
                                 "tmax_int","tmax_slope"))
tmin_coefs_lat_acari <- coefs_lat_acari %>% 
  select(tmin_slope, lat, weights, study, random_tmin_slope) %>% 
  mutate(study = as_factor(study))
  summary_tmin <- tibble(sum_effs$estimate[4],sum_effs$lower_ci[4],sum_effs$upper_ci[4],as.factor(43),0)
colnames(summary_tmin) <- colnames(tmin_coefs_lat_acari)

forest_tmin_lat_acari <- ggplot(tmin_coefs_lat_acari, aes(x = tmin_slope,
                                                         y = study))+
  geom_point(aes(color = as_factor(study),
                 size = weights))+
  geom_pointrange(aes(xmin = tmin_slope - random_tmin_slope,
                      xmax = tmin_slope + random_tmin_slope,
                      color = study))+
  geom_point(data = summary_tmin, aes(x = tmin_slope, y = study),
             color = "firebrick3",
             size = 2)+
  geom_pointrange(data = summary_tmin, aes(xmin = lat,
                                           xmax = weights),
                  color = "firebrick3",
                  size = 1.5)+
  labs(title = "Forest plot Acari slopes",
       x = "Tmin ~ lat (slope)",
       y = "Study")+
  geom_vline(xintercept = 0,
             color= "lightcoral",
             linetype = "dashed")+
  theme_cowplot()+
  theme(legend.position = "none")
forest_tmin_lat_acari
ggsave("forest_tmin_lat_acari.png",
       dpi = 300,
       width = 20,
       height = 15,
       units = "cm")
rectas_acari_lat <- sum_effs %>% 
  t() %>% 
  as_tibble() %>% 
  mutate_all(as.numeric) %>% 
colnames(rectas_acari_lat) <- sum_effs$parameter
forest_sum_fixeff_acari_lat <- ggplot(slopes_acari_lat, aes(estimate, parameter))+
  geom_point(aes(color = parameter), size= 5)+
  geom_pointrange(aes(xmin = lower_ci,
                      xmax = upper_ci,
                      color = parameter),
                  size = 1)+
  labs(x = "Slope (thermal traits ~ latitude)",
       y = "Thermal trait")+
  theme_cowplot()+
  scale_colour_manual(values = c("goldenrod","firebrick3","turquoise4"))+
  geom_vline(xintercept = 0,
             color= "lightcoral",
             linetype = "dashed")
  
forest_sum_fixeff_acari_lat
ggsave("forest_sum_fixeff_acari_lat.png",
       dpi = 300,
       width = 20,
       height = 15,
       units = "cm")
# 7. Orange example? ----------------
Orange
# let's generate a random covariate such as origin
source_region <- rep(c("Valencia","Murcia","Almeria", "Zaragoza","Alicante"), each = 7 )
orange_cov <- Orange %>% 
  mutate(region = source_region) %>% 
  glimpse()

f1 <- circumference ~ phi1 / (1 + exp(-(age - phi2)/phi3))
n1 <- nls(f1,
          data = orange_cov,
          start = list(phi1 = 200, phi2 = 700, phi3 = 350))
n2 <- nlme(f1,
           data = orange_cov,
           fixed = phi1 + phi2 + phi3 ~ as.factor(region),
           random = phi1 ~ 1,
           groups = ~ Tree,
           start = list(fixed = rep(c(200,700, 350), 5)),
           control = nlmeControl(pnlsTol = 10))
n2
