#- Pooled simulations dataset model fitting ----------------------------

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
library(ggdark)
#load data:
setwd("~/Dario Investigacion/IntRaPest/intrinsic_rates_pests/data")
IR_data_all <- read_csv("IR_data_all_clean.csv") %>% 
  rename(r = growth_rate,
         temp = temperature)
intrapest_simulations_rep <- read_csv("~/Dario Investigacion/IntRaPest/intrinsic_rates_pests/data/simulations_pooled.csv") %>% 
  select(-nrep) %>% 
  glimpse()

intrapest_simulations_single <- read_csv("simulations_less_dataset.csv")
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
# ...... a) without repetitions dataset --------------------------------------------

#ranges
explore_r1 <- intrapest_simulations_single %>% 
  select(r) %>% 
  as_vector()
range(explore_r1)
quantile(explore_r1)
outliers <- intrapest_simulations_single %>% 
  filter(abs(r) >=1) %>% 
  glimpse()
# ...... b) repeated dataset --------------------------------------------

#ranges
explore_r2 <- intrapest_simulations_rep %>% 
  select(r) %>% 
  as_vector()
range(explore_r2)
quantile(explore_r2)
outliers <- intrapest_simulations_single %>% 
  filter(abs(r) >=1) %>% 
  glimpse()

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


# 3. Model fitting summary effects --------------------------------------------------------

# ...... a) without repetitions dataset --------------------------------------------

# ...... b) repeated dataset --------------------------------------------

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

## without random-effects
gnls_fit3 <- gnls(r ~ briere1(a, temp = temp, Tmin, Tmax),
                  start = starts_nls2,
                  weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                    varExp(form = ~ temp)), # heteroscedasticity accounted for
                  data = IR_data_all,
                  control = gnlsControl(nlsTol = 100))
gnls_fit3 #fatality

# 4. Model fitting covariates --------------------------------------------------------

# ...... a) without repetitions dataset --------------------------------------------

# ...... b) repeated dataset --------------------------------------------

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



# 5. No temperature: --------------------------------------------------------
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

w# 6. ? ----------------


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
