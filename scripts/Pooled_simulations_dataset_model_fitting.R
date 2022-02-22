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

# ............. i) nlme() --------------------------------------------
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
                             Tmin=seq(8,15,by=1),
                             Tmax=seq(32,48,by=1)))

hemiptera_fit_brute <- nls2::nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                              data = hemiptera,
                              start = grid_br1,
                              algorithm = "brute-force",
                              trace = FALSE)
sum_hemiptera_nls2 <- summary(hemiptera_fit_brute)
starts_nls2_hemiptera <- coef(sum_hemiptera_nls2)[,1]

#### lepidopteralepidoptera
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





nlme_fit3 <- nlme(r ~ briere1(a,temp = temp,Tmin,Tmax),
                  start = starts_nls2,
                  fixed = a+Tmin+Tmax ~ 1, 
                  groups = ~ as.factor(id), #study level
                  weights = varComb(varFixed(value = ~vi), # inverse-variance weighting
                                    varExp(form = ~ temp)), # heteroscedasticity accounted for
                  data = acari, 
                  control = nlmeControl(msMaxIter = 100, #recommendation by the console
                                        maxIter =  100, #recommendation by the console
                                        pnlsTol = 1, #to achieve convergence
                                        sigma = 1 # to avoid within-studies variance modelling (for meta-analysis) 
                  )
)