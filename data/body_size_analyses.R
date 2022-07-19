# SCRIPT INFO --------------------
#     Authors: Dario San Segundo Molina, Sara Villen Perez, Ignacio Morales Castilla
#     Title: trait inferences meta-analyses
#     Aim: perform meta-analysis models to parameterised dataset of thermal traits
#     Date: March 2022
#
# 



# 1. Preparation ----------------------------------------------------------


library(tidyverse)
library(emmeans)
library(nlme)
library(here)
setwd(paste0(here(),"/data"))
thermal_traits_trans <- read_csv("thermal_traits_trans.csv")
thermal_traits_trans_order <- read_csv("thermal_traits_trans_order.csv")
thermal_traits_trans_fg <- read_csv("thermal_traits_trans_fg.csv")
thermal_traits_complete <- read_csv("thermal_traits_complete.csv")
thermal_traits_complete_order <- thermal_traits_complete %>% 
  filter(order %in% c("Lepidoptera", "Acari", "Hemiptera"))

species_list <- thermal_traits_complete %>% 
  group_by(id) %>% 
  summarise(spp = unique(paste(genus, species))) %>% 
  write_csv("species_list.csv")
#log-transform body sizes
thermal_traits_trans <- thermal_traits_trans %>% 
   mutate(body_length = log10(body_length))
thermal_traits_trans_order <- thermal_traits_trans_order %>% 
   mutate(body_length = log10(body_length))


# 2. Summary effects analysis ----------------------------------------------------------

# ~~ a) tmax ~ body size --------------------------------------------------------
# all raw 
tmax_bodysize_slope <- lme(tmax ~ body_length,
                           random = ~body_length|id,
                           weights = varFixed(~vi),
                           data = thermal_traits_complete,
                           na.action = na.omit,
                           control = lmeControl(sigma = 1, 
                                                msMaxIter = 1000))
tmax_bodysize_intercept <- lme(tmax ~ body_length,
                               random = ~ 1|id,
                               weights = varFixed(~vi),
                               data = thermal_traits_complete,
                               na.action = na.omit,
                               control = lmeControl(sigma = 1, 
                                                    msMaxIter = 100))
summary(tmax_bodysize_intercept)

# log(size)
log_tmax_bodysize_intercept <- lme(tmax ~ log10(body_length),
                                   random = ~ 1|id,
                                   weights = varFixed(~vi),
                                   data = thermal_traits_complete,
                                   na.action = na.omit,
                                   control = lmeControl(sigma = 1,
                                                msMaxIter = 100))

summary(log_tmax_bodysize_intercept) # significant and negative
# larger organisms have lower tmax.
intervals_log_tmax_bodysize_intercept <- intervals(log_tmax_bodysize_intercept, which = "fixed")
intervals_log_tmax_bodysize_intercept<- as.data.frame(intervals_log_tmax_bodysize_intercept$fixed) #probably give this as "slope per standard unit")

# trans(tmax)~log(size)
log_transtmax_bodysize_intercept <- lme(tmax ~ body_length,
                                   random = ~ 1|id,
                                   weights = varFixed(~vi),
                                   data = thermal_traits_trans,
                                   na.action = na.omit,
                                   control = lmeControl(sigma = 1,
                                                        msMaxIter = 100))

summary(log_transtmax_bodysize_intercept) # not significant

# ~~ b) tmin ~ body size --------------------------------------------------------
# all raw 
tmin_bodysize_slope <- lme(tmin ~ body_length,
                           random = ~body_length|id,
                           weights = varFixed(~vi),
                           data = thermal_traits_complete,
                           na.action = na.omit,
                           control = lmeControl(sigma = 1, 
                                                msMaxIter = 1000))
tmin_bodysize_intercept <- lme(tmin ~ body_length,
                               random = ~ 1|id,
                               weights = varFixed(~vi),
                               data = thermal_traits_complete,
                               na.action = na.omit,
                               control = lmeControl(sigma = 1, 
                                                    msMaxIter = 100))
summary(tmin_bodysize_intercept)

# log(size)
log_tmin_bodysize_intercept <- lme(tmin ~ log10(body_length),
                                   random = ~ 1|id,
                                   weights = varFixed(~vi),
                                   data = thermal_traits_complete,
                                   na.action = na.omit,
                                   control = lmeControl(sigma = 1,
                                                        msMaxIter = 100))

summary(log_tmin_bodysize_intercept) # not significant
intervals_log_tmin_bodysize_intercept <- intervals(log_tmin_bodysize_intercept, which = "fixed")
intervals_log_tmin_bodysize_intercept<- as.data.frame(intervals_log_tmin_bodysize_intercept$fixed) #probably give this as "slope per standard unit")

# trans(tmin)~log(size)
log_transtmin_bodysize_intercept <- lme(tmin ~ body_length,
                                        random = ~ 1|id,
                                        weights = varFixed(~vi),
                                        data = thermal_traits_trans,
                                        na.action = na.omit,
                                        control = lmeControl(sigma = 1,
                                                             msMaxIter = 100))

summary(log_transtmin_bodysize_intercept) # not significant

# ~~ c) topt ~ body size --------------------------------------------------------
# all raw 
topt_bodysize_slope <- lme(topt ~ body_length,
                           random = ~body_length|id,
                           weights = varFixed(~vi),
                           data = thermal_traits_complete,
                           na.action = na.omit,
                           control = lmeControl(sigma = 1, 
                                                msMaxIter = 1000))
topt_bodysize_intercept <- lme(topt ~ body_length,
                               random = ~ 1|id,
                               weights = varFixed(~vi),
                               data = thermal_traits_complete,
                               na.action = na.omit,
                               control = lmeControl(sigma = 1, 
                                                    msMaxIter = 100))
summary(topt_bodysize_intercept)

# log(size)
log_topt_bodysize_intercept <- lme(topt ~ log10(body_length),
                                   random = ~ 1|id,
                                   weights = varFixed(~vi),
                                   data = thermal_traits_complete,
                                   na.action = na.omit,
                                   control = lmeControl(sigma = 1,
                                                        msMaxIter = 100))

summary(log_topt_bodysize_intercept) # significant and negative
# larger organisms have lower topt.
intervals_log_topt_bodysize_intercept <- intervals(log_topt_bodysize_intercept, which = "fixed")
intervals_log_topt_bodysize_intercept<- as.data.frame(intervals_log_topt_bodysize_intercept$fixed) #probably give this as "slope per standard unit")

# trans(topt)~log(size)
log_transtopt_bodysize_intercept <- lme(topt ~ body_length,
                                        random = ~ 1|id,
                                        weights = varFixed(~vi),
                                        data = thermal_traits_trans,
                                        na.action = na.omit,
                                        control = lmeControl(sigma = 1,
                                                             msMaxIter = 100))

summary(log_transtopt_bodysize_intercept) # not significant

# 2. Subgroup analysis ----------------------------------------------------------
## we are using the option b) above. 
## in addition, we analyse first the pooled data set ad then by each order subsets
# ~~  a) order  ----------------------------------------------------------
# ~~~~  i. tmax  ----------------------------------------------------------
## pooled
tmax_order_bodysize_log <- lme(tmax ~ order*log10(body_length),
                               random = ~ 1|id,
                               weights = varFixed(~vi),
                               data = thermal_traits_complete_order,
                               na.action = na.omit,
                               control = lmeControl(sigma = 1,
                                                    msMaxIter = 200))
summary(tmax_order_bodysize_log)  #no differences in slope between orders

## lepidoptera
lepi <- thermal_traits_complete_order %>% 
  filter(order == "Lepidoptera")
tmax_body_size_lepi <- lme(tmax ~ log10(body_length),
                           random = ~ 1|id,
                           weights = varFixed(~vi),
                           data = lepi,
                           na.action = na.omit,
                           control = lmeControl(sigma = 1,
                                                msMaxIter = 200))
summary(tmax_body_size_lepi)  #no significant effect for lepidoptera.

## acari
acari <- thermal_traits_complete_order %>% 
  filter(order == "Acari")
tmax_body_size_acari <- lme(tmax ~ log10(body_length),
                           random = ~ 1|id,
                           weights = varFixed(~vi),
                           data = acari,
                           na.action = na.omit,
                           control = lmeControl(sigma = 1,
                                                msMaxIter = 200))
summary(tmax_body_size_acari)  #even worse (tetranychus problems... :/)

## hemiptera
hemiptera <- thermal_traits_complete_order %>% 
  filter(order == "Hemiptera")
tmax_body_size_hemiptera <- lme(tmax ~ log10(body_length),
                            random = ~ 1|id,
                            weights = varFixed(~vi),
                            data = hemiptera,
                            na.action = na.omit,
                            control = lmeControl(sigma = 1,
                                                 msMaxIter = 200))

summary(tmax_body_size_hemiptera)  #no significant

# ~~~~  ii. tmin  ----------------------------------------------------------
## pooled
tmin_order_bodysize_log <- lme(tmin ~ order*log10(body_length),
                               random = ~ 1|id,
                               weights = varFixed(~vi),
                               data = thermal_traits_complete_order,
                               na.action = na.omit,
                               control = lmeControl(sigma = 1,
                                                    msMaxIter = 200))
summary(tmin_order_bodysize_log)  
#difference of body size effect on ctmin among orders
coefs_tmin_order_bodysize_log <- as_tibble(coefficients(summary(tmin_order_bodysize_log))) %>% 
  mutate(estimate = -15.048259+Value)
emmeans(tmin_order_bodysize_log, list(pairwise ~ order), adjust = "tukey")
# differences between Hemiptera and Lepidoptera; Hemiptera with higher slope.

## lepidoptera
lepi <- thermal_traits_complete_order %>% 
  filter(order == "Lepidoptera")
tmin_body_size_lepi <- lme(tmin ~ log10(body_length),
                           random = ~ 1|id,
                           weights = varFixed(~vi),
                           data = lepi,
                           na.action = na.omit,
                           control = lmeControl(sigma = 1,
                                                msMaxIter = 200))
summary(tmin_body_size_lepi)  #no significant effect for lepidoptera.

## acari
acari <- thermal_traits_complete_order %>% 
  filter(order == "Acari")
tmin_body_size_acari <- lme(tmin ~ log10(body_length),
                            random = ~ 1|id,
                            weights = varFixed(~vi),
                            data = acari,
                            na.action = na.omit,
                            control = lmeControl(sigma = 1,
                                                 msMaxIter = 200))
summary(tmin_body_size_acari)  #even worse (tetranychus problems... :/)

## hemiptera
hemiptera <- thermal_traits_complete_order %>% 
  filter(order == "Hemiptera")
tmin_body_size_hemiptera <- lme(tmin ~ log10(body_length),
                                random = ~ 1|id,
                                weights = varFixed(~vi),
                                data = hemiptera,
                                na.action = na.omit,
                                control = lmeControl(sigma = 1,
                                                     msMaxIter = 200))

summary(tmin_body_size_hemiptera)  
#significant for hemiptera and positive; larger organisms have higher tmin; 
# counterintuitive

# ~~~~  iii. tbreadth  ----------------------------------------------------------
## pooled
thermal_breadth_order_bodysize_log <- lme(thermal_breadth ~ order*log10(body_length),
                                          random = ~ 1|id,
                                          weights = varFixed(~vi),
                                          data = thermal_traits_complete_order,
                                          na.action = na.omit,
                                          control = lmeControl(sigma = 1,
                                                               msMaxIter = 200))
summary(thermal_breadth_order_bodysize_log)  #no differences in slope between orders

## lepidoptera
lepi <- thermal_traits_complete_order %>% 
  filter(order == "Lepidoptera")
thermal_breadth_body_size_lepi <- lme(thermal_breadth ~ log10(body_length),
                           random = ~ 1|id,
                           weights = varFixed(~vi),
                           data = lepi,
                           na.action = na.omit,
                           control = lmeControl(sigma = 1,
                                                msMaxIter = 200))
summary(thermal_breadth_body_size_lepi)  #no significant effect for lepidoptera.

## acari
acari <- thermal_traits_complete_order %>% 
  filter(order == "Acari")
thermal_breadth_body_size_acari <- lme(thermal_breadth ~ log10(body_length),
                            random = ~ 1|id,
                            weights = varFixed(~vi),
                            data = acari,
                            na.action = na.omit,
                            control = lmeControl(sigma = 1,
                                                 msMaxIter = 200))
summary(thermal_breadth_body_size_acari)  #even worse (tetranychus problems... :/)

## hemiptera
hemiptera <- thermal_traits_complete_order %>% 
  filter(order == "Hemiptera")
thermal_breadth_body_size_hemiptera <- lme(thermal_breadth ~ log10(body_length),
                                random = ~ 1|id,
                                weights = varFixed(~vi),
                                data = hemiptera,
                                na.action = na.omit,
                                control = lmeControl(sigma = 1,
                                                     msMaxIter = 200))

summary(thermal_breadth_body_size_hemiptera)  #no significant

# ~~~~  iv. therm_window  ----------------------------------------------------------
## pooled
therm_window_order_bodysize_log <- lme(therm_window ~ order*log10(body_length),
                                          random = ~ 1|id,
                                          weights = varFixed(~vi),
                                          data = thermal_traits_complete_order,
                                          na.action = na.omit,
                                          control = lmeControl(sigma = 1,
                                                               msMaxIter = 200))
summary(therm_window_order_bodysize_log)  #no differences in slope between orders

## lepidoptera
lepi <- thermal_traits_complete_order %>% 
  filter(order == "Lepidoptera")
therm_window_body_size_lepi <- lme(therm_window ~ log10(body_length),
                                      random = ~ 1|id,
                                      weights = varFixed(~vi),
                                      data = lepi,
                                      na.action = na.omit,
                                      control = lmeControl(sigma = 1,
                                                           msMaxIter = 200))
summary(therm_window_body_size_lepi)  #no significant effect for lepidoptera.

## acari
acari <- thermal_traits_complete_order %>% 
  filter(order == "Acari")
therm_window_body_size_acari <- lme(therm_window ~ log10(body_length),
                                       random = ~ 1|id,
                                       weights = varFixed(~vi),
                                       data = acari,
                                       na.action = na.omit,
                                       control = lmeControl(sigma = 1,
                                                            msMaxIter = 200))
summary(therm_window_body_size_acari)  #even worse (tetranychus problems... :/)

## hemiptera
hemiptera <- thermal_traits_complete_order %>% 
  filter(order == "Hemiptera")
therm_window_body_size_hemiptera <- lme(therm_window ~ log10(body_length),
                                           random = ~ 1|id,
                                           weights = varFixed(~vi),
                                           data = hemiptera,
                                           na.action = na.omit,
                                           control = lmeControl(sigma = 1,
                                                                msMaxIter = 200))

summary(therm_window_body_size_hemiptera)  #no significant

# ~~~~  v. a_est  ----------------------------------------------------------
## pooled
a_est_order_bodysize_log <- lme(a_est ~ order*log10(body_length),
                                       random = ~ 1|id,
                                       weights = varFixed(~vi),
                                       data = thermal_traits_complete_order,
                                       na.action = na.omit,
                                       control = lmeControl(sigma = 1,
                                                            msMaxIter = 200))
summary(a_est_order_bodysize_log)  #no differences in slope between orders

## lepidoptera
lepi <- thermal_traits_complete_order %>% 
  filter(order == "Lepidoptera")
a_est_body_size_lepi <- lme(a_est ~ log10(body_length),
                                   random = ~ 1|id,
                                   weights = varFixed(~vi),
                                   data = lepi,
                                   na.action = na.omit,
                                   control = lmeControl(sigma = 1,
                                                        msMaxIter = 200))
summary(a_est_body_size_lepi)  #no significant effect for lepidoptera.

## acari
acari <- thermal_traits_complete_order %>% 
  filter(order == "Acari")
a_est_body_size_acari <- lme(a_est ~ log10(body_length),
                                    random = ~ 1|id,
                                    weights = varFixed(~vi),
                                    data = acari,
                                    na.action = na.omit,
                                    control = lmeControl(sigma = 1,
                                                         msMaxIter = 200))
summary(a_est_body_size_acari)  #even worse (tetranychus problems... :/)

## hemiptera
hemiptera <- thermal_traits_complete_order %>% 
  filter(order == "Hemiptera")
a_est_body_size_hemiptera <- lme(a_est ~ log10(body_length),
                                        random = ~ 1|id,
                                        weights = varFixed(~vi),
                                        data = hemiptera,
                                        na.action = na.omit,
                                        control = lmeControl(sigma = 1,
                                                             msMaxIter = 200))

summary(a_est_body_size_hemiptera)  #no significant

