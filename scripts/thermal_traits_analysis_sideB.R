# SCRIPT INFO --------------------
#     Authors: Dario San Segundo Molina, Sara Villen Perez, Ignacio Morales Castilla
#     Title: trait inferences meta-analyses
#     Aim: perform meta-analysis models to parameterised dataset of thermal traits
#     Date: March 2022
#
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
library(tidyverse)
library(emmeans)
library(nlme)
library(here)
setwd(paste0(here(),"/data"))
thermal_traits_trans <- read_csv("thermal_traits_trans.csv")
thermal_traits_trans_order <- read_csv("thermal_traits_trans_order.csv")
thermal_traits_trans_fg <- read_csv("thermal_traits_trans_fg.csv")
thermal_traits_complete <- read_csv("thermal_traits_complete.csv")
species_list <- thermal_traits_complete %>% 
  group_by(id) %>% 
  summarise(spp = unique(paste(genus, species))) %>% 
  write_csv("species_list.csv")
# Meta-analysis models ----
thermal_traits_trans <- thermal_traits_trans %>% 
  mutate(body_length = log10(body_length))
thermal_traits_trans_order <- thermal_traits_trans_order %>% 
  mutate(body_length = log10(body_length))
#tampoco es normal :/
#  a) tmax  ---- 
# ~~~~ i) random intercept  ---- 
thermal_traits_trans
tmax_intercept <- lme(tmax ~ 1,
                      random = ~1|id,
                      weights = varFixed(~vi),
                      data = thermal_traits_trans,
                      control = lmeControl(sigma = 1))
summary(tmax_intercept)
intervals_tmax_int <- bernr::bolker_ci(tmax_intercept,newdat = data.frame(x=1),pred_int = TRUE)
tmax_intercept_output <- intervals_tmax_int %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmax, estimate = .x)))

# ~~~~ ii) ~ lat  ---- 
## random slope & intercept
tmax_lat_slope <- lme(tmax ~ abs(lat),
                      random = ~abs(lat)|id,
                      weights = varFixed(~vi),
                      data = thermal_traits_trans,
                      control = lmeControl(sigma = 1))
summary(tmax_lat_slope)
intervals_tmax_lat_slope <- intervals(tmax_lat_slope, which = "fixed")
intervals_tmax_lat_slope <- as.data.frame(intervals_tmax_lat_slope$fixed) #probably give this as "slope per standard unit"
tmax_lat_slope_output <- intervals_tmax_lat_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmax, estimate = .x)))
# ~~~~ iii) ~ order  ---- 
## random slope & intercept
tmax_order <- lme(tmax ~ as_factor(order),
                  random = ~ as_factor(order)|id,
                  weights = varFixed(~vi),
                  data = thermal_traits_trans_order,
                  control = lmeControl(sigma = 1))
summary(tmax_order)
coefs_tmax_order <- as_tibble(coefficients(summary(tmax_order))) %>% 
  mutate(estimate = 33.419891+Value)
emmeans(tmax_order, list(pairwise ~ order), adjust = "tukey")

intervals_tmax_order <- bernr::bolker_ci(tmax_order,
                                         newdat = data.frame(order = c("Acari", "Hemiptera", "Lepidoptera")),
                                         pred_int = TRUE)
tmax_order_output <- intervals_tmax_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmax, estimate = .x))) %>% 
  mutate(order = c("Acari", "Hemiptera", "Lepidoptera")) %>% 
  print()

# ~~~~ iv) ~ feeding guild  ---- 
## random slope & intercept
tmax_feeding_guild <- lme(tmax ~ as_factor(feeding_guild),
                          random = ~ as_factor(feeding_guild)|id,
                          weights = varFixed(~vi),
                          data = thermal_traits_trans_fg,
                          control = lmeControl(sigma = 1))
summary(tmax_feeding_guild)
coefs_tmax_fg <- as_tibble(coefficients(summary(tmax_feeding_guild))) %>% 
  mutate(estimate = 33.419891+Value)
emmeans(tmax_feeding_guild, list(pairwise ~ feeding_guild), adjust = "tukey")

intervals_tmax_feeding_guild <- bernr::bolker_ci(tmax_feeding_guild,
                                                 newdat = data.frame(feeding_guild = c("borer", "chewer", "sucker")),
                                                 pred_int = TRUE)

tmax_feeding_guild_output <- intervals_tmax_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmax, estimate = .x))) %>% 
  mutate(feeding_guild = c("borer", "chewer", "sucker")) %>% 
  print()


# ~~~~ v) ~ year  ---- 
## random slope & intercept
tmax_year <- lme(tmax ~ Year,
                 random = ~ Year|id,
                 weights = varFixed(~vi),
                 data = thermal_traits_trans,
                 control = lmeControl(sigma = 1, msMaxIter = 100))
summary(tmax_year)
intervals_tmax_year <- intervals(tmax_year, which = "fixed")
intervals_tmax_year <- as.data.frame(intervals_tmax_year$fixed) #probably give this as "slope per standard unit"

tmax_year_output <- intervals_tmax_year %>% 
  select(lower, est., upper) %>% 
  slice(2) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmax, estimate = .x))) %>% 
  print()

# ~~~~ vi) ~ body size  ---- 
## random slope & intercept
tmax_bodysize_slope <- lme(tmax ~ body_length,
                           random = ~body_length|id,
                           weights = varFixed(~vi),
                          data = thermal_traits_trans,
                          control = lmeControl(sigma = 1,
                                               msMaxIter = 200))
summary(tmax_bodysize_slope)
intervals_tmax_bodysize_slope <- intervals(tmax_bodysize_slope, which = "fixed")
intervals_tmax_bodysize_slope <- as.data.frame(intervals_tmax_bodysize_slope$fixed) #probably give this as "slope per standard unit"
tmax_bodysize_slope_output <- intervals_tmax_bodysize_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmax, estimate = .x)))

#and without Lepidoptera? (might be exception to Bergmann's rule according to Huey & Kingsolver, 2008)
tmax_year <- lme(tmax ~ Year,
                 random = ~ Year|id,
                 weights = varFixed(~vi),
                 data = thermal_traits_trans,
                 control = lmeControl(sigma = 1, msMaxIter = 100))
summary(tmax_year)
intervals_tmax_year <- intervals(tmax_year, which = "fixed")
intervals_tmax_year <- as.data.frame(intervals_tmax_year$fixed) #probably give this as "slope per standard unit"

tmax_year_output <- intervals_tmax_year %>% 
  select(lower, est., upper) %>% 
  slice(2) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmax, estimate = .x))) %>% 
  print()

# ~~~~ vi) ~ body size + order  ---- 
## random intercept
tmax_bodysize_order <- lme(tmax ~ body_length + order + body_length*order,
                           random = ~1|id,
                           weights = varFixed(~vi),
                           data = thermal_traits_trans_order,
                           control = lmeControl(sigma = 1))
summary(tmax_bodysize_order)
# ~~~~~~~~~~ only lepidoptera  ---- 
lepi_trans <- thermal_traits_trans_order %>% 
  filter(order == "Lepidoptera")
tmax_bodysize_lepi <- lme(tmax ~ body_length,
                          random = ~body_length|id,
                          weights = varFixed(~vi),
                          data = lepi_trans,
                          control = lmeControl(sigma = 1))
summary(tmax_bodysize_lepi) #not significant
# ~~~~~~~~~~ only acari  ---- 
acari_trans <- thermal_traits_trans_order %>% 
  filter(order == "Acari")
tmax_bodysize_acari <- lme(tmax ~ body_length,
                          random = ~body_length|id,
                          weights = varFixed(~vi),
                          data = acari_trans,
                          control = lmeControl(sigma = 1,
                                               msMaxIter = 100))
summary(tmax_bodysize_acari) #not significant
# ~~~~~~~~~~ only hemiptera  ---- 
hemiptera_trans <- thermal_traits_trans_order %>% 
  filter(order == "Hemiptera")
tmax_bodysize_hemiptera<- lme(tmax ~ body_length,
                           random = ~body_length|id,
                           weights = varFixed(~vi),
                           data = hemiptera_trans,
                           control = lmeControl(sigma = 1,
                                                msMaxIter = 100))
summary(tmax_bodysize_hemiptera) #not significant



#  b) tmin  ---- 
# ~~~~ i) random intercept  ---- 

tmin_intercept <- lme(tmin ~ 1,
                      random = ~1|id,
                      weights = varFixed(~vi),
                      data = thermal_traits_trans,
                      control = lmeControl(sigma = 1))
summary(tmin_intercept)
intervals_tmin_int <- bernr::bolker_ci(tmin_intercept,newdat = data.frame(x=1),pred_int = TRUE)
tmin_intercept_output <- intervals_tmin_int %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmin, estimate = .x))) %>% 
  print()

# # ~~~~ ii) ~ lat  ---- 
## random slope & intercept
tmin_lat_slope <- lme(tmin ~ abs(lat),
                      random = ~abs(lat)|id,
                      weights = varFixed(~vi),
                      data = thermal_traits_trans,
                      control = lmeControl(sigma = 1))
summary(tmin_lat_slope)
intervals_tmin_lat_slope <- intervals(tmin_lat_slope, which = "fixed")
intervals_tmin_lat_slope <- as.data.frame(intervals_tmin_lat_slope$fixed) #probably give this as "slope per standard unit"
tmin_lat_slope_output <- intervals_tmin_lat_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmin, estimate = .x))) # and this for Tmin at lat=0?
# # ~~~~ iii) ~ order  ---- 
## random slope & intercept
tmin_order <- lme(tmin ~ as_factor(order),
                  random = ~ as_factor(order)|id,
                  weights = varFixed(~vi),
                  data = thermal_traits_trans_order,
                  control = lmeControl(sigma = 1))
summary(tmin_order)
emmeans(tmin_order, list(pairwise ~ order), adjust = "tukey")
intervals_tmin_order <- bernr::bolker_ci(tmin_order,
                                         newdat = data.frame(order = c("Acari", "Hemiptera", "Lepidoptera")),
                                         pred_int = TRUE)
tmin_order_output <- intervals_tmin_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmin, estimate = .x))) %>% 
  mutate(order = c("Acari", "Hemiptera", "Lepidoptera")) %>% 
  print()


# # ~~~~ iv) ~ feeding guild  ---- 
## random slope & intercept
tmin_feeding_guild <- lme(tmin ~ as_factor(feeding_guild),
                          random = ~ as_factor(feeding_guild)|id,
                          weights = varFixed(~vi),
                          data = thermal_traits_trans_fg,
                          control = lmeControl(sigma = 1))
summary(tmin_feeding_guild)
emmeans(tmin_feeding_guild, list(pairwise ~ feeding_guild), adjust = "tukey")
intervals_tmin_feeding_guild <- bernr::bolker_ci(tmin_feeding_guild,
                                         newdat = data.frame(feeding_guild = c("borer", "chewer", "sucker")),
                                         pred_int = TRUE)
tmin_feeding_guild_output <- intervals_tmin_feeding_guild %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmin, estimate = .x))) %>% 
  mutate(feeding_guild = c("borer", "chewer", "sucker")) %>% 
  print()


  # # ~~~~ v) ~ year  ---- 
## random slope & intercept
tmin_year <- lme(tmin ~ Year,
                 random = ~ Year|id,
                 weights = varFixed(~vi),
                 data = thermal_traits_trans,
                 control = lmeControl(sigma = 1,
                                      msMaxIter = 100))
summary(tmin_year)
intervals_tmin_year <- intervals(tmin_year, which = "fixed")
intervals_tmin_year <- as.data.frame(intervals_tmin_year$fixed) #probably give this as "slope per standard unit"

# ~~~~ vi) ~ body size  ---- 
## random slope & intercept
tmin_bodysize_slope <- lme(tmin ~ body_length,
                           random = ~body_length|id,
                           weights = varFixed(~vi),
                           data = thermal_traits_trans,
                           control = lmeControl(sigma = 1,
                                                maxIter = 200))
summary(tmin_bodysize_slope)
intervals_tmin_bodysize_slope <- intervals(tmin_bodysize_slope, which = "fixed")
intervals_tmin_bodysize_slope <- as.data.frame(intervals_tmin_bodysize_slope$fixed) #probably give this as "slope per standard unit"
tmin_bodysize_slope_output <- intervals_tmin_bodysize_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmin, estimate = .x)))

#and without Lepidoptera? (might be exception to Bergmann's rule according to Huey & Kingsolver, 2008)
# ~~~~ vi) ~ body size + order  ---- 
## random intercept
tmin_bodysize_order <- lme(tmin ~ body_length + order + body_length*order,
                           random = ~body_length*order|id,
                           weights = varFixed(~vi),
                           data = thermal_traits_trans_order,
                           control = lmeControl(sigma = 1,
                                                msMaxIter = 100))
summary(tmin_bodysize_order)
# ~~~~~~~~~~ only lepidoptera  ---- 
lepi_trans <- thermal_traits_trans_order %>% 
  filter(order == "Lepidoptera")
tmin_bodysize_lepi <- lme(tmin ~ body_length,
                          random = ~body_length|id,
                          weights = varFixed(~vi),
                          data = lepi_trans,
                          control = lmeControl(sigma = 1))
summary(tmin_bodysize_lepi) #not significant
# ~~~~~~~~~~ only acari  ---- 
acari_trans <- thermal_traits_trans_order %>% 
  filter(order == "Acari")
tmin_bodysize_acari <- lme(tmin ~ body_length,
                           random = ~body_length|id,
                           weights = varFixed(~vi),
                           data = acari_trans,
                           control = lmeControl(sigma = 1,
                                                msMaxIter = 200))
summary(tmin_bodysize_acari) #not significant
# ~~~~~~~~~~ only hemiptera  ---- 
hemiptera_trans <- thermal_traits_trans_order %>% 
  filter(order == "Hemiptera")
tmin_bodysize_hemiptera<- lme(tmin ~ body_length,
                              random = ~body_length|id,
                              weights = varFixed(~vi),
                              data = hemiptera_trans,
                              control = lmeControl(sigma = 1,
                                                   msMaxIter = 100))
summary(tmin_bodysize_hemiptera) #not significant




# c) topt  ---- 
# ~~~~ i) random intercept  ---- 

topt_intercept <- lme(topt ~ 1,
                      random = ~1|id,
                      weights = varFixed(~vi),
                      data = thermal_traits_trans,
                      control = lmeControl(sigma = 1))
summary(topt_intercept)
intervals_topt_int <- bernr::bolker_ci(topt_intercept,newdat = data.frame(x=1),pred_int = TRUE)
topt_intercept_output <- intervals_topt_int %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_topt, estimate = .x))) %>% 
  print()

# ~~~~ ii) ~ lat  ---- 
## random slope & intercept
topt_lat_slope <- lme(topt ~ abs(lat),
                      random = ~abs(lat)|id,
                      weights = varFixed(~vi),
                      data = thermal_traits_trans,
                      control = lmeControl(sigma = 1))
summary(topt_lat_slope)
intervals_topt_lat_slope <- intervals(topt_lat_slope, which = "fixed")
intervals_topt_lat_slope <- as.data.frame(intervals_topt_lat_slope$fixed) #probably give this as "slope per standard unit"
topt_lat_slope_output <- intervals_topt_lat_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_topt, estimate = .x))) # and this for topt at lat=0?
# ~~~~ iii) ~ order  ---- 
## random slope & intercept
topt_order <- lme(topt ~ as_factor(order),
                  random = ~ as_factor(order)|id,
                  weights = varFixed(~vi),
                  data = thermal_traits_trans_order,
                  control = lmeControl(sigma = 1))
summary(topt_order)
emmeans(topt_order, list(pairwise ~ order), adjust = "tukey")
intervals_topt_order <- bernr::bolker_ci(topt_order,
                                         newdat = data.frame(order = c("Acari", "Hemiptera", "Lepidoptera")),
                                         pred_int = TRUE)
topt_order_output <- intervals_topt_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_topt, estimate = .x))) %>% 
  mutate(order = c("Acari", "Hemiptera", "Lepidoptera")) %>% 
  print()


# ~~~~ iv) ~ feeding guild  ---- 
## random slope & intercept
topt_feeding_guild <- lme(topt ~ as_factor(feeding_guild),
                          random = ~ as_factor(feeding_guild)|id,
                          weights = varFixed(~vi),
                          data = thermal_traits_trans_fg,
                          control = lmeControl(sigma = 1))
summary(topt_feeding_guild)
emmeans(topt_feeding_guild, list(pairwise ~ feeding_guild), adjust = "tukey")
intervals_topt_feeding_guild <- bernr::bolker_ci(topt_feeding_guild,
                                                 newdat = data.frame(feeding_guild = c("borer", "chewer", "sucker")),
                                                 pred_int = TRUE)
topt_feeding_guild_output <- intervals_topt_feeding_guild %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_topt, estimate = .x))) %>% 
  mutate(feeding_guild = c("borer", "chewer", "sucker")) %>% 
  print()


# ~~~~ v) ~ year  ---- 
## random slope & intercept
topt_year <- lme(topt ~ Year,
                 random = ~ Year|id,
                 weights = varFixed(~vi),
                 data = thermal_traits_trans,
                 control = lmeControl(sigma = 1,
                                      msMaxIter = 100))
summary(topt_year)
intervals_topt_year <- intervals(topt_year, which = "fixed")
intervals_topt_year <- as.data.frame(intervals_topt_year$fixed) #probably give this as "slope per standard unit"
#  d) thermal_breadth  ---- 
# ~~~~ i) random intercept  ---- 
thermal_breadth_intercept <- lme(thermal_breadth ~ 1,
                                 weights = varFixed(~vi),
                                 random = ~1|id,
                                 data = thermal_traits_trans,
                                 control = lmeControl(sigma = 1))
summary(thermal_breadth_intercept)
intervals_thermal_breadth_int <- bernr::bolker_ci(thermal_breadth_intercept,
                                                  newdat = data.frame(x=1),
                                                  pred_int = TRUE)
thermal_breadth_intercept_output <- intervals_thermal_breadth_int %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_thermal_breadth, estimate = .x))) %>% 
  print()

# ~~~~ ii) ~ lat  ---- 
## random slope & intercept
thermal_breadth_lat_slope <- lme(thermal_breadth ~ abs(lat),
                      random = ~abs(lat)|id,
                      weights = varFixed(~vi),
                      data = thermal_traits_trans,
                      control = lmeControl(sigma = 1))
summary(thermal_breadth_lat_slope)
intervals_thermal_breadth_lat_slope <- intervals(thermal_breadth_lat_slope, which = "fixed")
intervals_thermal_breadth_lat_slope <- as.data.frame(intervals_thermal_breadth_lat_slope$fixed) #probably give this as "slope per standard unit"
thermal_breadth_lat_slope_output <- intervals_thermal_breadth_lat_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_thermal_breadth, estimate = .x))) # and this for thermal_breadth at lat=0?
# ~~~~ iii) ~ order  ---- 
## random slope & intercept
thermal_breadth_order <- lme(thermal_breadth ~ as_factor(order),
                             random = ~ as_factor(order)|id,
                             weights = varFixed(~vi),
                             data = thermal_traits_trans_order,
                             control = lmeControl(sigma = 1))
summary(thermal_breadth_order)
emmeans(thermal_breadth_order, list(pairwise ~ order), adjust = "tukey")
intervals_thermal_breadth_order <- bernr::bolker_ci(thermal_breadth_order,
                                         newdat = data.frame(order = c("Acari", "Hemiptera", "Lepidoptera")),
                                         pred_int = TRUE)
thermal_breadth_order_output <- intervals_thermal_breadth_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_thermal_breadth, estimate = .x))) %>% 
  mutate(order = c("Acari", "Hemiptera", "Lepidoptera")) %>% 
  print()


# ~~~~ iv) ~ feeding guild  ---- 
## random slope & intercept
thermal_breadth_feeding_guild <- lme(thermal_breadth ~ as_factor(feeding_guild),
                                     random = ~ as_factor(feeding_guild)|id,
                                     weights = varFixed(~vi),
                                     data = thermal_traits_trans_fg,
                                     control = lmeControl(sigma = 1))
summary(thermal_breadth_feeding_guild)
emmeans(thermal_breadth_feeding_guild, list(pairwise ~ feeding_guild), adjust = "tukey")
intervals_thermal_breadth_feeding_guild <- bernr::bolker_ci(thermal_breadth_feeding_guild,
                                                 newdat = data.frame(feeding_guild = c("borer", "chewer", "sucker")),
                                                 pred_int = TRUE)
thermal_breadth_feeding_guild_output <- intervals_thermal_breadth_feeding_guild %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_thermal_breadth, estimate = .x))) %>% 
  mutate(feeding_guild = c("borer", "chewer", "sucker")) %>% 
  print()


# ~~~~ v) ~ year  ---- 
## random slope & intercept
thermal_breadth_year <- lme(thermal_breadth ~ Year,
                            random = ~ Year|id,
                            weights = varFixed(~vi),
                            data = thermal_traits_trans,
                            control = lmeControl(sigma = 1,
                                                 msMaxIter = 100))
summary(thermal_breadth_year)
intervals_thermal_breadth_year <- intervals(thermal_breadth_year, which = "fixed")
intervals_thermal_breadth_year <- as.data.frame(intervals_thermal_breadth_year$fixed) #probably give this as "slope per standard unit"

# ~~~~ vi) ~ body size  ---- 
## random slope & intercept
thermal_breadth_bodysize_slope <- lme(thermal_breadth ~ body_length,
                                      random = ~body_length|id,
                                      weights = varFixed(~vi),
                                      data = thermal_traits_trans,
                                      control = lmeControl(sigma = 1))
summary(thermal_breadth_bodysize_slope)
intervals_thermal_breadth_bodysize_slope <- intervals(thermal_breadth_bodysize_slope, which = "fixed")
intervals_thermal_breadth_bodysize_slope <- as.data.frame(intervals_thermal_breadth_bodysize_slope$fixed) #probably give this as "slope per standard unit"
thermal_breadth_bodysize_slope_output <- intervals_thermal_breadth_bodysize_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_thermal_breadth, estimate = .x)))

#and without Lepidoptera? (might be exception to Bergmann's rule according to Huey & Kingsolver, 2008)
# ~~~~ vi) ~ body size + order  ---- 
## random intercept
thermal_breadth_bodysize_order <- lme(thermal_breadth ~ body_length + order + body_length*order,
                                      random = ~body_length*order|id,
                                      weights = varFixed(~vi),
                                      data = thermal_traits_trans_order,
                                      control = lmeControl(sigma = 1,
                                                           msMaxIter = 100))
summary(thermal_breadth_bodysize_order)
# ~~~~~~~~~~ only lepidoptera  ---- 
lepi_trans <- thermal_traits_trans_order %>% 
  filter(order == "Lepidoptera")
thermal_breadth_bodysize_lepi <- lme(thermal_breadth ~ body_length,
                                     random = ~body_length|id,
                                     weights = varFixed(~vi),
                                     data = lepi_trans,
                                     control = lmeControl(sigma = 1))
summary(thermal_breadth_bodysize_lepi) #not significant
# ~~~~~~~~~~ only acari  ---- 
acari_trans <- thermal_traits_trans_order %>% 
  filter(order == "Acari")
thermal_breadth_bodysize_acari <- lme(thermal_breadth ~ body_length,
                                      random = ~body_length|id,
                                      weights = varFixed(~vi),
                                      data = acari_trans,
                                      control = lmeControl(sigma = 1,
                                                           msMaxIter = 100))
summary(thermal_breadth_bodysize_acari) #not significant
# ~~~~~~~~~~ only hemiptera  ---- 
hemiptera_trans <- thermal_traits_trans_order %>% 
  filter(order == "Hemiptera")
thermal_breadth_bodysize_hemiptera<- lme(thermal_breadth ~ body_length,
                              random = ~body_length|id,
                              weights = varFixed(~vi),
                              data = hemiptera_trans,
                              control = lmeControl(sigma = 1,
                                                   msMaxIter = 100))
summary(thermal_breadth_bodysize_hemiptera) #not significant




#  e) thermal safety margin  ---- 
# ~~~~ i) random intercept  ---- 
tsm_intercept <- lme(tsm ~ 1,
                     weights = varFixed(~vi),
                     random = ~1|id,
                     data = thermal_traits_trans,
                     control = lmeControl(sigma = 1))
summary(tsm_intercept)
intervals_tsm_int <- bernr::bolker_ci(tsm_intercept,
                                                  newdat = data.frame(x=1),
                                                  pred_int = TRUE)
tsm_intercept_output <- intervals_tsm_int %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tsm, estimate = .x))) %>% 
  print()

# ~~~~ ii) ~ lat  ---- 
## random slope & intercept
tsm_lat_slope <- lme(tsm ~ abs(lat),
                                 random = ~abs(lat)|id,
                                 weights = varFixed(~vi),
                                 data = thermal_traits_trans,
                                 control = lmeControl(sigma = 1))
summary(tsm_lat_slope)
intervals_tsm_lat_slope <- intervals(tsm_lat_slope, which = "fixed")
intervals_tsm_lat_slope <- as.data.frame(intervals_tsm_lat_slope$fixed) #probably give this as "slope per standard unit"
tsm_lat_slope_output <- intervals_tsm_lat_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tsm, estimate = .x))) # and this for tsm at lat=0?
# ~~~~ iii) ~ order  ---- 
## random slope & intercept
tsm_order <- lme(tsm ~ as_factor(order),
                             random = ~ as_factor(order)|id,
                             weights = varFixed(~vi),
                             data = thermal_traits_trans_order,
                             control = lmeControl(sigma = 1))
summary(tsm_order)
emmeans(tsm_order, list(pairwise ~ order), adjust = "tukey")
intervals_tsm_order <- bernr::bolker_ci(tsm_order,
                                                    newdat = data.frame(order = c("Acari", "Hemiptera", "Lepidoptera")),
                                                    pred_int = TRUE)
tsm_order_output <- intervals_tsm_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tsm, estimate = .x))) %>% 
  mutate(order = c("Acari", "Hemiptera", "Lepidoptera")) %>% 
  print()


# ~~~~ iv) ~ feeding guild  ---- 
## random slope & intercept
tsm_feeding_guild <- lme(tsm ~ as_factor(feeding_guild),
                                     random = ~ as_factor(feeding_guild)|id,
                                     weights = varFixed(~vi),
                                     data = thermal_traits_trans_fg,
                                     control = lmeControl(sigma = 1))
summary(tsm_feeding_guild)
emmeans(tsm_feeding_guild, list(pairwise ~ feeding_guild), adjust = "tukey")
intervals_tsm_feeding_guild <- bernr::bolker_ci(tsm_feeding_guild,
                                                            newdat = data.frame(feeding_guild = c("borer", "chewer", "sucker")),
                                                            pred_int = TRUE)
tsm_feeding_guild_output <- intervals_tsm_feeding_guild %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_tsm, estimate = .x))) %>% 
  mutate(feeding_guild = c("borer", "chewer", "sucker")) %>% 
  print()


# ~~~~ v) ~ year  ---- 
## random slope & intercept
tsm_year <- lme(tsm ~ Year,
                random = ~ Year|id,
                weights = varFixed(~vi),
                data = thermal_traits_trans,
                control = lmeControl(sigma = 1,
                                                 msMaxIter = 100))
summary(tsm_year)
intervals_tsm_year <- intervals(tsm_year, which = "fixed")
intervals_tsm_year <- as.data.frame(intervals_tsm_year$fixed) #probably give this as "slope per standard unit"

#  f) thermal range   ---- 
# ~~~~ i) random intercept  ---- 
therm_range_intercept <- lme(therm_range ~ 1,
                             weights = varFixed(~vi),
                             random = ~1|id,
                             data = thermal_traits_trans,
                             control = lmeControl(sigma = 1))
summary(therm_range_intercept)
intervals_therm_range_int <- bernr::bolker_ci(therm_range_intercept,
                                      newdat = data.frame(x=1),
                                      pred_int = TRUE)
therm_range_intercept_output <- intervals_therm_range_int %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_therm_range, estimate = .x))) %>% 
  print()

# ~~~~ ii) ~ lat  ---- 
## random slope & intercept
therm_range_lat_slope <- lme(therm_range ~ abs(lat),
                             random = ~abs(lat)|id,
                             weights = varFixed(~vi),
                             data = thermal_traits_trans,
                             control = lmeControl(sigma = 1))
summary(therm_range_lat_slope)
intervals_therm_range_lat_slope <- intervals(therm_range_lat_slope, which = "fixed")
intervals_therm_range_lat_slope <- as.data.frame(intervals_therm_range_lat_slope$fixed) #probably give this as "slope per standard unit"
therm_range_lat_slope_output <- intervals_therm_range_lat_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_therm_range, estimate = .x))) # and this for therm_range at lat=0?
# ~~~~ iii) ~ order  ---- 
## random slope & intercept
therm_range_order <- lme(therm_range ~ as_factor(order),
                         random = ~ as_factor(order)|id,
                         weights = varFixed(~vi),
                         data = thermal_traits_trans_order,
                         control = lmeControl(sigma = 1))
summary(therm_range_order)
emmeans(therm_range_order, list(pairwise ~ order), adjust = "tukey")
intervals_therm_range_order <- bernr::bolker_ci(therm_range_order,
                                                newdat = data.frame(order = c("Acari", "Hemiptera", "Lepidoptera")),
                                                pred_int = TRUE)
therm_range_order_output <- intervals_therm_range_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_therm_range, estimate = .x))) %>% 
  mutate(order = c("Acari", "Hemiptera", "Lepidoptera")) %>% 
  print()


# ~~~~ iv) ~ feeding guild  ---- 
## random slope & intercept
therm_range_feeding_guild <- lme(therm_range ~ as_factor(feeding_guild),
                                 random = ~ as_factor(feeding_guild)|id,
                                 weights = varFixed(~vi),
                                 data = thermal_traits_trans_fg,
                                 control = lmeControl(sigma = 1))
summary(therm_range_feeding_guild)
emmeans(therm_range_feeding_guild, list(pairwise ~ feeding_guild), adjust = "tukey")
intervals_therm_range_feeding_guild <- bernr::bolker_ci(therm_range_feeding_guild,
                                                newdat = data.frame(feeding_guild = c("borer", "chewer", "sucker")),
                                                pred_int = TRUE)
therm_range_feeding_guild_output <- intervals_therm_range_feeding_guild %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_therm_range, estimate = .x))) %>% 
  mutate(feeding_guild = c("borer", "chewer", "sucker")) %>% 
  print()


# ~~~~ v) ~ year  ---- 
## random slope & intercept
therm_range_year <- lme(therm_range ~ Year,
                        random = ~ Year|id,
                        weights = varFixed(~vi),
                        data = thermal_traits_trans,
                        control = lmeControl(sigma = 1,
                                     msMaxIter = 100))
summary(therm_range_year)
intervals_therm_range_year <- intervals(therm_range_year, which = "fixed")
intervals_therm_range_year <- as.data.frame(intervals_therm_range_year$fixed) #probably give this as "slope per standard unit"

#  g) thermal window  ---- 
# ~~~~ i) random intercept  ---- 
therm_window_intercept <- lme(therm_window ~ 1,
                                 weights = varFixed(~vi),
                                 random = ~1|id,
                                 data = thermal_traits_trans,
                                 control = lmeControl(sigma = 1))
summary(therm_window_intercept)
intervals_therm_window_int <- bernr::bolker_ci(therm_window_intercept,
                                                  newdat = data.frame(x=1),
                                                  pred_int = TRUE)
therm_window_intercept_output <- intervals_therm_window_int %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_therm_window, estimate = .x))) %>% 
  print()

# ~~~~ ii) ~ lat  ---- 
## random slope & intercept
therm_window_lat_slope <- lme(therm_window ~ abs(lat),
                                 random = ~abs(lat)|id,
                                 weights = varFixed(~vi),
                                 data = thermal_traits_trans,
                                 control = lmeControl(sigma = 1))
summary(therm_window_lat_slope)
intervals_therm_window_lat_slope <- intervals(therm_window_lat_slope, which = "fixed")
intervals_therm_window_lat_slope <- as.data.frame(intervals_therm_window_lat_slope$fixed) #probably give this as "slope per standard unit"
therm_window_lat_slope_output <- intervals_therm_window_lat_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_therm_window, estimate = .x))) # and this for thermal_breadth at lat=0?
# ~~~~ iii) ~ order  ---- 
## random slope & intercept
therm_window_order <- lme(therm_window ~ as_factor(order),
                             random = ~ as_factor(order)|id,
                             weights = varFixed(~vi),
                             data = thermal_traits_trans_order,
                             control = lmeControl(sigma = 1))
summary(therm_window_order)
emmeans(therm_window_order, list(pairwise ~ order), adjust = "tukey")
intervals_therm_window_order <- bernr::bolker_ci(therm_window_order,
                                                    newdat = data.frame(order = c("Acari", "Hemiptera", "Lepidoptera")),
                                                    pred_int = TRUE)
therm_window_order_output <- intervals_therm_window_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_therm_window, estimate = .x))) %>% 
  mutate(order = c("Acari", "Hemiptera", "Lepidoptera")) %>% 
  print()


# ~~~~ iv) ~ feeding guild  ---- 
## random slope & intercept
therm_window_feeding_guild <- lme(therm_window ~ as_factor(feeding_guild),
                                     random = ~ as_factor(feeding_guild)|id,
                                     weights = varFixed(~vi),
                                     data = thermal_traits_trans_fg,
                                     control = lmeControl(sigma = 1))
summary(therm_window_feeding_guild)
emmeans(therm_window_feeding_guild, list(pairwise ~ feeding_guild), adjust = "tukey")
intervals_therm_window_feeding_guild <- bernr::bolker_ci(therm_window_feeding_guild,
                                                            newdat = data.frame(feeding_guild = c("borer", "chewer", "sucker")),
                                                            pred_int = TRUE)
therm_window_feeding_guild_output <- intervals_therm_window_feeding_guild %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_therm_window, estimate = .x))) %>% 
  mutate(feeding_guild = c("borer", "chewer", "sucker")) %>% 
  print()


# ~~~~ v) ~ year  ---- 
## random slope & intercept
therm_window_year <- lme(therm_window ~ Year,
                            random = ~ Year|id,
                            weights = varFixed(~vi),
                            data = thermal_traits_trans,
                            control = lmeControl(sigma = 1,
                                                 msMaxIter = 100))
summary(therm_window_year)
intervals_therm_window_year <- intervals(therm_window_year, which = "fixed")
intervals_therm_window_year <- as.data.frame(intervals_therm_window_year$fixed) #probably give this as "slope per standard unit"

# ~~~~ vi) ~ body size  ---- 
## random slope & intercept
therm_window_bodysize_slope <- lme(therm_window ~ body_length,
                                   random = ~body_length|id,
                                   weights = varFixed(~vi),
                                   data = thermal_traits_trans,
                                   control = lmeControl(sigma = 1))
summary(therm_window_bodysize_slope)
intervals_therm_window_bodysize_slope <- intervals(therm_window_bodysize_slope, which = "fixed")
intervals_therm_window_bodysize_slope <- as.data.frame(intervals_therm_window_bodysize_slope$fixed) #probably give this as "slope per standard unit"
therm_window_bodysize_slope_output <- intervals_therm_window_bodysize_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_therm_window, estimate = .x)))

#and without Lepidoptera? (might be exception to Bergmann's rule according to Huey & Kingsolver, 2008)
# ~~~~ vi) ~ body size + order  ---- 
## random intercept
therm_window_bodysize_order <- lme(therm_window ~ body_length + order + body_length*order,
                           random = ~body_length*order|id,
                           weights = varFixed(~vi),
                           data = thermal_traits_trans_order,
                           control = lmeControl(sigma = 1,
                                                msMaxIter = 100))
summary(therm_window_bodysize_order)
# ~~~~~~~~~~ only lepidoptera  ---- 
lepi_trans <- thermal_traits_trans_order %>% 
  filter(order == "Lepidoptera")
therm_window_bodysize_lepi <- lme(therm_window ~ body_length,
                          random = ~body_length|id,
                          weights = varFixed(~vi),
                          data = lepi_trans,
                          control = lmeControl(sigma = 1))
summary(therm_window_bodysize_lepi) #not significant
# ~~~~~~~~~~ only acari  ---- 
acari_trans <- thermal_traits_trans_order %>% 
  filter(order == "Acari")
therm_window_bodysize_acari <- lme(therm_window ~ body_length,
                           random = ~body_length|id,
                           weights = varFixed(~vi),
                           data = acari_trans,
                           control = lmeControl(sigma = 1,
                                                msMaxIter = 200))
summary(therm_window_bodysize_acari) #not significant
# ~~~~~~~~~~ only hemiptera  ---- 
hemiptera_trans <- thermal_traits_trans_order %>% 
  filter(order == "Hemiptera")
therm_window_bodysize_hemiptera<- lme(therm_window ~ body_length,
                              random = ~body_length|id,
                              weights = varFixed(~vi),
                              data = hemiptera_trans,
                              control = lmeControl(sigma = 1,
                                                   msMaxIter = 100))
summary(therm_window_bodysize_hemiptera) #not significant





#  h) r (proxy: a parameter)  ---- 
# ~~~~ i) random intercept  ---- 
a_est_intercept <- lme(a_est ~ 1,
                   weights = varFixed(~vi),
                   random = ~1|id,
                   data = thermal_traits_trans,
                   control = lmeControl(sigma = 1))
summary(a_est_intercept)


intervals_a_est_int <- bernr::bolker_ci(a_est_intercept,
                                               newdat = data.frame(x=1),
                                               pred_int = TRUE)
a_est_intercept_output <- intervals_a_est_int %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_a_est, estimate = .x))) %>% 
  print()

# ~~~~ ii) ~ lat  ---- 
## random slope & intercept
a_est_lat_slope <- lme(a_est ~ abs(lat),
                              random = ~abs(lat)|id,
                              weights = varFixed(~vi),
                              data = thermal_traits_trans,
                              control = lmeControl(sigma = 1))
summary(a_est_lat_slope)
intervals_a_est_lat_slope <- intervals(a_est_lat_slope, which = "fixed")
intervals_a_est_lat_slope <- as.data.frame(intervals_a_est_lat_slope$fixed) #probably give this as "slope per standard unit"
a_est_lat_slope_output <- intervals_a_est_lat_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_a_est, estimate = .x))) # and this for thermal_breadth at lat=0?
# ~~~~ iii) ~ order  ---- 
## random slope & intercept
a_est_order <- lme(a_est ~ as_factor(order),
                          random = ~ as_factor(order)|id,
                          weights = varFixed(~vi),
                          data = thermal_traits_trans_order,
                          control = lmeControl(sigma = 1))
summary(a_est_order)
emmeans(a_est_order, list(pairwise ~ order), adjust = "tukey")
intervals_a_est_order <- bernr::bolker_ci(a_est_order,
                                                 newdat = data.frame(order = c("Acari", "Hemiptera", "Lepidoptera")),
                                                 pred_int = TRUE)
a_est_order_output <- intervals_a_est_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_a_est, estimate = .x))) %>% 
  mutate(order = c("Acari", "Hemiptera", "Lepidoptera")) %>% 
  print()


# ~~~~ iv) ~ feeding guild  ---- 
## random slope & intercept
a_est_feeding_guild <- lme(a_est ~ as_factor(feeding_guild),
                                  random = ~ as_factor(feeding_guild)|id,
                                  weights = varFixed(~vi),
                                  data = thermal_traits_trans_fg,
                                  control = lmeControl(sigma = 1))
summary(a_est_feeding_guild)
emmeans(a_est_feeding_guild, list(pairwise ~ feeding_guild), adjust = "tukey")
intervals_a_est_feeding_guild <- bernr::bolker_ci(a_est_feeding_guild,
                                                         newdat = data.frame(feeding_guild = c("borer", "chewer", "sucker")),
                                                         pred_int = TRUE)
a_est_feeding_guild_output <- intervals_a_est_feeding_guild %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_a_est, estimate = .x))) %>% 
  mutate(feeding_guild = c("borer", "chewer", "sucker")) %>% 
  print()


# ~~~~ v) ~ year  ---- 
## random slope & intercept
a_est_year <- lme(a_est ~ Year,
                         random = ~ Year|id,
                         weights = varFixed(~vi),
                         data = thermal_traits_trans,
                         control = lmeControl(sigma = 1,
                                              msMaxIter = 100))
summary(a_est_year)
intervals_a_est_year <- intervals(a_est_year, which = "fixed")
intervals_a_est_year <- as.data.frame(intervals_a_est_year$fixed) #probably give this as "slope per standard unit"

# ~~~~ vi) ~ body size  ---- 
## random slope & intercept
a_est_bodysize_slope <- lme(a_est ~ body_length,
                                   random = ~body_length|id,
                                   weights = varFixed(~vi),
                                   data = thermal_traits_trans,
                                   control = lmeControl(sigma = 1, msMaxIter = 200))
summary(a_est_bodysize_slope)
intervals_a_est_bodysize_slope <- intervals(a_est_bodysize_slope, which = "fixed")
intervals_a_est_bodysize_slope <- as.data.frame(intervals_a_est_bodysize_slope$fixed) #probably give this as "slope per standard unit"
a_est_bodysize_slope_output <- intervals_a_est_bodysize_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_a_est, estimate = .x)))

#and without Lepidoptera? (might be exception to Bergmann's rule according to Huey & Kingsolver, 2008)
# ~~~~ vi) ~ body size + order  ---- 
## random intercept
a_est_bodysize_order <- lme(a_est ~ body_length + order + body_length*order,
                                   random = ~body_length*order|id,
                                   weights = varFixed(~vi),
                                   data = thermal_traits_trans_order,
                                   control = lmeControl(sigma = 1,
                                                        msMaxIter = 100))
summary(a_est_bodysize_order)
# ~~~~~~~~~~ only lepidoptera  ---- 
lepi_trans <- thermal_traits_trans_order %>% 
  filter(order == "Lepidoptera")
a_est_bodysize_lepi <- lme(a_est ~ body_length,
                                  random = ~body_length|id,
                                  weights = varFixed(~vi),
                                  data = lepi_trans,
                                  control = lmeControl(sigma = 1))
summary(a_est_bodysize_lepi) #not significant
# ~~~~~~~~~~ only acari  ---- 
acari_trans <- thermal_traits_trans_order %>% 
  filter(order == "Acari")
a_est_bodysize_acari <- lme(a_est ~ body_length,
                                   random = ~body_length|id,
                                   weights = varFixed(~vi),
                                   data = acari_trans,
                                   control = lmeControl(sigma = 1,
                                                        msMaxIter = 200))
summary(a_est_bodysize_acari) #not significant
# ~~~~~~~~~~ only hemiptera  ---- 
hemiptera_trans <- thermal_traits_trans_order %>% 
  filter(order == "Hemiptera")
a_est_bodysize_hemiptera<- lme(a_est ~ body_length,
                                      random = ~body_length|id,
                                      weights = varFixed(~vi),
                                      data = hemiptera_trans,
                                      control = lmeControl(sigma = 1,
                                                           msMaxIter = 100))
summary(a_est_bodysize_hemiptera) #not significant

# ~~~~ vii) ~ tmax ----
a_est_tmax<- lme(a_est ~ tmax,
                 random = ~tmax|id,
                 weights = varFixed(~vi),
                 data = thermal_traits_trans,
                 control = lmeControl(sigma = 1,
                                      msMaxIter = 100))
summary(a_est_tmax)

# maybe acari data is biaings, let's use lepidoptera
a_est_tmax_lepi<- lme(a_est ~ tmax,
                 random = ~tmax|id,
                 weights = varFixed(~vi),
                 data = lepi_trans,
                 control = lmeControl(sigma = 1,
                                      msMaxIter = 100))
summary(a_est_tmax_lepi)

#  i) t50_left  ---- 
# ~~~~ i) random intercept  ---- 
t50_left_intercept <- lme(t50_left ~ 1,
                     weights = varFixed(~vi),
                     random = ~1|id,
                     data = thermal_traits_trans,
                     control = lmeControl(sigma = 1))
summary(t50_left_intercept)
intervals_t50_left_int <- bernr::bolker_ci(t50_left_intercept,
                                      newdat = data.frame(x=1),
                                      pred_int = TRUE)
t50_left_intercept_output <- intervals_t50_left_int %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_t50_left, estimate = .x))) %>% 
  print()

# ~~~~ ii) ~ lat  ---- 
## random slope & intercept
t50_left_lat_slope <- lme(t50_left ~ abs(lat),
                     random = ~abs(lat)|id,
                     weights = varFixed(~vi),
                     data = thermal_traits_trans,
                     control = lmeControl(sigma = 1))
summary(t50_left_lat_slope)
intervals_t50_left_lat_slope <- intervals(t50_left_lat_slope, which = "fixed")
intervals_t50_left_lat_slope <- as.data.frame(intervals_t50_left_lat_slope$fixed) #probably give this as "slope per standard unit"
t50_left_lat_slope_output <- intervals_t50_left_lat_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_t50_left, estimate = .x))) # and this for t50_left at lat=0?
# ~~~~ iii) ~ order  ---- 
## random slope & intercept
t50_left_order <- lme(t50_left ~ as_factor(order),
                 random = ~ as_factor(order)|id,
                 weights = varFixed(~vi),
                 data = thermal_traits_trans_order,
                 control = lmeControl(sigma = 1))
summary(t50_left_order)
emmeans(t50_left_order, list(pairwise ~ order), adjust = "tukey")
intervals_t50_left_order <- bernr::bolker_ci(t50_left_order,
                                        newdat = data.frame(order = c("Acari", "Hemiptera", "Lepidoptera")),
                                        pred_int = TRUE)
t50_left_order_output <- intervals_t50_left_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_t50_left, estimate = .x))) %>% 
  mutate(order = c("Acari", "Hemiptera", "Lepidoptera")) %>% 
  print()


# ~~~~ iv) ~ feeding guild  ---- 
## random slope & intercept
t50_left_feeding_guild <- lme(t50_left ~ as_factor(feeding_guild),
                         random = ~ as_factor(feeding_guild)|id,
                         weights = varFixed(~vi),
                         data = thermal_traits_trans_fg,
                         control = lmeControl(sigma = 1))
summary(t50_left_feeding_guild)
emmeans(t50_left_feeding_guild, list(pairwise ~ feeding_guild), adjust = "tukey")
intervals_t50_left_feeding_guild <- bernr::bolker_ci(t50_left_feeding_guild,
                                                newdat = data.frame(feeding_guild = c("borer", "chewer", "sucker")),
                                                pred_int = TRUE)
t50_left_feeding_guild_output <- intervals_t50_left_feeding_guild %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_t50_left, estimate = .x))) %>% 
  mutate(feeding_guild = c("borer", "chewer", "sucker")) %>% 
  print()


# ~~~~ v) ~ year  ---- 
## random slope & intercept
t50_left_year <- lme(t50_left ~ Year,
                random = ~ Year|id,
                weights = varFixed(~vi),
                data = thermal_traits_trans,
                control = lmeControl(sigma = 1,
                                     msMaxIter = 100))
summary(t50_left_year)
intervals_t50_left_year <- intervals(t50_left_year, which = "fixed")
intervals_t50_left_year <- as.data.frame(intervals_t50_left_year$fixed) #probably give this as "slope per standard unit"

#  j) t50_right ---- 
# ~~~~ i) random intercept  ---- 
t50_right_intercept <- lme(t50_right ~ 1,
                          weights = varFixed(~vi),
                          random = ~1|id,
                          data = thermal_traits_trans,
                          control = lmeControl(sigma = 1))
summary(t50_right_intercept)
intervals_t50_right_int <- bernr::bolker_ci(t50_right_intercept,
                                           newdat = data.frame(x=1),
                                           pred_int = TRUE)
t50_right_intercept_output <- intervals_t50_right_int %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_t50_right, estimate = .x))) %>% 
  print()

# ~~~~ ii) ~ lat  ---- 
## random slope & intercept
t50_right_lat_slope <- lme(t50_right ~ abs(lat),
                          random = ~abs(lat)|id,
                          weights = varFixed(~vi),
                          data = thermal_traits_trans,
                          control = lmeControl(sigma = 1))
summary(t50_right_lat_slope)
intervals_t50_right_lat_slope <- intervals(t50_right_lat_slope, which = "fixed")
intervals_t50_right_lat_slope <- as.data.frame(intervals_t50_right_lat_slope$fixed) #probably give this as "slope per standard unit"
t50_right_lat_slope_output <- intervals_t50_right_lat_slope %>% 
  select(lower, est., upper) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_t50_right, estimate = .x))) # and this for t50_right at lat=0?
# ~~~~ iii) ~ order  ---- 
## random slope & intercept
t50_right_order <- lme(t50_right ~ as_factor(order),
                      random = ~ as_factor(order)|id,
                      weights = varFixed(~vi),
                      data = thermal_traits_trans_order,
                      control = lmeControl(sigma = 1))
summary(t50_right_order)
emmeans(t50_right_order, list(pairwise ~ order), adjust = "tukey")
intervals_t50_right_order <- bernr::bolker_ci(t50_right_order,
                                             newdat = data.frame(order = c("Acari", "Hemiptera", "Lepidoptera")),
                                             pred_int = TRUE)
t50_right_order_output <- intervals_t50_right_order %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_t50_right, estimate = .x))) %>% 
  mutate(order = c("Acari", "Hemiptera", "Lepidoptera")) %>% 
  print()


# ~~~~ iv) ~ feeding guild  ---- 
## random slope & intercept
t50_right_feeding_guild <- lme(t50_right ~ as_factor(feeding_guild),
                              random = ~ as_factor(feeding_guild)|id,
                              weights = varFixed(~vi),
                              data = thermal_traits_trans_fg,
                              control = lmeControl(sigma = 1))
summary(t50_right_feeding_guild)
emmeans(t50_right_feeding_guild, list(pairwise ~ feeding_guild), adjust = "tukey")
intervals_t50_right_feeding_guild <- bernr::bolker_ci(t50_right_feeding_guild,
                                                     newdat = data.frame(feeding_guild = c("borer", "chewer", "sucker")),
                                                     pred_int = TRUE)
t50_right_feeding_guild_output <- intervals_t50_right_feeding_guild %>% 
  select(pred, ci_l, ci_h, predint_l, predint_h) %>% 
  mutate(across(everything(), ~ backtransf(trans_var = bn_t50_right, estimate = .x))) %>% 
  mutate(feeding_guild = c("borer", "chewer", "sucker")) %>% 
  print()


# ~~~~ v) ~ year  ---- 
## random slope & intercept
t50_right_year <- lme(t50_right ~ Year,
                     random = ~ Year|id,
                     weights = varFixed(~vi),
                     data = thermal_traits_trans,
                     control = lmeControl(sigma = 1,
                                          msMaxIter = 100))
summary(t50_right_year)
intervals_t50_right_year <- intervals(t50_right_year, which = "fixed")
intervals_t50_right_year <- as.data.frame(intervals_t50_right_year$fixed) #probably give this as "slope per standard unit"
