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
thermal_traits_trans <- read_csv("thermal_traits_trans.csv")
thermal_traits_trans_order <- read_csv("thermal_traits_trans_order.csv")
thermal_traits_trans_fg <- read_csv("thermal_traits_trans_fg.csv")
thermal_traits_complete <- read_csv("thermal_traits_complete.csv")
species_list <- thermal_traits_complete %>% 
  group_by(id) %>% 
  summarise(spp = unique(paste(genus, species))) %>% 
  write_csv("species_list.csv")
# Meta-analysis models ----

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

