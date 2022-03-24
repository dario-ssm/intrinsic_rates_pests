# SCRIPT INFO --------------------
#     Authors: Dario San Segundo Molina, Sara Villen Perez, Ignacio Morales Castilla
#     Title: trait inferences meta-analyses
#     Aim: explore trait trends after model fitting from INTRAPEST database
#     Date: March 2022
#
# 
#_________________________ ####



# 1. Load Dataset & explore ranges ----
rm(list=ls())
#library(tidyverse)
library(tidyverse)
library(ggmap)
library(performance)
library(cowplot)
library(wesanderson)
library(nlme)
library(brms)
library(bestNormalize)
library(metafor)
library(emmeans)
#load paremeterised thermal traits dataset (coming from simulations_and_individual_gnlsfitting.R)
setwd(paste0(getwd(),"/data"))
thermal_traits_data_raw <- read_csv("parameters_individual_fitted.csv") %>% # or parameters_individual_fitted.csv (no dimitropoulou's way)
  rename(id=id, a_est = a_est, a_se = a_se, tmin = Tmin_est, tmin_se = Tmin_se,
         tmax = Tmax_est, tmax_se = Tmax_se, topt = Value, topt_se = Topt_se,
         start_a = starting_a, start_tmin = starting_Tmin, start_tmax = starting_Tmax )%>% 
  glimpse()

# ensemble dataset complete
ir_data_all <- read_csv("IR_data_all_clean.csv")
authors <- ir_data_all %>% 
  dplyr::select(Authors) %>% 
  separate(Authors,into = c("a","b"),sep = "," ) %>% 
  dplyr::select(1) %>% 
  rename(authors = a) %>% 
  print()

ir_data_complementary <- ir_data_all %>% 
  bind_cols(authors) %>% 
  dplyr::select(id, authors,Year, title, DOI, vi, order, family, genus, species, feeding_guild,
         lat, sd_median) %>% 
  group_by(id) %>% 
  mutate(vi = median(vi)) %>%
  summarise_all(unique) %>% 
  print()
counts <- thermal_traits_data_raw %>% 
  group_by(id) %>% 
  summarise(n = n_distinct(topt_se)) %>% 
  select(n) %>% 
  as_vector()
available_ids <- thermal_traits_data_raw %>% 
  group_by(id) %>% 
  summarise(n = n_distinct(topt_se)) %>% 
  select(id) %>% 
  as_vector()

ir_data_forloop <- tibble(id = NULL,
                          authors = NULL,
                          Year = NULL,
                          title = NULL,
                          DOI = NULL,
                          vi = NULL,
                          order = NULL,
                          family = NULL,
                          genus = NULL,
                          species = NULL,
                          feeding_guild = NULL,
                          lat = NULL,
                          sd_median = NULL)
for (i in 1:length(available_ids)){
  i_th_study <- ir_data_complementary %>% filter(id == available_ids[i])
  for (nrep in 1:counts[i]){
    ir_data_forloop <- ir_data_forloop %>% 
      bind_rows((i_th_study))
  }
}
ir_data_forloop
#let's check it out:
counts_forloop <- ir_data_forloop %>% 
  count(id) %>% 
  mutate(counts_checking = counts)
view(counts_forloop) #perfect!!!

# so we assemble the dataset
ir_data_complete <- thermal_traits_data_raw %>%
  select(-id) %>% 
  bind_cols(ir_data_forloop) %>% 
  filter(tmin != start_tmin &   # avoid convergence forced parameters
         tmax != start_tmax &   # avoid convergence forced parameters
         a_est != start_a) %>%  # avoid convergence forced parameters
  glimpse()
write_csv(ir_data_complete, "parameters_individual_fit_complemented.csv")
# ~~~~  Define fitting function for modelling -----------------------------------

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

# 2. Exploratory data analysis & manipulation ----
## apply filters 
ir_data_complete <- read_csv("parameters_individual_fit_complemented.csv")
range(ir_data_complete$tmin)
q_tmin <- quantile(ir_data_complete$tmin, probs = 0.05)
q_tmin_se <-quantile(ir_data_complete$tmin_se, probs = .95)
range(ir_data_complete$tmax)
q_tmax <- quantile(ir_data_complete$tmax,probs = .95)
q_tmax_se <- quantile(ir_data_complete$tmax_se,probs = .95)
range(ir_data_complete$topt)
q_topt <- quantile(ir_data_complete$topt,probs = .95)
q_topt_se <- quantile(ir_data_complete$topt_se,probs = .95)

ir_dataset_clean_estimates <- ir_data_complete %>% 
  filter(tmax < q_tmax &
           tmin > q_tmin &
           topt < q_topt) %>%
  glimpse()

ir_dataset_clean_se <- ir_data_complete %>% 
  filter(tmax_se < q_tmax_se &
           tmin_se < q_tmin_se &
           topt_se < q_topt_se) %>%
  glimpse()

ir_dataset_clean_se_order <- ir_dataset_clean_se %>% 
  filter(order == "Acari" |
           order == "Hemiptera" |
           order == "Lepidoptera") %>%
  glimpse()


ir_dataset_clean_se_fg <- ir_dataset_clean_se %>% 
  filter(feeding_guild == "borer" |
           feeding_guild == "chewer" |
           feeding_guild == "sucker") %>%
  glimpse()

#create map along coordinates (prepared to insert in a colored background with transparent ocean)
# map <- ggplot(data = thermal_traits_data, aes(x = lon, y = lat)) +
#   borders("world", colour = "transparent", fill = "white") +
#   geom_point(aes(colour = order), alpha = 0.35, size = 3)+theme_classic()+
#   theme_void()+
#   theme(legend.position = "bottom",panel.background = element_rect(fill="transparent",color=NA))
#   
# map
#ggsave("map_meta.png",height=15,width=25,units="cm",bg = "transparent")

# ~~~~ incorporate thermal breadth ----
#first create empty tibbles to overwrite later with a mutate:
thermal_breadth <- tibble(thermal_breadth = NULL)
thermal_tol_range <- tibble(thermal_tol_left = NULL,
                            thermal_tol_right = NULL)
# and do the loop for each row in the dataset
for(i in 1:length(ir_dataset_clean_se$tmin)){
  fit_10000 <- briere1(a = ir_dataset_clean_se$a_est[i], #simulate the curve with the parameters of the model
                         Tmin = ir_dataset_clean_se$tmin[i],
                         temp = seq(0,70, 0.0001), #more 0's, more continuous-like 
                         Tmax = ir_dataset_clean_se$tmax[i])
  tbl_fit <- tibble (temp = seq(0,70, 0.0001),fit = fit_10000)
  max_r_50 <-  (max(fit_10000, na.rm = TRUE)/2) #calculate maximum _r (but two values...)
  ##we need to partition the curve into 2 sides to obtain the two T_50 values
  half_left <- tbl_fit %>% filter(temp < ir_dataset_clean_se$topt[i])
  half_right <-tbl_fit %>% filter(temp >= ir_dataset_clean_se$topt[i])
  # and see temp_r50
  half_left_t50 <- half_left %>% 
    slice(max(which(half_left$fit <= max_r_50))) %>% 
    dplyr::select(temp) %>% 
    as_vector()
  half_right_t50 <- half_right %>% 
    slice(max(which(half_right$fit <= max_r_50))) %>% 
    dplyr::select(temp) %>% 
    as_vector()
  thermal_tol_range_i <- c(half_left_t50, half_right_t50) 
  thermal_breadth_i <- half_right_t50 - half_left_t50
  # and write it!
  thermal_breadth <- thermal_breadth %>% 
    bind_rows(thermal_breadth_i)
  thermal_tol_range <- thermal_tol_range %>% 
    bind_rows(thermal_tol_range_i)
  print(paste0(round(100*i/length(ir_dataset_clean_se$tmin), digits = 2), " %"))
}
thermal_breadth <- thermal_breadth %>% 
  rename(thermal_breadth = temp)
thermal_breadth_extremes <- thermal_tol_range %>% 
  rename(t50_left = temp...1, t50_right = temp...2)
save(thermal_breadth, file = "thermal_breadth.RData")
save(thermal_breadth_extremes, file = "thermal_breadth_extremes.RData")

## and obtain final dataset
#setwd("~/.../data")
#load(thermal_breadth.RData)
#load(thermal_breadth_extremes.RData)

thermal_traits_complete <- ir_dataset_clean_se %>% 
  bind_cols(thermal_breadth, thermal_breadth_extremes) %>% 
  mutate(tsm = tmax-topt)
write_csv(thermal_traits_complete, "thermal_traits_complete.csv")

# ~~~~ model assumptions (N, homosced.) ----
# thermal_traits_complete <- read_csv("thermal_traits_complete.csv") %>% 
#   mutate(tsm = tmax-topt)


## tmin
hist(thermal_traits_complete$tmin, breaks = 30)
shapiro.test(thermal_traits_complete$tmin)
## topt
hist(thermal_traits_complete$topt, breaks = 30)
shapiro.test(thermal_traits_complete$topt)
## tmax
hist(thermal_traits_complete$tmax, breaks = 30)
shapiro.test(thermal_traits_complete$tmax)
## thermal breadth
hist(thermal_traits_complete$thermal_breadth, breaks = 30)
shapiro.test(thermal_traits_complete$thermal_breadth)
## thermal breadth
hist(thermal_traits_complete$thermal_breadth, breaks = 30)
shapiro.test(thermal_traits_complete$thermal_breadth)
## a
hist(thermal_traits_complete$a_est, breaks = 30)
shapiro.test(thermal_traits_complete$a_est)

#we need to transform...

# ~~~~ variable transformations ----
## bestnormalize

tmax_scan_normalizer <- bestNormalize(thermal_traits_complete$tmax)
# OrderNorm transformation (ORQ) is automatically applied
hist(tmax_scan_normalizer$x.t, breaks = 30)
## as an example, we select more common transformation suggested by this function:
arsin <- arcsinh_x(thermal_traits_complete$tmax)
bc1 <- boxcox(thermal_traits_complete$tmax)
log1 <- bestNormalize::log_x(thermal_traits_complete$tmax)
sqrt1 <- bestNormalize::sqrt_x(thermal_traits_complete$tmax)
par(mfrow=c(2,3))
hist(arsin$x.t, 30)
hist(bc1$x.t, 30)
hist(log1$x.t, 30)
hist(sqrt1$x.t, 30)
hist(log(thermal_traits_complete$tmax),30)
hist(tmax_scan_normalizer$x.t,30) #
## so we transform every trait
bn_tmax <- bestNormalize(thermal_traits_complete$tmax)
bn_tmin <- bestNormalize(thermal_traits_complete$tmin)
bn_topt <- bestNormalize(thermal_traits_complete$topt)
bn_thermal_breadth <- bestNormalize(thermal_traits_complete$thermal_breadth)
bn_a <- bestNormalize(thermal_traits_complete$a_est)
bn_tsm <- bestNormalize(thermal_traits_complete$tsm)


#then we can model anything but we need to back-transform the estimate with the following structure:
backtransf <- function(trans_var, estimate){
  # transformation 
  bt_param = predict(trans_var,estimate, inverse = TRUE) 
}

# ~~~~ transformed dataset ----
thermal_traits_trans <- thermal_traits_complete %>% 
  mutate(tmax = bn_tmax$x.t,
         tmin = bn_tmin$x.t,
         topt = bn_topt$x.t,
         thermal_breadth = bn_thermal_breadth$x.t,
         a_est = bn_a$x.t,
         tsm = bn_tsm$x.t)
thermal_traits_trans_order <- thermal_traits_trans %>% 
  filter(order == "Acari" |
           order == "Hemiptera" |
           order == "Lepidoptera") %>%
  glimpse()


thermal_traits_trans_fg <- thermal_traits_trans %>% 
  filter(feeding_guild == "borer" |
           feeding_guild == "chewer" |
           feeding_guild == "sucker") %>%
  glimpse()

# 3. Plots  ----
# since they are also exploratory and for visualization purposes, we used the non-transformed dataset
# although analyses were carried out using the ORQ-transformed one.
#_ _ 3.1 Thermal traits across latitude  ----
## tmin
Tmin_to_lat_plot <- ggplot(thermal_traits_complete,aes(abs(lat),tmin))+
  geom_point(alpha=0.032,
             color="skyblue4")+
  geom_smooth(method="lm",color="skyblue4",fill="skyblue2")+
  labs(title= "Tmin ~ latitude",
       x= "latitude",
       y= "Tmin")+
  theme_classic()
Tmin_to_lat_plot 

tmin_to_lat_order_plot <- ggplot(thermal_traits_complete,aes(abs(lat),tmin))+
  geom_point(alpha=0.032,
             color="skyblue4",)+
  geom_smooth(method="lm",color="skyblue4",fill="skyblue2")+
  labs(title= "Tmin ~ latitude",
       x= "latitude",
       y= "Tmin")+
  facet_wrap(.~order)+
  theme_light()
tmin_to_lat_order_plot

## tmax
Tmax_to_lat_plot <- ggplot(thermal_traits_complete,aes(abs(lat),tmax))+
  geom_point(alpha=0.032,
             color="red4")+
  geom_smooth(method="lm",color="red4",fill="red3")+
  labs(title= "Tmax ~ latitude",
       x= "latitude",
       y= "Tmax")+
  theme_classic()
Tmax_to_lat_plot
tmax_to_lat_order_plot <- ggplot(thermal_traits_complete,aes(abs(lat),tmax))+
  geom_point(alpha=0.032,
             color="red4")+
  geom_smooth(method="lm",color="red4",fill="red3")+
  labs(title= "Tmax ~ latitude",
       x= "latitude",
       y= "Tmax")+
  facet_wrap(.~order)+
  theme_light()
tmax_to_lat_order_plot

## topt
Topt_to_lat_plot <- ggplot(thermal_traits_complete,aes(abs(lat),topt))+
  geom_point(alpha=0.032,
             color="mediumorchid4")+ 
  geom_smooth(method="lm",color="mediumorchid4",fill="mediumorchid3")+
  labs(title= "Topt ~ latitude",
       x= "latitude",
       y= "Topt")+
  theme_classic()
Topt_to_lat_plot
topt_to_lat_order_plot <- ggplot(thermal_traits_complete,aes(abs(lat),topt))+
  geom_point(alpha=0.032,
             color="mediumorchid4")+
  geom_smooth(method="lm",color="mediumorchid4",fill="mediumorchid3")+
  labs(title= "Topt ~ latitude",
       x= "latitude",
       y= "Topt")+
  facet_wrap(.~order)+
  theme_light()
topt_to_lat_order_plot

lms_plot_grid <- plot_grid(Tmin_to_lat_plot,Topt_to_lat_plot,Tmax_to_lat_plot,
                           nrow = 1,labels = c("A","B","C"))
lms_plot_grid
all_loess_combined <- ggplot(ir_dataset_clean_se)+
  geom_point(aes(x=abs(lat),y=tmin),color="skyblue4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=tmin),color="skyblue4",fill="skyblue2")+
  geom_point(aes(x=abs(lat),y=tmax),color="red4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=tmax),color="red4",fill="red3")+
  geom_point(aes(x=abs(lat),y=topt),color="mediumorchid4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=topt),color="mediumorchid4",fill="mediumorchid3")+
  labs(title= "Thermal traits ~ latitude",x= "latitude",y= "Thermal traits (?C)")+
  theme_classic()
all_loess_combined  

all_lms_combined <- ggplot(ir_dataset_clean_se)+
  geom_point(aes(x=abs(lat),y=tmin),color="skyblue4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=tmin),color="skyblue4",fill="skyblue2", method = "lm")+
  geom_point(aes(x=abs(lat),y=tmax),color="red4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=tmax),color="red4",fill="red3", method = "lm")+
  geom_point(aes(x=abs(lat),y=topt),color="mediumorchid4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=topt),color="mediumorchid4",fill="mediumorchid3", method = "lm")+
  labs(title= "Thermal traits ~ latitude",x= "latitude",y= "Thermal traits (?C)")+
  theme_classic()
all_lms_combined  

# 4. Meta-analysis models ----
#tampoco es normal :/
# _ _ a) tmax  ---- 
# _ _ _ _i) random intercept  ---- 
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

# _ _ _ _ii) ~ lat  ---- 
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
# _ _ _ _iii) ~ order  ---- 
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

# _ _ _ _iv) ~ feeding guild  ---- 
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


# _ _ _ _v) ~ year  ---- 
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

# _ _ b) tmin  ---- 
# _ _ _ _i) random intercept  ---- 

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

# _ _ _ _ii) ~ lat  ---- 
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
  mutate(across(everything(), ~ backtransf(trans_var = bn_tmin, estimate = .x))) # and this for Tmin at lat=0º
# _ _ _ _iii) ~ order  ---- 
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


# _ _ _ _iv) ~ feeding guild  ---- 
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


  # _ _ _ _v) ~ year  ---- 
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

# _ _ c) topt  ---- 
# _ _ _ _i) random intercept  ---- 

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

# _ _ _ _ii) ~ lat  ---- 
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
  mutate(across(everything(), ~ backtransf(trans_var = bn_topt, estimate = .x))) # and this for topt at lat=0º
# _ _ _ _iii) ~ order  ---- 
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


# _ _ _ _iv) ~ feeding guild  ---- 
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


# _ _ _ _v) ~ year  ---- 
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
# _ _ d) thermal_breadth  ---- 
# _ _ _ _i) random intercept  ---- 
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

# _ _ _ _ii) ~ lat  ---- 
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
  mutate(across(everything(), ~ backtransf(trans_var = bn_thermal_breadth, estimate = .x))) # and this for thermal_breadth at lat=0º
# _ _ _ _iii) ~ order  ---- 
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


# _ _ _ _iv) ~ feeding guild  ---- 
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


# _ _ _ _v) ~ year  ---- 
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

# _ _ e) thermal safety margin  ---- 
# _ _ _ _i) random intercept  ---- 
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

# _ _ _ _ii) ~ lat  ---- 
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
  mutate(across(everything(), ~ backtransf(trans_var = bn_tsm, estimate = .x))) # and this for tsm at lat=0º
# _ _ _ _iii) ~ order  ---- 
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


# _ _ _ _iv) ~ feeding guild  ---- 
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


# _ _ _ _v) ~ year  ---- 
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
