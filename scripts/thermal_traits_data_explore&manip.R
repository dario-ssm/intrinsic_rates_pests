# SCRIPT INFO --------------------
#     Authors: Dario San Segundo Molina, Sara Villen Perez, Ignacio Morales Castilla
#     Title: thermal traits pararemeterised data exploration, manipulation & variable transformation
#     Aim: explore trait trends after model fitting from INTRAPEST database
#     Date: March 2022
#
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# 1. Load Dataset & explore ranges ----
rm(list=ls())
#library(tidyverse)
library(tidyverse)
library(performance)
library(nlme)
library(bestNormalize)
library(networkD3)

#shortcut for later analyses in thermal_traits_analysis_sideB script 
## AVOID to run it otherwise
thermal_traits_complete <- read_csv("thermal_traits_complete.csv")
bn_tmax <- bestNormalize(thermal_traits_complete$tmax)
bn_tmin <- bestNormalize(thermal_traits_complete$tmin)
bn_topt <- bestNormalize(thermal_traits_complete$topt)
bn_thermal_breadth <- bestNormalize(thermal_traits_complete$thermal_breadth)
bn_a <- bestNormalize(thermal_traits_complete$a_est)
bn_tsm <- bestNormalize(thermal_traits_complete$tsm)
bn_therm_range <- bestNormalize(thermal_traits_complete$therm_range)
bn_therm_window <- bestNormalize(thermal_traits_complete$therm_window)
#then we can model anything but we need to back-transform the estimate with the following structure:
backtransf <- function(trans_var, estimate){
  # transformation 
  bt_param = predict(trans_var,estimate, inverse = TRUE) 
}
## End of shortcut


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

# 2. Exploratory data analysis & manipulation ----
## see relationships order vs. feeding guild
struir_data_complete
counts_fg <- thermal_traits_complete %>%
  group_by(id) %>% 
  summarise(feeding_guild = unique(feeding_guild)) %>% 
  count(feeding_guild)
structures <- thermal_traits_complete %>%
  select(id, order, feeding_guild, genus, species) %>% 
  group_by(id) %>% 
  summarise_all(unique) %>% 
  count(feeding_guild,order) %>% 
  print()

borers_composition <- structures %>% 
  filter(feeding_guild == "borer") %>% 
  mutate(perc= round(100*n/sum(n),1)) %>% 
  print()
chewers_composition <- structures %>% 
  filter(feeding_guild == "chewer") %>% 
  mutate(perc= round(100*n/sum(n),1)) %>% 
  print()
suckers_composition <- structures %>% 
  filter(feeding_guild == "sucker") %>% 
  mutate(perc= round(100*n/sum(n),1)) %>% 
  print()

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
counts_by_order <- ir_dataset_clean_se_order %>% 
  group_by(id) %>% 
  summarise(order = unique(order)) %>% 
  count(order)
#create map along coordinates (prepared to insert in a colored background with transparent ocean)
# map <- ggplot(data = thermal_traits_data, aes(x = lon, y = lat)) +
#   borders("world", colour = "transparent", fill = "white") +
#   geom_point(aes(colour = order), alpha = 0.35, size = 3)+theme_classic()+
#   theme_void()+
#   theme(legend.position = "bottom",panel.background = element_rect(fill="transparent",color=NA))
#   
# map
#ggsave("map_meta.png",height=15,width=25,units="cm",bg = "transparent")

# ~~~~ a) incorporate thermal breadth ----
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

thermal_traits_complete <- thermal_traits_complete %>% 
  bind_cols(thermal_breadth, thermal_breadth_extremes) %>% 
  mutate(tsm = tmax-topt,
         therm_range = tmax-tmin,
         therm_window = topt-tmin)
write_csv(thermal_traits_complete, "thermal_traits_complete.csv")
thermal_traits_order <- thermal_traits_complete %>% 
  filter(order == "Acari" |
           order == "Hemiptera" |
           order == "Lepidoptera")
thermal_traits_fg <- thermal_traits_complete %>% 
  filter(feeding_guild == "borer" |
           feeding_guild == "chewer" |
           feeding_guild == "sucker")
# ~~~~ b) model assumptions (N, homosced.) ----
# thermal_traits_complete <-read_csv("thermal_traits_complete.csv")
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

# ~~~~ c) variable transformations ----
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
bn_therm_range <- bestNormalize(thermal_traits_complete$therm_range)
bn_therm_window <- bestNormalize(thermal_traits_complete$therm_window)
#then we can model anything but we need to back-transform the estimate with the following structure:
backtransf <- function(trans_var, estimate){
  # transformation 
  bt_param = predict(trans_var,estimate, inverse = TRUE) 
}

# ~~~~ d) transformed dataset ----
thermal_traits_trans <- thermal_traits_complete %>% 
  mutate(tmax = bn_tmax$x.t,
         tmin = bn_tmin$x.t,
         topt = bn_topt$x.t,
         thermal_breadth = bn_thermal_breadth$x.t,
         a_est = bn_a$x.t,
         tsm = bn_tsm$x.t,
         therm_range = bn_therm_range$x.t,
         therm_window = bn_therm_window$x.t) 
write_csv(thermal_traits_trans, "thermal_traits_trans.csv")

# subset for most represented taxa
thermal_traits_trans_order <- thermal_traits_trans %>% 
  filter(order == "Acari" |
           order == "Hemiptera" |
           order == "Lepidoptera") %>%
  glimpse()
write_csv(thermal_traits_trans_order, "thermal_traits_trans_order.csv")

# subset for most represented feeding guilds
thermal_traits_trans_fg <- thermal_traits_trans %>% 
  filter(feeding_guild == "borer" |
           feeding_guild == "chewer" |
           feeding_guild == "sucker") %>%
  glimpse()
write_csv(thermal_traits_trans_fg, "thermal_traits_trans_fg.csv")
