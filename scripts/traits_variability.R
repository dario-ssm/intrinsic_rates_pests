# SCRIPT INFO --------------------
#     Authors: Dar?o San Segundo Molina, Sara Vill?n P?rez, Ignacio Morales Castilla
#     Title: trait inferences: Variability assessment
#     Aim: explore traits variability to infer ecological patterns
#     Date: March 2022
#
# Side A (from Dimitropoulo's method)
#_________________________ ####

# 1. Load Datasets ----
rm(list=ls())
#library(tidyverse)
library(tidyverse)
library(performance)
library(cowplot)
library(nlme)
library(ggthemes)
library(car)
setwd(paste0(getwd(),"/data"))
# ~~~~ a) traits inferred from our simulations ----
thermal_traits_simulated <- read_csv("parameters_individual_fit_complemented.csv") %>% 
  glimpse()
q_tmin_se <-quantile(thermal_traits_simulated$tmin_se, probs = .95)
q_tmax_se <- quantile(thermal_traits_simulated$tmax_se,probs = .95)
q_topt_se <- quantile(thermal_traits_simulated$topt_se,probs = .95)
thermal_traits_clean <- thermal_traits_simulated %>% 
  filter(tmax_se < q_tmax_se &
           tmin_se < q_tmin_se &
           topt_se < q_topt_se) %>%
  glimpse()

thermal_traits_clean_order <- thermal_traits_clean %>% 
  filter(order == "Acari" |
           order == "Hemiptera" |
           order == "Lepidoptera") %>%
  glimpse()

thermal_traits_clean_fg <- thermal_traits_clean %>% 
  filter(feeding_guild == "borer" |
           feeding_guild == "chewer" |
           feeding_guild == "sucker") %>%
  glimpse()

# ~~~~ b) pooled simulated dataset  ----

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

# ~~~~ c) Equations (Brière-1)  --------------------------------------------

briere1 <- function(a, temp, Tmin, Tmax){
  a*temp*(temp-Tmin)*(Tmax-temp)^(1/2)
}
# Since Topt is not a direct parameter of the model, it can be derived from Tmin and Tmax
# according to Marchioro & Foerster, 2011:
Topt <- function(Tmin,Tmax,m){
  Topt=((2*m*Tmax+(m+1)*Tmin)+sqrt(4*(m^2)*(Tmax^2)+((m+1)^2)*(Tmin^2)-4*(m^2)*Tmin*Tmax))/(4*m+2)
  return(Topt)
}

# 2. Graphical evaluation ----

# ~~~~ a) thermal traits simulated ----
# ~~~~~~~~ (1) traits variability ~ order~----
boxplot_traits_ord<- ggplot(data=thermal_traits_clean_order, aes(x=order))+
  geom_boxplot(aes(y=tmin),
               fill= "turquoise4",
               color = "lightgrey")+
  geom_boxplot(aes(y=tmax),
               fill="indianred3",
               color = "lightgrey")+
  labs(title = "Thermal traits across taxa",
       x = "Order", 
       y = "Thermal traits (ºC)")+
  ggthemes::theme_few()+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))+
  theme(legend.position = "none")
boxplot_traits_ord
ggsave("boxplot_traits_order.png", dpi = 300, width = 20, height = 20, units = "cm")
iq_range_tmin_order <- quantile(thermal_traits_clean_order$tmin, .75) - quantile(thermal_traits_clean_order$tmin, .25)
iq_range_tmax_order <- quantile(thermal_traits_clean_order$tmax, .75) - quantile(thermal_traits_clean_order$tmax, .25)

violin_traits_ord<- ggplot(data=thermal_traits_clean_order, aes(x=order))+
  geom_violin(aes(y=tmin),
              fill= "turquoise4",
              color = "lightgrey")+
  geom_violin(aes(y=tmax),
              fill="indianred3",
              color = "lightgrey")+
  labs(title = "Thermal traits across taxa",
       x = "Order", 
       y = "Thermal traits (ºC)")+
  ggthemes::theme_few()+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))+
  theme(legend.position = "none")
violin_traits_ord
ggsave("violin_traits_order.png", dpi = 300, width = 20, height = 20, units = "cm")
# ~~~~~~~~ (2) traits variability ~ feeding guild ~----
boxplot_traits_fg <- ggplot(data=thermal_traits_clean_fg, aes(x=feeding_guild))+
  geom_boxplot(aes(y=tmin),
               fill= "turquoise4",
               color = "lightgrey")+
  geom_boxplot(aes(y=tmax),
               fill="indianred3",
               color = "lightgrey")+
  labs(title = "Thermal traits across taxa",
       x = "Feeding guild", 
       y = "Thermal traits (ºC)")+
  ggthemes::theme_few()+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))+
  theme(legend.position = "none")
boxplot_traits_fg
ggsave("boxplot_traits_feeding_guild.png", dpi = 300, width = 20, height = 20, units = "cm")
iq_range_tmin_fg <- quantile(thermal_traits_clean_fg$tmin, .75) - quantile(thermal_traits_clean_fg$tmin, .25)
iq_range_tmax_fg <- quantile(thermal_traits_clean_fg$tmax, .75) - quantile(thermal_traits_clean_fg$tmax, .25)

violin_traits_fg <- ggplot(data=thermal_traits_clean_fg, aes(x=feeding_guild))+
  geom_violin(aes(y=tmin),
              fill= "turquoise4",
              color = "lightgrey")+
  geom_violin(aes(y=tmax),
              fill="indianred3",
              color = "lightgrey")+
  labs(title = "Thermal traits across feeding guilds",
       x = "Feeding guild", 
       y = "Thermal traits (ºC)")+
  ggthemes::theme_few()+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))+
  theme(legend.position = "none")
violin_traits_fg
ggsave("violin_traits_fg.png", dpi = 300, width = 20, height = 20, units = "cm")

# ~~~~~~~~ (3) traits variability ~ lat ~----
all_lms_combined <- ggplot(thermal_traits_clean)+
  geom_point(aes(x=abs(lat),y=tmin),color="skyblue4",alpha=0.1, position = position_jitter(width = 3))+
  geom_smooth(aes(x=abs(lat),y=tmin),color="skyblue4",fill="skyblue2", method = "lm")+
  geom_point(aes(x=abs(lat),y=tmax),color="red4",alpha=0.1, position = position_jitter(width = 3))+
  geom_smooth(aes(x=abs(lat),y=tmax),color="red4",fill="red3", method = "lm")+
  geom_point(aes(x=abs(lat),y=topt),color="mediumorchid4",alpha=0.02, position = position_jitter(width = 3))+
  geom_smooth(aes(x=abs(lat),y=topt),color="mediumorchid4",fill="mediumorchid3", method = "lm", alpha = 0.01)+
  labs(title= "Thermal traits ~ latitude",x= "latitude",y= "Thermal traits (ºC)")+
  theme_few()
all_lms_combined  
ggsave("lms_traits_lat.png", dpi = 300, width = 10, height = 20, units = "cm")

# ~~~~~~~~ (4) traits variability ~ 1 ----

iq_range_tmin_all <- quantile(thermal_traits_clean$tmin, .75) - quantile(thermal_traits_clean$tmin, .25)
iq_range_tmax_all <- quantile(thermal_traits_clean$tmax, .75) - quantile(thermal_traits_clean$tmax, .25)

violin_traits_all <- ggplot(data=thermal_traits_clean,aes(x=1))+
  geom_violin(aes(y=tmin),
              fill= "turquoise4",
              color = "lightgrey")+
  geom_violin(aes(y=tmax),
              fill="indianred3",
              color = "lightgrey")+
  labs(title = "Thermal traits across taxa",
       x = element_blank(),
       y = "Thermal traits (ºC)")+
  ggthemes::theme_few()+
  theme(axis.text.x = element_blank())+
  theme(legend.position = "none")
violin_traits_all
ggsave("violin_traits_all.png", dpi = 300, width = 10, height = 20, units = "cm")

# ~~~~ b) Pooled simulated dataset ----
# note that for this approach we need to infer parameters with nonlinear regression 
# (package nlme) and then use fitted values to study heterogeneity

# ~~~~~~~~ (1) traits variability ~ order~----
load(file = "nlme_fit_sims_order.RData")
sum_nlme_fit_sims_order <-summary(nlme_fit_sims_order)
sum_nlme_order_df <- as_tibble(sum_nlme_fit_sims_order$tTable) 
nlme_fit_sims_order
## let's ensemble the fixeff
#a
a_order_23 <- sum_nlme_order_df %>% 
  slice(2:3) %>% 
  mutate(estimate = 0.000209+Value)
a_order_1 <- sum_nlme_order_df %>% 
  slice(1) %>% 
  mutate(estimate = 0.000209)
a_order <- a_order_1 %>% 
  bind_rows(a_order_23) %>% 
  mutate(parameter = rep("a", 3))
#tmin
tmin_order_23 <- sum_nlme_order_df %>% 
  slice(5:6) %>% 
  mutate(estimate = 11.8 +Value)
tmin_order_1 <- sum_nlme_order_df %>% 
  slice(4) %>% 
  mutate(estimate = 11.8 )
tmin_order <- tmin_order_1 %>% 
  bind_rows(tmin_order_23) %>% 
  mutate(parameter = rep("tmin", 3))
#tmax
tmax_order_23 <- sum_nlme_order_df %>% 
  slice(8:9) %>% 
  mutate(estimate = 36.7 +Value)
tmax_order_1 <- sum_nlme_order_df %>% 
  slice(7) %>% 
  mutate(estimate = 36.7 )
tmax_order <- tmax_order_1 %>% 
  bind_rows(tmax_order_23) %>% 
  mutate(parameter = rep("tmax", 3))
## assembly:
nlme_order_estimates <- a_order %>% 
  bind_rows(tmin_order, tmax_order) %>% 
  select(-Value, - DF) %>% 
  relocate(5,4,1,2,3) %>% 
  print()
## doooooon't know what to dooo

# ~~~~~~~~ (2) traits variability ~ 1~----
simulated_intrapest
load(file = "nlme_fit_sims.RData")
sum_nlme_fit_sims <- summary(nlme_fit_sims)
relative_taus <- as_tibble(as.numeric(VarCorr(nlme_fit_sims)[1:3,2])) %>% 
  mutate(estimate = fixed.effects(nlme_fit_sims)) %>% 
  mutate(rel_tau = value/estimate) %>% 
  mutate(parameter = c("a", "tmin", "tmax")) %>% 
  print()
#prediction
# boxplot_rel_tau <- ggplot(relative_taus, aes(x = 1, y = estimate))+
#   geom_point(aes(color = parameter))+
#   geom_pointrange(aes(ymin = estimate - rel_tau,
#                       ymax = estimate + rel_tau,
#                       color = parameter))+
#   facet_wrap(.~parameter)+
#   theme_few()
# boxplot_rel_tau
# 
plot_taus <- ggplot(relative_taus, aes(x = parameter, y = rel_tau))+
  geom_point(aes(color = parameter))+
  geom_line()+
  theme_few()

# 3. Statistical Tests ----
# ~~~~ a) thermal traits simulated ----
# ~~~~~~~~  (1) F-test (see Herrando-Pérez, 2019) ----
variances_tmin <- thermal_traits_clean %>% 
  select(tmin, vi, id) %>% 
  mutate(parameter = rep("tmin")) %>% 
  rename(estimate = tmin)
variances_tmax <- thermal_traits_clean %>% 
  select(tmax, vi, id) %>% 
  mutate(parameter = rep("tmax")) %>% 
  rename(estimate = tmax)
variances_traits = variances_tmin %>% 
  bind_rows(variances_tmax) %>% 
  glimpse()
#test variance differences
var.test(estimate ~ parameter, variances_traits, 
         alternative = "greater") #que tmin tiene más varianza que tmax
# ~~~~~~~~  (2) Levene's test (see Hoffmann, 2013) ----

leveneTest(y = estimate ~ parameter, data = variances_traits) # also significantly different :)

# ~~~~~~~~ (3) lme modelling output ----
lme_4vars <- lme(tmin ~ 1,
                 random = ~1|id,
                 weights = varFixed(~vi),
                 control = lmeControl(sigma = 1),
                 data = thermal_traits_clean)
lme_4vars2 <- lme(tmax ~ 1,
                 random = ~1|id,
                 weights = varFixed(~vi),
                 control = lmeControl(sigma = 1),
                 data = thermal_traits_clean)


tau_rel_tmin <- as.numeric(VarCorr(lme_4vars)[1,2])/lme_4vars$coefficients$fixed
tau_rel_tmax <- as.numeric(VarCorr(lme_4vars2)[1,2])/lme_4vars2$coefficients$fixed
# similar...