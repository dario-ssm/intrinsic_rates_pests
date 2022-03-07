# SCRIPT INFO --------------------
#     Authors: Dar?o San Segundo Molina, Sara Vill?n P?rez, Ignacio Morales Castilla
#     Title: trait inferences
#     Aim: explore trait trends after model fitting from INTRAPEST database
#     Date: October 2021
#
# Side A (from Dimitropoulo's method)
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
thermal_traits_data_raw <- read_csv("parameters_individual_fitted.csv") %>% # or parameters_individual_fitted.csv (no dimitropoulou's way)
  rename(id=id, a_est = a_est, a_se = a_se, tmin = Tmin_est, tmin_se = Tmin_se,
         tmax = Tmax_est, tmax_se = Tmax_se, topt = Value, topt_se = Topt_se,
         start_a = starting_a, start_tmin = starting_Tmin, start_tmax = starting_Tmax )%>% 
  glimpse()

# ensemble dataset complete
ir_data_all <- read_csv("C:/Users/dario-ssm/Documents/Dario Investigacion/IntRaPest/intrinsic_rates_pests/data/IR_data_all_clean.csv")
authors <- ir_data_all %>% 
  select(Authors) %>% 
  separate(Authors,into = c("a","b"),sep = "," ) %>% 
  select(1) %>% 
  rename(authors = a) %>% 
  print()

ir_data_complementary <- ir_data_all %>% 
  bind_cols(authors) %>% 
  select(id, authors,Year, title, DOI, vi, order, family, genus, species, feeding_guild,
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
#write_csv(ir_data_complete, "parameters_individual_fit_complemented.csv")

# 2. Exploratory data analysis ----
## apply filters 
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

# 3. Linear regressions ----
# _ _ 3.1. Tmin to lat  ---- 
tmin_lat <- lm(tmin~abs(lat),
               data = ir_dataset_clean_se)
tmin_lat_sum <- summary(tmin_lat)
performance::check_model(tmin_lat) #quite good

# _ _ 3.2. Tmax to lat  ---- 
tmax_lat <- lm(tmax~abs(lat),
               data = ir_dataset_clean_se)
tmax_lat_sum <- summary(tmax_lat)
performance::check_model(tmax_lat) 

# _ _ 3.3. Topt to lat  ----
topt_lat <- lm(topt~abs(lat),
               data = ir_dataset_clean_se)
topt_lat_sum <- summary(topt_lat)
performance::check_model(topt_lat) 

# _ _ 3.4. Tmax to Tmin  ----
tmax_tmin <- lm(tmax~tmin, data = ir_dataset_clean_estimates)
tmax_tmin_sum <- summary(tmax_tmin)
check_model(tmax_tmin)
tmax_tmin_sum
# _ _ 3.5. a to lat  ----
a_lat <- lm(a_est~abs(lat),
            data = ir_dataset_clean_se)
a_lat_sum <- summary(a_lat)
performance::check_model(a_lat)


# 4. Plots  ----
#_ _ 4.1 Thermal traits across latitude  ----
## tmin
Tmin_to_lat_plot <- ggplot(ir_dataset_clean_se,aes(abs(lat),tmin))+
  geom_point(alpha=0.1,
             color="skyblue4",
             position = position_jitter(width = 3))+
  geom_smooth(method="lm",color="skyblue4",fill="skyblue2")+
  labs(title= "Tmin ~ latitude",
       x= "latitude",
       y= "Tmin")+
  theme_classic()
Tmin_to_lat_plot #should we remove those whose start = estimate?

tmin_to_lat_order_plot <- ggplot(ir_dataset_clean_se_order,aes(abs(lat),tmin))+
  geom_point(alpha=0.1,
             color="skyblue4",
             position = position_jitter(width = 3))+
  geom_smooth(method="lm",color="skyblue4",fill="skyblue2")+
  labs(title= "Tmin ~ latitude",
       x= "latitude",
       y= "Tmin")+
  facet_wrap(.~order)+
  theme_light()
tmin_to_lat_order_plot

## tmax
Tmax_to_lat_plot <- ggplot(ir_dataset_clean_se,aes(abs(lat),tmax))+
  geom_point(alpha=0.1,
             color="red4",
             position = position_jitter(width = 3))+
  geom_smooth(method="lm",color="red4",fill="red3")+
  labs(title= "Tmax ~ latitude",
       x= "latitude",
       y= "Tmax")+
  theme_classic()
Tmax_to_lat_plot
tmax_to_lat_order_plot <- ggplot(ir_dataset_clean_se_order,aes(abs(lat),tmax))+
  geom_point(alpha=0.1,
             color="red4",
             position = position_jitter(width = 3))+
  geom_smooth(method="lm",color="red4",fill="red3")+
  labs(title= "Tmax ~ latitude",
       x= "latitude",
       y= "Tmax")+
  facet_wrap(.~order)+
  theme_light()
tmax_to_lat_order_plot

## topt
Topt_to_lat_plot <- ggplot(ir_dataset_clean_se,aes(abs(lat),topt))+
  geom_point(alpha=0.1,
             color="mediumorchid4",
             position = position_jitter(width = 3))+ 
  geom_smooth(method="lm",color="mediumorchid4",fill="mediumorchid3")+
  labs(title= "Topt ~ latitude",
       x= "latitude",
       y= "Topt")+
  theme_classic()
Topt_to_lat_plot
topt_to_lat_order_plot <- ggplot(ir_dataset_clean_se_order,aes(abs(lat),topt))+
  geom_point(alpha=0.1,
             color="mediumorchid4",
             position = position_jitter(width = 3))+
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
  geom_point(aes(x=abs(lat),y=tmin),color="skyblue4",alpha=0.1, position = position_jitter(width = 3))+
  geom_smooth(aes(x=abs(lat),y=tmin),color="skyblue4",fill="skyblue2")+
  geom_point(aes(x=abs(lat),y=tmax),color="red4",alpha=0.1, position = position_jitter(width = 3))+
  geom_smooth(aes(x=abs(lat),y=tmax),color="red4",fill="red3")+
  geom_point(aes(x=abs(lat),y=topt),color="mediumorchid4",alpha=0.1, position = position_jitter(width = 3))+
  geom_smooth(aes(x=abs(lat),y=topt),color="mediumorchid4",fill="mediumorchid3")+
  labs(title= "Thermal traits ~ latitude",x= "latitude",y= "Thermal traits (ºC)")+
  theme_classic()
all_loess_combined  

all_lms_combined <- ggplot(ir_dataset_clean_se)+
  geom_point(aes(x=abs(lat),y=tmin),color="skyblue4",alpha=0.1, position = position_jitter(width = 3))+
  geom_smooth(aes(x=abs(lat),y=tmin),color="skyblue4",fill="skyblue2", method = "lm")+
  geom_point(aes(x=abs(lat),y=tmax),color="red4",alpha=0.1, position = position_jitter(width = 3))+
  geom_smooth(aes(x=abs(lat),y=tmax),color="red4",fill="red3", method = "lm")+
  geom_point(aes(x=abs(lat),y=topt),color="mediumorchid4",alpha=0.1, position = position_jitter(width = 3))+
  geom_smooth(aes(x=abs(lat),y=topt),color="mediumorchid4",fill="mediumorchid3", method = "lm")+
  labs(title= "Thermal traits ~ latitude",x= "latitude",y= "Thermal traits (ºC)")+
  theme_classic()
all_lms_combined  


all_lms_combined


#_ _ 4.2 Thermal traits across taxa  ----
#prepare dataframe
counts <- ir_dataset_clean_se %>%
  count(order)%>%
  select(n)%>%
  as_vector()

mean_by_order <- ir_dataset_clean_se %>%
  select(tmin,tmax,topt,tmin_se,tmax_se,topt_se,order) %>% 
  group_by(order) %>%
  summarise(across(c(tmin,tmax,topt,tmin_se,tmax_se,topt_se),
                   ~ mean(.x, na.rm = TRUE)))

## first we need subsets for each trait to plot pointrange correctly
#Tmax
lower_tmax <- tapply(ir_dataset_clean_se$tmax,
                     ir_dataset_clean_se$order, min)
upper_tmax <- tapply(ir_dataset_clean_se$tmax,
                     ir_dataset_clean_se$order, max)
tmax_ranges <- tibble(order = mean_by_order$order,
                      Tmax =mean_by_order$tmax,
                      n=counts,
                      lower=lower_tmax,
                      upper=upper_tmax)

#Tmin
lower_tmin <- tapply(ir_dataset_clean_se$tmin,
                     ir_dataset_clean_se$order, min)
upper_tmin <- tapply(ir_dataset_clean_se$tmin,
                     ir_dataset_clean_se$order, max)
tmin_ranges <- tibble(order = mean_by_order$order,
                      Tmin =mean_by_order$tmin,
                      n=counts,
                      lower=lower_tmin,
                      upper=upper_tmin)
#Topt
lower_topt <- tapply(ir_dataset_clean_se$topt,
                     ir_dataset_clean_se$order, min)
upper_topt <- tapply(ir_dataset_clean_se$topt,
                     ir_dataset_clean_se$order, max)
topt_ranges <- tibble(order = mean_by_order$order,
                      Topt =mean_by_order$topt,
                      n=counts,
                      lower=lower_topt,
                      upper=upper_topt)

#now we plot
CombinedAcross_order_traits <- ggplot()+
  geom_point(data=tmax_ranges ,aes(x=order,y=Tmax,size=n),
             color="indianred3",
             position = position_nudge(x=-0.25))+
  geom_pointrange(data=tmax_ranges ,aes(x=order,y=Tmax,ymin=lower,ymax=upper),
                  color="indianred3",
                  position = position_nudge(x=-0.25))+
  geom_point(data=topt_ranges ,aes(x=order,y=Topt,size=n),
             color="mediumorchid4")+
  geom_pointrange(data=topt_ranges ,aes(x=order,y=Topt,ymin=lower,ymax=upper),
                  color="mediumorchid4")+
  geom_point(data=tmin_ranges ,aes(x=order,y=Tmin,size=n),
             color="turquoise4",
             position = position_nudge(x=0.25))+
  geom_pointrange(data=tmin_ranges ,aes(x=order,y=Tmin,ymin=lower,ymax=upper),
                  color="turquoise4",
                  position = position_nudge(x=0.25))+
  labs(title = "Thermal traits across taxa",
       subtitle= "Circles' size are proportional to sample size",
       x = "Order",
       y = "Thermal traits (?C)")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))+
  geom_text(data=tmin_ranges,aes(x=order,
                                 y=47,
                                 label = paste("n=",n, sep = ""),
                                 fontface=3),
            color="lightslategrey",
            parse = FALSE)+
  theme(legend.position = "none")
CombinedAcross_order_traits

ggsave("thermal_traits_acrossTaxa.png",dpi=300,
       width = 20,height = 20,units="cm")

# boxplots
boxplot_traits_ord <- ggplot(data=thermal_traits_data, aes(x=order))+
  geom_boxplot(aes(y=Tmin_est_nls),
               fill= "turquoise4")+
  geom_boxplot(aes(y=Tmax_est_nls),
               fill="indianred3")+
  labs(title = "Thermal traits across taxa",
       x = "Order",
       y = "Thermal traits (?C)")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))+
  theme(legend.position = "none")

boxplot_traits_ord

#pointranges with mean se
tmax_se_avg <- tapply(ir_dataset_clean_se$tmax_se,
                      ir_dataset_clean_se$order,
                      mean)
tmax_se_pointrange <- tibble(order = mean_by_order$order,
                             Tmax =mean_by_order$tmax,
                             n=counts,
                             tmax_se_avg)
tmin_se_avg <- tapply(ir_dataset_clean_se$tmin_se,
                      ir_dataset_clean_se$order,
                      mean)
tmin_se_pointrange <- tibble(order = mean_by_order$order,
                             Tmin =mean_by_order$tmin,
                             n=counts,
                             tmin_se_avg)
topt_se_avg <- tapply(ir_dataset_clean_se$topt_se,
                      ir_dataset_clean_se$order,
                      mean)
topt_se_pointrange <- tibble(order = mean_by_order$order,
                             Topt =mean_by_order$topt,
                             n=counts,
                             topt_se_avg)

traits_across_orders_se <- ggplot()+
  geom_point(data=tmax_se_pointrange ,aes(x=order,y=Tmax,size=n),
             color="indianred3",
             position = position_nudge(x=-0.25))+
  geom_pointrange(data=tmax_se_pointrange ,aes(x=order,y=Tmax,
                                               ymin=Tmax-tmax_se_avg,
                                               ymax=Tmax+tmax_se_avg),
                  color="indianred3",
                  position = position_nudge(x=-0.25))+
  geom_point(data=topt_se_pointrange ,aes(x=order,y=Topt,size=n),
             color="mediumorchid4")+
  geom_pointrange(data=topt_se_pointrange ,aes(x=order,y=Topt,
                                               ymin=Topt-topt_se_avg,
                                               ymax=Topt+topt_se_avg),
                  color="mediumorchid4")+
  geom_point(data=tmin_se_pointrange ,aes(x=order,y=Tmin,size=n),
             color="turquoise4",
             position = position_nudge(x=0.25))+
  geom_pointrange(data=tmin_se_pointrange ,aes(x=order,y=Tmin,
                                               ymin=Tmin-tmin_se_avg,
                                               ymax=Tmin+tmin_se_avg),
                  color="turquoise4",
                  position = position_nudge(x=0.25))+
  labs(title = "Thermal traits across taxa",
       subtitle= "Error bars represent average standard error within each group. 
Circles' size are proportional to sample size",
       x = "Order",
       y = "Thermal traits (?C)")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))+
  geom_text(data=tmax_se_pointrange,aes(x=order,
                                        y=40,
                                        label = paste("n=",n, sep = ""),
                                        fontface=3),
            color="lightslategrey",
            parse = FALSE)+
  theme(legend.position = "none")
traits_across_orders_se
ggsave("traits_across_orders_se.png", dpi = 300,
       width = 20,height = 20,units="cm")

## ANOVAs ##
tmin_4anova <- thermal_traits_data %>%
  select(order,feeding_guild,Tmin_est_nls,Tmin_se_nls)
tmin_anova_order <- lm(Tmin_est_nls~order,data=tmin_4anova)
summary(tmin_anova_order) #no differences
tmin_anova_guild <- lm(Tmin_est_nls~feeding_guild,data=tmin_4anova)
summary(tmin_anova_guild) #no differences

tmax_4anova <- thermal_traits_data %>%
  select(order,feeding_guild,Tmax_est_nls,Tmax_se_nls)
tmax_anova_order <- lm(Tmax_est_nls~order,data=tmax_4anova)
summary(tmax_anova_order) #no differences
tmax_anova_guild <- lm(Tmax_est_nls~feeding_guild,data=tmax_4anova)
summary(tmax_anova_guild) #only chewers are different from the others, with lower tmax

# 5. Meta-analysis models ----
# _ _ a) tmax  ---- 
# _ _ _ _i) random intercept  ---- 

tmax_intercept <- lme(tmax ~ 1,
                      random = ~1|id,
                      weights = varFixed(~vi),
                      data = ir_dataset_clean_se,
                      control = lmeControl(sigma = 1))
summary(tmax_intercept)
performance::check_model(tmax_intercept)
# _ _ _ _ii) ~ lat  ---- 
## random slope & intercept
tmax_lat_slope <- lme(tmax ~ abs(lat),
                      random = ~abs(lat)|id,
                      weights = varFixed(~vi),
                      data = ir_dataset_clean_se,
                      control = lmeControl(sigma = 1))
summary(tmax_lat_slope)
# _ _ _ _iii) ~ order  ---- 
## random slope & intercept
tmax_order <- lme(tmax ~ as_factor(order),
                  random = ~ as_factor(order)|id,
                  weights = varFixed(~vi),
                  data = ir_dataset_clean_se_order,
                  control = lmeControl(sigma = 1))
summary(tmax_order)
coefs_tmax_order <- as_tibble(coefficients(summary(tmax_order))) %>% 
  mutate(estimate = 33.419891+Value)
performance::check_model(tmin_order)
# _ _ _ _iv) ~ feeding guild  ---- 
## random slope & intercept
tmax_feeding_guild <- lme(tmax ~ as_factor(feeding_guild),
                          random = ~ as_factor(feeding_guild)|id,
                          weights = varFixed(~vi),
                          data = ir_dataset_clean_se_fg,
                          control = lmeControl(sigma = 1))
summary(tmax_feeding_guild)
coefs_tmax_fg <- as_tibble(coefficients(summary(tmax_feeding_guild))) %>% 
  mutate(estimate = 33.419891+Value)
performance::check_model(tmin_order)

# _ _ _ _v) ~ year  ---- 
## random slope & intercept
tmax_year <- lme(tmax ~ Year,
                 random = ~ Year|id,
                 weights = varFixed(~vi),
                 data = ir_dataset_clean_se,
                 control = lmeControl(sigma = 1))
summary(tmax_year)
performance::check_model(tmax_year)
# _ _ b) tmin  ---- 
# _ _ _ _i) random intercept  ---- 

tmin_intercept <- lme(tmin ~ 1,
                      random = ~1|id,
                      weights = varFixed(~vi),
                      data = ir_dataset_clean_se,
                      control = lmeControl(sigma = 1))
summary(tmin_intercept)
performance::check_model(tmin_intercept)

# _ _ _ _ii) ~ lat  ---- 
## random slope & intercept
tmin_lat_slope <- lme(tmin ~ abs(lat),
                      random = ~abs(lat)|id,
                      weights = varFixed(~vi),
                      data = ir_dataset_clean_se,
                      control = lmeControl(sigma = 1))
summary(tmin_lat_slope)
performance::check_model(tmin_intercept)

# _ _ _ _iii) ~ order  ---- 
## random slope & intercept
tmin_order <- lme(tmin ~ as_factor(order),
                  random = ~ as_factor(order)|id,
                  weights = varFixed(~vi),
                  data = ir_dataset_clean_se_order,
                  control = lmeControl(sigma = 1))
summary(tmin_order)
anova(tmin_order)

cld(lsmeans(tmin_order, specs=~as_factor(order)))
coefs_tmin_order <- as_tibble(coefficients(summary(tmin_order))) %>% 
  mutate(estimate = 10.814+Value)
performance::check_model(tmin_order)

# _ _ _ _iv) ~ feeding guild  ---- 
## random slope & intercept
tmin_feeding_guild <- lme(tmin ~ as_factor(feeding_guild),
                          random = ~ as_factor(feeding_guild)|id,
                          weights = varFixed(~vi),
                          data = ir_dataset_clean_se_fg,
                          control = lmeControl(sigma = 1))
summary(tmin_feeding_guild)
anova(tmin_feeding_guild)
coefs_tmin_feeding_guild <- as_tibble(coefficients(summary(tmin_feeding_guild))) %>% 
  mutate(estimate = 12.369563  +Value)
  # _ _ _ _v) ~ year  ---- 
## random slope & intercept
tmin_year <- lme(tmin ~ Year,
                 random = ~ Year|id,
                 weights = varFixed(~vi),
                 data = ir_dataset_clean_se,
                 control = lmeControl(sigma = 1))
summary(tmin_year)


# _ _ c) topt  ---- 
# _ _ _ _i) random intercept  ---- 

topt_intercept <- lme(topt ~ 1,
                      random = ~1|id,
                      weights = varFixed(~vi),
                      data = ir_dataset_clean_estimates,
                      control = lmeControl(sigma = 1))
summary(topt_intercept)
# _ _ _ _ii) ~ lat  ---- 
## random slope & intercept
topt_lat_slope <- lme(topt ~ abs(lat),
                      random = ~abs(lat)|id,
                      weights = varFixed(~vi),
                      data = ir_dataset_clean_se,
                      control = lmeControl(sigma = 1))
summary(topt_lat_slope)
# _ _ _ _iii) ~ order  ---- 
## random slope & intercept
topt_order <- lme(topt ~ as_factor(order),
                  random = ~ as_factor(order)|as_factor(id),
                  weights = varFixed(~vi),
                  data = ir_dataset_clean_se_order,
                  control = lmeControl(sigma = 1))
summary(topt_order)
coefs_topt_order <- as_tibble(coefficients(summary(topt_order))) %>% 
  mutate(estimate = 28.200073   +Value)
# _ _ _ _iv) ~ feeding guild  ---- 
## random slope & intercept
topt_feeding_guild <- lme(topt ~ as_factor(feeding_guild),
                          random = ~ as_factor(feeding_guild)|id,
                          weights = varFixed(~vi),
                          data = ir_dataset_clean_se_fg,
                          control = lmeControl(sigma = 1))
summary(topt_feeding_guild)
coefs_topt_order <- as_tibble(coefficients(summary(topt_feeding_guild))) %>% 
  mutate(estimate = 27.684007    +Value)
# _ _ _ _v) ~ year  ---- 
## random slope & intercept
topt_year <- lme(topt ~ Year,
                 random = ~ Year|id,
                 weights = varFixed(~vi),
                 data = ir_dataset_clean_se,
                 control = lmeControl(sigma = 1))
summary(topt_year)

