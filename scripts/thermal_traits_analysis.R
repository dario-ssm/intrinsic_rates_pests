#### SCRIPT INFO ####     
#     Authors: Darío San Segundo Molina, Sara Villén Pérez, Ignacio Morales Castilla
#     Title: trait inferences
#     Aim: explore trait trends after model fitting from INTRAPEST database
#     Date: October 2021
#_________________________ ####



# 1. Load Dataset ----
rm(list=ls())
#library(tidyverse)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)
library(ggmap)
library(performance)
library(cowplot)
library(wesanderson)
thermal_traits_data_raw <- read_csv("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/dataset_params_metaanalysis.csv") %>%
  glimpse()
## apply filters 
#NOT NECESSARY once the Analyses_nls&bootstrap script will have run
#acari_raw <- thermal_traits_data_raw %>% #ensemble all acari
 # filter(order == "Acari>Prostigmata" | 
  #         order == "Acari>Trombidiformes")%>%
  #mutate(order = "Acari")
#thermal_traits_data <- thermal_traits_data_raw %>%
 # filter(order != "Acari>Prostigmata" &
  #         order != "Acari>Trombidiformes")%>%
  # bind_rows(acari_raw)
thermal_traits_data <- thermal_traits_data_raw %>%
  filter(Tmin_est_nls >0 &
         Tmin_se_nls <10 &
         Tmax_se_nls <10 &
         Tmax_est_nls < 45) %>%
  glimpse() #15 studies discarded

Tmin_est <- thermal_traits_data$Tmin_est_nls
Tmin_se <- thermal_traits_data$Tmin_se_nls
Tmax_est <- thermal_traits_data$Tmax_est_nls
Tmax_se <- thermal_traits_data$Tmax_se_nls
Topt_est <- thermal_traits_data$Topt_est
Topt_se <- thermal_traits_data$Topt_se_delta
lat <- thermal_traits_data$lat
lon <- thermal_traits_data$lon

# 2. Exploratory data analysis ----
# _ _ a) Tmin ----
summary(Tmin_est)
hist(Tmin_est) #quite normal
shapiro.test(Tmin_est) #normal
# _ _ b) Tmax ----
summary(Tmax_est)
hist(Tmax_est)   #quite normal
shapiro.test(Tmax_est) #normal
# _ _ c) lon, lat ----
summary(lat)
hist(lat)    
#create map along coordinates (prepared to insert in a colored background with transparent ocean)
map <- ggplot(data = thermal_traits_data, aes(x = lon, y = lat)) +
  borders("world", colour = "transparent", fill = "white") +
  geom_point(aes(colour = order), alpha = 0.35, size = 3)+theme_classic()+
  theme_void()+
  theme(legend.position = "bottom",panel.background = element_rect(fill="transparent",color=NA))
  
map
#ggsave("map_meta.png",height=15,width=25,units="cm",bg = "transparent")
# _ _ d) Topt ----
summary(Topt_est)
hist(Topt_est) #quite normal
shapiro.test(Topt_est) #normal

# 3. Linear regressions ----
# _ _ 3.1. Tmin to lat  ---- 
Tmin_lat <- lm(Tmin_est_nls~abs(lat), data = thermal_traits_data)
check_model(Tmin_lat) #quite good
Tmin_lat_sum <- summary(Tmin_lat)
Tmin_lat_sum

# _ _ 3.2. Tmax to lat  ---- 
Tmax_lat <- lm(Tmax_est_nls~abs(lat), data = thermal_traits_data)
check_model(Tmax_lat) #quite good
Tmax_lat_sum <- summary(Tmax_lat)
Tmax_lat_sum

# _ _ 3.3. Topt to lat  ----
Topt_lat <- lm(Topt_est~abs(lat), data = thermal_traits_data)
check_model(Topt_lat)
Topt_lat_sum <- summary(Topt_lat)
Topt_lat_sum

# _ _ 3.4. Tmax to Tmin  ----
Tmax_Tmin <- lm(Tmax_est_nls~Tmin_est_nls, data = thermal_traits_data)
check_model(Tmax_Tmin)
Tmax_Tmin_sum <- summary(Tmax_Tmin)
Tmax_Tmin_sum

# 4. Plots  ----
#_ _ 4.1 Thermal traits across latitude  ----

Tmin_to_lat_plot <- ggplot(thermal_traits_data,aes(abs(lat),Tmin_est_nls))+
  geom_point(alpha=0.5,color="skyblue4")+
  geom_smooth(method="lm",color="skyblue4",fill="skyblue2")+
  labs(title= "Tmin ~ latitude",
       x= "latitude",
       y= "Tmin")+
  theme_classic()
Tmin_to_lat_plot

Tmax_to_lat_plot <- ggplot(thermal_traits_data,aes(abs(lat),Tmax_est_nls))+
  geom_point(alpha=0.5,color="red3")+
  geom_smooth(method="lm",color="red4",fill="red3")+
  labs(title= "Tmax ~ latitude",
       x= "latitude",
       y= "Tmax")+
  theme_classic()
Tmax_to_lat_plot

Topt_to_lat_plot <- ggplot(thermal_traits_data,aes(abs(lat),Topt_est))+
  geom_point(alpha=0.5,color="mediumorchid4")+
  geom_smooth(method="lm",color="mediumorchid4",fill="mediumorchid3")+
  labs(title= "Topt ~ latitude",
       x= "latitude",
       y= "Topt")+
  theme_classic()
Topt_to_lat_plot

lms_plot_grid <- plot_grid(Tmin_to_lat_plot,Topt_to_lat_plot,Tmax_to_lat_plot,
                            nrow = 1,labels = c("A","B","C"))

lms_plot_grid

all_loess_combined <- ggplot(thermal_traits_data)+
  geom_point(aes(x=abs(lat),y=Tmin_est_nls),color="skyblue4",alpha=0.5)+
  geom_smooth(aes(x=abs(lat),y=Tmin_est_nls),color="skyblue4",fill="skyblue2")+
  geom_point(aes(x=abs(lat),y=Tmax_est_nls),color="red4",alpha=0.5)+
  geom_smooth(aes(x=abs(lat),y=Tmax_est_nls),color="red4",fill="red3")+
  geom_point(aes(x=abs(lat),y=Topt_est),color="mediumorchid4",alpha=0.5)+
  geom_smooth(aes(x=abs(lat),y=Topt_est),color="mediumorchid4",fill="mediumorchid3")+
  labs(title= "Thermal traits ~ latitude",x= "latitude",y= "Thermal traits (ºC)")+
  theme_classic()
all_loess_combined  

all_lms_combined <- ggplot(thermal_traits_data)+
  geom_point(aes(x=abs(lat),y=Tmin_est_nls),color="skyblue4",alpha=0.5)+
  geom_smooth(aes(x=abs(lat),y=Tmin_est_nls),
              method="lm",
              color="skyblue4",
              fill="skyblue2")+
  geom_point(aes(x=abs(lat),y=Tmax_est_nls),color="red4",alpha=0.5)+
  geom_smooth(aes(x=abs(lat),y=Tmax_est_nls),
              method="lm",
              color="red4",
              fill="red3")+
  geom_point(aes(x=abs(lat),y=Topt_est),color="mediumorchid4",alpha=0.5)+
  geom_smooth(aes(x=abs(lat),y=Topt_est),
              method="lm",
              color="mediumorchid4",
              fill="mediumorchid3")+
  labs(title= "Thermal traits ~ latitude",x= "latitude",y= "Thermal traits (ºC)")+
  theme_classic()

all_lms_combined

#by order
Tmin_to_lat_plot_order <- ggplot(thermal_traits_data,aes(abs(lat),Tmin_est_nls))+
  geom_point(alpha=0.5,aes(color=order))+
  geom_smooth(method="lm",aes(color=order,fill=order))+
  facet_wrap(.~order)+
  labs(title= "Tmin ~ latitude",
       x= "latitude",
       y= "Tmin")+
  theme_light()+
  theme(legend.position = "none")
Tmin_to_lat_plot_order

Tmax_to_lat_plot_order <- ggplot(thermal_traits_data,aes(abs(lat),Tmax_est_nls))+
  geom_point(alpha=0.5,aes(color=order))+
  geom_smooth(method="lm",aes(color=order,fill=order))+
  facet_wrap(.~order)+
  labs(title= "Tmax ~ latitude",
       x= "latitude",
       y= "Tmax")+
  theme_light()+
  theme(legend.position = "none")
Tmax_to_lat_plot_order

Topt_to_lat_plot_order <- ggplot(thermal_traits_data,aes(abs(lat),Topt_est))+
  geom_point(alpha=0.5,aes(color=order))+
  geom_smooth(method="lm",aes(color=order,fill=order))+
  facet_wrap(.~order)+
  labs(title= "Topt ~ latitude",
       x= "latitude",
       y= "Topt")+
  theme_light()+
  theme(legend.position = "none")
Topt_to_lat_plot_order
#_ _ 4.2 Thermal traits across taxa  ----
#prepare dataframe
counts <- thermal_traits_data %>%
  count(order)%>%
  select(n)%>%
  as_vector()

acari <- thermal_traits_data %>%
  filter(order == "Acari>Prostigmata" | 
           order == "Acari>Trombidiformes")%>%
  mutate(order = "Acari")%>%
  select(Tmin_est_nls,Tmax_est_nls,Topt_est,Tmin_se_nls,Tmax_se_nls,Topt_se_delta,
         order)

mean_by_order <- thermal_traits_data %>%
  select(Tmin_est_nls,Tmax_est_nls,Topt_est,Tmin_se_nls,Tmax_se_nls,Topt_se_delta,
         order) %>%
  filter(order != "Acari>Prostigmata" &
           order != "Acari>Trombidiformes")%>%
  bind_rows(acari)%>%
  group_by(order) %>%
  summarise(across(c(Tmin_est_nls,Tmax_est_nls,
                     Topt_est,Tmin_se_nls,Tmax_se_nls,
                     Topt_se_delta), ~ mean(.x, na.rm = TRUE)))

## first we need subsets for each trait to plot pointrange correctly
#Tmax
lower_tmax <- tapply(thermal_traits_data$Tmax_est_nls,
                    thermal_traits_data$order, min)
upper_tmax <- tapply(thermal_traits_data$Tmax_est_nls,
                     thermal_traits_data$order, max)
tmax_ranges <- tibble(order = mean_by_order$order,
                      Tmax =mean_by_order$Tmax_est_nls,
                      n=counts,
                      lower=lower_tmax,
                      upper=upper_tmax)

#Tmin
lower_tmin <- tapply(thermal_traits_data$Tmin_est_nls,
                     thermal_traits_data$order, min)
upper_tmin <- tapply(thermal_traits_data$Tmin_est_nls,
                     thermal_traits_data$order, max)
tmin_ranges <- tibble(order = mean_by_order$order,
                      Tmin =mean_by_order$Tmin_est_nls,
                      n=counts,
                      lower=lower_tmin,
                      upper=upper_tmin)
#Topt
lower_topt <- tapply(thermal_traits_data$Topt_est,
                     thermal_traits_data$order, min)
upper_topt <- tapply(thermal_traits_data$Topt_est,
                     thermal_traits_data$order, max)
topt_ranges <- tibble(order = mean_by_order$order,
                      Topt =mean_by_order$Topt_est,
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
       y = "Thermal traits (ºC)")+
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
       y = "Thermal traits (ºC)")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))+
  theme(legend.position = "none")

boxplot_traits_ord

#pointranges with mean se
tmax_se_avg <- tapply(thermal_traits_data$Tmax_se_nls,
                      thermal_traits_data$order,
                      mean)
tmax_se_pointrange <- tibble(order = mean_by_order$order,
                        Tmax =mean_by_order$Tmax_est_nls,
                        n=counts,
                        tmax_se_avg)
tmin_se_avg <- tapply(thermal_traits_data$Tmin_se_nls,
                      thermal_traits_data$order,
                      mean)
tmin_se_pointrange <- tibble(order = mean_by_order$order,
                             Tmin =mean_by_order$Tmin_est_nls,
                             n=counts,
                             tmin_se_avg)
topt_se_avg <- tapply(thermal_traits_data$Topt_se_delta,
                      thermal_traits_data$order,
                      mean)
topt_se_pointrange <- tibble(order = mean_by_order$order,
                             Topt =mean_by_order$Topt_est,
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
       y = "Thermal traits (ºC)")+
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
