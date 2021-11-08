#  
#     Script: Meta-analysis using metafor package
#     Aim: estimate response of pests (intrinsic rate of increase) to temperature
#          across categorical variables (taxa, feeding guilds, host plant, etc)
#    This is the second script for this purpose. Here, we will use intrinsic rates directly as
#    effect sizes,since sd is available for each observation in the subset.

#### 1. Load packages ####
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(metafor)
library(fitdistrplus)
library(performance)
library(cowplot)
library(purrr)
#### 2. Explore dataset ####
####.... a) Data preparation ####

IR_data_metafor_raw<-read_csv("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/IR_4meta4.csv")
acari <- IR_data_metafor_raw %>%
  filter(order == "Acari>Prostigmata" | 
           order == "Acari>Trombidiformes")%>%
  mutate(order = "Acari")
IR_data_metafor_raw2 <- IR_data_metafor_raw %>%
  filter(order != "Acari>Prostigmata" &
           order != "Acari>Trombidiformes")%>%
  bind_rows(acari)

miners <- IR_data_metafor_raw %>%
  filter(feeding_guild == "leafminer" | 
           order == "miner")%>%
  mutate(order = "miners")
IR_data_metafor_raw2 <- IR_data_metafor_raw2 %>%
  filter(feeding_guild != "leafminer" | 
           order != "miner") %>%
  bind_rows(miners)


# insert id
IR_id <- IR_data_metafor_raw2 %>%
  distinct(`Article Title`)%>%
  mutate(id=row_number())
IR_data_metafor <- inner_join(IR_data_metafor_raw2,IR_id,by="Article Title")
#select interest variables
eff_size_prep <- IR_data_metafor %>%
  dplyr::select(id,ir_mean,ir_se,ir_n,temp,order,feeding_guild,lon,lat,Authors,Year)%>%
  mutate(vi=(ir_se*sqrt(ir_n))^2)%>% #transform se into variance
  mutate(Authors= word(Authors,1,sep = ","))%>% #use only first author to simplify
  glimpse()
####.... b) Data exploration ####
## effect size distribution
ir <- eff_size_prep$ir_mean
hist(ir) #not normal
hist(log(ir))
boxplot(log(ir))
boxplot(eff_size_prep$vi)
summary(eff_size_prep$vi) #hay que eliminar outliers
plot(eff_size_prep$vi)

eff_size_clean <- eff_size_prep %>%
  filter(vi<=.01)


#### 3. Intrinsic Rate as effect size ####
# We can use intrinsic rates of increase as effet size since it is already standardized
# with requirement of normality assumptions. The following analysis is performed without
# any transformation, but a log-transform for this approach might be more appropriate

#### .. 3.1. Mixed-Effects meta-regression ####
# simply not included because it accounts for temperature rather than for across-study
#   variation.

#### .. 3.2. MultiVariate RE-meta-regression temp####
re_mv_temp <- rma.mv(yi = ir_mean,
                     V = vi, #not sure if V is the sampling variance for each eff-size
                     random = ~1|id,
                     mods = ~temp,
                     data = eff_size_clean)
summary(re_mv_temp)
# temp is significant
forest(re_mv_temp) #not useful
regplot(re_mv_temp,
        xlim=c(5,50),
        xlab="Temperature (ºC)",lcol = "darkred",
        col = alpha(colour = "slategray",.3),shade = "lightcoral",
        pch = 16,alpha=.3)
funnel(re_mv_temp)

#### .. 3.3. MultiVariate Mixed-Efects order####
me_mv_order <- rma.mv(yi = ir_mean,
                     V = vi, #not sure if V is the sampling variance for each eff-size
                     random = ~1|id,
                     mods = ~order,
                     data = eff_size_clean)
summary(me_mv_order)
# non-significant Qbet (p=0.0636); significant between Lepidoptera and chewer
# AIC = 29589
# BIC = 29618
forest(me_mv_order) #not useful
funnel(me_mv_order)

#### .. 3.4. Mix.Eff. MV metareg temp-order####
#### ...... a) temp+order####

me_mv_temp_add_ord <- rma.mv(yi = ir_mean,
                         V = vi,
                         random = ~1|id,
                         mods = ~temp+order,
                         data = eff_size_clean)
summary(me_mv_temp_add_ord)
# temp is significant
# differences of intrinsic rates between Acari and Coleoptera and Acari and Lepidoptera (more numerous)
# AIC = 28844.9954   
# BIC = 28877.6440   
forest(me_mv_temp_add_ord) #not useful
funnel(me_mv_temp_add_ord)

#### ...... b) temp*order####

me_mv_temp_int_ord <- rma.mv(yi = ir_mean,
                             V = vi,
                             random = ~1|id,
                             mods = ~temp*order,
                             data = eff_size_clean)
summary(me_mv_temp_int_ord)
# temp is significant
# differences between main groups (Acari, Hemiptera, Lepidoptera)
# temperature affects each order differently
# much better than for previous:
# AIC = 4800.5094      
# BIC = 4854.5964    
forest(me_mv_temp_int_ord) #not useful
funnel(me_mv_temp_int_ord)

#### .. 3.5. Mix.Eff. MV metareg temp-feeding guild####
#### ...... a) temp+feeding_guild####
me_mv_temp_add_fg <- rma.mv(yi = ir_mean,
                             V = vi,
                             random = ~1|id,
                             mods = ~temp+feeding_guild,
                             data = eff_size_clean)
summary(me_mv_temp_add_fg)
# temp is significant
# no differences among feeding guilds (statistically significant)
# AIC = 28840.6308      
# BIC = 28866.0743         
forest(me_mv_temp_add_fg) #not useful
funnel(me_mv_temp_add_fg)

#### ...... b) temp*feeding_guild####

me_mv_temp_int_fg <- rma.mv(yi = ir_mean,
                             V = vi,
                             random = ~1|id,
                             mods = ~temp*feeding_guild,
                             data = eff_size_clean)
summary(me_mv_temp_int_fg)
# temp is significant
# differences between main groups (Acari, Hemiptera, Lepidoptera)
# temperature affects each feeding guild differently
# similar outcome:
# AIC = 28345.1073         
# BIC = 28384.9317       
forest(me_mv_temp_int_fg) #not useful
funnel(me_mv_temp_int_fg)

#### .. 3.6. Mix.Eff. MV order x feeding guild####
me_mv_order_int_fg <- rma.mv(yi = ir_mean,
                            V = vi,
                            random = ~1|id,
                            mods = ~order*feeding_guild,
                            data = eff_size_clean)
summary(me_mv_order_int_fg) #not working
forest(me_mv_order_int_fg) #not useful
funnel(me_mv_order_int_fg)

#### .. 3.7. Aggregate by study ####

# first we need to compute an unique study-level effet size.
# for that aim, we compute the index:
#          estimate = mean(ir_mean)/mean(temp) for each study
#          variance = mean(vi) for each study

#### ...... prev) compute summary effect size by study ####

agg_eff_sizes_var<- eff_size_clean %>%
  group_by(id, order, feeding_guild,Authors,Year,lat)%>%
  summarise(vi = mean(vi))

agg_eff_sizes_estimate<- eff_size_clean %>%
  group_by(id, order, feeding_guild,Authors,Year,lat)%>%
  summarise(ir_mean = mean(ir_mean)/mean(temp))


agg_eff_sizes<- tibble(agg_eff_sizes_estimate$id,
                       agg_eff_sizes_estimate$order,
                       agg_eff_sizes_estimate$feeding_guild,
                       agg_eff_sizes_estimate$Authors,
                       agg_eff_sizes_estimate$Year,
                       agg_eff_sizes_estimate$lat,
                       agg_eff_sizes_estimate$ir_mean,
                       agg_eff_sizes_var$vi)
colnames(agg_eff_sizes) <- c(colnames(agg_eff_sizes_estimate),"vi")
agg_eff_sizes



length(agg_eff_sizes_estimate$ir_mean)


#### .. 3.8. Subgsetting for forest plots ####
#### ...... a) order ####
# lepidoptera:
re_mv_lepi <- rma.mv(yi = ir_mean,
                      V = vi, #not sure if V is the sampling variance for each eff-size
                      subset=(order=="Lepidoptera"),
                      random = ~1|id,
                      mods = ~temp,
                      data = eff_size_clean)
sum_lepi <- summary(re_mv_lepi)
coefs_lepi <- as.data.frame(sum_lepi$beta)
ci_lower_lepi <- as.data.frame(sum_lepi$ci.lb)
ci_upper_lepi <- as.data.frame(sum_lepi$ci.ub)
# acari:
re_mv_acari <- rma.mv(yi = ir_mean,
                     V = vi, #not sure if V is the sampling variance for each eff-size
                     subset=(order=="Acari"),
                     random = ~1|id,
                     mods = ~temp,
                     data = eff_size_clean)
sum_acari <- summary(re_mv_acari)
coefs_acari <- as.data.frame(sum_acari$beta)
ci_lower_acari <- as.data.frame(sum_acari$ci.lb)
ci_upper_acari <- as.data.frame(sum_acari$ci.ub)
# hemiptera:
re_mv_hemiptera <- rma.mv(yi = ir_mean,
                      V = vi, #not sure if V is the sampling variance for each eff-size
                      subset=(order=="Hemiptera"),
                      mods = ~temp,
                      random = ~1|id,
                      data = eff_size_clean)
sum_hemip <- summary(re_mv_hemiptera)
coefs_hemip <- as.data.frame(sum_hemip$beta)
ci_lower_hemip <- as.data.frame(sum_hemip$ci.lb)
ci_upper_hemip <- as.data.frame(sum_hemip$ci.ub)
# collect model estimates and corresponding variances
estimates <- c(coefs_lepi[2,],coefs_hemip[2,],coefs_acari[2,])
ci_lower <- c(ci_lower_lepi[2,], ci_lower_hemip[2,], ci_lower_acari[2,])
ci_upper <- c(ci_upper_lepi[2,], ci_upper_hemip[2,], ci_upper_acari[2,])
### create vector with labels
labels <- c("Lepidoptera", "Hemiptera", "Acari")

visual_re_mv_subset <- tibble(labels,estimates,ci_lower,ci_upper)
# see sample size
eff_size_clean %>% 
  filter(!is.na(order)) %>% 
  group_by(order) %>% 
  count()

subsets_plot_metareg <- ggplot(data=visual_re_mv_subset,aes(x=estimates,y=labels))+
  geom_point(aes(color=labels),size=3)+
  geom_pointrange(aes(y=labels,
                      x=estimates,
                      xmin=ci_lower,
                      xmax=ci_upper,
                      color=labels),
                  size=1)+
  labs(x="Slopes", 
       y= "Order",
       title = "Multivariate Random-Effects",
       subtitle = "Moderator = Temperature (meta-regression)")+
  theme_classic()+
  theme(legend.position = "none")+
  annotate(geom= "text", x = .0068, y = "Lepidoptera", 
           label = "n = 74", color = "slategrey",size=3.15)+
  annotate(geom= "text", x = .007, y = "Hemiptera", 
           label = "n = 78", color = "slategrey",size=3.15)+
  annotate(geom= "text", x = .0109, y = "Acari", 
           label = "n = 98", color = "slategrey",size=3.15)

subsets_plot_metareg

#### ...... b) feeding_guild ####
# suckers:
re_mv_sucker <- rma.mv(yi = ir_mean,
                     V = vi, #not sure if V is the sampling variance for each eff-size
                     subset=(feeding_guild=="sucker"),
                     random = ~1|id,
                     mods = ~temp,
                     data = eff_size_clean)
sum_sucker <- summary(re_mv_sucker)
coefs_sucker <- as.data.frame(sum_sucker$beta)
ci_lower_sucker <- as.data.frame(sum_sucker$ci.lb)
ci_upper_sucker <- as.data.frame(sum_sucker$ci.ub)
# chewers:
re_mv_chewer <- rma.mv(yi = ir_mean,
                      V = vi, #not sure if V is the sampling variance for each eff-size
                      subset=(feeding_guild=="chewer"),
                      random = ~1|id,
                      mods = ~temp,
                      data = eff_size_clean)
sum_chewer <- summary(re_mv_chewer)
coefs_chewer <- as.data.frame(sum_chewer$beta)
ci_lower_chewer <- as.data.frame(sum_chewer$ci.lb)
ci_upper_chewer <- as.data.frame(sum_chewer$ci.ub)
# miners:
re_mv_miner <- rma.mv(yi = ir_mean,
                          V = vi, #not sure if V is the sampling variance for each eff-size
                          subset=(feeding_guild=="miner"),
                          mods = ~temp,
                          random = ~1|id,
                          data = eff_size_clean)
sum_miner <- summary(re_mv_miner)
coefs_miner <- as.data.frame(sum_miner$beta)
ci_lower_miner <- as.data.frame(sum_miner$ci.lb)
ci_upper_miner <- as.data.frame(sum_miner$ci.ub)

#suckers:
re_mv_borer <- rma.mv(yi = ir_mean,
                       V = vi, #not sure if V is the sampling variance for each eff-size
                       subset=(feeding_guild=="borer"),
                       random = ~1|id,
                       mods = ~temp,
                       data = eff_size_clean)
sum_borer <- summary(re_mv_borer)
coefs_borer <- as.data.frame(sum_borer$beta)
ci_lower_borer <- as.data.frame(sum_borer$ci.lb)
ci_upper_borer <- as.data.frame(sum_borer$ci.ub)

# collect model estimates and corresponding variances
estimates <- c(coefs_sucker[2,],coefs_chewer[2,],coefs_miner[2,],coefs_borer[2,])
ci_lower <- c(ci_lower_sucker[2,], ci_lower_chewer[2,], ci_lower_miner[2,],ci_lower_borer[2,])
ci_upper <- c(ci_upper_sucker[2,], ci_upper_chewer[2,], ci_upper_miner[2,],ci_upper_borer[2,])
### create vector with labels
labels <- c("Suckers", "Chewers", "Miners", "Borers")

visual_re_mv_subset_fg <- tibble(labels,estimates,ci_lower,ci_upper)
# see sample size
eff_size_clean %>% 
  filter(!is.na(feeding_guild)) %>% 
  group_by(feeding_guild) %>% 
  count()


subsets_plot_metareg_fg <- ggplot(data=visual_re_mv_subset_fg,aes(x=estimates,y=labels))+
  geom_point(aes(color=labels),size=3)+
  geom_pointrange(aes(y=labels,
                      x=estimates,
                      xmin=ci_lower,
                      xmax=ci_upper,
                      color=labels),
                  size=1)+
  labs(x="Slopes", 
       y= "Feeding Guild",
       title = "Multivariate Random-Effects",
       subtitle = "Moderator = Temperature (meta-regression)")+
  theme_classic()+
  theme(legend.position = "none")+
  annotate(geom= "text", x = .0025, y = "Suckers", 
           label = "n = 186 ", color = "slategrey",size=3.15)+
  annotate(geom= "text", x = .008, y = "Miners", 
           label = "n = 25 ", color = "slategrey",size=3.15)+
  annotate(geom= "text", x = .007, y = "Chewers", 
           label = "n = 44 ", color = "slategrey",size=3.15)+
  annotate(geom= "text", x = .009, y = "Borers", 
           label = "n = 41 ", color = "slategrey",size=3.15)

subsets_plot_metareg_fg

subsets_plot_metareg_fg <- ggplot(data=visual_re_mv_subset_fg,aes(x=estimates,y=labels))+
  geom_point(aes(color=labels),size=3)+
  geom_pointrange(aes(y=labels,
                      x=estimates,
                      xmin=ci_lower,
                      xmax=ci_upper,
                      color=labels),
                  size=1)+
  labs(x="Slopes", 
       y= "Feeding Guild",
       title = "Multivariate Random-Effects",
       subtitle = "Moderator = Temperature (meta-regression)")+
  theme_classic()+
  theme(legend.position = "none")

subsets_plot_metareg_grid <-ggplot(data=visual_re_mv_subset_fg,aes(x=estimates,y=labels))+
  geom_point(aes(color=labels),size=3)+
  geom_pointrange(aes(y=labels,
                      x=estimates,
                      xmin=ci_lower,
                      xmax=ci_upper,
                      color=labels),
                  size=1)+
  labs(x="Slopes", 
       y= "Feeding Guild",
       title = "Multivariate Random-Effects",
       subtitle = "Moderator = Temperature (meta-regression)")+
  theme_classic()+
  theme(legend.position = "none")


rma

#### 4. Thermal Traits inference ####
#### .... 4.1.- Prepare data ####
#load data
thermal_traits_data_raw <- read_csv("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/dataset_params_metaanalysis.csv") %>%
  dplyr::select(-7) %>%
  glimpse()

## apply filters 
acari_raw <- thermal_traits_data_raw %>% #ensemble all acari
  filter(order == "Acari>Prostigmata" | 
           order == "Acari>Trombidiformes")%>%
  mutate(order = "Acari")

thermal_traits_data <- thermal_traits_data_raw %>%
  filter(order != "Acari>Prostigmata" &
           order != "Acari>Trombidiformes")%>%
  bind_rows(acari_raw)%>%
  filter(Tmin_est_nls >0 &
           Tmin_se_nls <10 &
           Tmax_se_nls <10 &
           Tmax_est_nls < 45) %>%  #all the authors whose papers report variance
  filter(Authors == "de Campos" |
           Authors == "Li" |
           Authors == "Schlesener" |
           Authors == "del Pino" |
           Authors == "Tian" |
           Authors == "Islam" |
           Authors == "Govindan" |
           Authors == "Choudhary" |
           Authors == "Saeed" |
           Authors == "Chiu" |
           Authors == "Lacerda" |
           Authors == "Barbosa" |
           Authors == "Ahn" |
           Authors == "Yazdanpanah" |
           Authors == "Soh" |
           Authors == "Qin" |
           Authors == "Sanchez-Ramos" |
           Authors == "Liao" |
           Authors == "Geng" |
           Authors == "Ngowi" |
           Authors == "Guo" |
           Authors == "Jiang" |
           Authors == "Chi" |
           Authors == "Basirat" |
           Authors == "Chen" |
           Authors == "Khadioli" |
           Authors == "Martinez" |
           Authors == "Pilkington" |
           Authors == "Fand" |
           Authors == "Kumar" |
           Authors == "Golizadeh" |
           Authors == "Amiri" |
           Authors == "Zahiri" |
           Authors == "Aysal" |
           Authors == "Kuo" |
           Authors == "Toapanta" |
           Authors == "Premachandra" |
           Authors == "Satar" |
           Authors == "Hoddle" |
           Authors == "Turak" |
           Authors == "Dreyer" |
           Authors == "Lin" |
           Authors == "Ullah" |
           Authors == "Hasanvand" |
           Authors == "Bahirai" |
           Authors == "Fidelis" |
           Authors == "Xie" |
           Authors == "Gotoh" |
           Authors == "Vangansbeke" |
           Authors == "Karami-Jamour" |
           Authors == "Kasap"
         ) %>%
  glimpse() #subset from those which have an standard error in primary studies for sample size

#define variables
Tmin_est <- thermal_traits_data$Tmin_est_nls
Tmin_se <- thermal_traits_data$Tmin_se_nls
Tmax_est <- thermal_traits_data$Tmax_est_nls
Tmax_se <- thermal_traits_data$Tmax_se_nls
Topt_est <- thermal_traits_data$Topt_est
Topt_se <- thermal_traits_data$Topt_se_delta
lat <- thermal_traits_data$lat
lon <- thermal_traits_data$lon

#### .... 4.2.- Check model assumptions ####
#normality
#### ........ a) Tmin ####
summary(Tmin_est)
hist(Tmin_est) #quite normal
shapiro.test(Tmin_est) #normal
#### ........ b) Tmax ####
summary(Tmax_est)
hist(Tmax_est)   #quite normal
shapiro.test(Tmax_est) #normal
#### ........ c) lon, lat ####
summary(lat)
hist(lat)    
#create map along coordinates
map <- ggplot(data = thermal_traits_data, aes(x = lon, y = lat)) +
  borders("world", colour = "gray90", fill = "gray90") +
  geom_point(aes(colour = order), alpha = 0.35, size = 3)+theme_classic()+
  theme(legend.position = "bottom")
map
#### ........ d) Topt ####
summary(Topt_est)
hist(Topt_est) #quite normal
shapiro.test(Topt_est) #normal

## CONCLUSION: we can use thermal traits as an effect size metric
## nota: falta incluir en el dataset los valores promedio de varianzas
##       para pesar los estudios en este tipo de metaanálisis

#### .... 4.3. Analyses ####
#### ........a) traits ~ lat ####
## Tmin
Tmin_lat <- lm(Tmin_est_nls~abs(lat), data = thermal_traits_data)
performance::check_model(Tmin_lat) #quite good
Tmin_lat_sum <- summary(Tmin_lat)
Tmin_lat_sum
## Tmax
Tmax_lat <- lm(Tmax_est_nls~abs(lat), data = thermal_traits_data)
performance::check_model(Tmax_lat) #quite good
Tmax_lat_sum <- summary(Tmax_lat)
Tmax_lat_sum
## Topt
Topt_lat <- lm(Topt_est~abs(lat), data = thermal_traits_data)
check_model(Topt_lat)
Topt_lat_sum <- summary(Topt_lat)
Topt_lat_sum
## Tmax ~ Tmin
Tmax_Tmin <- lm(Tmax_est_nls~Tmin_est_nls, data = thermal_traits_data)
performance::check_model(Tmax_Tmin)
Tmax_Tmin_sum <- summary(Tmax_Tmin)
Tmax_Tmin_sum

#### ........b) plots ####

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


#### ............. Traits ~ order ####
# Tmax
tmax_order <- lm(Tmax_est_nls~order, data = thermal_traits_data)
check_model(tmax_order) #quite good
tmax_order <- summary(tmax_order)

# Tmin
tmin_order <- lm(Tmin_est_nls~order, data = thermal_traits_data)
check_model(tmin_order) #quite good
sum_tmin_order <- summary(tmin_order)


#### ............. Traits ~ feeding guild ####
# Tmax
tmax_fg <- lm(Tmax_est_nls~feeding_guild, data = thermal_traits_data)
sum_tmax_fg <- summary(tmax_fg)
anova_tmax_fg <- aov(tmax_fg) # sí hay diferencias
TukeyHSD(anova_tmax_fg)
# Tmin
tmin_fg <- lm(Tmin_est_nls~feeding_guild, data = thermal_traits_data)
plot(tmin_fg)
sum_tmin_fg <- summary(tmin_fg)
anova(tmin_fg) # no hay
anova_tmin_fg <- aov(tmin_fg) 
TukeyHSD(anova_tmin_fg)

#### ........ c) variance across traits ####

#prepare dataframe
counts <- thermal_traits_data %>%
  count(order)%>%
  dplyr::select(n)%>%
  as_vector()
mean_by_order <- thermal_traits_data %>%
  dplyr::select(Tmin_est_nls,Tmax_est_nls,Topt_est,Tmin_se_nls,Tmax_se_nls,Topt_se_delta,
         order) %>%
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
       subtitle= "Circles' size are proportional to sample size.
Error bars represent the range of estimates",
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
       subtitle= "Error bars represent average SE for parameter estimate within each group. 
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









#### EXTRA: metafor examples ####
dat.konstantopoulos2011
#### .... 1. mixed-effects ####
res1 <- rma(dat.konstantopoulos2011,yi,vi,
    mods =~school)

forest(res1)
funnel(res1)

#### .... 2. mv ####
dat.hasselblad1998
res1 <- rma.mv(dat.hasselblad1998,xi,vi,
            mods =~school,
            random = ~1|study)

forest(res1)
funnel(res1)








#### .. b) plot regression model####
regplot(ME_metareg_temp,
        xlim=c(5,50),
        xlab="Temperature (ºC)",lcol = "darkred",
        col = alpha(colour = "slategray",.3),shade = "lightcoral",
        pch = 16,alpha=.3)