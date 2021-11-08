#### SCRIPT INFO ####     
#     Authors: Darío San Segundo Molina, Sara Villén Pérez, Ignacio Morales Castilla
#     Title: Exploratory analysis from INTRAPEST database 
#     Aim: explore first relationships in the database, as well as literature biases and trends
#     Date: October 2021
#_________________________ ####
 
#### 1. Data loading and preparation ####
rm(list=ls())
library(ggplot2)
library(readr)
library(dplyr)
library(wesanderson)
library(svglite)
library(performance)
library(qqplotr)
library(moments)

# After first filter, now we have different categories
filtered2 <- read_csv("/Users/Ecologia/Documents/USUARIOS/DARÍO/filtered_review.csv", 
                      col_names=TRUE) %>%
  filter(Filter == "yes")%>%
  glimpse()
#See counts of combinations of categories:
nn <-filtered2 %>%
  count(Approach, Subapproach,Lab_response)
# See grouping categories:
filtered2 %>%
  count(Approach)%>%
  View()
filtered2 %>%
  count(Approach, Subapproach)%>%
  View()
# Extract  a csv with each filtered dataset
intrinsic_growth <- filtered2 %>%
  filter(Approach == "experimental" &
           Subapproach == "laboratory" &
           Lab_response == "rm" )%>%
  write_csv(path = "/Users/Ecologia/Documents/USUARIOS/DARÍO/intrinsic_growth2.csv")

#### 2. Variable counting ####
intrinsic_classified <- read_csv("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/intrinsic_growth2.csv",
                                 col_names = TRUE) %>%
  filter(Filter_3 == "yes" & n==1) # pass filter
# 1s are for first row of each study)
#by order
intrinsic_classified %>%
  count(order)%>%
  View()
#by family
intrinsic_classified %>%
  count(family)%>%
  View()
#by feeding_guilds
intrinsic_classified %>%
  count(feeding_guild)%>%
  View()
#by host plant family
intrinsic_classified %>%
  count(h_p_family)%>%
  View()
#by diet class
intrinsic_classified %>%
  count(diet_family)%>%
  View()
#by both order and family
intrinsic_classified %>%
  count(order,family)%>%
  View()

#by both family and feeding guild
intrinsic_classified %>%
  count(family,feeding_guild)%>%
  View()
intrinsic_classified %>%
  count(order,feeding_guild)%>%
  View()

#### 3. Exploratory INTRAPEST database analysis ####
#### _ _ _ 3.1. Tibble ensembling ####
IR_data <-read_csv("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/intrinsic_growth2_upd.csv") %>%
  filter(Filter_3 == "yes") %>% #select only those which have passed the filter 3 in the source csv
  mutate(growth_rate = replace(growth_rate, growth_rate <= 0,0))%>% #transform negative values into 0
  filter(as.numeric(temperature)<50)%>% # exclude possible mistranscription in the dataset ensambling process
  as_tibble()%>% #convert to tibble object
  glimpse() # similar to str()
#since some numerical variables have been read as characters, we change them to dbl
IR_data$temperature <- as.numeric(IR_data$temperature)
IR_data$growth_rate <- as.numeric(IR_data$growth_rate)
IR_data$error <- as.numeric(IR_data$error)
IR_data$RH <- as.numeric(IR_data$RH)
IR_data$lat <- as.numeric(IR_data$lat)
IR_data$n_1 <- as.numeric(IR_data$n_1)
#now we count how many points by paper we have
IR_data %>% 
  count(`Article Title`)%>%
  View()
#### _ _ _ 3.2. Subset: n and  s^2 ####
IR_suitable4metaanalysis <- IR_data %>% 
  filter(error != "no" &  #exclude lines without error and/or sample
         error != "NA" &
          n_1 != "no" &
           n_1 != "NA")%>%
  mutate(title =`Article Title`)%>% #create a column without spacing to facilitate requiring
  select(-`Article Title`)
#Group all acari into one category
  acari <- IR_data %>%
  filter(order == "Acari>Prostigmata" | #select both categories of acari
           order == "Acari>Trombidiformes")%>%
  mutate(order = "Acari") #change the content into just "Acari"
IR_data <-IR_data %>%
  filter(order != "Acari>Prostigmata" & #exclude acari rows
           order != "Acari>Trombidiformes")%>% #exclude acari rows
  bind_rows(acari) #bind "Acari" rows

#select all columns except original Article Title
#now we write the csv which is suitable to posterior meta-analysis techniques
write_csv(IR_suitable4metaanalysis,"IR_metaanalysis_suitable.csv")

#which is the resulting sample size after subsetting?
counts4meta <- IR_suitable4metaanalysis %>%
  count(title) #see how many rows for each title
length(unique(counts4meta$title)) #see number of different titles
#### _ _ _ 3.3. Check model assumptions ####
# linearity,normality,regressions
#### _ _ _ _ _ _ a) histograms ####
hist(IR_data$temperature)
hist(IR_data$growth_rate) #poisson with lambda = 1 ?
shapiro.test(IR_data$temperature) #temperature not-normal? by Saphiro test
shapiro.test(IR_data$growth_rate) #r not-normal? by Saphiro test
ks.test(IR_data$temperature,"pnorm")#temperature not-normal? by Kolmogorov-Smirnov
ks.test(IR_data$growth_rate,"pnorm") #r not normal by Kolmogorov-Smirnov
#### _ _ _ _ _ _ b) skewness and kurtosis ####
#skewness and kurtosis, package moments
skewness(IR_data$temperature,na.rm=TRUE) # just slightly left-skewed
skewness (IR_data$growth_rate,na.rm=TRUE) #highly right-skewed
kurtosis(IR_data$temperature,na.rm=TRUE) #a bit flattened
kurtosis(IR_data$growth_rate,na.rm=TRUE) #strongly peaked
#### _ _ _ _ _ _ c) check_model function ####
check_model(lm1) # package performance

#### _ _ _ _ _ _ d) Transformations ####
# since intrinsic rate of increase in our dataset is not normally distributed, 
# one option to work with it is to use a logarithmic transformation, according to
# Koricheva et al. (2013), and given that r_m is log-normally distributed (Dillingham
# et al. 2016)
hist(IR_data$growth_rate)
hist(sqrt(IR_data$growth_rate))
hist(log(IR_data$growth_rate))
hist(IR_data$growth_rate^(1/3))
#maybe a sqrt is a better option
shapiro.test(sqrt(IR_data$growth_rate)) #not normal
shapiro.test(IR_data$growth_rate^(1/3)) #not normal
shapiro.test(log(IR_data$growth_rate)) #NA
ks.test(log_rm_data$rm,"pnorm") 

#homocedasticity
fligner.test(growth_rate~temperature,data = IR_data) #no variance homogeneity found

#### _ _ _ 3.4. Scatterplots ####
#### _ _ _ _ _ a) all pooled ####
scatter_IR_grouped <- ggplot(data=IR_data, aes(x=temperature, y=growth_rate))+
  geom_point(alpha=0.2,aes(color=order))+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette="Dark2")
#ggsave(filename = "scatter_rm_all.png",units = "cm" ,height = 20, width = 15)
#### _ _ _ _ _ b) order ####
# one plot
scatter_IR_grouped_se <- ggplot(data=IR_data, aes(x=temperature, y=growth_rate,color=order,fill = order))+
  geom_point(alpha=0.4)+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette="Dark2")
#ggsave(filename = "scatter_rm_all.png",units = "cm" ,height = 20, width = 15)

# grided (group acari first)
scatter_IR_grided_scalesFixed <- ggplot(data=IR_data, aes(x=temperature, y=growth_rate,color=order,fill = order))+
  geom_point(alpha=0.2)+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  facet_wrap(.~order,nrow = 1)+
  theme_bw()+
  theme(legend.position = "bottom")
#ggsave(filename = "scatter_IR_grided.png",units = "cm",height = 12, width = 25)
scatter_IR_grided_scalesfree <- ggplot(data=IR_data, aes(x=temperature, y=growth_rate,color=order,fill = order))+
  geom_point(alpha=0.2)+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  facet_wrap(.~order,nrow = 1,scales = "free")+
  theme_bw()+
  theme(legend.position = "bottom")
## let's see hemiptera by family
hemiptera <- IR_data %>%
  filter(order == "Hemiptera")%>%
  glimpse()
hemiptera_curve_grid<- ggplot(hemiptera, aes(x=temperature,y=growth_rate,color=family,fill=family))+
  geom_point(alpha=0.25)+
  labs(title = "Intrinsic rate of increase values for Hemiptera",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  facet_wrap(.~family,nrow=2)+
  theme_bw()+
  theme(legend.position = "none")
# ggsave(filename = "hemiptera_family_grided.png",units = "cm" ,height = 20, width = 15)
# let's see lepidoptera
lepidoptera <- IR_data %>%
  filter(order=="Lepidoptera")%>%
  glimpse()
lepidoptera_curve_gridfam<- ggplot(lepidoptera, aes(x=temperature,y=growth_rate,color=family,fill=family))+
  geom_point(alpha=0.25)+
  labs(title = "Intrinsic rate of increase values for Lepidoptera",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  facet_wrap(.~family)+
  theme_bw()+
  theme(legend.position = "none")
#ggsave(filename = "lepidoptera_by_fam.png",units = "cm" ,height = 20, width = 15)

#### _ _ _ _ _ c) feeding guild ####
# first we group miners and leafminers into "miners"
mines <- IR_data %>%
  filter(feeding_guild == "leafminer" | 
           feeding_guild == "miner")%>%
  mutate(feeding_guild = "miner")
IR_data <-IR_data %>%
  filter(feeding_guild != "miner" &
           feeding_guild != "leafminer")%>%
  bind_rows(mines)
#subset with interest variables
rm_feeding <- IR_data %>%
  select(temperature,growth_rate,lat,feeding_guild,order)
#feeding_guild
scatter_rmfeed_grided_by_feedinguild <- ggplot(data=rm_feeding, aes(x=temperature, y=growth_rate))+
  geom_point(aes(color=feeding_guild),alpha=0.4)+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth(aes(color=feeding_guild))+
  facet_wrap(.~feeding_guild,nrow = 1,scales = "free")+
  theme_bw()+
  theme(legend.position = "bottom")
#ggsave(filename = "scatter_feeding2.png",units = "cm" ,height = 12, width = 25)

#### _ _ _ _ _ d) latitude ####
scatter_IR_grouped_lat <- ggplot(data=IR_data, aes(x=temperature, y=growth_rate))+
  geom_point(alpha=0.5,size=3,aes(color=abs(lat)),position = position_jitter(w = 2, h = 0))+
     labs(subtitle = "Intrinsic rate of increase values and latitude",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()+
  scale_colour_gradient(low = "lightcoral",high = "cyan4")
#ggsave(filename = "scatter_latitude.png",units = "cm" ,height = 20, width = 15)

#latitude effects by order
scatter_IR_grouped_lat_gridorder <- ggplot(data=IR_data, aes(x=temperature, y=growth_rate))+
  geom_point(alpha=0.3,size=3,aes(color=abs(lat)),position = position_jitter(w = 2, h = 0))+
  labs(subtitle = "Intrinsic rate of increase values and latitude",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  facet_wrap(.~order,nrow = 2)+
  theme_bw()+
  scale_colour_gradient(low = "lightcoral",high = "cyan4")+
  theme(legend.position = "bottom")
#ggsave(filename = "scatter_latitude_gridord",units = "cm" ,height = 20, width = 15)

# let's discretize lat variable into 4 groups to easily visualize trends
IR_data_disLat <- IR_data %>%
  mutate(lat=abs(lat)) %>% #absolute value
  mutate(lat= cut_number(lat,n=4, #make 4 equally distributed categories
                         labels=c("Tropical","warm","mild","cold temperate")))
  #drop_na(lat)
write_csv(IR_data_disLat,file="intrinsic_dislat.csv")
# grid por orden, smooths para categorías de latitud
scatter_grid_disLat2 <- ggplot(data=IR_data_disLat,aes(x=temperature, y=growth_rate))+
  geom_point(alpha=0.3,size=3,aes(color=order),position=position_jitter(w=2,h=0))+
  labs(subtitle = "Intrinsic rate of increase values and latitude",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()+
  facet_wrap(.~lat,nrow = 1,scales = "free")+
  theme(legend.position = "bottom")
#ggsave(filename = "smooths_lat_gridbylat2.png",units = "cm" ,height = 12, width = 25)

#### 4. Possibilities for outliers ####
# Some aphid data have uncommon high r values
IR_data_without_outlierAphids <- IR_data%>%
filter(`UT (Unique WOS ID)` != "WOS:000071638600063"&
         `UT (Unique WOS ID)`!= "WOS:A1995RM24300002" &
         growth_rate <= 0.6)
#Repetir análisis anteriores con el dataset sobreescrito
#### CONCLUSIONS ####
# 1) we conclude that it might be useful to include lat as a covariable.
# 2) data are not NORMAL, so the regression should follow generalized methods
# 3) data are not LINEAR, so we will have to use a non-linear approach (ex. Brière-1)
# 4) data have random variables such as Study (var. "Article Title"). Latitude as a covariate
# 5) data's variances are not homogeneous.
# 6) MODEL: nonlinear meta-regression with mixed-effects and corrections for variance, non-normality

