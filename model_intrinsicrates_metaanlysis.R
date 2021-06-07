#  
#     Script: DATA ANALYSIS FOR SYSTEMATIC REVIEW:
#     Aim: estimate response of pests (intrinsic rate of increase) to temperature
#          across categorical variables (taxa, feeding guilds, host plant, etc)


#### 1. Load Dataset ####
rm(list=ls())
#library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(svglite)
#load data:
IR_data <-read_csv("/Users/Ecologia/Documents/USUARIOS/DARÍO/intrinsic_growth2.csv") %>%
  filter(Filter_3 == "yes") %>% #select only those which have passed the filter 3 in the source csv
  mutate(growth_rate = replace(growth_rate, growth_rate <= 0,0))%>% #we assign 0 to negative values of intrinsic rates (biological nonsense)
  filter(as.numeric(temperature)<50) %>% #exclude possible missleading points
  select(Authors,`Article Title`,Filter_3,Approach,Subapproach,temperature,growth_rate,
         error,n_1,order,family,genus,species,feeding_guild,h_p_family,
         diet_family,RH,daylength,lat,lon)%>% #subset variables
glimpse() # similar as str()
#since some numerical variables are read as characters, we change them to dbl
IR_data$temperature <- as.numeric(IR_data$temperature)
IR_data$growth_rate <- as.numeric(IR_data$growth_rate)
IR_data$error <- as.numeric(IR_data$n_1)
IR_data$RH <- as.numeric(IR_data$RH)
IR_data$lat <- as.numeric(IR_data$lat)

#### 2. Define fitting function ####
#load training artificial data
train_data <- read.table("train.csv",sep=";",dec=",",header=TRUE) #inventados siguiendo la curva aprox.
colnames(train_data) <- c("temp","r")
train_data$temp <- as.numeric(train_data$temp)
train_data$r <- as.numeric(train_data$r)
View(train_data)
ggplot(train_data,aes(x=temp,y=r))+
  geom_point()+
  stat_smooth(method= "loess")+
  theme_classic()

# We use a Brière-1 model (Briere et al., 1999)
briere1 <- function(a, temp, Tmin, Tmax){
  a*temp*(temp-Tmin)*(Tmax-temp)^(1/2)
}


#### 3. Train data ####

fitted_br1 <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                   data = train_data,
                   start = list(a = 0.0002,
                                Tmin =7.5,
                                Tmax= 33.5))
# problemas de ajuste. Ver starting con plots sucesivos para escoger el 
# orden de magnitud del parámetro a
par(mfrow=c(2,2))
temp <- train_data$temp
growth <- train_data$r
plot(temp,growth)
plot(temp,briere1(a=0.0002,Tmin=7.5,Tmax=33.5,temp)) #ese valor de la a va mejor
plot(temp,briere1(a=0.00002,Tmin=7.5,Tmax=33.5,temp))
plot(temp,briere1(a=0.2,Tmin=7.5,Tmax=33.5,temp))
#repetimos el análisis
fitted_br1 <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                  data = train_data,
                  start = list(a = 0.00017,
                               Tmin =8,
                               Tmax= 39)) #cojo 39 porque con Tmax<39 no sale
coef(fitted_br1)
AIC(fitted_br1)
BIC(fitted_br1)
#plot
temp <- train_data$temp
growth <- train_data$r
pred <- tibble(temp,growth,predict(fitted_br1))
colnames(pred) <- c("temp","growth","fit")
View(pred)
briere_plot_train <- ggplot(pred,aes(temp,fit))+
  geom_point(aes(x=temp,y=growth))+
  geom_line(color="lightcoral",size=2)+
  theme_bw()+
  geom_smooth(aes(x=temp,y=growth))+
  labs(y= "Intrinsic rate of increase",
       x="Temperature (ºC)",
       title="Brière-1 fit",
       subtitle="train data")+
  annotate(geom = "text", x = 10, y = 0.25, 
           label = "a = 0.00007026,
Tmin = -3.26ºC
Tmax = 38.88ºC", hjust = 0, vjust = 1, size = 4)
briere_plot_train
#voy probando con Tmax aumentando 1ºC cada vez hasta que salga que converge

#### 4. Database examples: probar con datos de ácaros de nuestra revisión ####
#### _ _ a) Acari ####
acari <- IR_data %>%
  filter(order == "Acari>Prostigmata" | 
           order == "Acari>Trombidiformes")%>%
  mutate(order = "Acari")%>%
  glimpse()
acari_test <- tibble(temp=acari$temperature,r=acari$growth_rate)
fitted_br1 <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                  data = acari_test,
                  start = list(a = 0.00017,
                               Tmin =4,
                               Tmax= 38))
coef(fitted_br1)
AIC(fitted_br1)
BIC(fitted_br1)
#plot
temp <- acari_test$temp
growth <- acari_test$r
pred <- tibble(temp,growth,predict(fitted_br1))
colnames(pred) <- c("temp","growth","fit")
View(pred)
briere_plot_acari <- ggplot(pred,aes(temp,fit))+
  geom_point(aes(x=temp,y=growth))+
  theme_bw()+
  geom_smooth(aes(x=temp,y=growth),fill="lightblue")+
  geom_line(color="lightcoral",size=2)+
  labs(y= "Intrinsic rate of increase",
       x="Temperature (ºC)",
       title="Brière-1 fit",
       subtitle="Order = Acari")+
  annotate(geom = "text", x = 12, y = 0.15, 
           label = "a = 0.000116,
Tmin = 8.14ºC
Tmax = 41.08ºC", hjust = 0, vjust = 1, size = 4)
briere_plot_acari

#### _ _ b) Lepidoptera ####
# a ver con otro subset...
lepidoptera <- IR_data %>%
  filter(order == "Lepidoptera")%>%
  glimpse()
lepidoptera_test <- tibble(temp=lepidoptera$temperature,r=lepidoptera$growth_rate)
# ver la loess para estimar a ojo valores Tmin y Tmax
ggplot(lepidoptera_test,aes(x=temp,y=r))+
  geom_point()+
  stat_smooth(method= "loess")+
  theme_classic()
# see starting values a
temp <- lepidoptera_test$temp
growth <- lepidoptera_test$r
plot(temp,growth)
plot(temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp))
plot(temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp)) 
plot(temp,briere1(a=0.00008,Tmin=8,Tmax=35,temp))#nos quedamos con este, más próximo al loess

fitted_br1 <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                  data = lepidoptera_test,
                  start = list(a = 0.00008,
                               Tmin =-5, #hasta aquí no cortaba... ver el loess
                               Tmax= 38)) # a partir de 38

predict(fitted_br1)
summary(fitted_br1)
AIC(fitted_br1)
BIC(fitted_br1)
coef(fitted_br1)
#vamos a dibujarlo
pred <- na.exclude(tibble(temp,growth))
pred <- pred%>%
  mutate(fit=predict(fitted_br1))
colnames(pred) <- c("temp","growth","fit")
View(pred)
briere_plot_lepidoptera <- ggplot(pred,aes(temp,fit))+
  geom_point(aes(x=temp,y=growth))+
  theme_bw()+
  geom_smooth(aes(x=temp,y=growth),fill="lightblue")+
  geom_line(color="lightcoral",size=2)+
  labs(y= "Intrinsic rate of increase",
       x="Temperature (ºC)",
       title="Brière-1 fit",
       subtitle="Order = Lepidoptera")+
  annotate(geom = "text", x = 10, y = 0.2, 
           label = "a = 0.0000609,
Tmin = 6.42ºC
Tmax = 37.52ºC", hjust = 0, vjust = 1, size = 4)
briere_plot_lepidoptera
#### _ _ c) Diptera ####

#otro subset...
diptera <- IR_data %>%
  filter(order == "Diptera")%>%
  glimpse()
diptera_test <- tibble(temp=diptera$temperature,r=diptera$growth_rate)
ggplot(diptera_test,aes(x=temp,y=r))+
  geom_point()+
  stat_smooth(method= "loess")+
  theme_classic()

# see starting values

temp <- diptera_test$temp
growth <- diptera_test$r
par(mfrow=c(2,2))
plot(temp,growth)
plot(temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp))
plot(temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp)) #nos quedamos con este otra vez, más próximo al loess
plot(temp,briere1(a=0.00015,Tmin=8,Tmax=35,temp))


fitted_br1 <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                  data = diptera_test,
                  start = list(a = 0.0001,
                               Tmin =8,
                               Tmax= 35),
                  trace = FALSE)
#aquí sí sale!
summary(fitted_br1)
AIC(fitted_br1)
BIC(fitted_br1)
coef(fitted_br1)
#vamos a dibujarlo
pred <- tibble(temp,growth,predict(fitted_br1))
colnames(pred) <- c("temp","growth","fit")
View(pred)
briere_plot_diptera <- ggplot(pred,aes(temp,fit))+
  geom_point(aes(x=temp,y=growth))+
  theme_bw()+
  geom_smooth(aes(x=temp,y=growth),fill="lightblue")+
  geom_line(color="lightcoral",size=2)+
  labs(y= "Intrinsic rate of increase",
       x="Temperature (ºC)",
       title="Brière-1 fit",
       subtitle="Order = Diptera")+
  annotate(geom = "text", x = 10, y = 0.15, 
label = "a = 0.000075,
Tmin = 7.33ºC
Tmax = 38.12ºC", hjust = 0, vjust = 1, size = 4)
briere_plot_diptera
#### _ _ d) Hemiptera ####
hemiptera <- IR_data %>%
  filter(order == "Hemiptera")%>%
  glimpse()
hemiptera_test <- tibble(temp=hemiptera$temperature,r=hemiptera$growth_rate)
ggplot(hemiptera_test,aes(x=temp,y=r))+
  geom_point()+
  stat_smooth(method= "loess")+
  theme_classic()
# see starting values
temp <- hemiptera_test$temp
growth <- hemiptera_test$r
par(mfrow=c(2,2))
plot(temp,growth)
plot(temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp))
plot(temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp)) 
plot(temp,briere1(a=0.00015,Tmin=8,Tmax=35,temp))#nos quedamos con este más próximo al loess


fitted_br1 <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                  data = hemiptera_test,
                  start = list(a = 0.00015,
                               Tmin =8,
                               Tmax= 40),
                  trace = FALSE)
#no converge

#let's use only aphids
aphids <- hemiptera%>%
  filter(family=="Aphididae")%>%
  glimpse()
aphids_test <- tibble(temp=aphids$temperature,r=aphids$growth_rate)
ggplot(aphids_test,aes(x=temp,y=r))+
  geom_point()+
  stat_smooth(method= "loess")+
  theme_classic()
# see starting values
temp <- aphids_test$temp
growth <- aphids_test$r
par(mfrow=c(2,2))
plot(temp,growth)
plot(temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp))#nos quedamos con este más próximo al loess
plot(temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp)) 
plot(temp,briere1(a=0.00015,Tmin=8,Tmax=35,temp))


fitted_br1 <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                  data = aphids_test,
                  start = list(a = 0.0002,
                               Tmin =8,
                               Tmax= 39),
                  trace = FALSE)
## no hay forma... eliminamos outliers
aphids_clean <- aphids%>%
  filter(growth_rate <0.6)
  glimpse()
aphids_test2 <- tibble(temp=aphids_clean$temperature,r=aphids_clean$growth_rate)
ggplot(aphids_test2,aes(x=temp,y=r))+
  geom_point()+
  stat_smooth(method= "loess")+
  theme_classic()
temp <- aphids_test2$temp
growth <- aphids_test2$r
par(mfrow=c(2,2))
plot(temp,growth)
plot(temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp))#nos quedamos con este más próximo al loess
plot(temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp)) 
plot(temp,briere1(a=0.00015,Tmin=8,Tmax=35,temp))
 fitted_br1 <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                   data = aphids_test2,
                   start = list(a = 0.0002,
                                Tmin =8,
                                Tmax= 35),
                   trace = FALSE)  
 #NO HAY FORMA
summary(fitted_br1)
AIC(fitted_br1)
BIC(fitted_br1)
coef(fitted_br1)
#aquí tampoco sale
summary(fitted_br1)

#### _ _ e) Coleoptera ####
coleoptera <- IR_data %>%
  filter(order == "Coleoptera")%>%
  glimpse()
coleoptera_test <- tibble(temp=coleoptera$temperature,r=coleoptera$growth_rate)
ggplot(coleoptera_test,aes(x=temp,y=r))+
  geom_point()+
  stat_smooth(method= "loess")+
  theme_classic()
# see starting values
temp <- coleoptera_test$temp
growth <- coleoptera_test$r
par(mfrow=c(2,2))
plot(temp,growth)
plot(temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp))
plot(temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp)) 
plot(temp,briere1(a=0.00006,Tmin=8,Tmax=35,temp))#nos quedamos con este más próximo al loess

fitted_br1 <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                  data = hemiptera_test,
                  start = list(a = 0.00006,
                               Tmin =12,
                               Tmax= 40),
                  trace = FALSE)
#Error in nls(r ~ briere1(a, temp, Tmin, Tmax), 
#### 5. All plots combined ####
library(cowplot)
plot_grid(briere_plot_acari,briere_plot_diptera,briere_plot_lepidoptera,
          nrow = 1,labels = c("A","B","C"))
