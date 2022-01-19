#### SCRIPT INFO ####     
#     Authors: Darío San Segundo Molina, Sara Villén Pérez, Ignacio Morales Castilla
#     Title: Brière-1 Modeling
#     Aim: fit models to each study to obtain ecologically-informative parameters (traits)
#     Date: October 2021
#_________________________ ####


#### 1. Load Dataset ####
rm(list=ls())
#library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(svglite)
library(nlstools)
library(tidyr)
library(nls2)
library(msm)
library(magrittr)
library(cowplot)
library(car)
library(nlme)
library(brms)
#load data: (different values for each temperature within a single study have been summarised to obtain one-treatment <- one response)
IR_data <-read_csv("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/Intrinsic_metaanalysis/IR_metaanalysis_suitable.csv") %>%
  filter(Filter_3 == "yes") %>% #select only those which have passed the filter 3 in the source csv
  #mutate(growth_rate = replace(growth_rate, growth_rate <= 0,0))%>% #we assign 0 to negative values of intrinsic rates (biological nonsense?)
  filter(as.numeric(temperature)<50) %>% #exclude possible missleading points
  dplyr::select(Authors,title,Filter_3,Approach,Subapproach,temperature,growth_rate,
         error,n_1,order,family,genus,species,feeding_guild,h_p_family,
         diet_family,RH,daylength,lat,lon)%>%
  glimpse()
#### 2. Define fitting function ####
# We use a Brière-1 model (Briere et al., 1999).This equation represents a trade-off 
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
#Try 
Topt(Tmin=8.14,Tmax=41.08,m=2) #ok

#### 3. Database examples: test with acari data from the database ####
#### _ _ a) Acari ####
acari <- IR_data %>%
  filter(order == "Acari>Prostigmata" | 
           order == "Acari>Trombidiformes")%>%
  mutate(order = "Acari")%>%
  glimpse()
acari_test <- tibble(temp=acari$temperature,r=acari$growth_rate)
#let's see exploratory data
temp <- acari_test$temp
growth <- acari_test$r
par(mfrow=c(2,2))
plot(temp,growth)
#now let's see predicted form different starting values for the function
plot(temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp))
plot(temp,briere1(a=0.0002,Tmin=8,Tmax=40,temp))
plot(temp,briere1(a=0.00015,Tmin=8,Tmax=40,temp))
# we make a grid to scan among a range of possible starting values
grid_br1_acari <- expand.grid(list(a=seq(0.00015,0.00025,by=0.00001),
                               Tmin=seq(5,15,by=0.5),
                               Tmax=seq(33,43,by=0.5)))
# and try to converge the model along the grid
fitted_br1_acari_grid <- nls2(r ~ briere1(a,temp,Tmin,Tmax),
                        data = acari_test,
                        start = grid_br1_acari,
                        algorithm = "brute-force",
                        trace=TRUE)
sum_fitted_br1_acari_grid <- summary(fitted_br1_acari_grid)
#now we use the estimates as starting values
fitted_br1_acari <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                        data = acari_test,
                        start = sum_fitted_br1_acari_grid$coefficients[,1])
model_acari_sum <-summary(fitted_br1_acari) #summary
coefs_acari <- data.frame(model_acari_sum$coefficients[,1:2],row.names = NULL) #call coeffs
#Topt is computed aaccording to Marchioro & Foerster (2011); formula above
Topt_acari <- Topt(Tmin=model_acari_sum$coefficients[2,1], #compute Topt
                   Tmax=model_acari_sum$coefficients[3,1],
                   m=2) #en Brière-1, m=2 siempre
#Topt standard error is calculated via deltamethod (see Ritz & Streibig, 2008)
Topt_acari_se <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                         coef(fitted_br1_acari), vcov(fitted_br1_acari))
# now we save the parameters into a data.frame
params_br1_acari <- data.frame(coefs_acari[1,],
                               coefs_acari[2,],
                               coefs_acari[3,],
                               Topt_acari,
                               Topt_acari_se,
                               "all")
colnames(params_br1_acari) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "which_signif")

params_br1_acari
#plot
temp <- acari_test$temp
growth <- acari_test$r
pred <- tibble(temp,growth,predict(fitted_br1_acari))
colnames(pred) <- c("temp","growth","fit")
briere_plot_acari <- ggplot(pred,aes(temp,fit))+
  geom_point(aes(x=temp,y=growth))+
  theme_bw()+
  geom_smooth(aes(x=temp,y=growth),fill="lightblue")+
  geom_line(color="lightcoral",size=2)+
  labs(y= "Intrinsic rate of increase",
       x="Temperature (ºC)",
       title="Brière-1 fit",
       subtitle="Order = Acari")+
  annotate(geom = "text", x = 10, y = 0.4, 
           label = "a = 0.000116,
Tmin = 8.14ºC
Topt = 33.79ºC
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
grid_br1_lepi <- expand.grid(list(a=seq(0.00005,0.00015, by = 0.00001),
                                   Tmin=seq(0,15,by=0.5),
                                   Tmax=seq(30,40,by=0.5)))
# and try to converge the model along the grid
fitted_br1_lepi_grid <- nls2(r ~ briere1(a,temp,Tmin,Tmax),
                              data = lepidoptera_test,
                              start = grid_br1_lepi,
                              algorithm = "brute-force",
                              trace=TRUE)
sum_fitted_br1_lepi_grid <- summary(fitted_br1_lepi_grid)

fitted_br1_lepidoptera <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                              data = lepidoptera_test,
                              start = sum_fitted_br1_lepi_grid$coefficients[,1])
model_lepidoptera_sum <- summary(fitted_br1_lepidoptera)
coefs_lepi <- data.frame(model_lepidoptera_sum$coefficients[,1:2],row.names = NULL)
Topt_lepidoptera <- Topt(Tmin=model_lepidoptera_sum$coefficients[2,1],
                   Tmax=model_lepidoptera_sum$coefficients[3,1],
                   m=2) #en Brière-1, m=2 siempre
Topt_lepidoptera_se <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                             coef(fitted_br1_lepidoptera), vcov(fitted_br1_lepidoptera))

params_br1_lepidoptera <- data.frame(coefs_lepi[1,],
                                     coefs_lepi[2,],
                                     coefs_lepi[3,],
                                     Topt_lepidoptera,
                                     Topt_lepidoptera_se,
                                     "all")
colnames(params_br1_lepidoptera) <- c("a_est_nls","a_se_nls",
                                "Tmin_est_nls","Tmin_se_nls",
                                "Tmax_est_nls","Tmax_se_nls",
                                "Topt_est","Topt_se_delta",
                                "which_signif")

params_br1_lepidoptera

#vamos a dibujarlo
pred <- na.exclude(tibble(temp,growth))
pred <- pred%>%
  mutate(fit=predict(fitted_br1_lepidoptera))
colnames(pred) <- c("temp","growth","fit")
#View(pred)
briere_plot_lepidoptera <- ggplot(pred,aes(temp,fit))+
  geom_point(aes(x=temp,y=growth))+
  theme_bw()+
  geom_smooth(aes(x=temp,y=growth),fill="lightblue")+
  geom_line(color="lightcoral",size=2)+
  labs(y= "Intrinsic rate of increase",
       x="Temperature (ºC)",
       title="Brière-1 fit",
       subtitle="Order = Lepidoptera")+
  annotate(geom = "text", x = 8, y = 0.25, 
           label = "a = 0.0000609,
Tmin = 6.42ºC
Topt = 30.73ºC
Tmax = 37.52ºC", hjust = 0, vjust = 1, size = 4)
briere_plot_lepidoptera

Parameters_lepidoptera <- data.frame(model_lepidoptera_sum$parameters[,1:2],
                                     Boot_fitbr1$estiboot[,1],
                                     Boot_fitbr1$estiboot[,2])
colnames(Parameters_acari) <- c("nls_estimate","nls_SE","Bootstrap_estimate","Bootstrap_SE")
Parameters_acari

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


fitted_br1_diptera <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                          data = diptera_test,
                          start = list(a = 0.0001,
                                       Tmin =8,
                                       Tmax= 35),
                          trace = FALSE)
#aquí sí sale!
#vamos a dibujarlo
pred <- tibble(temp,growth,predict(fitted_br1_diptera))
colnames(pred) <- c("temp","growth","fit")
briere_plot_diptera <- ggplot(pred,aes(temp,fit))+
  geom_point(aes(x=temp,y=growth))+
  theme_bw()+
  geom_smooth(aes(x=temp,y=growth),fill="lightblue")+
  geom_line(color="lightcoral",size=2)+
  labs(y= "Intrinsic rate of increase",
       x="Temperature (ºC)",
       title="Brière-1 fit",
       subtitle="Order = Diptera")+
  annotate(geom = "text", x = 8, y = 0.18, 
           label = "a = 0.000075,
Tmin = 7.33ºC
Topt =31.33ºC
Tmax = 38.12ºC", hjust = 0, vjust = 1, size = 4)
briere_plot_diptera

predict(fitted_br1_diptera)
model_diptera_sum <- summary(fitted_br1_diptera)
coefs_diptera <- data.frame(model_diptera_sum$coefficients[,1:2],row.names = NULL)
Topt_diptera <- Topt(Tmin=model_diptera_sum$coefficients[2,1],
                         Tmax=model_diptera_sum$coefficients[3,1],
                         m=2) #en Brière-1, m=2 siempre
Topt_diptera_se <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                                   coef(fitted_br1_diptera), vcov(fitted_br1_diptera))

params_br1_diptera <- data.frame(coefs_diptera[1,],
                                     coefs_diptera[2,],
                                     coefs_diptera[3,],
                                     Topt_diptera,
                                     Topt_diptera_se,
                                     "all")
colnames(params_br1_diptera) <- c("a_est_nls","a_se_nls",
                                      "Tmin_est_nls","Tmin_se_nls",
                                      "Tmax_est_nls","Tmax_se_nls",
                                      "Topt_est","Topt_se_delta",
                                      "which_signif")

params_br1_diptera

#### _ _ d) Hemiptera ####
hemiptera <- IR_data %>%
  filter(order == "Hemiptera")%>%
  filter(growth_rate!=0)%>%
  glimpse()
hemiptera_test <- tibble(temp=hemiptera$temperature,r=hemiptera$growth_rate)
ggplot(hemiptera_test,aes(x=temp,y=r))+
  geom_point()+
  stat_smooth(method= "loess")+
  theme_bw()
# see starting values
temp <- hemiptera_test$temp
growth <- hemiptera_test$r
par(mfrow=c(2,2))
plot(temp,growth)
plot(temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp))
plot(temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp)) 
plot(temp,briere1(a=0.00017,Tmin=8,Tmax=35,temp))#nos quedamos con este más próximo al loess
grid_br1_hemiptera <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(30,45,by=0.5)))
fitted_br1_hemiptera <- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                     data = hemiptera_test,
                     start = grid_br1_hemiptera,
                     algorithm = "brute-force",
                     trace = TRUE)
summary(fitted_br1_hemiptera)
#usamos esos como starting vals
fitted_br1_hemiptera <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                  data = hemiptera_test,
                  start = list(a = 0.00007,
                               Tmin = 0,
                               Tmax= 40)) #me salía error number iterations exceeded 50, 
# pero con este argumento sigue sacando valores
summary(fitted_br1) #solo se puede sacar el Tmax y a(40.51ºC)

#### _ _ _ d.1- Hemiptera> aphids ####
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
lmaph<-lm(r~temp,data=aphids_test)
plot(lmaph)
fitted_br1_aphids <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                         data = aphids_test,
                         start = list(a = 0.0002,
                                      Tmin =-20,
                                      Tmax= 34),
                         trace = FALSE)
summary(fitted_br1_aphids)
briere_plot_aphids <- ggplot(pred,aes(temp,fit))+
  geom_point(aes(x=temp,y=growth))+
  theme_bw()+
  geom_smooth(aes(x=temp,y=growth),fill="lightblue")+
  geom_line(color="lightcoral",size=2)+
  labs(y= "Intrinsic rate of increase",
       x="Temperature (ºC)",
       title="Brière-1 fit",
       subtitle="Order = Hemiptera, Family = Aphididae")+
  annotate(geom = "text", x = 10, y = 0.65, 
           label = "a = 0.000175,
Tmin = ns
Topt =26ºC
Tmax = 32.42ºC", hjust = 0, vjust = 1, size = 4)
briere_plot_aphids
model_aphids_sum <- summary(fitted_br1_aphids)
coefs_aphids <- data.frame(model_aphids_sum$coefficients[,1:2],row.names = NULL)
Topt_aphids <- Topt(Tmin=model_aphids_sum$coefficients[2,1],
                     Tmax=model_aphids_sum$coefficients[3,1],
                     m=2) #en Brière-1, m=2 siempre
Topt_aphids_se <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                               coef(fitted_br1_aphids), vcov(fitted_br1_aphids))

params_br1_aphids <- data.frame(coefs_aphids[1,],
                                 coefs_aphids[2,],
                                 coefs_aphids[3,],
                                 Topt_aphids,
                                 Topt_aphids_se,
                                 "a,Tmax")
colnames(params_br1_aphids) <- c("a_est_nls","a_se_nls",
                                  "Tmin_est_nls","Tmin_se_nls",
                                  "Tmax_est_nls","Tmax_se_nls",
                                  "Topt_est","Topt_se_delta",
                                  "which_signif")

params_br1_aphids



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
grid_br1_coleoptera <- expand.grid(list(a=seq(0.00001,0.0001,by=0.00001),
                                        Tmin=seq(5,15,by=0.5),
                                        Tmax=seq(25,45,by=0.5)))
fitted_br1_coleoptera <- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                              data = coleoptera_test,
                              start = grid_br1_coleoptera,
                              algorithm = "brute-force",
                              trace = TRUE)
summary(fitted_br1_coleoptera)
fitted_br1_coleoptera<- nls(r ~ briere1(a,temp,Tmin,Tmax),
                         data = coleoptera_test,
                         start = list(a = 0.00004,
                                      Tmin =7.5,
                                      Tmax= 36),
                         trace = FALSE)
model_coleoptera_sum <- summary(fitted_br1_coleoptera)
coefs_coleoptera <- data.frame(model_coleoptera_sum$coefficients[,1:2],row.names = NULL)
Topt_coleoptera <- Topt(Tmin=model_coleoptera_sum$coefficients[2,1],
                     Tmax=model_coleoptera_sum$coefficients[3,1],
                     m=2) #en Brière-1, m=2 siempre
Topt_coleoptera_se <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                               coef(fitted_br1_coleoptera), vcov(fitted_br1_coleoptera))

params_br1_coleoptera <- data.frame(coefs_coleoptera[1,],
                                 coefs_coleoptera[2,],
                                 coefs_coleoptera[3,],
                                 Topt_coleoptera,
                                 Topt_coleoptera_se,
                                 "a,Tmax")
colnames(params_br1_coleoptera) <- c("a_est_nls","a_se_nls",
                                  "Tmin_est_nls","Tmin_se_nls",
                                  "Tmax_est_nls","Tmax_se_nls",
                                  "Topt_est","Topt_se_delta",
                                  "which_signif")

params_br1_coleoptera

parameters_pooled <- rbind(params_br1_acari,params_br1_aphids,
                           params_br1_lepidoptera,params_br1_diptera,
                           params_br1_coleoptera)
parameters_pooled <- cbind(parameters_pooled,c("acari","aphids",
                                               "lepidoptera","diptera",
                                               "coleoptera"))
colnames(parameters_pooled)[10] <- "order"

####_ _ f) All plots combined ####
plot_grid(briere_plot_acari,briere_plot_diptera,briere_plot_lepidoptera, briere_plot_aphids,
          nrow = 2,labels = c("A","B","C","D"))

#### 4. Pooled data plots ####
IR_test <- IR_data_all %>%
  select(r=growth_rate,temp=temperature,order,lat,hostplant=h_p_family,title)%>%
  glimpse()
#### _ _ a) r ~ temp ####
#explorar el scatterplot agrupado
all_together <- ggplot(IR_test,aes(x=temp,y=r))+
  geom_point(color="turquoise4",alpha=0.2)+
  geom_smooth(color="lightcoral",fill="lightcoral")+
  theme_classic()
all_together

grid_br1_all <- expand.grid(list(a=seq(0.00001,0.0001,by=0.00001),
                                        Tmin=seq(-5,10,by=0.5),
                                        Tmax=seq(32,46,by=0.5)))
fitted_br1_all<- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                              data = IR_test,
                              start = grid_br1_all,
                              algorithm = "brute-force",
                              trace = TRUE)
sum_all_grid <- summary(fitted_br1_all)
fitt_br1_pooled <- nls(r~briere1(a,temp,Tmin,Tmax),
                       data= IR_test,
                       start=sum_all_grid$parameters[,1])

sum_all <-summary(fitt_br1_pooled)
coefs_pooled <-data.frame(coef(fitt_br1_pooled))
Topt_pooled<- Topt(Tmin=coefs_pooled[2,],
                   Tmax=coefs_pooled[3,],
                   m=2)

#plot
pred <- na.exclude(tibble(temp,growth))
pred <- pred%>%
  mutate(fit=predict(fitt_br1_pooled))
colnames(pred) <- c("temp","growth","fit")
View(pred)
briere_plot_pooled <- ggplot(pred,aes(temp,fit))+
  geom_point(aes(x=temp,y=growth),alpha=0.2,color="turquoise4",
             position=position_jitter(width=1))+
  theme_bw()+
  geom_smooth(aes(x=temp,y=growth),color="turquoise4",fill="turquoise3")+
  geom_line(color="lightcoral",size=1.2)+
  labs(y= "Intrinsic rate of increase",
       x="Temperature (ºC)",
       title="Brière-1 fit",
       subtitle="Pooled across taxa")+
  annotate(geom = "text", x = 10, y = 0.6, 
           label = "a = 0.000046,
Tmin = -6.77ºC (ns)
Topt = 32.18ºC (ns)
Tmax = 40.99ºC", hjust = 0, vjust = 1, size = 4)
briere_plot_pooled

#### 5. nlsTools & nlme packages ####
library(nlstools)
acari # acari data
#intervalos de confianza para un parámetro
confint2(fitted_br1_acari, "Tmin", level = 0.95, method = c("asymptotic", "profile"))

#resampling bootstrap
Boot_fitbr1 <- nlsBoot(fitted_br1_acari, niter = 999)
Boot_fitbr1$coefboot
Boot_fitbr1$bootCI
Boot_fitbr1$estiboot #not working?
plot(Boot_fitbr1)
#intregions
fitbr1Cont<-nlsContourRSS (fitted_br1_acari, lseq = 100, exp = 2)
plot(fitbr1Cont)

#residuals: check NORMALITY and HOMOSCEDASTICITY
fitbr1_resid <- nlsResiduals(fitted_br1_acari) #neither normal nor homoscedastic
car::leveneTest(acari_test$r,acari_test$temp) #homoscedastic though via test
plot(fitbr1_resid,which=0) #todos
plot(fitbr1_resid,which=5) #histogram of residuals
plot(fitbr1_resid,which=6) #qqplot
test.nlsResiduals(fitbr1_resid) #test de normalidad Shapiro-Wilk

#preview to scan data in searching for starting values
dev.off()
preview (formula=r ~ briere1(a,temp,Tmin,Tmax),
         data = acari_test,
         start = list(a = 0.00014, #el que mejor me ajusta el RSS
                      Tmin =8.6,
                      Tmax= 39),variable = 1)

# if model assumptions are not accomplished, parameter estimation is reliable but
# confidence intervals and statistical inference based on these estimates might not be
# reliable. 

# Dealing with assumptions violations in nlmer package
fitted_br1_var <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                       data = acari_test,
                       start = sum_fitted_br1_acari_grid$coefficients[,1],
                       weights = varPower())
summary(fitted_br1_var)

#transformations box-cox
library(nlrwr)
BC_fitted_br1 <- boxcox.nls(fitted_br1)
summary(BC_fitted_br1)

#### 6 individual studies regression with Bootstrap ####
## first: assign a number to each unique article as ID
IR_data_ID <- IR_data %>%
  distinct(title)%>%
  mutate(id=row_number()) #one id per distinct paper in a new dataframe
IR_data_all <- inner_join(IR_data,IR_data_ID,by='title')  #merge both dataframes
startVals_list <- list() #to store the starting values now on

#### _ _ Study 1 ####
IR_data_ID1 <- IR_data_all %>%
  filter(id==1)
par(mfrow=c(2,2))
plot(IR_data_ID1$temperature,IR_data_ID1$growth_rate)
plot(IR_data_ID1$temperature,briere1(a=0.0002,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))
plot(IR_data_ID1$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))# <-
plot(IR_data_ID1$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))
grid_br1_1 <- expand.grid(list(a=seq(0.00001,0.00025,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(28,40,by=0.5)))
fitted_br1_1 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID1,
                     start = grid_br1_1,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_br1_1_grid <- summary(fitted_br1_1)
startVals_1 <- c(sum_br1_1_grid$coefficients[,1])
# use those estimates as starting values
fitted_br1_1 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID1,
                     start = startVals_1,
                     trace = FALSE)
nlme_br1_1 <- nlme(model= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    start = startVals_1,
                    fixed = ~IR_data_all$order,
                    random ~ 1|id,
                    data = IR_data_all)

sum_br1_1 <- summary(fitted_br1_1)
startVals_1 <- c(sum_br1_1$coefficients[,1])
startVals_list[[1]] <- startVals_1 #llenar lista starting values

boot_br1_1 <- nlsBoot(fitted_br1_1, niter = 999)
coefs_1 <- data.frame(sum_br1_1$coefficients[,1:2],row.names = NULL)
boot_1 <- data.frame(boot_br1_1$estiboot)
Topt_est_1 <- Topt(Tmin=coefs_1[2,1],
                   Tmax=coefs_1[3,1],
                   m=2)
Topt_se_1 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                         coef(fitted_br1_1), vcov(fitted_br1_1))

params_br1_1 <- data.frame(coefs_1[1,],
                           coefs_1[2,],
                           coefs_1[3,],
                           Topt_est_1,
                           Topt_se_1,
                           boot_1[1,],
                           boot_1[2,],
                           boot_1[3,],
                           startVals_1[1],
                           startVals_1[2],
                           startVals_1[3],
                           "a,Tmax")
colnames(params_br1_1) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "starting_a","starting_Tmin",
                            "starting_Tmax","which_signif")

params_br1_1

#### _ _ Study 2 ####
IR_data_ID2 <- IR_data_all %>%
  filter(id==2)
par(mfrow=c(2,2))
plot(IR_data_ID2$temperature,IR_data_ID2$growth_rate)
plot(IR_data_ID2$temperature,briere1(a=0.00003,Tmin=8,Tmax=35,temp=IR_data_ID2$temperature))# <-
plot(IR_data_ID2$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID2$temperature))
plot(IR_data_ID2$temperature,briere1(a=0.00001,Tmin=8,Tmax=35,temp=IR_data_ID2$temperature))
grid_br1_2 <- expand.grid(list(a=seq(0.00001,0.0001,by=0.00001),
                               Tmin=seq(5,15,by=0.5),
                               Tmax=seq(28,40,by=0.5)))
fitted_br1_2 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID2,
                     start = grid_br1_2,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_br1_2 <- summary(fitted_br1_2)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_2 <- c(sum_br1_1$coefficients[,1])
startVals_list[[2]] <- startVals_2 #llenar lista starting values

# probar con esos parámetros como iniciales:


fitted_br1_2 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID2,
                    start = startVals_2)
sum_br1_2 <- summary(fitted_br1_2)
boot_br1_2 <- nlsBoot(fitted_br1_2, niter = 999)
coefs_2 <- data.frame(sum_br1_2$coefficients[,1:2],row.names = NULL)
boot_2 <- data.frame(boot_br1_2$estiboot)
Topt_est_2 <- Topt(Tmin=coefs_2[2,1],
                   Tmax=coefs_2[3,1],
                   m=2)
Topt_se_2 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                         coef(fitted_br1_2), vcov(fitted_br1_2))

params_br1_2 <- data.frame(coefs_2[1,],
                           coefs_2[2,],
                           coefs_2[3,],
                           Topt_est_2,
                           Topt_se_2,
                           boot_2[1,],
                           boot_2[2,],
                           boot_2[3,],
                           startVals_1[1],
                           startVals_1[2],
                           startVals_1[3],
                           "all")
colnames(params_br1_2) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "starting_a","starting_Tmin",
                            "starting_Tmax","which_signif")

params_br1_2

#### _ _ Study 3 ####
IR_data_ID3 <- IR_data_all %>%
  filter(id==3)
par(mfrow=c(2,2))
plot(IR_data_ID1$temperature,IR_data_ID1$growth_rate)
plot(IR_data_ID1$temperature,briere1(a=0.0002,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))
plot(IR_data_ID1$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))# <-
plot(IR_data_ID1$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))
grid_br1_3 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                               Tmin=seq(0,10,by=0.5),
                               Tmax=seq(28,38,by=0.5)))
fitted_br1_3 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID3,
                     start = grid_br1_3,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_br1_3 <- summary(fitted_br1_3)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_3 <- c(sum_br1_3$coefficients[,1])
startVals_list[[3]] <- startVals_3 #llenar lista starting values

# probar con esos parámetros como iniciales:
fitted_br1_3 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID3,
                    start = startVals_3)

sum_br1_3 <- summary(fitted_br1_3)
boot_br1_3 <- nlsBoot(fitted_br1_3, niter = 999)
coefs_3 <- data.frame(sum_br1_3$coefficients[,1:2],row.names = NULL)
boot_3 <- data.frame(boot_br1_3$estiboot)
Topt_est_3 <- Topt(Tmin=coefs_3[2,1],
                   Tmax=coefs_3[3,1],
                   m=2)
Topt_se_3 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                         coef(fitted_br1_3), vcov(fitted_br1_3))

params_br1_3 <- data.frame(coefs_3[1,],
                           coefs_3[2,],
                           coefs_3[3,],
                           Topt_est_3,
                           Topt_se_3,
                           boot_3[1,],
                           boot_3[2,],
                           boot_3[3,],
                           startVals_1[1],
                           startVals_1[2],
                           startVals_1[3],
                           "Tmax")
colnames(params_br1_3) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "starting_a","starting_Tmin",
                            "starting_Tmax","which_signif")

params_br1_3


#### _ _ Study 4 ####
IR_data_ID4 <- IR_data_all %>%
  filter(id==4)
par(mfrow=c(2,2))
plot(IR_data_ID4$temperature,IR_data_ID4$growth_rate)
plot(IR_data_ID4$temperature,briere1(a=0.00003,Tmin=8,Tmax=35,temp=IR_data_ID4$temperature))
plot(IR_data_ID4$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID4$temperature))# <-
plot(IR_data_ID4$temperature,briere1(a=0.00001,Tmin=8,Tmax=35,temp=IR_data_ID4$temperature))
grid_br1_4 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                               Tmin=seq(5,15,by=0.5),
                               Tmax=seq(28,38,by=0.5)))
fitted_br1_4 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID4,
                     start = grid_br1_4,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_br1_4 <- summary(fitted_br1_4)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_4 <- c(sum_br1_4$coefficients[,1])
startVals_list[[4]] <- startVals_4 #llenar lista starting values
fitted_br1_4 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID4,
                    start = startVals_4)

sum_br1_4 <- summary(fitted_br1_4)
boot_br1_4 <- nlsBoot(fitted_br1_4, niter = 999)
coefs_4 <- data.frame(sum_br1_4$coefficients[,1:2],row.names = NULL)
boot_4 <- data.frame(boot_br1_4$estiboot)
Topt_est_4 <- Topt(Tmin=coefs_4[2,1],
                   Tmax=coefs_4[3,1],
                   m=2)
Topt_se_4 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                         coef(fitted_br1_4), vcov(fitted_br1_4))

params_br1_4 <- data.frame(coefs_4[1,],
                           coefs_4[2,],
                           coefs_4[3,],
                           Topt_est_4,
                           Topt_se_4,
                           boot_4[1,],
                           boot_4[2,],
                           boot_4[3,],
                           startVals_1[1],
                           startVals_1[2],
                           startVals_1[3],
                           "none")
colnames(params_br1_4) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "starting_a","starting_Tmin",
                            "starting_Tmax","which_signif")
params_br1_4

#### _ _ Study 5 ####
IR_data_ID5 <- IR_data_all %>%
  filter(id==5)
par(mfrow=c(2,2))
plot(IR_data_ID5$temperature,IR_data_ID5$growth_rate)
plot(IR_data_ID5$temperature,briere1(a=0.00023,Tmin=8,Tmax=35,temp=IR_data_ID5$temperature))# <-
plot(IR_data_ID5$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID5$temperature))
plot(IR_data_ID5$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID5$temperature))
grid_br1_5 <- expand.grid(list(a=seq(0.0002,0.0003,by=0.00001),
                               Tmin=seq(5,15,by=0.5),
                               Tmax=seq(28,38,by=0.5)))
fitted_br1_5 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID5,
                     start = grid_br1_5,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_br1_5 <- summary(fitted_br1_5)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_5 <- c(sum_br1_5$coefficients[,1])
startVals_list[[5]] <- startVals_5 #llenar lista starting values

#con esos valores introducimos nuevos iniciales
#fitted_br1_5 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
 #                 data = IR_data_ID5,
  #              start = startVals_5)

sum_br1_5 <- summary(fitted_br1_5)
boot_br1_5 <- nlsBoot(fitted_br1_5, niter = 999)
coefs_5 <- data.frame(sum_br1_5$coefficients[,1:2],row.names = NULL)
boot_5 <- data.frame(boot_br1_5$estiboot)
Topt_est_5 <- Topt(Tmin=coefs_5[2,1],
                   Tmax=coefs_5[3,1],
                   m=2)
Topt_se_5 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                         coef(fitted_br1_5), vcov(fitted_br1_5))

params_br1_5 <- data.frame(coefs_5[1,],
                           coefs_5[2,],
                           coefs_5[3,],
                           Topt_est_5,
                           Topt_se_5,
                           boot_5[1,],
                           boot_5[2,],
                           boot_5[3,],
                           startVals_1[1],
                           startVals_1[2],
                           startVals_1[3],
                           "all")
colnames(params_br1_5) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "starting_a","starting_Tmin",
                            "starting_Tmax","which_signif")
params_br1_5

#### _ _ Study 6 ####
IR_data_ID6 <- IR_data_all %>%
  filter(id==6)
par(mfrow=c(2,2))
plot(IR_data_ID6$temperature,IR_data_ID6$growth_rate)
plot(IR_data_ID6$temperature,briere1(a=0.00023,Tmin=8,Tmax=35,temp=IR_data_ID6$temperature))
plot(IR_data_ID6$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID6$temperature))
plot(IR_data_ID6$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID6$temperature))# <-
grid_br1_6 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(28,45,by=0.5)))
fitted_br1_6 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID6,
                     start = grid_br1_6,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_br1_6 <- summary(fitted_br1_6)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_6 <- c(sum_br1_6$coefficients[,1])
startVals_list[[6]] <- startVals_6 #llenar lista starting values
#con esos valores introducimos nuevos iniciales
fitted_br1_6 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID6,
                    start = startVals_6)

sum_br1_6 <- summary(fitted_br1_6)
boot_br1_6 <- nlsBoot(fitted_br1_6, niter = 999)
coefs_6 <- data.frame(sum_br1_6$coefficients[,1:2],row.names = NULL)
boot_6 <- data.frame(boot_br1_6$estiboot)
Topt_est_6 <- Topt(Tmin=coefs_6[2,1],
                   Tmax=coefs_6[3,1],
                   m=2)
Topt_se_6 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                         coef(fitted_br1_6), vcov(fitted_br1_6))

params_br1_6 <- data.frame(coefs_6[1,],
                           coefs_6[2,],
                           coefs_6[3,],
                           Topt_est_6,
                           Topt_se_6,
                           boot_6[1,],
                           boot_6[2,],
                           boot_6[3,],
                           startVals_6[1],
                           startVals_6[2],
                           startVals_6[3],
                           "Tmin,Tmax")
colnames(params_br1_6) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "starting_a","starting_Tmin",
                            "starting_Tmax","which_signif")
params_br1_6
#### _ _ Study 7 ####
IR_data_ID7 <- IR_data_all %>%
  filter(id==7)
par(mfrow=c(2,2))
plot(IR_data_ID7$temperature,IR_data_ID7$growth_rate)
plot(IR_data_ID7$temperature,briere1(a=0.00023,Tmin=8,Tmax=35,temp=IR_data_ID7$temperature))
plot(IR_data_ID7$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID7$temperature))
plot(IR_data_ID7$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID7$temperature))# <-
grid_br1_7 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                               Tmin=seq(5,15,by=0.5),
                               Tmax=seq(28,42,by=0.5)))
fitted_br1_7 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID7,
                     start = grid_br1_7,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_br1_7 <- summary(fitted_br1_7)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_7 <- c(sum_br1_7$coefficients[,1])
startVals_list[[7]] <- startVals_7 #llenar lista starting values
fitted_br1_7 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID7,
                    start = startVals_7)

sum_br1_7 <- summary(fitted_br1_7)
boot_br1_7 <- nlsBoot(fitted_br1_7, niter = 999)
coefs_7 <- data.frame(sum_br1_7$coefficients[,1:2],row.names = NULL)
boot_7 <- data.frame(boot_br1_7$estiboot)
Topt_est_7 <- Topt(Tmin=coefs_7[2,1],
                   Tmax=coefs_7[3,1],
                   m=2)
Topt_se_7 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                         coef(fitted_br1_7), vcov(fitted_br1_7))

params_br1_7 <- data.frame(coefs_7[1,],
                           coefs_7[2,],
                           coefs_7[3,],
                           Topt_est_7,
                           Topt_se_7,
                           boot_7[1,],
                           boot_7[2,],
                           boot_7[3,],
                           startVals_7[1],
                           startVals_7[2],
                           startVals_7[3],
                           "Tmin,Tmax")
colnames(params_br1_7) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "starting_a","starting_Tmin",
                            "starting_Tmax","which_signif")
params_br1_7
#### _ _ Study 8 ####
IR_data_ID8 <- IR_data_all %>%
  filter(id==8)
par(mfrow=c(2,2))
plot(IR_data_ID8$temperature,IR_data_ID8$growth_rate)
plot(IR_data_ID8$temperature,briere1(a=0.00023,Tmin=8,Tmax=35,temp=IR_data_ID8$temperature))
plot(IR_data_ID8$temperature,briere1(a=0.00025,Tmin=8,Tmax=35,temp=IR_data_ID8$temperature))
plot(IR_data_ID8$temperature,briere1(a=0.00027,Tmin=8,Tmax=35,temp=IR_data_ID8$temperature))# <-
grid_br1_8 <- expand.grid(list(a=seq(0.0001,0.00025,by=0.00001),
                               Tmin=seq(-5,10,by=0.5),
                               Tmax=seq(28,36,by=0.5)))

fitted_br1_8 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID8,
                     start = grid_br1_8,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_br1_8 <- summary(fitted_br1_8)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_8 <- c(sum_br1_8$coefficients[,1])
startVals_list[[8]] <- startVals_8 #llenar lista starting values
fitted_br1_8 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                   data = IR_data_ID8,
                  start = startVals_8)
sum_br1_8 <- summary(fitted_br1_8)
boot_br1_8 <- nlsBoot(fitted_br1_8, niter = 999)
coefs_8 <- data.frame(sum_br1_8$coefficients[,1:2],row.names = NULL)
boot_8 <- data.frame(boot_br1_8$estiboot)
Topt_est_8 <- Topt(Tmin=coefs_8[2,1],
                   Tmax=coefs_8[3,1],
                   m=2)
Topt_se_8 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                         coef(fitted_br1_8), vcov(fitted_br1_8))

params_br1_8 <- data.frame(coefs_8[1,],
                           coefs_8[2,],
                           coefs_8[3,],
                           Topt_est_8,
                           Topt_se_8,
                           NA, NA,
                           NA,NA,
                           NA,NA,
                           startVals_8[1],
                           startVals_8[2],
                           startVals_8[3],
                           "Tmax")
colnames(params_br1_8) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "starting_a","starting_Tmin",
                            "starting_Tmax","which_signif")
params_br1_8

#### _ _ Study 9 ####
IR_data_ID9 <- IR_data_all %>%
  filter(id==9)
par(mfrow=c(2,2))
plot(IR_data_ID9$temperature,IR_data_ID9$growth_rate)
plot(IR_data_ID9$temperature,briere1(a=0.00023,Tmin=8,Tmax=35,temp=IR_data_ID9$temperature))
plot(IR_data_ID9$temperature,briere1(a=0.00008,Tmin=8,Tmax=35,temp=IR_data_ID9$temperature))# <-
plot(IR_data_ID9$temperature,briere1(a=0.00007,Tmin=8,Tmax=35,temp=IR_data_ID9$temperature))
grid_br1_9 <- expand.grid(list(a=seq(0.00006,0.00012,by=0.00001),
                               Tmin=seq(5,20,by=0.5),
                               Tmax=seq(28,45,by=0.5)))
fitted_br1_9 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID9,
                     start = grid_br1_9,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_br1_9 <- summary(fitted_br1_9)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_9 <- c(sum_br1_9$coefficients[,1])
startVals_list[[9]] <- startVals_9 
fitted_br1_9 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID9,
                    start = list(a=0.00008,Tmin=15,Tmax=33.5))

sum_br1_9 <- summary(fitted_br1_9)
boot_br1_9 <- nlsBoot(fitted_br1_9, niter = 999)
coefs_9 <- data.frame(sum_br1_9$coefficients[,1:2],row.names = NULL)
boot_9 <- data.frame(boot_br1_9$estiboot)
Topt_est_9 <- Topt(Tmin=coefs_9[2,1],
                   Tmax=coefs_9[3,1],
                   m=2)
Topt_se_9 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                         coef(fitted_br1_9), vcov(fitted_br1_9))

params_br1_9 <- data.frame(coefs_9[1,],
                           coefs_9[2,],
                           coefs_9[3,],
                           Topt_est_9,
                           Topt_se_9,
                           boot_9[1,],
                           boot_9[2,],
                           boot_9[3,],
                           startVals_9[1],
                           startVals_9[2],
                           startVals_9[3],
                           "all")
colnames(params_br1_9) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "starting_a","starting_Tmin",
                            "starting_Tmax","which_signif")
params_br1_9

#### _ _ Study 10 ####
IR_data_ID10 <- IR_data_all %>%
  filter(id==10)
par(mfrow=c(2,2))
plot(IR_data_ID10$temperature,IR_data_ID10$growth_rate)
plot(IR_data_ID10$temperature,briere1(a=0.00023,Tmin=8,Tmax=35,temp=IR_data_ID10$temperature))
plot(IR_data_ID10$temperature,briere1(a=0.00015,Tmin=18,Tmax=35.5,temp=IR_data_ID10$temperature))
plot(IR_data_ID10$temperature,briere1(a=0.0002,Tmin=18,Tmax=35,temp=IR_data_ID10$temperature))# <-
grid_br1_10 <- expand.grid(list(a=seq(0.00015,0.0003,by=0.00001),
                                Tmin=seq(12,22,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_10 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID10,
                      start = grid_br1_10,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_10 <- summary(fitted_br1_10)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_10 <- c(sum_br1_10$coefficients[,1])
startVals_list[[10]] <- startVals_10 
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_10 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID10,
                     start = startVals_10)

sum_br1_10 <- summary(fitted_br1_10)
boot_br1_10 <- nlsBoot(fitted_br1_10, niter = 999)
coefs_10 <- data.frame(sum_br1_10$coefficients[,1:2],row.names = NULL)
boot_10 <- data.frame(boot_br1_10$estiboot)
Topt_est_10 <- Topt(Tmin=coefs_10[2,1],
                    Tmax=coefs_10[3,1],
                    m=2)
Topt_se_10 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_10), vcov(fitted_br1_10))

params_br1_10 <- data.frame(coefs_10[1,],
                            coefs_10[2,],
                            coefs_10[3,],
                            Topt_est_10,
                            Topt_se_10,
                            boot_10[1,],
                            boot_10[2,],
                            boot_10[3,],
                            startVals_10[1],
                            startVals_10[2],
                            startVals_10[3],
                            "Tmax")
colnames(params_br1_10) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_10

#### _ _ Study 11 ####
IR_data_ID11 <- IR_data_all %>%
  filter(id==11)
par(mfrow=c(2,2))
plot(IR_data_ID11$temperature,IR_data_ID11$growth_rate)
plot(IR_data_ID11$temperature,briere1(a=0.00023,Tmin=13,Tmax=33,temp=IR_data_ID11$temperature))
plot(IR_data_ID11$temperature,briere1(a=0.00020,Tmin=13,Tmax=31,temp=IR_data_ID11$temperature))
plot(IR_data_ID11$temperature,briere1(a=0.00021,Tmin=12,Tmax=31,temp=IR_data_ID11$temperature))# <-
grid_br1_11 <- expand.grid(list(a=seq(0.00018,0.00025,by=0.00001),
                                Tmin=seq(5,16,by=0.5),
                                Tmax=seq(28,35,by=0.5)))
fitted_br1_11 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID11,
                      start = grid_br1_11,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_11 <- summary(fitted_br1_11)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_11 <- c(sum_br1_11$coefficients[,1])
startVals_list[[11]] <- startVals_11 
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_11 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID11,
                     start = startVals_11)

sum_br1_11 <- summary(fitted_br1_11)
boot_br1_11 <- nlsBoot(fitted_br1_11, niter = 999)
coefs_11 <- data.frame(sum_br1_11$coefficients[,1:2],row.names = NULL)
boot_11 <- data.frame(boot_br1_11$estiboot)
Topt_est_11 <- Topt(Tmin=coefs_11[2,1],
                    Tmax=coefs_11[3,1],
                    m=2)
Topt_se_11 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_11), vcov(fitted_br1_11))

params_br1_11 <- data.frame(coefs_11[1,],
                            coefs_11[2,],
                            coefs_11[3,],
                            Topt_est_11,
                            Topt_se_11,
                            boot_11[1,],
                            boot_11[2,],
                            boot_11[3,],
                            startVals_11[1],
                            startVals_11[2],
                            startVals_11[3],
                            "all")
colnames(params_br1_11) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_11
#### _ _ Study 12 ####
IR_data_ID12 <- IR_data_all %>%
  filter(id==12)
par(mfrow=c(2,2))
plot(IR_data_ID12$temperature,IR_data_ID12$growth_rate)
plot(IR_data_ID12$temperature,briere1(a=0.00023,Tmin=13,Tmax=33,temp=IR_data_ID12$temperature))
plot(IR_data_ID12$temperature,briere1(a=0.00025,Tmin=13,Tmax=35,temp=IR_data_ID12$temperature))
plot(IR_data_ID12$temperature,briere1(a=0.00021,Tmin=13,Tmax=36,temp=IR_data_ID12$temperature))# <-
grid_br1_12 <- expand.grid(list(a=seq(0.00018,0.00025,by=0.00001),
                                Tmin=seq(10,16,by=0.5),
                                Tmax=seq(33,39,by=0.5)))
fitted_br1_12 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID12,
                      start = grid_br1_12,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_12 <- summary(fitted_br1_12)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_12 <- c(sum_br1_12$coefficients[,1])
startVals_list[[12]] <- startVals_12 
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste

fitted_br1_12 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID12,
                     start = startVals_12)

sum_br1_12 <- summary(fitted_br1_12)
sum_br1_12 #todos significativos
boot_br1_12 <- nlsBoot(fitted_br1_12, niter = 999)
coefs_12 <- data.frame(sum_br1_12$coefficients[,1:2],row.names = NULL)
boot_12 <- data.frame(boot_br1_12$estiboot)
Topt_est_12 <- Topt(Tmin=coefs_12[2,1],
                    Tmax=coefs_12[3,1],
                    m=2)
Topt_se_12 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_12), vcov(fitted_br1_12))

params_br1_12 <- data.frame(coefs_12[1,],
                            coefs_12[2,],
                            coefs_12[3,],
                            Topt_est_12,
                            Topt_se_12,
                            boot_12[1,],
                            boot_12[2,],
                            boot_12[3,],
                            startVals_12[1],
                            startVals_12[2],
                            startVals_12[3],
                            "all")
colnames(params_br1_12) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_12

#### _ _ Study 13 ####
IR_data_ID13 <- IR_data_all %>%
  filter(id==13)
par(mfrow=c(2,2))
plot(IR_data_ID13$temperature,IR_data_ID13$growth_rate)
plot(IR_data_ID13$temperature,briere1(a=0.00023,Tmin=13,Tmax=33,temp=IR_data_ID13$temperature))
plot(IR_data_ID13$temperature,briere1(a=0.00020,Tmin=8,Tmax=36,temp=IR_data_ID13$temperature))
plot(IR_data_ID13$temperature,briere1(a=0.00015,Tmin=5,Tmax=35,temp=IR_data_ID13$temperature))# <-
grid_br1_13 <- expand.grid(list(a=seq(0.00012,0.00025,by=0.00001),
                                Tmin=seq(0,15,by=0.5),
                                Tmax=seq(30,39,by=0.5)))
fitted_br1_13 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID13,
                      start = grid_br1_13,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_13 <- summary(fitted_br1_13)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_13 <- c(sum_br1_13$coefficients[,1])
startVals_list[[13]] <- startVals_13 
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_13 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID13,
                     start = startVals_13)

sum_br1_13 <- summary(fitted_br1_13)
sum_br1_13 
boot_br1_13 <- nlsBoot(fitted_br1_13, niter = 999)
coefs_13 <- data.frame(sum_br1_13$coefficients[,1:2],row.names = NULL)
boot_13 <- data.frame(boot_br1_13$estiboot)
Topt_est_13 <- Topt(Tmin=coefs_13[2,1],
                    Tmax=coefs_13[3,1],
                    m=2)
Topt_se_13 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_13), vcov(fitted_br1_13))

params_br1_13 <- data.frame(coefs_13[1,],
                            coefs_13[2,],
                            coefs_13[3,],
                            Topt_est_13,
                            Topt_se_13,
                            boot_13[1,],
                            boot_13[2,],
                            boot_13[3,],
                            startVals_13[1],
                            startVals_13[2],
                            startVals_13[3],
                            "a,Tmax")
colnames(params_br1_13) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_13

#### _ _ Study 14 ####
IR_data_ID14 <- IR_data_all %>%
  filter(id==14)
par(mfrow=c(2,2))
plot(IR_data_ID14$temperature,IR_data_ID14$growth_rate)
plot(IR_data_ID14$temperature,briere1(a=0.0002,Tmin=10,Tmax=34,temp=IR_data_ID14$temperature))
plot(IR_data_ID14$temperature,briere1(a=0.00015,Tmin=10,Tmax=34,temp=IR_data_ID14$temperature))
plot(IR_data_ID14$temperature,briere1(a=0.00009,Tmin=10,Tmax=34,temp=IR_data_ID14$temperature))# <-
grid_br1_14 <- expand.grid(list(a=seq(0.00006,0.00018,by=0.00001),
                                Tmin=seq(5,15,by=0.5),
                                Tmax=seq(30,39,by=0.5)))
fitted_br1_14 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID14,
                      start = grid_br1_14,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_14 <- summary(fitted_br1_14)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_14 <- c(sum_br1_14$coefficients[,1])
startVals_list[[14]] <- startVals_14 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_14 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID14,
                     start = startVals_14)

sum_br1_14 <- summary(fitted_br1_14)
sum_br1_14 #Tmin no significativo, a no significativo
boot_br1_14 <- nlsBoot(fitted_br1_14, niter = 999) #sale mal...
coefs_14 <- data.frame(sum_br1_14$coefficients[,1:2],row.names = NULL)
boot_14 <- data.frame(boot_br1_14$estiboot)
Topt_est_14 <- Topt(Tmin=coefs_14[2,1],
                    Tmax=coefs_14[3,1],
                    m=2)
Topt_se_14 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_14), vcov(fitted_br1_14))

params_br1_14 <- data.frame(coefs_14[1,],
                            coefs_14[2,],
                            coefs_14[3,],
                            Topt_est_14,
                            Topt_se_14,
                            NA,NA,
                            NA,NA,
                            NA,NA,
                            startVals_14[1],
                            startVals_14[2],
                            startVals_14[3],
                            "Tmax")
colnames(params_br1_14) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_14

#### _ _ Study 15 ####
IR_data_ID15 <- IR_data_all %>%
  filter(id==15)
par(mfrow=c(2,2))
plot(IR_data_ID15$temperature,IR_data_ID15$growth_rate)
plot(IR_data_ID15$temperature,briere1(a=0.00023,Tmin=13,Tmax=33,temp=IR_data_ID15$temperature))
plot(IR_data_ID15$temperature,briere1(a=0.00022,Tmin=13,Tmax=30,temp=IR_data_ID15$temperature))
plot(IR_data_ID15$temperature,briere1(a=0.00035,Tmin=13,Tmax=30,temp=IR_data_ID15$temperature))# <-
grid_br1_15 <- expand.grid(list(a=seq(0.00015,0.0003,by=0.00001),
                                Tmin=seq(8,18,by=0.5),
                                Tmax=seq(28,35,by=0.5)))
fitted_br1_15 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID15,
                      start = grid_br1_15,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_15 <- summary(fitted_br1_15)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_15 <- c(sum_br1_15$coefficients[,1])
startVals_list[[15]] <- startVals_15 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_15 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID15,
                     start = startVals_15)

sum_br1_15 <- summary(fitted_br1_15)
sum_br1_15 #Tmin no significativo y a tampoco
boot_br1_15 <- nlsBoot(fitted_br1_15, niter = 999)
coefs_15 <- data.frame(sum_br1_15$coefficients[,1:2],row.names = NULL)
boot_15 <- data.frame(boot_br1_15$estiboot)
Topt_est_15 <- Topt(Tmin=coefs_15[2,1],
                    Tmax=coefs_15[3,1],
                    m=2)
Topt_se_15 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_15), vcov(fitted_br1_15))

params_br1_15 <- data.frame(coefs_15[1,],
                            coefs_15[2,],
                            coefs_15[3,],
                            Topt_est_15,
                            Topt_se_15,
                            boot_15[1,],
                            boot_15[2,],
                            boot_15[3,],
                            startVals_15[1],
                            startVals_15[2],
                            startVals_15[3],
                            "Tmax")
colnames(params_br1_15) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_15

#### _ _ Study 16 ####
IR_data_ID16 <- IR_data_all %>%
  filter(id==16)
par(mfrow=c(2,2))
plot(IR_data_ID16$temperature,IR_data_ID16$growth_rate)
plot(IR_data_ID16$temperature,briere1(a=0.00023,Tmin=13,Tmax=33,temp=IR_data_ID16$temperature))
plot(IR_data_ID16$temperature,briere1(a=0.00025,Tmin=13,Tmax=35,temp=IR_data_ID16$temperature))
plot(IR_data_ID16$temperature,briere1(a=0.00020,Tmin=13,Tmax=35,temp=IR_data_ID16$temperature))# <-
grid_br1_16 <- expand.grid(list(a=seq(0.00016,0.00028,by=0.00001),
                                Tmin=seq(8,18,by=0.5),
                                Tmax=seq(32,40,by=0.5)))
fitted_br1_16 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID16,
                      start = grid_br1_16,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_16 <- summary(fitted_br1_16)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_16 <- c(sum_br1_16$coefficients[,1])
startVals_list[[16]] <- startVals_16 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_16 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID16,
                     start = startVals_16)

sum_br1_16 <- summary(fitted_br1_16)
sum_br1_16 
boot_br1_16 <- nlsBoot(fitted_br1_16, niter = 999)
coefs_16 <- data.frame(sum_br1_16$coefficients[,1:2],row.names = NULL)
boot_16 <- data.frame(boot_br1_16$estiboot)
Topt_est_16 <- Topt(Tmin=coefs_16[2,1],
                    Tmax=coefs_16[3,1],
                    m=2)
Topt_se_16 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_16), vcov(fitted_br1_16))

params_br1_16 <- data.frame(coefs_16[1,],
                            coefs_16[2,],
                            coefs_16[3,],
                            Topt_est_16,
                            Topt_se_16,
                            boot_16[1,],
                            boot_16[2,],
                            boot_16[3,],
                            startVals_16[1],
                            startVals_16[2],
                            startVals_16[3],
                            "all")
colnames(params_br1_16) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_16

#### _ _ Study 17 ####
IR_data_ID17 <- IR_data_all %>%
  filter(id==17)
par(mfrow=c(2,2))
plot(IR_data_ID17$temperature,IR_data_ID17$growth_rate)
plot(IR_data_ID17$temperature,briere1(a=0.0001,Tmin=8,Tmax=33,temp=IR_data_ID17$temperature))
plot(IR_data_ID17$temperature,briere1(a=0.0001,Tmin=15,Tmax=30,temp=IR_data_ID17$temperature))
plot(IR_data_ID17$temperature,briere1(a=0.00023,Tmin=17,Tmax=30,temp=IR_data_ID17$temperature))
grid_br1_17 <- expand.grid(list(a=seq(0.00018,0.00035,by=0.00001),
                                Tmin=seq(12,20,by=0.5),
                                Tmax=seq(28,35,by=0.5)))
fitted_br1_17 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID17,
                      start = grid_br1_17,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_17 <- summary(fitted_br1_17)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_17 <- c(sum_br1_17$coefficients[,1])
startVals_list[[17]] <- startVals_17 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_17 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID17,
                     start = startVals_17) # NO SALE

sum_br1_17 <- summary(fitted_br1_17)
sum_br1_17 #Tmin no significativo y a tampoco
boot_br1_17 <- nlsBoot(fitted_br1_17, niter = 999)
coefs_17 <- data.frame(sum_br1_17$coefficients[,1:2],row.names = NULL)
boot_17 <- data.frame(boot_br1_17$estiboot)
Topt_est_17 <- Topt(Tmin=coefs_17[2,1],
                    Tmax=coefs_17[3,1],
                    m=2)
Topt_se_17 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_17), vcov(fitted_br1_17))

params_br1_17 <- data.frame(coefs_17[1,],
                            coefs_17[2,],
                            coefs_17[3,],
                            Topt_est_17,
                            Topt_se_17,
                            boot_17[1,],
                            boot_17[2,],
                            boot_17[3,],
                            startVals_17[1],
                            startVals_17[2],
                            startVals_17[3],
                            "Tmax")
colnames(params_br1_17) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_17

#### _ _ Study 18 ####
IR_data_ID18 <- IR_data_all %>%
  filter(id==18)
par(mfrow=c(2,2))
plot(IR_data_ID18$temperature,IR_data_ID18$growth_rate)
plot(IR_data_ID18$temperature,briere1(a=0.00005,Tmin=13,Tmax=33,temp=IR_data_ID18$temperature))
plot(IR_data_ID18$temperature,briere1(a=0.00005,Tmin=13,Tmax=30,temp=IR_data_ID18$temperature))
plot(IR_data_ID18$temperature,briere1(a=0.00005,Tmin=8,Tmax=30,temp=IR_data_ID18$temperature))# <-
grid_br1_18 <- expand.grid(list(a=seq(0.00003,0.0001,by=0.00001),
                                Tmin=seq(6,18,by=0.5),
                                Tmax=seq(28,35,by=0.5)))
fitted_br1_18 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID18,
                      start = grid_br1_18,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_18 <- summary(fitted_br1_18)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_18 <- c(sum_br1_18$coefficients[,1])
startVals_list[[18]] <- startVals_18 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_18 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID18,
                     start = startVals_18)#no hay forma

sum_br1_18 <- summary(fitted_br1_18)
sum_br1_18 
boot_br1_18 <- nlsBoot(fitted_br1_18, niter = 999)
coefs_18 <- data.frame(sum_br1_18$coefficients[,1:2],row.names = NULL)
boot_18 <- data.frame(boot_br1_18$estiboot)
Topt_est_18 <- Topt(Tmin=coefs_18[2,1],
                    Tmax=coefs_18[3,1],
                    m=2)
Topt_se_18 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_18), vcov(fitted_br1_18))

params_br1_18 <- data.frame(coefs_18[1,],
                            coefs_18[2,],
                            coefs_18[3,],
                            Topt_est_18,
                            Topt_se_18,
                            boot_18[1,],
                            boot_18[2,],
                            boot_18[3,],
                            startVals_18[1],
                            startVals_18[2],
                            startVals_18[3],
                            "Tmax")
colnames(params_br1_18) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_18

#### _ _ Study 19 ####
IR_data_ID19 <- IR_data_all %>%
  filter(id==19)
par(mfrow=c(2,2))
plot(IR_data_ID19$temperature,IR_data_ID19$growth_rate)
plot(IR_data_ID19$temperature,briere1(a=0.00015,Tmin=15,Tmax=35,temp=IR_data_ID19$temperature))
plot(IR_data_ID19$temperature,briere1(a=0.00017,Tmin=15,Tmax=35,temp=IR_data_ID19$temperature))
plot(IR_data_ID19$temperature,briere1(a=0.0002,Tmin=17,Tmax=35,temp=IR_data_ID19$temperature))# <-
grid_br1_19 <- expand.grid(list(a=seq(0.00012,0.00030,by=0.00001),
                                Tmin=seq(11,19,by=0.5),
                                Tmax=seq(30,38,by=0.5)))
fitted_br1_19 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID19,
                      start = grid_br1_19,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_19 <- summary(fitted_br1_19)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_19 <- c(sum_br1_19$coefficients[,1])
startVals_list[[19]] <- startVals_19 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_19 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID19,
                     start = startVals_19)

sum_br1_19 <- summary(fitted_br1_19)
sum_br1_19 
boot_br1_19 <- nlsBoot(fitted_br1_19, niter = 999)
coefs_19 <- data.frame(sum_br1_19$coefficients[,1:2],row.names = NULL)
boot_19 <- data.frame(boot_br1_19$estiboot)
Topt_est_19 <- Topt(Tmin=coefs_19[2,1],
                    Tmax=coefs_19[3,1],
                    m=2)
Topt_se_19 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_19), vcov(fitted_br1_19))

params_br1_19 <- data.frame(coefs_19[1,],
                            coefs_19[2,],
                            coefs_19[3,],
                            Topt_est_19,
                            Topt_se_19,
                            boot_19[1,],
                            boot_19[2,],
                            boot_19[3,],
                            startVals_19[1],
                            startVals_19[2],
                            startVals_19[3],
                            "all")
colnames(params_br1_19) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_19

#### _ _ Study 20 ####
IR_data_ID20 <- IR_data_all %>%
  filter(id==20)
par(mfrow=c(2,2))
plot(IR_data_ID20$temperature,IR_data_ID20$growth_rate)
plot(IR_data_ID20$temperature,briere1(a=0.00002,Tmin=16,Tmax=31.5,temp=IR_data_ID20$temperature))
plot(IR_data_ID20$temperature,briere1(a=0.00009,Tmin=14,Tmax=30.5,temp=IR_data_ID20$temperature))
plot(IR_data_ID20$temperature,briere1(a=0.00009,Tmin=15,Tmax=31,temp=IR_data_ID20$temperature))# <-
grid_br1_20 <- expand.grid(list(a=seq(0.00001,0.00015,by=0.00001),
                                Tmin=seq(11,19,by=0.5),
                                Tmax=seq(28,35,by=0.5)))
fitted_br1_20 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID20,
                      start = grid_br1_20,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_20 <- summary(fitted_br1_20)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_20 <- c(sum_br1_20$coefficients[,1])
startVals_list[[20]] <- startVals_20 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_20 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID20,
                     start = startVals_20)#no hay forma

sum_br1_20 <- summary(fitted_br1_20)
sum_br1_20 
boot_br1_20 <- nlsBoot(fitted_br1_20, niter = 999)
coefs_20 <- data.frame(sum_br1_20$coefficients[,1:2],row.names = NULL)
boot_20 <- data.frame(boot_br1_20$estiboot)
Topt_est_20 <- Topt(Tmin=coefs_20[2,1],
                    Tmax=coefs_20[3,1],
                    m=2)
Topt_se_20 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_20), vcov(fitted_br1_20))

params_br1_20 <- data.frame(coefs_20[1,],
                            coefs_20[2,],
                            coefs_20[3,],
                            Topt_est_20,
                            Topt_se_20,
                            boot_20[1,],
                            boot_20[2,],
                            boot_20[3,],
                            startVals_20[1],
                            startVals_20[2],
                            startVals_20[3],
                            "Tmax")
colnames(params_br1_20) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_20

#### _ _ Study 21 ####
IR_data_ID21 <- IR_data_all %>%
  filter(id==21)
par(mfrow=c(2,2))
plot(IR_data_ID21$temperature,IR_data_ID21$growth_rate)
plot(IR_data_ID21$temperature,briere1(a=0.0001,Tmin=15,Tmax=35,temp=IR_data_ID21$temperature))
plot(IR_data_ID21$temperature,briere1(a=0.00009,Tmin=12,Tmax=35,temp=IR_data_ID21$temperature))
plot(IR_data_ID21$temperature,briere1(a=0.00009,Tmin=13,Tmax=35,temp=IR_data_ID21$temperature))# <-
grid_br1_21 <- expand.grid(list(a=seq(0.00007,0.00013,by=0.00001),
                                Tmin=seq(10,19,by=0.5),
                                Tmax=seq(30,37,by=0.5)))
fitted_br1_21 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID21,
                      start = grid_br1_21,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_21 <- summary(fitted_br1_21)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_21 <- c(sum_br1_21$coefficients[,1])
startVals_list[[21]] <- startVals_21 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_21 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID21,
                     start = startVals_21)#no hay forma

sum_br1_21 <- summary(fitted_br1_21)
sum_br1_21 
boot_br1_21 <- nlsBoot(fitted_br1_21, niter = 999)
coefs_21 <- data.frame(sum_br1_21$coefficients[,1:2],row.names = NULL)
boot_21 <- data.frame(boot_br1_21$estiboot)
Topt_est_21 <- Topt(Tmin=coefs_21[2,1],
                    Tmax=coefs_21[3,1],
                    m=2)
Topt_se_21 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_21), vcov(fitted_br1_21))

params_br1_21 <- data.frame(coefs_21[1,],
                            coefs_21[2,],
                            coefs_21[3,],
                            Topt_est_21,
                            Topt_se_21,
                            boot_21[1,],
                            boot_21[2,],
                            boot_21[3,],
                            startVals_21[1],
                            startVals_21[2],
                            startVals_21[3],
                            "all")
colnames(params_br1_21) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_21

#### _ _ Study 22 ####
IR_data_ID22 <- IR_data_all %>%
  filter(id==22)
par(mfrow=c(2,2))
plot(IR_data_ID22$temperature,IR_data_ID22$growth_rate)
plot(IR_data_ID22$temperature,briere1(a=0.0001,Tmin=15,Tmax=35,temp=IR_data_ID22$temperature))
plot(IR_data_ID22$temperature,briere1(a=0.00009,Tmin=12,Tmax=40,temp=IR_data_ID22$temperature))
plot(IR_data_ID22$temperature,briere1(a=0.00004,Tmin=13,Tmax=42,temp=IR_data_ID22$temperature))# <-
grid_br1_22 <- expand.grid(list(a=seq(0.00001,0.00015,by=0.00001),
                                Tmin=seq(5,19,by=0.5),
                                Tmax=seq(35,45,by=0.5)))
fitted_br1_22 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID22,
                      start = grid_br1_22,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_22 <- summary(fitted_br1_22)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_22 <- c(sum_br1_22$coefficients[,1])
startVals_list[[22]] <- startVals_22 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste

# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_22 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID22,
                     start = startVals_22) #NO sale
sum_br1_22 <- summary(fitted_br1_22)
sum_br1_22 
boot_br1_22 <- nlsBoot(fitted_br1_22, niter = 999)
coefs_22 <- data.frame(sum_br1_22$coefficients[,1:2],row.names = NULL)
boot_22 <- data.frame(boot_br1_22$estiboot)
Topt_est_22 <- Topt(Tmin=coefs_22[2,1],
                    Tmax=coefs_22[3,1],
                    m=2)
Topt_se_22 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_22), vcov(fitted_br1_22))

params_br1_22 <- data.frame(coefs_22[1,],
                            coefs_22[2,],
                            coefs_22[3,],
                            Topt_est_22,
                            Topt_se_22,
                            boot_22[1,],
                            boot_22[2,],
                            boot_22[3,],
                            startVals_22[1],
                            startVals_22[2],
                            startVals_22[3],
                            "Tmax")
colnames(params_br1_22) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_22

#### _ _ Study 23 ####
IR_data_ID23 <- IR_data_all %>%
  filter(id==23)
par(mfrow=c(2,2))
plot(IR_data_ID23$temperature,IR_data_ID23$growth_rate)
plot(IR_data_ID23$temperature,briere1(a=0.0001,Tmin=20,Tmax=33,temp=IR_data_ID23$temperature))
plot(IR_data_ID23$temperature,briere1(a=0.00015,Tmin=20,Tmax=32,temp=IR_data_ID23$temperature))
plot(IR_data_ID23$temperature,briere1(a=0.00019,Tmin=20,Tmax=31,temp=IR_data_ID23$temperature))# <-
grid_br1_23 <- expand.grid(list(a=seq(0.0001,0.00025,by=0.00001),
                                Tmin=seq(15,25,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_23 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID23,
                      start = grid_br1_23,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_23 <- summary(fitted_br1_23)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_23 <- c(sum_br1_23$coefficients[,1])
startVals_list[[23]] <- startVals_23 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_23 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID23,
                     start = list(a=0.00021,Tmin=21.5,Tmax=33))
sum_br1_23 <- summary(fitted_br1_23)
sum_br1_23 
boot_br1_23 <- nlsBoot(fitted_br1_23, niter = 999)
coefs_23 <- data.frame(sum_br1_23$coefficients[,1:2],row.names = NULL)
boot_23 <- data.frame(boot_br1_23$estiboot)
Topt_est_23 <- Topt(Tmin=coefs_23[2,1],
                    Tmax=coefs_23[3,1],
                    m=2)
Topt_se_23 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_23), vcov(fitted_br1_23))

params_br1_23 <- data.frame(coefs_23[1,],
                            coefs_23[2,],
                            coefs_23[3,],
                            Topt_est_23,
                            Topt_se_23,
                            boot_23[1,],
                            boot_23[2,],
                            boot_23[3,],
                            startVals_23[1],
                            startVals_23[2],
                            startVals_23[3],
                            "Tmax")
colnames(params_br1_23) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_23

#### _ _ Study 24 ####
IR_data_ID24 <- IR_data_all %>%
  filter(id==24)
par(mfrow=c(2,2))
plot(IR_data_ID24$temperature,IR_data_ID24$growth_rate)
plot(IR_data_ID24$temperature,briere1(a=0.0001,Tmin=20,Tmax=33,temp=IR_data_ID24$temperature))
plot(IR_data_ID24$temperature,briere1(a=0.0003,Tmin=8,Tmax=30,temp=IR_data_ID24$temperature))
plot(IR_data_ID24$temperature,briere1(a=0.00025,Tmin=8,Tmax=30,temp=IR_data_ID24$temperature))# <-
grid_br1_24 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(0,10,by=0.5),
                                Tmax=seq(25,35,by=0.5)))
fitted_br1_24 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID24,
                      start = grid_br1_24,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_24 <- summary(fitted_br1_24)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_24 <- c(sum_br1_24$coefficients[,1])
startVals_list[[24]] <- startVals_24 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_24 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID24,
                     start = startVals_24) #NO sale
sum_br1_24 <- summary(fitted_br1_24)
sum_br1_24 
boot_br1_24 <- nlsBoot(fitted_br1_24, niter = 999)
coefs_24 <- data.frame(sum_br1_24$coefficients[,1:2],row.names = NULL)
boot_24 <- data.frame(boot_br1_24$estiboot)
Topt_est_24 <- Topt(Tmin=coefs_24[2,1],
                    Tmax=coefs_24[3,1],
                    m=2)
Topt_se_24 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_24), vcov(fitted_br1_24))

params_br1_24 <- data.frame(coefs_24[1,],
                            coefs_24[2,],
                            coefs_24[3,],
                            Topt_est_24,
                            Topt_se_24,
                            boot_24[1,],
                            boot_24[2,],
                            boot_24[3,],
                            startVals_24[1],
                            startVals_24[2],
                            startVals_24[3],
                            "Tmax")
colnames(params_br1_24) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_24

#### _ _ Study 25 ####
IR_data_ID25 <- IR_data_all %>%
  filter(id==25)
par(mfrow=c(2,2))
plot(IR_data_ID25$temperature,IR_data_ID25$growth_rate)
plot(IR_data_ID25$temperature,briere1(a=0.00002,Tmin=5,Tmax=30,temp=IR_data_ID25$temperature))
plot(IR_data_ID25$temperature,briere1(a=0.00002,Tmin=5,Tmax=25,temp=IR_data_ID25$temperature))
plot(IR_data_ID25$temperature,briere1(a=0.000015,Tmin=5,Tmax=18,temp=IR_data_ID25$temperature))# <-
grid_br1_25 <- expand.grid(list(a=seq(0.00001,0.00003,by=0.000001),
                                Tmin=seq(-5,10,by=0.5),
                                Tmax=seq(15,35,by=0.5)))
fitted_br1_25 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID25,
                      start = grid_br1_25,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_25 <- summary(fitted_br1_25)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_25 <- c(sum_br1_25$coefficients[,1])
startVals_list[[25]] <- startVals_25 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_25 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID25,
                     start = list(a=0.00003,Tmin=0,Tmax=20)) #NO sale
# no sale params_br1_25

#### _ _ Study 26 ####
IR_data_ID26 <- IR_data_all %>%
  filter(id==26)
par(mfrow=c(2,2))
plot(IR_data_ID26$temperature,IR_data_ID26$growth_rate)
plot(IR_data_ID26$temperature,briere1(a=0.0001,Tmin=20,Tmax=33,temp=IR_data_ID26$temperature))
plot(IR_data_ID26$temperature,briere1(a=0.0002,Tmin=12,Tmax=36,temp=IR_data_ID26$temperature))
plot(IR_data_ID26$temperature,briere1(a=0.00015,Tmin=13,Tmax=35.5,temp=IR_data_ID26$temperature))# <-
grid_br1_26 <- expand.grid(list(a=seq(0.00001,0.0001,by=0.00001),
                                Tmin=seq(0,12,by=0.5),
                                Tmax=seq(30,40,by=0.5)))
fitted_br1_26 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID26,
                      start = grid_br1_26,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_26 <- summary(fitted_br1_26)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_26 <- c(sum_br1_26$coefficients[,1])
startVals_list[[26]] <- startVals_26 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_26 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID26,
                   start = startVals_26)
sum_br1_26 <- summary(fitted_br1_26)
sum_br1_26 
boot_br1_26 <- nlsBoot(fitted_br1_26, niter = 999)
coefs_26 <- data.frame(sum_br1_26$coefficients[,1:2],row.names = NULL)
boot_26 <- data.frame(boot_br1_26$estiboot)
Topt_est_26 <- Topt(Tmin=coefs_26[2,1],
                    Tmax=coefs_26[3,1],
                    m=2)
Topt_se_26 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_26), vcov(fitted_br1_26))

params_br1_26 <- data.frame(coefs_26[1,],
                            coefs_26[2,],
                            coefs_26[3,],
                            Topt_est_26,
                            Topt_se_26,
                            boot_26[1,],
                            boot_26[2,],
                            boot_26[3,],
                            startVals_26[1],
                            startVals_26[2],
                            startVals_26[3],
                            "Tmax")
colnames(params_br1_26) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_26

#### _ _ Study 27 ####
IR_data_ID27 <- IR_data_all %>%
  filter(id==27)
par(mfrow=c(2,2))
plot(IR_data_ID27$temperature,IR_data_ID27$growth_rate)
#no

#### _ _ Study 28 ####
IR_data_ID28 <- IR_data_all %>%
  filter(id==28)
par(mfrow=c(2,2))
plot(IR_data_ID28$temperature,IR_data_ID28$growth_rate)
plot(IR_data_ID28$temperature,briere1(a=0.0001,Tmin=15,Tmax=35,temp=IR_data_ID28$temperature))
plot(IR_data_ID28$temperature,briere1(a=0.0001,Tmin=14,Tmax=35,temp=IR_data_ID28$temperature))
grid_br1_28 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(30,40,by=0.5)))
fitted_br1_28 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID28,
                      start = grid_br1_28,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_28 <- summary(fitted_br1_28)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_28 <- c(sum_br1_28$coefficients[,1])
startVals_list[[28]] <- startVals_28 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_28 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID28,
                     start = startVals_28)
sum_br1_28 <- summary(fitted_br1_28)
sum_br1_28 
boot_br1_28 <- nlsBoot(fitted_br1_28, niter = 999)
coefs_28 <- data.frame(sum_br1_28$coefficients[,1:2],row.names = NULL)
boot_28 <- data.frame(boot_br1_28$estiboot)
Topt_est_28 <- Topt(Tmin=coefs_28[2,1],
                    Tmax=coefs_28[3,1],
                    m=2)
Topt_se_28 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_28), vcov(fitted_br1_28))

params_br1_28 <- data.frame(coefs_28[1,],
                            coefs_28[2,],
                            coefs_28[3,],
                            Topt_est_28,
                            Topt_se_28,
                            boot_28[1,],
                            boot_28[2,],
                            boot_28[3,],
                            startVals_28[1],
                            startVals_28[2],
                            startVals_28[3],
                            "Tmin,Tmax")
colnames(params_br1_28) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_28

#### _ _ Study 29 ####
IR_data_ID29 <- IR_data_all %>%
  filter(id==29)
par(mfrow=c(2,2))
plot(IR_data_ID29$temperature,IR_data_ID29$growth_rate)
#no

#### _ _ Study 30 #### 
IR_data_ID30 <- IR_data_all %>%
  filter(id==30)
par(mfrow=c(2,2))
plot(IR_data_ID30$temperature,IR_data_ID30$growth_rate)

#### _ _ Study 31 ####
IR_data_ID31 <- IR_data_all %>%
  filter(id==31)
par(mfrow=c(2,2))
plot(IR_data_ID31$temperature,IR_data_ID31$growth_rate)
plot(IR_data_ID31$temperature,briere1(a=0.0001,Tmin=15,Tmax=35,temp=IR_data_ID31$temperature))
plot(IR_data_ID31$temperature,briere1(a=0.0002,Tmin=15,Tmax=32,temp=IR_data_ID31$temperature))
plot(IR_data_ID31$temperature,briere1(a=0.00023,Tmin=12,Tmax=32,temp=IR_data_ID31$temperature)) # < --
grid_br1_31 <- expand.grid(list(a=seq(0.0001,0.00015,by=0.00001),
                                Tmin=seq(5,15,by=0.5),
                                Tmax=seq(28,36,by=0.5)))
fitted_br1_31 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID31,
                      start = grid_br1_31,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_31 <- summary(fitted_br1_31)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_31 <- c(sum_br1_31$coefficients[,1])
startVals_list[[31]] <- startVals_31 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_31 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID31,
                     start = startVals_31)
sum_br1_31 <- summary(fitted_br1_31)
sum_br1_31 
boot_br1_31 <- nlsBoot(fitted_br1_31, niter = 999)
coefs_31 <- data.frame(sum_br1_31$coefficients[,1:2],row.names = NULL)
boot_31 <- data.frame(boot_br1_31$estiboot)
Topt_est_31 <- Topt(Tmin=coefs_31[2,1],
                    Tmax=coefs_31[3,1],
                    m=2)
Topt_se_31 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_31), vcov(fitted_br1_31))

params_br1_31 <- data.frame(coefs_31[1,],
                            coefs_31[2,],
                            coefs_31[3,],
                            Topt_est_31,
                            Topt_se_31,
                            boot_31[1,],
                            boot_31[2,],
                            boot_31[3,],
                            startVals_31[1],
                            startVals_31[2],
                            startVals_31[3],
                            "Tmax")
colnames(params_br1_31) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_31

#### _ _ Study 32 ####
IR_data_ID32 <- IR_data_all %>%
  filter(id==32)
par(mfrow=c(2,2))
plot(IR_data_ID32$temperature,IR_data_ID32$growth_rate)
plot(IR_data_ID32$temperature,briere1(a=0.0001,Tmin=15,Tmax=35,temp=IR_data_ID32$temperature))
plot(IR_data_ID32$temperature,briere1(a=0.00016,Tmin=14,Tmax=37,temp=IR_data_ID32$temperature))
plot(IR_data_ID32$temperature,briere1(a=0.00013,Tmin=15,Tmax=38,temp=IR_data_ID32$temperature)) # < --
grid_br1_32 <- expand.grid(list(a=seq(0.0001,0.00015,by=0.00001),
                                Tmin=seq(14,20,by=0.5),
                                Tmax=seq(35,42,by=0.5)))
fitted_br1_32 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID32,
                      start = grid_br1_32,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_32 <- summary(fitted_br1_32)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_32 <- c(sum_br1_32$coefficients[,1])
startVals_list[[32]] <- startVals_32 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_32 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID32,
                     start = startVals_32) 
sum_br1_32 <- summary(fitted_br1_32)
sum_br1_32 
boot_br1_32 <- nlsBoot(fitted_br1_32, niter = 999)
coefs_32 <- data.frame(sum_br1_32$coefficients[,1:2],row.names = NULL)
boot_32 <- data.frame(boot_br1_32$estiboot)
Topt_est_32 <- Topt(Tmin=coefs_32[2,1],
                    Tmax=coefs_32[3,1],
                    m=2)
Topt_se_32 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_32), vcov(fitted_br1_32))
params_br1_32 <- data.frame(coefs_32[1,],
                            coefs_32[2,],
                            coefs_32[3,],
                            Topt_est_32,
                            Topt_se_32,
                            boot_32[1,],
                            boot_32[2,],
                            boot_32[3,],
                            startVals_32[1],
                            startVals_32[2],
                            startVals_32[3],
                            "all")
colnames(params_br1_32) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_32

#### _ _ Study 33 ####
IR_data_ID33 <- IR_data_all %>%
  filter(id==33)
par(mfrow=c(2,2))
plot(IR_data_ID33$temperature,IR_data_ID33$growth_rate)
plot(IR_data_ID33$temperature,briere1(a=0.0001,Tmin=15,Tmax=35,temp=IR_data_ID33$temperature))
plot(IR_data_ID33$temperature,briere1(a=0.0001,Tmin=18,Tmax=28,temp=IR_data_ID33$temperature))
plot(IR_data_ID33$temperature,briere1(a=0.00013,Tmin=15,Tmax=26,temp=IR_data_ID33$temperature)) # < --
grid_br1_33 <- expand.grid(list(a=seq(0.0001,0.0002,by=0.00001),
                                Tmin=seq(12,20,by=0.5),
                                Tmax=seq(23,30,by=0.5)))
fitted_br1_33 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID33,
                      start = grid_br1_33,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_33 <- summary(fitted_br1_33)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_33 <- c(sum_br1_33$coefficients[,1])
startVals_list[[33]] <- startVals_33 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_33 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID33,
                     start = startVals_33)
sum_br1_33 <- summary(fitted_br1_33)
sum_br1_33 
boot_br1_33 <- nlsBoot(fitted_br1_33, niter = 999)
coefs_33 <- data.frame(sum_br1_33$coefficients[,1:2],row.names = NULL)
boot_33 <- data.frame(boot_br1_33$estiboot)
Topt_est_33 <- Topt(Tmin=coefs_33[2,1],
                    Tmax=coefs_33[3,1],
                    m=2)
Topt_se_33 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_33), vcov(fitted_br1_33))


params_br1_33 <- data.frame(coefs_33[1,],
                            coefs_33[2,],
                            coefs_33[3,],
                            Topt_est_33,
                            Topt_se_33,
                            boot_33[1,],
                            boot_33[2,],
                            boot_33[3,],
                            startVals_33[1],
                            startVals_33[2],
                            startVals_33[3],
                            "all")
colnames(params_br1_33) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")


#### _ _ Study 34 ####
IR_data_ID34 <- IR_data_all %>%
  filter(id==34)
par(mfrow=c(2,2))
plot(IR_data_ID34$temperature,IR_data_ID34$growth_rate)
plot(IR_data_ID34$temperature,briere1(a=0.0001,Tmin=19,Tmax=31,temp=IR_data_ID34$temperature))
plot(IR_data_ID33$temperature,briere1(a=0.0002,Tmin=19,Tmax=30.5,temp=IR_data_ID34$temperature))
plot(IR_data_ID33$temperature,briere1(a=0.0002,Tmin=19,Tmax=30,temp=IR_data_ID34$temperature)) # < --
grid_br1_34 <- expand.grid(list(a=seq(0.0002,0.0003,by=0.00001),
                                Tmin=seq(13,23,by=0.5),
                                Tmax=seq(25,35,by=0.5)))
fitted_br1_34 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID34,
                      start = grid_br1_34,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_34 <- summary(fitted_br1_34)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_34 <- c(sum_br1_34$coefficients[,1])
startVals_list[[34]] <- startVals_34 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_34 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID34,
                     start = startVals_34)
sum_br1_34 <- summary(fitted_br1_34)
sum_br1_34 
boot_br1_34 <- nlsBoot(fitted_br1_34, niter = 999)
coefs_34 <- data.frame(sum_br1_34$coefficients[,1:2],row.names = NULL)
boot_34 <- data.frame(boot_br1_34$estiboot)
Topt_est_34 <- Topt(Tmin=coefs_34[2,1],
                    Tmax=coefs_34[3,1],
                    m=2)
Topt_se_34 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_34), vcov(fitted_br1_34))


params_br1_34 <- data.frame(coefs_34[1,],
                            coefs_34[2,],
                            coefs_34[3,],
                            Topt_est_34,
                            Topt_se_34,
                            boot_34[1,],
                            boot_34[2,],
                            boot_34[3,],
                            startVals_34[1],
                            startVals_34[2],
                            startVals_34[3],
                            "Tmin,Tmax")
colnames(params_br1_34) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_34

#### _ _ Study 35 ####
IR_data_ID35 <- IR_data_all %>%
  filter(id==35)
par(mfrow=c(2,2))
plot(IR_data_ID35$temperature,IR_data_ID35$growth_rate)
plot(IR_data_ID35$temperature,briere1(a=0.0001,Tmin=19,Tmax=31,temp=IR_data_ID35$temperature))
plot(IR_data_ID35$temperature,briere1(a=0.0002,Tmin=16,Tmax=37,temp=IR_data_ID35$temperature))
plot(IR_data_ID35$temperature,briere1(a=0.0001,Tmin=16,Tmax=37,temp=IR_data_ID35$temperature)) # < --
grid_br1_35 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(13,23,by=0.5),
                                Tmax=seq(34,40,by=0.5)))
fitted_br1_35 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID35,
                      start = grid_br1_35,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_35 <- summary(fitted_br1_35)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_35 <- c(sum_br1_35$coefficients[,1])
startVals_list[[35]] <- startVals_35 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_35 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID35,
                     start = startVals_35)
sum_br1_35 <- summary(fitted_br1_35)
sum_br1_35 
boot_br1_35 <- nlsBoot(fitted_br1_35, niter = 999)
coefs_35 <- data.frame(sum_br1_35$coefficients[,1:2],row.names = NULL)
boot_35 <- data.frame(boot_br1_35$estiboot)
Topt_est_35 <- Topt(Tmin=coefs_35[2,1],
                    Tmax=coefs_35[3,1],
                    m=2)
Topt_se_35 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_35), vcov(fitted_br1_35))


params_br1_35 <- data.frame(coefs_35[1,],
                            coefs_35[2,],
                            coefs_35[3,],
                            Topt_est_35,
                            Topt_se_35,
                            boot_35[1,],
                            boot_35[2,],
                            boot_35[3,],
                            startVals_35[1],
                            startVals_35[2],
                            startVals_35[3],
                            "Tmax")
colnames(params_br1_35) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_35

#### _ _ Study 36 ####
IR_data_ID36 <- IR_data_all %>%
  filter(id==36)
par(mfrow=c(2,2))
plot(IR_data_ID36$temperature,IR_data_ID36$growth_rate)
plot(IR_data_ID36$temperature,briere1(a=0.0001,Tmin=19,Tmax=31,temp=IR_data_ID36$temperature))
plot(IR_data_ID36$temperature,briere1(a=0.0002,Tmin=8,Tmax=32,temp=IR_data_ID36$temperature))
plot(IR_data_ID36$temperature,briere1(a=0.00007,Tmin=5,Tmax=31,temp=IR_data_ID36$temperature))
grid_br1_36 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(0,15,by=0.5),
                                Tmax=seq(25,35,by=0.5)))
fitted_br1_36 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID36,
                      start = grid_br1_36,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_36 <- summary(fitted_br1_36)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_36 <- c(sum_br1_36$coefficients[,1])
startVals_list[[36]] <- startVals_36 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_36 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID36,
                     start = startVals_36)
sum_br1_36 <- summary(fitted_br1_36)
sum_br1_36 
boot_br1_36 <- nlsBoot(fitted_br1_36, niter = 999)
coefs_36 <- data.frame(sum_br1_36$coefficients[,1:2],row.names = NULL)
boot_36 <- data.frame(boot_br1_36$estiboot)
Topt_est_36 <- Topt(Tmin=coefs_36[2,1],
                    Tmax=coefs_36[3,1],
                    m=2)
Topt_se_36 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_36), vcov(fitted_br1_36))


params_br1_36 <- data.frame(coefs_36[1,],
                            coefs_36[2,],
                            coefs_36[3,],
                            Topt_est_36,
                            Topt_se_36,
                            boot_36[1,],
                            boot_36[2,],
                            boot_36[3,],
                            startVals_36[1],
                            startVals_36[2],
                            startVals_36[3],
                            "a,Tmax")
colnames(params_br1_36) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_36

#### _ _ Study 37 ####
IR_data_ID37 <- IR_data_all %>%
  filter(id==37)
par(mfrow=c(2,2))
plot(IR_data_ID37$temperature,IR_data_ID37$growth_rate)
plot(IR_data_ID37$temperature,briere1(a=0.0001,Tmin=10,Tmax=31,temp=IR_data_ID37$temperature))
plot(IR_data_ID37$temperature,briere1(a=0.0001,Tmin=10,Tmax=42,temp=IR_data_ID37$temperature))
plot(IR_data_ID37$temperature,briere1(a=0.00009,Tmin=10,Tmax=42,temp=IR_data_ID37$temperature))
grid_br1_37 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(0,15,by=0.5),
                                Tmax=seq(35,48,by=0.5)))
fitted_br1_37 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID37,
                      start = grid_br1_37,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_37 <- summary(fitted_br1_37)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_37 <- c(sum_br1_37$coefficients[,1])
startVals_list[[37]] <- startVals_37 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_37 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID37,
                     start = startVals_37) #no
sum_br1_37 <- summary(fitted_br1_37)
sum_br1_37 
boot_br1_37 <- nlsBoot(fitted_br1_37, niter = 999)
coefs_37 <- data.frame(sum_br1_37$coefficients[,1:2],row.names = NULL)
boot_37 <- data.frame(boot_br1_37$estiboot)
Topt_est_37 <- Topt(Tmin=coefs_37[2,1],
                    Tmax=coefs_37[3,1],
                    m=2)
Topt_se_37 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_37), vcov(fitted_br1_37))


params_br1_37 <- data.frame(coefs_37[1,],
                            coefs_37[2,],
                            coefs_37[3,],
                            Topt_est_37,
                            Topt_se_37,
                            boot_37[1,],
                            boot_37[2,],
                            boot_37[3,],
                            startVals_37[1],
                            startVals_37[2],
                            startVals_37[3],
                            "none")
colnames(params_br1_37) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_37

#### _ _ Study 38 ####
IR_data_ID38 <- IR_data_all %>%
  filter(id==38)
par(mfrow=c(2,2))
plot(IR_data_ID38$temperature,IR_data_ID38$growth_rate)
plot(IR_data_ID38$temperature,briere1(a=0.0001,Tmin=10,Tmax=31,temp=IR_data_ID38$temperature))
plot(IR_data_ID38$temperature,briere1(a=0.0001,Tmin=13,Tmax=26,temp=IR_data_ID38$temperature))
plot(IR_data_ID38$temperature,briere1(a=0.0001,Tmin=13,Tmax=25,temp=IR_data_ID38$temperature))
grid_br1_38 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(8,18,by=0.5),
                                Tmax=seq(20,32,by=0.5)))
fitted_br1_38 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID38,
                      start = grid_br1_38,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_38 <- summary(fitted_br1_38)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_38 <- c(sum_br1_38$coefficients[,1])
startVals_list[[38]] <- startVals_38 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_38 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID38,
                     start = startVals_38)
sum_br1_38 <- summary(fitted_br1_38)
sum_br1_38 
boot_br1_38 <- nlsBoot(fitted_br1_38, niter = 999)
coefs_38 <- data.frame(sum_br1_38$coefficients[,1:2],row.names = NULL)
boot_38 <- data.frame(boot_br1_38$estiboot)
Topt_est_38 <- Topt(Tmin=coefs_38[2,1],
                    Tmax=coefs_38[3,1],
                    m=2)
Topt_se_38 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_38), vcov(fitted_br1_38))


params_br1_38 <- data.frame(coefs_38[1,],
                            coefs_38[2,],
                            coefs_38[3,],
                            Topt_est_38,
                            Topt_se_38,
                            boot_38[1,],
                            boot_38[2,],
                            boot_38[3,],
                            startVals_38[1],
                            startVals_38[2],
                            startVals_38[3],
                            "all")
colnames(params_br1_38) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_38

#### _ _ Study 39 ####
IR_data_ID39 <- IR_data_all %>%
  filter(id==39)
par(mfrow=c(2,2))
plot(IR_data_ID39$temperature,IR_data_ID39$growth_rate)
plot(IR_data_ID39$temperature,briere1(a=0.0001,Tmin=8,Tmax=27,temp=IR_data_ID39$temperature))
plot(IR_data_ID39$temperature,briere1(a=0.0001,Tmin=8,Tmax=26,temp=IR_data_ID39$temperature))
plot(IR_data_ID39$temperature,briere1(a=0.0001,Tmin=8,Tmax=25.5,temp=IR_data_ID39$temperature))
grid_br1_39 <- expand.grid(list(a=seq(0.00001,0.00007,by=0.00001),
                                Tmin=seq(0,15,by=0.5),
                                Tmax=seq(20,30,by=0.5)))
fitted_br1_39 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID39,
                      start = grid_br1_39,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_39 <- summary(fitted_br1_39)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_39 <- c(sum_br1_39$coefficients[,1])
startVals_list[[39]] <- startVals_39 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_39 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID39,
                     start = startVals_39)
sum_br1_39 <- summary(fitted_br1_39)
sum_br1_39 
boot_br1_39 <- nlsBoot(fitted_br1_39, niter = 999)
coefs_39 <- data.frame(sum_br1_39$coefficients[,1:2],row.names = NULL)
boot_39 <- data.frame(boot_br1_39$estiboot)
Topt_est_39 <- Topt(Tmin=coefs_39[2,1],
                    Tmax=coefs_39[3,1],
                    m=2)
Topt_se_39 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_39), vcov(fitted_br1_39))


params_br1_39 <- data.frame(coefs_39[1,],
                            coefs_39[2,],
                            coefs_39[3,],
                            Topt_est_39,
                            Topt_se_39,
                            boot_39[1,],
                            boot_39[2,],
                            boot_39[3,],
                            startVals_39[1],
                            startVals_39[2],
                            startVals_39[3],
                            "a,Tmax")
colnames(params_br1_39) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_39

#### _ _ Study 40 ####
IR_data_ID40 <- IR_data_all %>%
  filter(id==40)
par(mfrow=c(2,2))
plot(IR_data_ID40$temperature,IR_data_ID40$growth_rate)
plot(IR_data_ID40$temperature,briere1(a=0.0001,Tmin=9,Tmax=33,temp=IR_data_ID40$temperature))
plot(IR_data_ID40$temperature,briere1(a=0.0002,Tmin=9,Tmax=33,temp=IR_data_ID40$temperature))
grid_br1_40 <- expand.grid(list(a=seq(0.00015,0.00025,by=0.00001),
                                Tmin=seq(5,15,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_40 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID40,
                      start = grid_br1_40,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_40 <- summary(fitted_br1_40)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_40 <- c(sum_br1_40$coefficients[,1])
startVals_list[[40]] <- startVals_40 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_40 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID40,
                     start = startVals_40)
sum_br1_40 <- summary(fitted_br1_40)
sum_br1_40 
boot_br1_40 <- nlsBoot(fitted_br1_40, niter = 999)
coefs_40 <- data.frame(sum_br1_40$coefficients[,1:2],row.names = NULL)
boot_40 <- data.frame(boot_br1_40$estiboot)
Topt_est_40 <- Topt(Tmin=coefs_40[2,1],
                    Tmax=coefs_40[3,1],
                    m=2)
Topt_se_40 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_40), vcov(fitted_br1_40))


params_br1_40 <- data.frame(coefs_40[1,],
                            coefs_40[2,],
                            coefs_40[3,],
                            Topt_est_40,
                            Topt_se_40,
                            boot_40[1,],
                            boot_40[2,],
                            boot_40[3,],
                            startVals_40[1],
                            startVals_40[2],
                            startVals_40[3],
                            "all")
colnames(params_br1_40) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_40

#### _ _ Study 41 ####
IR_data_ID41 <- IR_data_all %>%
  filter(id==41)
par(mfrow=c(2,2))
plot(IR_data_ID41$temperature,IR_data_ID41$growth_rate)
plot(IR_data_ID41$temperature,briere1(a=0.0001,Tmin=13,Tmax=33,temp=IR_data_ID41$temperature))
plot(IR_data_ID41$temperature,briere1(a=0.00012,Tmin=15,Tmax=34,temp=IR_data_ID41$temperature))
grid_br1_41 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(8,18,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_41 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID41,
                      start = grid_br1_41,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_41 <- summary(fitted_br1_41)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_41 <- c(sum_br1_41$coefficients[,1])
startVals_list[[41]] <- startVals_41 # viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_41 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID41,
                     start = startVals_41)
sum_br1_41 <- summary(fitted_br1_41)
sum_br1_41 
boot_br1_41 <- nlsBoot(fitted_br1_41, niter = 999)
coefs_41 <- data.frame(sum_br1_41$coefficients[,1:2],row.names = NULL)
boot_41 <- data.frame(boot_br1_41$estiboot)
Topt_est_41 <- Topt(Tmin=coefs_41[2,1],
                    Tmax=coefs_41[3,1],
                    m=2)
Topt_se_41 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_41), vcov(fitted_br1_41))


params_br1_41 <- data.frame(coefs_41[1,],
                            coefs_41[2,],
                            coefs_41[3,],
                            Topt_est_41,
                            Topt_se_41,
                            boot_41[1,],
                            boot_41[2,],
                            boot_41[3,],
                            startVals_41[1],
                            startVals_41[2],
                            startVals_41[3],
                            "Tmax")
colnames(params_br1_41) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_41

#### _ _ Study 42 #### NOOO
IR_data_ID42 <- IR_data_all %>%
  filter(id==42)
par(mfrow=c(2,2))
plot(IR_data_ID42$temperature,IR_data_ID42$growth_rate)

#### _ _ Study 43 ####
IR_data_ID43 <- IR_data_all %>%
  filter(id==43)
par(mfrow=c(2,2))
plot(IR_data_ID43$temperature,IR_data_ID43$growth_rate)
plot(IR_data_ID43$temperature,briere1(a=0.0002,Tmin=14,Tmax=38,temp=IR_data_ID43$temperature))
plot(IR_data_ID43$temperature,briere1(a=0.00017,Tmin=14,Tmax=38,temp=IR_data_ID43$temperature))
grid_br1_43 <- expand.grid(list(a=seq(0.0001,0.0002,by=0.00001),
                                Tmin=seq(5,15,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_43 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID43,
                      start = grid_br1_43,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_43 <- summary(fitted_br1_43)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_43 <- c(sum_br1_43$coefficients[,1])
startVals_list[[43]] <- startVals_43 #llenar lista starting values
fitted_br1_43 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID43,
                     start = startVals_43)
sum_br1_43 <- summary(fitted_br1_43)
sum_br1_43 
boot_br1_43 <- nlsBoot(fitted_br1_43, niter = 999)
coefs_43 <- data.frame(sum_br1_43$coefficients[,1:2],row.names = NULL)
boot_43 <- data.frame(boot_br1_43$estiboot)
Topt_est_43 <- Topt(Tmin=coefs_43[2,1],
                    Tmax=coefs_43[3,1],
                    m=2)
Topt_se_43 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_43), vcov(fitted_br1_43))

params_br1_43 <- data.frame(coefs_43[1,],
                            coefs_43[2,],
                            coefs_43[3,],
                            Topt_est_43,
                            Topt_se_43,
                            boot_43[1,],
                            boot_43[2,],
                            boot_43[3,],
                            startVals_43[1],
                            startVals_43[2],
                            startVals_43[3],
                            "all")
colnames(params_br1_43) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_43

#### _ _ Study 44 ####
IR_data_ID44 <- IR_data_all %>%
  filter(id==44)
par(mfrow=c(2,2))
plot(IR_data_ID44$temperature,IR_data_ID44$growth_rate)
plot(IR_data_ID44$temperature,briere1(a=0.0001,Tmin=5,Tmax=30,temp=IR_data_ID44$temperature))
plot(IR_data_ID44$temperature,briere1(a=0.00001,Tmin=0,Tmax=30,temp=IR_data_ID44$temperature))
grid_br1_44 <- expand.grid(list(a=seq(0.00001,0.00002,by=0.000001),
                                Tmin=seq(-5,10,by=0.5),
                                Tmax=seq(22,32,by=0.5)))
fitted_br1_44 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID44,
                      start = grid_br1_44,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_44 <- summary(fitted_br1_44)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_44 <- c(sum_br1_44$coefficients[,1])
startVals_list[[44]] <- startVals_44 #llenar lista starting values
fitted_br1_44 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID44,
                     start = startVals_44) #No sale
## no sale, eliminamos


#### _ _ Study 45 ####
IR_data_ID45 <- IR_data_all %>%
  filter(id==45)
par(mfrow=c(2,2))
plot(IR_data_ID45$temperature,IR_data_ID45$growth_rate)
plot(IR_data_ID45$temperature,briere1(a=0.0002,Tmin=16,Tmax=35,temp=IR_data_ID45$temperature))
plot(IR_data_ID45$temperature,briere1(a=0.00017,Tmin=16,Tmax=35,temp=IR_data_ID45$temperature))
plot(IR_data_ID45$temperature,briere1(a=0.00013,Tmin=16,Tmax=35,temp=IR_data_ID45$temperature))
grid_br1_45 <- expand.grid(list(a=seq(0.0001,0.0002,by=0.00001),
                                Tmin=seq(12,22,by=0.5),
                                Tmax=seq(30,40,by=0.5)))
fitted_br1_45 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID45,
                      start = grid_br1_45,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_45 <- summary(fitted_br1_45)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_45 <- c(sum_br1_45$coefficients[,1])
startVals_list[[45]] <- startVals_45 #llenar lista starting values
fitted_br1_45 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID45,
                     start = startVals_45)
sum_br1_45 <- summary(fitted_br1_45)
sum_br1_45 
boot_br1_45 <- nlsBoot(fitted_br1_45, niter = 999)
coefs_45 <- data.frame(sum_br1_45$coefficients[,1:2],row.names = NULL)
boot_45 <- data.frame(boot_br1_45$estiboot)
Topt_est_45 <- Topt(Tmin=coefs_45[2,1],
                    Tmax=coefs_45[3,1],
                    m=2)
Topt_se_45 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_45), vcov(fitted_br1_45))

params_br1_45 <- data.frame(coefs_45[1,],
                            coefs_45[2,],
                            coefs_45[3,],
                            Topt_est_45,
                            Topt_se_45,
                            boot_45[1,],
                            boot_45[2,],
                            boot_45[3,],
                            startVals_45[1],
                            startVals_45[2],
                            startVals_45[3],
                            "all")
colnames(params_br1_45) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_45


#### _ _ Study 47 ####
IR_data_ID46 <- IR_data_all %>%
  filter(id==46)
par(mfrow=c(2,2))
plot(IR_data_ID46$temperature,IR_data_ID46$growth_rate)
plot(IR_data_ID46$temperature,briere1(a=0.0002,Tmin=16,Tmax=35,temp=IR_data_ID46$temperature))
plot(IR_data_ID46$temperature,briere1(a=0.00002,Tmin=16,Tmax=35.5,temp=IR_data_ID46$temperature))
plot(IR_data_ID46$temperature,briere1(a=0.00001,Tmin=16,Tmax=35.5,temp=IR_data_ID46$temperature))
grid_br1_46 <- expand.grid(list(a=seq(0.00001,0.00002,by=0.00001),
                                Tmin=seq(12,22,by=0.5),
                                Tmax=seq(30,40,by=0.5)))
fitted_br1_46 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID46,
                      start = grid_br1_46,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_46 <- summary(fitted_br1_46)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_46 <- c(sum_br1_46$coefficients[,1])
startVals_list[[46]] <- startVals_46 #llenar lista starting values
fitted_br1_46 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID46,
                     start = startVals_46)
sum_br1_46 <- summary(fitted_br1_46)
sum_br1_46 
boot_br1_46 <- nlsBoot(fitted_br1_46, niter = 999)
coefs_46 <- data.frame(sum_br1_46$coefficients[,1:2],row.names = NULL)
boot_46 <- data.frame(boot_br1_46$estiboot)
Topt_est_46 <- Topt(Tmin=coefs_46[2,1],
                    Tmax=coefs_46[3,1],
                    m=2)
Topt_se_46 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_46), vcov(fitted_br1_46))

params_br1_46 <- data.frame(coefs_46[1,],
                            coefs_46[2,],
                            coefs_46[3,],
                            Topt_est_46,
                            Topt_se_46,
                            boot_46[1,],
                            boot_46[2,],
                            boot_46[3,],
                            startVals_46[1],
                            startVals_46[2],
                            startVals_46[3],
                            "all")
colnames(params_br1_46) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_46

#### _ _ Study 47 ####
IR_data_ID47 <- IR_data_all %>%
  filter(id==47)
par(mfrow=c(2,2))
plot(IR_data_ID47$temperature,IR_data_ID47$growth_rate)
plot(IR_data_ID47$temperature,briere1(a=0.0002,Tmin=16,Tmax=38,temp=IR_data_ID47$temperature))
plot(IR_data_ID47$temperature,briere1(a=0.0001,Tmin=16,Tmax=38,temp=IR_data_ID47$temperature))
plot(IR_data_ID47$temperature,briere1(a=0.00007,Tmin=16,Tmax=38,temp=IR_data_ID47$temperature))
grid_br1_47 <- expand.grid(list(a=seq(0.00001,0.0001,by=0.00001),
                                Tmin=seq(12,22,by=0.5),
                                Tmax=seq(30,40,by=0.5)))
fitted_br1_47 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID47,
                      start = grid_br1_47,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_47 <- summary(fitted_br1_47)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_47 <- c(sum_br1_47$coefficients[,1])
startVals_list[[47]] <- startVals_47 #llenar lista starting values
fitted_br1_47 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID47,
                     start = startVals_47)
sum_br1_47 <- summary(fitted_br1_47)
sum_br1_47 
boot_br1_47 <- nlsBoot(fitted_br1_47, niter = 999)
coefs_47 <- data.frame(sum_br1_47$coefficients[,1:2],row.names = NULL)
boot_47 <- data.frame(boot_br1_47$estiboot)
Topt_est_47 <- Topt(Tmin=coefs_47[2,1],
                    Tmax=coefs_47[3,1],
                    m=2)
Topt_se_47 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_47), vcov(fitted_br1_47))

params_br1_47 <- data.frame(coefs_47[1,],
                            coefs_47[2,],
                            coefs_47[3,],
                            Topt_est_47,
                            Topt_se_47,
                            boot_47[1,],
                            boot_47[2,],
                            boot_47[3,],
                            startVals_47[1],
                            startVals_47[2],
                            startVals_47[3],
                            "Tmin,Tmax")
colnames(params_br1_47) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_47

#### _ _ Study 48 ####
IR_data_ID48 <- IR_data_all %>%
  filter(id==48)
par(mfrow=c(2,2))
plot(IR_data_ID48$temperature,IR_data_ID48$growth_rate)
plot(IR_data_ID48$temperature,briere1(a=0.00002,Tmin=16,Tmax=29,temp=IR_data_ID48$temperature))
plot(IR_data_ID48$temperature,briere1(a=0.0001,Tmin=16,Tmax=27,temp=IR_data_ID48$temperature))
plot(IR_data_ID48$temperature,briere1(a=0.00011,Tmin=16,Tmax=27,temp=IR_data_ID48$temperature))
grid_br1_48 <- expand.grid(list(a=seq(0.0001,0.0002,by=0.00001),
                                Tmin=seq(12,22,by=0.5),
                                Tmax=seq(22,32,by=0.5)))
fitted_br1_48 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID48,
                      start = grid_br1_48,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_48 <- summary(fitted_br1_48)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_48 <- c(sum_br1_48$coefficients[,1])
startVals_list[[48]] <- startVals_48 #llenar lista starting values
fitted_br1_48 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID48,
                     start = startVals_48)
sum_br1_48 <- summary(fitted_br1_48)
sum_br1_48 
boot_br1_48 <- nlsBoot(fitted_br1_48, niter = 999)
coefs_48 <- data.frame(sum_br1_48$coefficients[,1:2],row.names = NULL)
boot_48 <- data.frame(boot_br1_48$estiboot)
Topt_est_48 <- Topt(Tmin=coefs_48[2,1],
                    Tmax=coefs_48[3,1],
                    m=2)
Topt_se_48 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_48), vcov(fitted_br1_48))

params_br1_48 <- data.frame(coefs_48[1,],
                            coefs_48[2,],
                            coefs_48[3,],
                            Topt_est_48,
                            Topt_se_48,
                            boot_48[1,],
                            boot_48[2,],
                            boot_48[3,],
                            startVals_48[1],
                            startVals_48[2],
                            startVals_48[3],
                            "Tmax")
colnames(params_br1_48) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_48

#### _ _ Study 49 ####
IR_data_ID49 <- IR_data_all %>%
  filter(id==49)
par(mfrow=c(2,2))
plot(IR_data_ID49$temperature,IR_data_ID49$growth_rate)
plot(IR_data_ID49$temperature,briere1(a=0.0001,Tmin=16,Tmax=42,temp=IR_data_ID49$temperature))
plot(IR_data_ID49$temperature,briere1(a=0.0002,Tmin=16,Tmax=42,temp=IR_data_ID49$temperature))
plot(IR_data_ID49$temperature,briere1(a=0.0001,Tmin=10,Tmax=42,temp=IR_data_ID49$temperature))
grid_br1_49 <- expand.grid(list(a=seq(0.00008,0.0002,by=0.00001),
                                Tmin=seq(5,18,by=0.5),
                                Tmax=seq(34,48,by=0.5)))
fitted_br1_49 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID49,
                      start = grid_br1_49,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_49 <- summary(fitted_br1_49)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_49 <- c(sum_br1_49$coefficients[,1])
startVals_list[[49]] <- startVals_49 #llenar lista starting values
fitted_br1_49 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID49,
                     start = startVals_49)
sum_br1_49 <- summary(fitted_br1_49)
sum_br1_49 
boot_br1_49 <- nlsBoot(fitted_br1_49, niter = 999)
coefs_49 <- data.frame(sum_br1_49$coefficients[,1:2],row.names = NULL)
boot_49 <- data.frame(boot_br1_49$estiboot)
Topt_est_49 <- Topt(Tmin=coefs_49[2,1],
                    Tmax=coefs_49[3,1],
                    m=2)
Topt_se_49 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_49), vcov(fitted_br1_49))

params_br1_49 <- data.frame(coefs_49[1,],
                            coefs_49[2,],
                            coefs_49[3,],
                            Topt_est_49,
                            Topt_se_49,
                            boot_49[1,],
                            boot_49[2,],
                            boot_49[3,],
                            startVals_49[1],
                            startVals_49[2],
                            startVals_49[3],
                            "a,Tmax")
colnames(params_br1_49) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_49

#### _ _ Study 50 ####
IR_data_ID50 <- IR_data_all %>%
  filter(id==50)
par(mfrow=c(2,2))
plot(IR_data_ID50$temperature,IR_data_ID50$growth_rate)
plot(IR_data_ID50$temperature,briere1(a=0.0001,Tmin=14,Tmax=32,temp=IR_data_ID50$temperature))
plot(IR_data_ID50$temperature,briere1(a=0.0004,Tmin=14,Tmax=32,temp=IR_data_ID50$temperature))
plot(IR_data_ID50$temperature,briere1(a=0.0002,Tmin=14,Tmax=32,temp=IR_data_ID50$temperature))
grid_br1_50 <- expand.grid(list(a=seq(0.0001,0.0009,by=0.0001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_50 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID50,
                      start = grid_br1_50,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_50 <- summary(fitted_br1_50)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_50 <- c(sum_br1_50$coefficients[,1])
startVals_list[[50]] <- startVals_50 #llenar lista starting values
fitted_br1_50 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID50,
                     start = startVals_50)
sum_br1_50 <- summary(fitted_br1_50)
sum_br1_50 
boot_br1_50 <- nlsBoot(fitted_br1_50, niter = 999)
coefs_50 <- data.frame(sum_br1_50$coefficients[,1:2],row.names = NULL)
boot_50 <- data.frame(boot_br1_50$estiboot)
Topt_est_50 <- Topt(Tmin=coefs_50[2,1],
                    Tmax=coefs_50[3,1],
                    m=2)
Topt_se_50 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_50), vcov(fitted_br1_50))

params_br1_50 <- data.frame(coefs_50[1,],
                            coefs_50[2,],
                            coefs_50[3,],
                            Topt_est_50,
                            Topt_se_50,
                            boot_50[1,],
                            boot_50[2,],
                            boot_50[3,],
                            startVals_50[1],
                            startVals_50[2],
                            startVals_50[3],
                            "all")
colnames(params_br1_50) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_50

#### _ _ Study 51 ####
IR_data_ID51 <- IR_data_all %>%
  filter(id==51)
par(mfrow=c(2,2))
plot(IR_data_ID51$temperature,IR_data_ID51$growth_rate)
plot(IR_data_ID51$temperature,briere1(a=0.0001,Tmin=22,Tmax=32.5,temp=IR_data_ID51$temperature))
plot(IR_data_ID51$temperature,briere1(a=0.00015,Tmin=22,Tmax=32.5,temp=IR_data_ID51$temperature))
plot(IR_data_ID51$temperature,briere1(a=0.0002,Tmin=22,Tmax=32.5,temp=IR_data_ID51$temperature))
grid_br1_51 <- expand.grid(list(a=seq(0.0001,0.00025,by=0.00001),
                                Tmin=seq(15,25,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_51 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID51,
                      start = grid_br1_51,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_51 <- summary(fitted_br1_51)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_51 <- c(sum_br1_51$coefficients[,1])
startVals_list[[51]] <- startVals_51 #llenar lista starting values
fitted_br1_51 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID51,
                     start = startVals_51)
sum_br1_51 <- summary(fitted_br1_51)
sum_br1_51 
boot_br1_51 <- nlsBoot(fitted_br1_51, niter = 999)
coefs_51 <- data.frame(sum_br1_51$coefficients[,1:2],row.names = NULL)
boot_51 <- data.frame(boot_br1_51$estiboot)
Topt_est_51 <- Topt(Tmin=coefs_51[2,1],
                    Tmax=coefs_51[3,1],
                    m=2)
Topt_se_51 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_51), vcov(fitted_br1_51))

params_br1_51 <- data.frame(coefs_51[1,],
                            coefs_51[2,],
                            coefs_51[3,],
                            Topt_est_51,
                            Topt_se_51,
                            boot_51[1,],
                            boot_51[2,],
                            boot_51[3,],
                            startVals_51[1],
                            startVals_51[2],
                            startVals_51[3],
                            "all")
colnames(params_br1_51) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_51

#### _ _ Study 52 ####
IR_data_ID52 <- IR_data_all %>%
  filter(id==52)
par(mfrow=c(2,2))
plot(IR_data_ID52$temperature,IR_data_ID52$growth_rate)
plot(IR_data_ID52$temperature,briere1(a=0.0001,Tmin=16,Tmax=36,temp=IR_data_ID52$temperature))
plot(IR_data_ID52$temperature,briere1(a=0.00025,Tmin=16,Tmax=36,temp=IR_data_ID52$temperature))
plot(IR_data_ID52$temperature,briere1(a=0.00026,Tmin=18,Tmax=33,temp=IR_data_ID52$temperature))
grid_br1_52 <- expand.grid(list(a=seq(0.0002,0.0003,by=0.00001),
                                Tmin=seq(13,23,by=0.5),
                                Tmax=seq(30,43,by=0.5)))
fitted_br1_52 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID52,
                      start = grid_br1_52,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_52 <- summary(fitted_br1_52)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_52 <- c(sum_br1_52$coefficients[,1])
startVals_list[[52]] <- startVals_52 #llenar lista starting values
fitted_br1_52 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID52,
                     start = startVals_52)
sum_br1_52 <- summary(fitted_br1_52)
sum_br1_52 
boot_br1_52 <- nlsBoot(fitted_br1_52, niter = 999)
coefs_52 <- data.frame(sum_br1_52$coefficients[,1:2],row.names = NULL)
boot_52 <- data.frame(boot_br1_52$estiboot)
Topt_est_52 <- Topt(Tmin=coefs_52[2,1],
                    Tmax=coefs_52[3,1],
                    m=2)
Topt_se_52 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_52), vcov(fitted_br1_52))

params_br1_52 <- data.frame(coefs_52[1,],
                            coefs_52[2,],
                            coefs_52[3,],
                            Topt_est_52,
                            Topt_se_52,
                            boot_52[1,],
                            boot_52[2,],
                            boot_52[3,],
                            startVals_52[1],
                            startVals_52[2],
                            startVals_52[3],
                            "all")
colnames(params_br1_52) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_52

#### _ _ Study 53 ####
IR_data_ID53 <- IR_data_all %>%
  filter(id==53)
par(mfrow=c(2,2))
plot(IR_data_ID53$temperature,IR_data_ID53$growth_rate)
plot(IR_data_ID53$temperature,briere1(a=0.0002,Tmin=14,Tmax=37,temp=IR_data_ID53$temperature))
plot(IR_data_ID53$temperature,briere1(a=0.0004,Tmin=14,Tmax=37,temp=IR_data_ID53$temperature))
plot(IR_data_ID53$temperature,briere1(a=0.00035,Tmin=14,Tmax=36.5,temp=IR_data_ID53$temperature))
grid_br1_53 <- expand.grid(list(a=seq(0.0002,0.0004,by=0.00001),
                                Tmin=seq(6,20,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_53 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID53,
                      start = grid_br1_53,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_53 <- summary(fitted_br1_53)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_53 <- c(sum_br1_53$coefficients[,1])
startVals_list[[53]] <- startVals_53 #llenar lista starting values
fitted_br1_53 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID53,
                     start = startVals_53)
sum_br1_53 <- summary(fitted_br1_53)
sum_br1_53 
boot_br1_53 <- nlsBoot(fitted_br1_53, niter = 999)
coefs_53 <- data.frame(sum_br1_53$coefficients[,1:2],row.names = NULL)
boot_53 <- data.frame(boot_br1_53$estiboot)
Topt_est_53 <- Topt(Tmin=coefs_53[2,1],
                    Tmax=coefs_53[3,1],
                    m=2)
Topt_se_53 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_53), vcov(fitted_br1_53))

params_br1_53 <- data.frame(coefs_53[1,],
                            coefs_53[2,],
                            coefs_53[3,],
                            Topt_est_53,
                            Topt_se_53,
                            boot_53[1,],
                            boot_53[2,],
                            boot_53[3,],
                            startVals_53[1],
                            startVals_53[2],
                            startVals_53[3],
                            "all")
colnames(params_br1_53) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_53

#### _ _ Study 54 ####
IR_data_ID54 <- IR_data_all %>%
  filter(id==54)
par(mfrow=c(2,2))
plot(IR_data_ID54$temperature,IR_data_ID54$growth_rate)
plot(IR_data_ID54$temperature,briere1(a=0.0001,Tmin=15,Tmax=32,temp=IR_data_ID54$temperature))
plot(IR_data_ID54$temperature,briere1(a=0.0003,Tmin=15,Tmax=32,temp=IR_data_ID54$temperature))
plot(IR_data_ID54$temperature,briere1(a=0.0003,Tmin=15,Tmax=31.5,temp=IR_data_ID54$temperature))
grid_br1_54 <- expand.grid(list(a=seq(0.0002,0.0004,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_54 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID54,
                      start = grid_br1_54,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_54 <- summary(fitted_br1_54)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_54 <- c(sum_br1_54$coefficients[,1])
startVals_list[[54]] <- startVals_54 #llenar lista starting values
fitted_br1_54 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID54,
                     start = startVals_54)
sum_br1_54 <- summary(fitted_br1_54)
sum_br1_54 
boot_br1_54 <- nlsBoot(fitted_br1_54, niter = 999)
coefs_54 <- data.frame(sum_br1_54$coefficients[,1:2],row.names = NULL)
boot_54 <- data.frame(boot_br1_54$estiboot)
Topt_est_54 <- Topt(Tmin=coefs_54[2,1],
                    Tmax=coefs_54[3,1],
                    m=2)
Topt_se_54 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_54), vcov(fitted_br1_54))

params_br1_54 <- data.frame(coefs_54[1,],
                            coefs_54[2,],
                            coefs_54[3,],
                            Topt_est_54,
                            Topt_se_54,
                            boot_54[1,],
                            boot_54[2,],
                            boot_54[3,],
                            startVals_54[1],
                            startVals_54[2],
                            startVals_54[3],
                            "Tmin,Tmax")
colnames(params_br1_54) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_54

#### _ _ Study 55 ####
IR_data_ID55 <- IR_data_all %>%
  filter(id==55)
par(mfrow=c(2,2))
plot(IR_data_ID55$temperature,IR_data_ID55$growth_rate)
plot(IR_data_ID55$temperature,briere1(a=0.0001,Tmin=17,Tmax=34,temp=IR_data_ID55$temperature))
plot(IR_data_ID55$temperature,briere1(a=0.00012,Tmin=17,Tmax=35,temp=IR_data_ID55$temperature))
grid_br1_55 <- expand.grid(list(a=seq(0.00002,0.0002,by=0.00001),
                                Tmin=seq(12,22,by=0.5),
                                Tmax=seq(30,40,by=0.5)))
fitted_br1_55 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID55,
                      start = grid_br1_55,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_55 <- summary(fitted_br1_55)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_55 <- c(sum_br1_55$coefficients[,1])
startVals_list[[55]] <- startVals_55 #llenar lista starting values
fitted_br1_55 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID55,
                     start = startVals_55)
sum_br1_55 <- summary(fitted_br1_55)
sum_br1_55 
boot_br1_55 <- nlsBoot(fitted_br1_55, niter = 999)
coefs_55 <- data.frame(sum_br1_55$coefficients[,1:2],row.names = NULL)
boot_55 <- data.frame(boot_br1_55$estiboot)
Topt_est_55 <- Topt(Tmin=coefs_55[2,1],
                    Tmax=coefs_55[3,1],
                    m=2)
Topt_se_55 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_55), vcov(fitted_br1_55))

params_br1_55 <- data.frame(coefs_55[1,],
                            coefs_55[2,],
                            coefs_55[3,],
                            Topt_est_55,
                            Topt_se_55,
                            boot_55[1,],
                            boot_55[2,],
                            boot_55[3,],
                            startVals_55[1],
                            startVals_55[2],
                            startVals_55[3],
                            "none")
colnames(params_br1_55) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_55

#### _ _ Study 56 ####
IR_data_ID56 <- IR_data_all %>%
  filter(id==56)
par(mfrow=c(2,2))
plot(IR_data_ID56$temperature,IR_data_ID56$growth_rate)
plot(IR_data_ID56$temperature,briere1(a=0.0001,Tmin=17,Tmax=38,temp=IR_data_ID56$temperature))
plot(IR_data_ID56$temperature,briere1(a=0.0003,Tmin=17,Tmax=39,temp=IR_data_ID56$temperature))
plot(IR_data_ID56$temperature,briere1(a=0.00025,Tmin=15,Tmax=39,temp=IR_data_ID56$temperature))
grid_br1_56 <- expand.grid(list(a=seq(0.0001,0.0004,by=0.00001),
                                Tmin=seq(6,20,by=0.5),
                                Tmax=seq(32,46,by=0.5)))
fitted_br1_56 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID56,
                      start = grid_br1_56,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_56 <- summary(fitted_br1_56)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_56 <- c(sum_br1_56$coefficients[,1])
startVals_list[[56]] <- startVals_56 #llenar lista starting values
fitted_br1_56 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID56,
                     start = startVals_56)
sum_br1_56 <- summary(fitted_br1_56)
sum_br1_56 
boot_br1_56 <- nlsBoot(fitted_br1_56, niter = 999)
coefs_56 <- data.frame(sum_br1_56$coefficients[,1:2],row.names = NULL)
boot_56 <- data.frame(boot_br1_56$estiboot)
Topt_est_56 <- Topt(Tmin=coefs_56[2,1],
                    Tmax=coefs_56[3,1],
                    m=2)
Topt_se_56 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_56), vcov(fitted_br1_56))

params_br1_56 <- data.frame(coefs_56[1,],
                            coefs_56[2,],
                            coefs_56[3,],
                            Topt_est_56,
                            Topt_se_56,
                            boot_56[1,],
                            boot_56[2,],
                            boot_56[3,],
                            startVals_56[1],
                            startVals_56[2],
                            startVals_56[3],
                            "Tmax")
colnames(params_br1_56) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_56

#### 7. List of starting values ####
startVals_list
startVals_df <- startVals_list %>%
  as_tibble(.name_repair = "minimal")
colnames(startVals_df) <- seq(1,109,1) #check it out (probably 109 <- 56)
rownames(startVals_df) <- c("a","Tmin","Tmax")
startVals_df_long <- startVals_df %>%
 t()%>%
  as_tibble()%>%
  write_csv(file = "startVals.csv")

#### 8. Ensemble data frame with parameter estimates####
list_pattern <- ls(pattern = "params_br1_")#list objects with that word
df_list <- mget(list_pattern) #combines it into a df
ids_raw <- str_sub(list_pattern,12,15) #extract ids
ids <-as.numeric(ids_raw) #to sort them
parameters <- df_list %>%
  bind_rows() %>%
  bind_cols(ids) %>%
  mutate(id=...19) %>%
  select(-...19)%>%
  group_by_all()%>%
  summarise(id=sort(id))

write_csv(file="parameters_briere1_by_study.csv")

#### 9. Merge parameters df with info####
parameters <- read_csv("parameters_briere1_by_study.csv")
colnames(IR_data_all)
IR_data_all$id <- as.numeric(IR_data_all$id)
variables <- c("order","family","feeding_guild",
               "daylength","lat","lon","id") #columns of interest
data4params <- IR_data_all %>%
  select(order,family,feeding_guild,daylength,
         lat,lon,id) %>%
  set_colnames(variables)%>% #change variable names
  group_by_all() %>%
  summarise(id=unique(id)) %>%
  arrange(id)
#select only those with estimated parameters
filter_params <- unique(parameters$id)
filter_params_vec <- as.vector(filter_params)
dataset_params_metaanalysis <- data4params %>%
  filter(id %in% filter_params)%>% #apply a filter with numbers contain in that vector
  bind_cols(parameters) %>%
  glimpse()
write_csv(dataset_params_metaanalysis,file="dataset_params_metaanalysis.csv")  

#### 10. simulate random data points for each study ####
# firstly we convert standard error into variance
IR_data_sd <- IR_data %>%
  mutate(stdev=error*sqrt(n_1))
IR_data_ID <- IR_data_sd %>%
  distinct(title)%>%
  mutate(id=row_number()) #one id per distinct paper in a new dataframe
IR_data_all <- inner_join(IR_data_sd,IR_data_ID,by='title')
#### _ _ Study 1####
IR_data_ID1 <- IR_data_all %>%
  filter(id==1) %>%
  glimpse()
# first we assume normal distribution
## simulate data points
sim_10_ID1 <- rnorm(IR_data_ID1$n_1[1],mean = IR_data_ID1$growth_rate[1],sd =IR_data_ID1$stdev[1])
sim_15_ID1 <- rnorm(IR_data_ID1$n_1[2],mean = IR_data_ID1$growth_rate[2],sd =IR_data_ID1$stdev[2])
sim_20_ID1 <- rnorm(IR_data_ID1$n_1[3],mean = IR_data_ID1$growth_rate[3],sd =IR_data_ID1$stdev[3])
sim_25_ID1 <- rnorm(IR_data_ID1$n_1[4],mean = IR_data_ID1$growth_rate[4],sd =IR_data_ID1$stdev[4])
sim_30_ID1 <- rnorm(IR_data_ID1$n_1[5],mean = IR_data_ID1$growth_rate[5],sd =IR_data_ID1$stdev[5])
sim_33_ID1 <- rnorm(IR_data_ID1$n_1[6],mean = IR_data_ID1$growth_rate[6],sd =IR_data_ID1$stdev[6])
## make a tibble for each data point
IR_data_ID1_sim10 <- tibble(temp=rep(10,length(sim_10_ID1)),r=sim_10_ID1) 
IR_data_ID1_sim15 <- tibble(temp=rep(15,length(sim_15_ID1)),r=sim_15_ID1) 
IR_data_ID1_sim20 <- tibble(temp=rep(20,length(sim_20_ID1)),r=sim_20_ID1) 
IR_data_ID1_sim25 <- tibble(temp=rep(25,length(sim_25_ID1)),r=sim_25_ID1) 
IR_data_ID1_sim30 <- tibble(temp=rep(30,length(sim_30_ID1)),r=sim_30_ID1) 
IR_data_ID1_sim33 <- tibble(temp=rep(33,length(sim_33_ID1)),r=sim_33_ID1) 
## merge all of them
IR_data_ID1_sim <- bind_rows(IR_data_ID1_sim10,IR_data_ID1_sim15,
                              IR_data_ID1_sim20,IR_data_ID1_sim25,
                              IR_data_ID1_sim30,IR_data_ID1_sim33)

##let's explore data
par(mfrow=c(2,2))
plot(IR_data_ID1_sim$temp,IR_data_ID_1_sim$r)
plot(IR_data_ID1_sim$temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp=IR_data_ID_1_sim$temp))
plot(IR_data_ID1_sim$temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID_1_sim$temp))
plot(IR_data_ID1_sim$temp,briere1(a=0.00014,Tmin=8,Tmax=35,temp=IR_data_ID_1_sim$temp))
## now let's search starting values across a grid
grid_br1_1 <- expand.grid(list(a=seq(0.00010,0.00020,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(28,40,by=0.5)))
fitted_br1_1 <- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                     data = IR_data_ID_1_sim,
                     start = grid_br1_1,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_grid_1 <- summary(fitted_br1_1) #save the summary of this first scan
starVals_ID1_sim <- sum_grid_1$parameters[,1] #these are the starting values for
# the next model fitting function

## fit the model
fitted_br1_1sim <- nls(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
    data = IR_data_ID_1_sim,
    start = starVals_1_sim,
    trace = FALSE)
sum_br1_1sim <- summary(fitted_br1_1sim) #save summary

## alternative for heterogeneity of variances (nlme::gnls())
fitted_br1_1sim

fitted_br1_ID1sim_varPower <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                   data = IR_data_ID1_sim,
                                   start =starVals_ID1_sim,
                                   weights = varPower())
fitted_br1_ID1sim_varExp <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                 data = IR_data_ID1_sim,
                                 start =starVals_ID1_sim,
                                 weights = varExp())
fitted_br1_ID1sim_varConstPower <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                        data = IR_data_ID1_sim,
                                        start =starVals_ID1_sim,
                                        weights = varConstPower())

anova(fitted_br1_ID1sim_varPower,fitted_br1_ID1sim_varExp,fitted_br1_ID1sim_varConstPower)

## compare original vs. simulated 
#(we first need to run the Study 1 code in paragraph 6 of the script)
compare <- tibble(parameters = c("a","Tmin","Tmax"),
                  source_model=sum_br1_1$coefficients[,1],
                  signif_source = c("*","ns","***"),
                  simul_model=sum_br1_1sim$coefficients[,1],
                  signif_simul = c("***","*","***"),) %>%
  View()

#### _ _ Study 2####
IR_data_ID2 <- IR_data_all %>%
  filter(id==2) %>%
  dplyr::select(temperature)
  glimpse()
# first we assume normal distribution
## simulate data points
temps_ID2 <- IR_data_ID2$temperature %>%
  print()
sim_16_ID2 <- rnorm(IR_data_ID2$n_1[1],mean = IR_data_ID2$growth_rate[1],sd =IR_data_ID2$stdev[1])
sim_20_ID2 <- rnorm(IR_data_ID2$n_1[2],mean = IR_data_ID2$growth_rate[2],sd =IR_data_ID2$stdev[2])
sim_24_ID2 <- rnorm(IR_data_ID2$n_1[3],mean = IR_data_ID2$growth_rate[3],sd =IR_data_ID2$stdev[3])
sim_28_ID2 <- rnorm(IR_data_ID2$n_1[4],mean = IR_data_ID2$growth_rate[4],sd =IR_data_ID2$stdev[4])
sim_32_ID2 <- rnorm(IR_data_ID2$n_1[5],mean = IR_data_ID2$growth_rate[5],sd =IR_data_ID2$stdev[5])
## make a tibble for each data point
IR_data_ID2_sim16 <- tibble(temp=rep(16,length(sim_16_ID2)),r=sim_16_ID2) 
IR_data_ID2_sim20 <- tibble(temp=rep(20,length(sim_20_ID2)),r=sim_20_ID2) 
IR_data_ID2_sim24 <- tibble(temp=rep(24,length(sim_24_ID2)),r=sim_24_ID2) 
IR_data_ID2_sim28 <- tibble(temp=rep(28,length(sim_28_ID2)),r=sim_28_ID2) 
IR_data_ID2_sim32 <- tibble(temp=rep(32,length(sim_32_ID2)),r=sim_32_ID2) 
## merge all of them
IR_data_ID2_sim <- bind_rows(IR_data_ID2_sim16,IR_data_ID2_sim20,
                             IR_data_ID2_sim24,IR_data_ID2_sim28,
                             IR_data_ID2_sim32)

##let's explore data
par(mfrow=c(2,2))
plot(IR_data_ID2_sim$temp,IR_data_ID2_sim$r)
plot(IR_data_ID2_sim$temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp=IR_data_ID2_sim$temp))
plot(IR_data_ID2_sim$temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID2_sim$temp))
plot(IR_data_ID2_sim$temp,briere1(a=0.00014,Tmin=8,Tmax=35,temp=IR_data_ID2_sim$temp))
## now let's search starting values across a grid
grid_br1_ID2 <- expand.grid(list(a=seq(0.00010,0.00020,by=0.00001),
                               Tmin=seq(5,20,by=0.5),
                               Tmax=seq(32,42,by=0.5)))
fitted_br1_ID2 <- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                     data = IR_data_ID2_sim,
                     start = grid_br1_ID2,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_grid_ID2 <- summary(fitted_br1_ID2) #save the summary of this first scan
starVals_ID2_sim <- sum_grid_ID2$parameters[,1] #these are the starting values for
# the next model fitting function

## fit the model
fitted_br1_ID2sim <- nls(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                       data = IR_data_ID2_sim,
                       start = starVals_ID2_sim,
                       trace = FALSE)
sum_br1_ID2sim <- summary(fitted_br1_ID2sim) #save summary
## compare original vs. simulated 
#(we first need to run the Study ID2 code in paragraph 6 of the script)
compare <- tibble(parameters = c("a","Tmin","Tmax"),
                  source_model=sum_br1_2$coefficients[,1],
                  signif_source = c("*","*","**"),
                  simul_model=sum_br1_ID2sim$coefficients[,1],
                  signif_simul = c("***","***","***"),
                  bootstrap = boot_2$Estimate) %>%
  View()

fitted_br1_ID2sim_varPower <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                   data = IR_data_ID2_sim,
                                   start =starVals_ID2_sim,
                                   weights = varPower())
fitted_br1_ID2sim_varExp <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                 data = IR_data_ID2_sim,
                                 start =starVals_ID2_sim,
                                 weights = varExp())
fitted_br1_ID2sim_varConstPower <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                        data = IR_data_ID2_sim,
                                        start =starVals_ID2_sim,
                                        weights = varConstPower())
fitted_br1_ID2sim_varComb <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                  data = IR_data_ID2_sim,
                                  start =starVals_ID2_sim,
                                  weights = varComb())

anova(fitted_br1_ID2sim_varPower,fitted_br1_ID2sim_varExp)

#### _ _ Study 3####
IR_data_ID3 <- IR_data_all %>%
  filter(id==3) %>%
  glimpse()
# first we assume normal distribution
## simulate data points
temps_ID3 <- IR_data_ID3$temperature %>%
  print()
sim_13_ID3 <- rnorm(IR_data_ID3$n_1[1],mean = IR_data_ID3$growth_rate[1],sd =IR_data_ID3$stdev[1])
sim_18_ID3 <- rnorm(IR_data_ID3$n_1[2],mean = IR_data_ID3$growth_rate[2],sd =IR_data_ID3$stdev[2])
sim_23_ID3 <- rnorm(IR_data_ID3$n_1[3],mean = IR_data_ID3$growth_rate[3],sd =IR_data_ID3$stdev[3])
sim_25_ID3 <- rnorm(IR_data_ID3$n_1[4],mean = IR_data_ID3$growth_rate[4],sd =IR_data_ID3$stdev[4])
sim_28_ID3 <- rnorm(IR_data_ID3$n_1[5],mean = IR_data_ID3$growth_rate[5],sd =IR_data_ID3$stdev[5])
## make a tibble for each data point
IR_data_ID3_sim13 <- tibble(temp=rep(13,length(sim_13_ID3)),r=sim_13_ID3) 
IR_data_ID3_sim18 <- tibble(temp=rep(18,length(sim_18_ID3)),r=sim_18_ID3) 
IR_data_ID3_sim23 <- tibble(temp=rep(23,length(sim_23_ID3)),r=sim_23_ID3) 
IR_data_ID3_sim25 <- tibble(temp=rep(25,length(sim_25_ID3)),r=sim_25_ID3) 
IR_data_ID3_sim28 <- tibble(temp=rep(32,length(sim_28_ID3)),r=sim_28_ID3) 
## merge all of them
IR_data_ID3_sim <- bind_rows(IR_data_ID3_sim13,IR_data_ID3_sim18,
                             IR_data_ID3_sim23,IR_data_ID3_sim25,
                             IR_data_ID3_sim28)

##let's explore data
par(mfrow=c(2,2))
plot(IR_data_ID3_sim$temp,IR_data_ID3_sim$r)
plot(IR_data_ID3_sim$temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp=IR_data_ID3_sim$temp))
plot(IR_data_ID3_sim$temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID3_sim$temp))
plot(IR_data_ID3_sim$temp,briere1(a=0.00014,Tmin=8,Tmax=35,temp=IR_data_ID3_sim$temp))
## now let's search starting values across a grid
grid_br1_ID3 <- expand.grid(list(a=seq(0.00015,0.00025,by=0.00001),
                                 Tmin=seq(5,20,by=0.5),
                                 Tmax=seq(28,42,by=0.5)))
fitted_br1_ID3 <- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                       data = IR_data_ID3_sim,
                       start = grid_br1_ID3,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_grid_ID3 <- summary(fitted_br1_ID3) #save the summary of this first scan
starVals_ID3_sim <- sum_grid_ID3$parameters[,1] #these are the starting values for
# the next model fitting function

## fit the model
fitted_br1_ID3sim <- nls(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                         data = IR_data_ID3_sim,
                         start = starVals_ID3_sim,
                         trace = FALSE) # does not converge; we'll use the previous grid as final
sum_br1_ID3sim <- summary(fitted_br1_ID3) #save summary

fitted_br1_ID3sim_varPower <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                            data = IR_data_ID3_sim,
                            start =starVals_ID3_sim,
                            weights = varPower())
fitted_br1_ID3sim_varExp <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                   data = IR_data_ID3_sim,
                                   start =starVals_ID3_sim,
                                   weights = varExp())
fitted_br1_ID3sim_varConstPower <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                   data = IR_data_ID3_sim,
                                   start =starVals_ID3_sim,
                                   weights = varConstPower())
fitted_br1_ID3sim_varComb <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                   data = IR_data_ID3_sim,
                                   start =starVals_ID3_sim,
                                   weights = varComb())

anova(fitted_br1_ID3sim_varPower,fitted_br1_ID3sim_varExp,fitted_br1_ID3sim_varConstPower)
## compare original vs. simulated 
#(we first need to run the Study ID3 code in paragraph 6 of the script)
compare_ID3 <- tibble(parameters = c("a","Tmin","Tmax"),
                  source_model=sum_br1_3$coefficients[,1],
                  signif_source = c("*","*","**"),
                  simul_model=sum_br1_ID3sim$coefficients[,1],
                  signif_simul = c("*","**","***"),
                  bootstrap = boot_3$Estimate) %>%
  View()
### QUITAR PUNTO SI LA VARIANZA ES TAN GRANDE (ver el plot de simulación)

#### _ _ Study 4####
IR_data_ID4 <- IR_data_all %>%
  filter(id==4) %>%
  glimpse()
# first we assume normal distribution
## simulate data points
temps_ID4 <- IR_data_ID4$temperature %>%
  print()
sim_15_ID4 <- rnorm(IR_data_ID4$n_1[1],mean = IR_data_ID4$growth_rate[1],sd =IR_data_ID4$stdev[1])
sim_20_ID4 <- rnorm(IR_data_ID4$n_1[2],mean = IR_data_ID4$growth_rate[2],sd =IR_data_ID4$stdev[2])
sim_25_ID4 <- rnorm(IR_data_ID4$n_1[3],mean = IR_data_ID4$growth_rate[3],sd =IR_data_ID4$stdev[3])
sim_30_ID4 <- rnorm(IR_data_ID4$n_1[4],mean = IR_data_ID4$growth_rate[4],sd =IR_data_ID4$stdev[4])
## make a tibble for each data point
IR_data_ID4_sim15 <- tibble(temp=rep(15,length(sim_15_ID4)),r=sim_15_ID4) 
IR_data_ID4_sim20 <- tibble(temp=rep(20,length(sim_20_ID4)),r=sim_20_ID4) 
IR_data_ID4_sim25 <- tibble(temp=rep(25,length(sim_25_ID4)),r=sim_25_ID4) 
IR_data_ID4_sim30 <- tibble(temp=rep(30,length(sim_30_ID4)),r=sim_30_ID4) 
## merge all of them
IR_data_ID4_sim <- bind_rows(IR_data_ID4_sim15,IR_data_ID4_sim20,
                             IR_data_ID4_sim25,IR_data_ID4_sim30)

##let's explore data
par(mfrow=c(2,2))
plot(IR_data_ID4_sim$temp,IR_data_ID4_sim$r)
plot(IR_data_ID4_sim$temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp=IR_data_ID4_sim$temp))
plot(IR_data_ID4_sim$temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID4_sim$temp))
plot(IR_data_ID4_sim$temp,briere1(a=0.00014,Tmin=8,Tmax=35,temp=IR_data_ID4_sim$temp))
## now let's search starting values across a grid
grid_br1_ID4 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                 Tmin=seq(8,20,by=0.5),
                                 Tmax=seq(28,38,by=0.5)))
fitted_br1_ID4 <- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                       data = IR_data_ID4_sim,
                       start = grid_br1_ID4,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_grid_ID4 <- summary(fitted_br1_ID4) #save the summary of this first scan
starVals_ID4_sim <- sum_grid_ID4$parameters[,1] #these are the starting values for
# the next model fitting function

## fit the model
fitted_br1_ID4sim <- nls(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                         data = IR_data_ID4_sim,
                         start = starVals_ID4_sim,
                         trace = FALSE) 
sum_br1_ID4sim <- summary(fitted_br1_ID4sim) #save summary
## alternative for heterogeneity of variances (nlme::gnls())
fitted_br1_ID4sim_varPower <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                   data = IR_data_ID4_sim,
                                   start =starVals_ID4_sim,
                                   weights = varPower())
fitted_br1_ID4sim_varExp <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                 data = IR_data_ID4_sim,
                                 start =starVals_ID4_sim,
                                 weights = varExp())
fitted_br1_ID4sim_varConstPower <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                        data = IR_data_ID4_sim,
                                        start =starVals_ID4_sim,
                                        weights = varConstPower())

anova(fitted_br1_ID4sim_varPower,fitted_br1_ID4sim_varExp,fitted_br1_ID4sim_varConstPower)

## compare original vs. simulated 
#(we first need to run the Study ID4 code in paragraph 6 of the script)
compare_ID4 <- tibble(parameters = c("a","Tmin","Tmax"),
                      source_model=sum_br1_4$coefficients[,1],
                      signif_source = c("ns","ns","ns"),
                      simul_model=sum_br1_ID4sim$coefficients[,1],
                      signif_simul = c("***","***","***"),
                      simul_varPower = sum_br1_ID4sim_var$coefficients,
                      signif_sim_varPower= c("***","***","***"),
                      bootstrap = boot_4$Estimate) %>%
  View()


#### _ _ Study 5####
IR_data_ID5 <- IR_data_all %>%
  filter(id==5) %>%
  glimpse()
# first we assume normal distribution
## simulate data points
temps_ID5 <- IR_data_ID5$temperature %>%
  print()
sim_13_ID5 <- rnorm(IR_data_ID5$n_1[1],mean = IR_data_ID5$growth_rate[1],sd =IR_data_ID5$stdev[1])
sim_18_ID5 <- rnorm(IR_data_ID5$n_1[2],mean = IR_data_ID5$growth_rate[2],sd =IR_data_ID5$stdev[2])
sim_23_ID5 <- rnorm(IR_data_ID5$n_1[3],mean = IR_data_ID5$growth_rate[3],sd =IR_data_ID5$stdev[3])
sim_28_ID5 <- rnorm(IR_data_ID5$n_1[4],mean = IR_data_ID5$growth_rate[4],sd =IR_data_ID5$stdev[4])
## make a tibble for each data point
IR_data_ID5_sim13 <- tibble(temp=rep(13,length(sim_13_ID5)),r=sim_13_ID5) 
IR_data_ID5_sim18 <- tibble(temp=rep(18,length(sim_18_ID5)),r=sim_18_ID5) 
IR_data_ID5_sim23 <- tibble(temp=rep(23,length(sim_23_ID5)),r=sim_23_ID5) 
IR_data_ID5_sim28 <- tibble(temp=rep(28,length(sim_28_ID5)),r=sim_28_ID5) 
## merge all of them
IR_data_ID5_sim <- bind_rows(IR_data_ID5_sim13,IR_data_ID5_sim18,
                             IR_data_ID5_sim23,IR_data_ID5_sim28)

##let's explore data
par(mfrow=c(2,2))
plot(IR_data_ID5_sim$temp,IR_data_ID5_sim$r)
plot(IR_data_ID5_sim$temp,briere1(a=0.0002,Tmin=8,Tmax=35,temp=IR_data_ID5_sim$temp))
plot(IR_data_ID5_sim$temp,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID5_sim$temp))
plot(IR_data_ID5_sim$temp,briere1(a=0.00014,Tmin=8,Tmax=35,temp=IR_data_ID5_sim$temp))
## now let's search starting values across a grid
grid_br1_ID5 <- expand.grid(list(a=seq(0.00015,0.00025,by=0.00001),
                                 Tmin=seq(5,15,by=0.5),
                                 Tmax=seq(28,38,by=0.5)))
fitted_br1_ID5 <- nls2(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                       data = IR_data_ID5_sim,
                       start = grid_br1_ID5,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_grid_ID5 <- summary(fitted_br1_ID5) #save the summary of this first scan
starVals_ID5_sim <- sum_grid_ID5$parameters[,1] #these are the starting values for
# the next model fitting function

## fit the model
fitted_br1_ID5sim <- nls(formula= r ~ briere1(a,temp = temp,Tmin,Tmax),
                         data = IR_data_ID5_sim,
                         start = starVals_ID5_sim,
                         trace = FALSE) 
sum_br1_ID5sim <- summary(fitted_br1_ID5sim) #save summary
## alternative for heterogeneity of variances (nlme::gnls())
fitted_br1_ID5sim_varPower <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                   data = IR_data_ID5_sim,
                                   start =starVals_ID5_sim,
                                   weights = varPower())
fitted_br1_ID5sim_varExp <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                 data = IR_data_ID5_sim,
                                 start =starVals_ID5_sim,
                                 weights = varExp())
fitted_br1_ID5sim_varConstPower <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                                        data = IR_data_ID5_sim,
                                        start =starVals_ID5_sim,
                                        weights = varConstPower())

anova(fitted_br1_ID5sim_varPower,fitted_br1_ID5sim_varExp,fitted_br1_ID5sim_varConstPower)

## compare original vs. simulated 
#(we first need to run the Study ID4 code in paragraph 6 of the script)
compare_ID4 <- tibble(parameters = c("a","Tmin","Tmax"),
                      source_model=sum_br1_4$coefficients[,1],
                      signif_source = c("ns","ns","ns"),
                      simul_model=sum_br1_ID4sim$coefficients[,1],
                      signif_simul = c("***","***","***"),
                      simul_varPower = sum_br1_ID4sim_var$coefficients,
                      signif_sim_varPower= c("***","***","***"),
                      bootstrap = boot_4$Estimate) %>%
  View()












#### 11. Evaluate distribution at each temperature #### 
# should we try log-norm distribution assumption? (rlnorm() function)
hemiptera_first <- hemiptera %>%
  filter(temperature==20)
hist(hemiptera_first$growth_rate)
dens1 <- density(hemiptera_first$growth_rate)
lnorm_hemi<- hemiptera_first %>%
  mutate(r_trans = log(growth_rate))
hist(lnorm_hemi$r_trans,breaks = 8)
dens2 <- density(lnorm_hemi$r_trans)
par(mfrow=c(1,2))
plot(dens1)
plot(dens2)
shapiro.test(hemiptera_first$growth_rate)
shapiro.test(lnorm_hemi$r_trans)

acari_first <- acari %>%
  filter(temperature==30) #mejora especialmente en temperaturas altas
hist(acari_first$growth_rate)
dens1 <- density(acari_first$growth_rate)
lnorm_acari<- acari_first %>%
  mutate(r_trans = log(growth_rate))
hist(lnorm_acari$r_trans,breaks = 8)
dens2 <- density(lnorm_acari$r_trans)
par(mfrow=c(1,2))
plot(dens1)
plot(dens2)
shapiro.test(acari_first$growth_rate)
shapiro.test(which(!is.na(lnorm_acari$r_trans)))



#lepi
lepidoptera_first <- lepidoptera %>%
  filter(temperature==30)
hist(lepidoptera_first$growth_rate)
dens1 <- density(lepidoptera_first$growth_rate)
lnorm_lepi<- lepidoptera_first %>%
  mutate(r_trans = log(growth_rate))
hist(lnorm_lepi$r_trans,breaks = 8)
dens2 <- density(lnorm_lepi$r_trans)
par(mfrow=c(1,2))
plot(dens1)
plot(dens2)
shapiro.test(lepidoptera_first$growth_rate)
shapiro.test(lnorm_lepi$r_trans)

#pooled
all_first <- IR_data %>%
  filter(temperature==35) #at 15, 20, 25, 30 it usually becomes more normal after log-transform
hist(all_first$growth_rate)
dens1 <- density(all_first$growth_rate)
lnorm_all<- all_first %>%
  mutate(r_trans = log(growth_rate))
hist(lnorm_all$r_trans,breaks = 8)
dens2 <- density(lnorm_all$r_trans)
par(mfrow=c(1,2))
plot(dens1)
plot(dens2)
shapiro.test(all_first$growth_rate)
shapiro.test(which(!is.na(lnorm_all$r_trans)))

#### 12. Evaluate homoscedasticity ####
car::leveneTest(IR_data_ID2_sim$r,IR_data_ID2_sim$temp) #homoscedastic though via test

ID2sim_resid <- nlsResiduals(fitted_br1_ID2sim) #neither normal nor homoscedastic
plot(ID2sim_resid,which=0) #todos

prueba_r1 <- rnorm(25,0.1,0.1)
prueba_r2 <- rnorm(25,0.15,0.1)
prueba_r3 <- rnorm(25,0.25,0.6)
prueba_r4 <- rnorm(25,0.08,1.5)
prueba_r <- c(prueba_r1,prueba_r2,prueba_r3,prueba_r4)
prueba_t <- c(rep(15,25),rep(20,25),rep(25,25),rep(30,25))
prueba <- tibble(temp=prueba_t,r=prueba_r)
plot(prueba)

#### 13. try nlme() ####
plot(IR_data_all$temperature,IR_data_all$growth_rate)
grid_all <- expand.grid(list(a=seq(0.00010,0.00020,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(28,45,by=0.5)))
fitted_br1_start <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_all,
                     start = grid_all,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_starting <- summary(fitted_br1_start)
sum_starting

starts = c(sum_starting$coefficients[,1],temperature=20)[c(1,4,2,3)]

nlme_br1_1 <- nlme(model= growth_rate ~ briere1(a,temperature,Tmin,Tmax),
                   start = starts,
                   fixed = list(a~1,temperature~1,Tmin~1,Tmax~1),
                   random ~ 1|id,
                   data = IR_data_all)
#Error in chol.default((value + t(value))/2) : 
# the leading minor of order 1 is not positive definite

#let's try to remove 0's in growth_rate function
IR_data_all_mod4nlme <- IR_data_all %>%
  filter(growth_rate >=0.01)
length(unique(IR_data_all_mod4nlme$id)) #it only removes points, not complete studies
nlme_br1_1 <- nlme(model= growth_rate ~ briere1(a,temperature,Tmin,Tmax),
                   start = starts,
                   fixed = list(a~1,temperature~1,Tmin~1,Tmax~1),
                   random ~ 1|id,
                   data = IR_data_all_mod4nlme)
#same error :(
nlme_br1_1 <- nlme(growth_rate ~ briere1(a,temperature,Tmin,Tmax),
                   data =IR_data_all,
                   start = sum_starting$coefficients[,1],
                   fixed = a+Tmin+Tmax ~1,
                   random ~ 1|id)


acari <- IR_data_all %>%
  filter(order=="Acari>Prostigmata")
nlme_br1_acari <- nlme(model= growth_rate ~ briere1(a,temperature,Tmin,Tmax),
                   start = c(a=0.0002,Tmin=8,Tmax=40),
                   fixed = a+Tmin+Tmax ~1,
                   random =a+Tmin+Tmax~1|id,
                   data = acari)
nlme_br1_acari
summary(nlme_br1_acari)

nlme_br1_all <- nlme(model= growth_rate ~ briere1(a,temperature,Tmin,Tmax),
                       start = c(a=0.0001,Tmin=10,Tmax=40),
                       fixed = a+Tmin+Tmax ~1,
                       random =a+Tmin+Tmax~1|id,
                       data = IR_data_all,
                       control = nlmeControl(pnlsTol = 1, msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
#either with tolerance = 1 (def = 1e-6) or pnlsTol =1(def. = 1e-3)
nlme_br1_all <- nlme(model= growth_rate ~ briere1(a,temperature,Tmin,Tmax),
                     start = c(a=0.0001,Tmin=10,Tmax=40),
                     fixed = a+Tmin+Tmax ~1,
                     random =a+Tmin+Tmax~1|id,
                     data = IR_data_all,
                     control = nlmeControl(pnlsTol = 1,msVerbose = TRUE))

summary(nlme_br1_all)
#probably errors and problems come from high correlations between a and Tmin
nlme_br1_acari <- nlme(model= growth_rate ~ briere1(a,temperature,Tmin,Tmax),
                     start = c(a=0.0002,Tmin=8,Tmax=40),
                     fixed = a+Tmin+Tmax ~1,
                     random =a+Tmin+Tmax~1|id,
                     data = acari,
                     control = nlmeControl(pnlsTol = , msVerbose = TRUE))#to avoid error of singularity in backsolve at level 0; block 1
#maxIter problem now; not solved with maxIter increase;
# it is solved with better starting values (changing Tmin into 8 rather than 6)
lepi <- IR_data_all %>%
  filter(order=="Lepidoptera")

nlme_br1_lepi <- nlme(model= growth_rate ~ briere1(a,temperature,Tmin,Tmax),
                       start = c(a=0.0000609,Tmin=6.42,Tmax=37.52),
                       fixed = a+Tmin+Tmax ~1,
                       random =a+Tmin+Tmax~1|id,
                       data = lepi)
                       control = nlmeControl(pnlsTol = , msVerbose = TRUE))
nlme_br1_lepi

#### 14. try brm() ####
syntax_example <- brm(formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
            data = kidney, family = lognormal(),
            prior = c(set_prior("normal(0,5)", class = "b"), #
                      set_prior("cauchy(0,2)", class = "sd"),
                      set_prior("lkj(2)", class = "cor")),
            warmup = 1000,
            iter = 2000, 
            chains = 4,
            control = list(adapt_delta = 0.95))

brm_br1_all <- brm(formula = growth_rate ~ temperature + (1|id), #model + random effects
                   data = IR_data_all,
                   family = "lognormal",
                   nonlinear = a + Tmin + Tmax ~ 1,
                   prior = c(c(set_prior("normal(4.187e-05, 2e-05)", nlpar = "a"),
                               set_prior("normal(-5.028, 10)", nlpar = "Tmin"),
                               set_prior("normal(45.45,5.78)", nlpar = "Tmax"))),
                   warmup = 1000,
                   iter = 2000, 
                   chains = 4,
                   control = list(adapt_delta = 0.95))

#example meta-analysis https://www.barelysignificant.com/slides/RGUG2019#97
prior4 <- c(
  prior(normal(0, 1), coef = intercept),
  prior(cauchy(0, 1), class = sd)
)
mod4 <- brm(
  yi | se(sqrt(vi) ) ~ 0 + intercept + (1|study) + (1|experiment),
  data = d,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 2000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

### from https://rpubs.com/aforren1/orange-nonlinear 

f1 <- circumference ~ phi1 / (1 + exp(-(age - phi2)/phi3)) #briere1 is already defined in our example
n1 <- nls(f1,
          data = Orange, 
          start = list(phi1 = 200, phi2 = 700, phi3 = 350)) #n1 is our estimate with nls (fitt_br1_pooled)

prior_1 <- c(set_prior("normal(200, 50)", nlpar = "phi1"),
             set_prior("normal(700, 50)", nlpar = "phi2"),
             set_prior("normal(350, 50)", nlpar = "phi3"))
n1_b <- brm(f1, 
            data = Orange,
            nonlinear = phi1 + phi2 + phi3 ~ 1,
            prior = prior_1, 
            chains = 3)
# translate it
f1 <- growth_rate ~ a*temperature*(temperature-Tmin)*sqrt(Tmax-temperature) #briere1 is already defined in our example
n1 <- fitt_br1_pooled
sum_all
prior_1 <- c(set_prior("normal(4.187e-05, 2e-05)", nlpar = "a"),
             set_prior("normal(-5.028, 10)", nlpar = "Tmin"),
             set_prior("normal(45.45,5.78)", nlpar = "Tmax"))
n1_b <- brm(f1,
            data = IR_data_all,
            nonlinear = a + Tmin + Tmax ~ 1,
            prior = prior_1, 
            chains = 4)

prior_summary(prior_1)
