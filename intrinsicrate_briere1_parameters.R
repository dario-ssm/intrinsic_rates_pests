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
library(nlstools)
library(tidyr)
library(nls2)
library(msm)
library(magrittr)
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

# We use a Brière-1 model (Briere et al., 1999)
briere1 <- function(a, temp, Tmin, Tmax){
  a*temp*(temp-Tmin)*(Tmax-temp)^(1/2)
}

Topt <- function(Tmin,Tmax,m){
  Topt=((2*m*Tmax+(m+1)*Tmin)+sqrt(4*(m^2)*(Tmax^2)+((m+1)^2)*(Tmin^2)-4*(m^2)*Tmin*Tmax))/(4*m+2)
  return(Topt)
}
Topt(Tmin=8.14,Tmax=41.08,m=2) #sale 33.79ºC (ver ejemplo en acari)

#### 3. Train data ####
#load training artificial data
train_data <- read.table("train.csv",sep=";",dec=",",header=TRUE) #inventados siguiendo la curva aprox.
colnames(train_data) <- c("temp","r")
train_data$temp <- as.numeric(train_data$temp)
train_data$r <- as.numeric(train_data$r)
#View(train_data)
ggplot(train_data,aes(x=temp,y=r))+
  geom_point()+
  stat_smooth(method= "loess")+
  theme_classic()
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
summary(fitted_br1)
coef(fitted_br1)
AIC(fitted_br1)
BIC(fitted_br1)
#plot
temp <- train_data$temp
growth <- train_data$r
pred <- tibble(temp,growth,predict(fitted_br1))
colnames(pred) <- c("temp","growth","fit")
pred %>% 
  filter(fit==max(fit))
#View(pred)
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
Tmin = ns,
Topt = 30.5ºC
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
#vamos a ver initials (hasta que RSS sea mínimo, primero a, luego valores)
preview(r~briere1(a,temp,Tmin,Tmax),
        data=acari_test,
        start=list(a=0.00011,
                   Tmin=4,
                   Tmax=39))

fitted_br1_acari <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                        data = acari_test,
                        start = list(a = 0.00011,
                                     Tmin =4,
                                     Tmax= 39))
model_acari_sum <-summary(fitted_br1_acari)
coef(fitted_br1_acari)
Topt_acari <- Topt(Tmin=model_acari_sum$coefficients[2,1],
                   Tmax=model_acari_sum$coefficients[3,1],
                   m=2) #en Brière-1, m=2 siempre
Topt_acari
AIC(fitted_br1_acari)
BIC(fitted_br1_acari)
#sacamos error con bootstrapping
Boot_fitbr1 <- nlsBoot(fitted_br1_acari, niter = 999)
Boot_fitbr1$coefboot
Boot_fitbr1$bootCI
Boot_fitbr1$estiboot #interesante esto también
plot(Boot_fitbr1)
Parameters_acari <- data.frame(model_acari_sum$parameters[,1:2],
                               Boot_fitbr1$estiboot[,1],
                               Boot_fitbr1$estiboot[,2])
colnames(Parameters_acari) <- c("nls_estimate","nls_SE","Bootstrap_estimate","Bootstrap_SE")
Parameters_acari
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

fitted_br1_lepidoptera <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                              data = lepidoptera_test,
                              start = list(a = 0.00008,
                                           Tmin =-5, #hasta aquí no cortaba... ver el loess
                                           Tmax= 38)) # a partir de 38

predict(fitted_br1_lepidoptera)
model_lepidoptera_sum <- summary(fitted_br1_lepidoptera)
AIC(fitted_br1_lepidoptera)
BIC(fitted_br1_lepidoptera)
coefs_lepi <-data.frame(coef(fitted_br1_lepidoptera))
Topt_lepidoptera <- Topt(Tmin=coefs_lepi[2,],
                         Tmax=coefs_lepi[3,],
                         m=2)
Topt_lepidoptera
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
summary(fitted_br1_diptera)
AIC(fitted_br1_diptera)
BIC(fitted_br1_diptera)
coef(fitted_br1_diptera)
coefs_diptera <-data.frame(coef(fitted_br1_diptera))
Topt_diptera <- Topt(Tmin=coefs_diptera[2,],
                     Tmax=coefs_diptera[3,],
                     m=2)
Topt_diptera
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


fitted_br1 <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                  data = hemiptera_test,
                  start = list(a = 0.00017,
                               Tmin =8,
                               Tmax= 40),
                  trace = FALSE,
                  nls.control(warnOnly=TRUE)) #me salía error number iterations exceeded 50, 
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
AIC(fitted_br1)
BIC(fitted_br1)
coef(fitted_br1)
#vamos a dibujarlo
pred <- tibble(temp,growth,predict(fitted_br1))
colnames(pred) <- c("temp","growth","fit")
pred %>%
  filter(fit == max(fit)) #26ºC
#View(pred)
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



#getInitials?
#function   for the model

#model_initial<-function(mCall,LHS,data){
xy<-sortedXyData(mCall[['temp']],LHS,data)
fit<-lm(xy[,'y']~xy[,'x'])
coefs<-coef(fit)
a<-coefs[1]    
Tmin<-coefs[2]
Tmax<-coefs[3]
value<-c(a,Tmin,Tmax)

names(value)<-mCall[c('a','Tmin','Tmax')]
value
}

Self_starter<-selfStart(briere1,model_initial,c('a','Tmin','Tmax'))

print(getInitial(r ~ Self_starter(temp,a,TO,Tl), data=aphids_test), digits = 3)

## no hay forma... eliminamos outliers
aphids_clean <- aphids%>%
  filter(growth_rate <0.4)%>%
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
                  start = c(a = 0.0002,
                            Tmin =8,
                            Tmax= 39),
                  trace = FALSE)  
#quitando los ceros y los outliers, sale, pero solo el Tmax (34.12ºC)
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
grid_br1_coleoptera <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                        Tmin=seq(5,15,by=0.5),
                                        Tmax=seq(25,45,by=0.5)))
fitted_br1_coleoptera <- nls2(formula= r ~ briere1(a,temp,Tmin,Tmax),
                              data = coleoptera_test,
                              start = grid_br1_coleoptera,
                              algorithm = "brute-force",
                              trace = TRUE)
summary(fitted_br1_coleoptera)
# Tmin = 11; Tmax=36,
Topt_coleoptera <- Topt(Tmin=11,Tmax=36,m=2)

#Error in nls(r ~ briere1(a, temp, Tmin, Tmax), 
#### 5. All plots combined ####
library(cowplot)
plot_grid(briere_plot_acari,briere_plot_diptera,briere_plot_lepidoptera, briere_plot_aphids,
          nrow = 2,labels = c("A","B","C","D"))

#### 6. More complex models ####
#### _ _ 6.1. Lepidoptera ####
lepidoptera_lat <- IR_data%>%
  filter(order=="Lepidoptera")%>%
  select(order,r=growth_rate,temp=temperature,lat)%>%
  mutate(lat=abs(lat))%>%
  mutate(order=)
glimpse()


fitted_br1_lat <- nls(r ~ briere1(a,temp,Tmin,Tmax),
                      data = lepidoptera_lat,
                      start = list(a = 0.00008,
                                   Tmin =-5, #hasta aquí no cortaba... ver el loess
                                   Tmax= 38))
#### 7. Pooled ####
IR_test <- IR_data %>%
  select(r=growth_rate,temp=temperature,order,lat,hostplant=h_p_family,title=`Article Title`)%>%
  glimpse()
#### _ _ a) r ~ temp ####
#explorar el scatterplot agrupado
all_together <- ggplot(IR_test,aes(x=temp,y=r))+
  geom_point(color="turquoise4",alpha=0.2)+
  geom_smooth(color="lightcoral",fill="lightcoral")+
  theme_classic()
par(mfrow=c(2,2))
growth <- IR_test$r
temp <- IR_test$temp
plot(temp,growth)
plot(temp,briere1(a=0.0002,Tmin=8,Tmax=40,temp))
plot(temp,briere1(a=0.0001,Tmin=8,Tmax=40,temp))# <- este
plot(temp,briere1(a=0.00015,Tmin=8,Tmax=40,temp))

preview(formula=r~briere1(a,temp,Tmin,Tmax),
        data=IR_test,
        start=list(a=0.0001,Tmin=8,Tmax=40))

fitt_br1_pooled <- nls(r~briere1(a,temp,Tmin,Tmax),
                       data= IR_test,
                       start=list(a=0.0001,Tmin=8,Tmax=40))

summary(fitt_br1_pooled)
coefs_pooled <-data.frame(coef(fitt_br1_pooled))
Topt_pooled<- Topt(Tmin=coefs_pooled[2,],
                   Tmax=coefs_pooled[3,],
                   m=2)

AIC(fitt_br1_pooled)
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

#### _ _ b) r ~ temp+order ####
fitt_br1_pooled_order <- nls(r~briere1(a,temp,Tmin,Tmax)+order,
                             data= IR_test,
                             start=list(a=0.001,Tmin=8,Tmax=40))
temp_sized <- seq(from = min(temp),to = max(temp))
seq(1,10,by = 1)
#### 8. nlstools package ####
library(nlstools)
fitted_br1 # lo llamamos aquí
#intervalos de confianza para un parámetro
confint2(fitted_br1, "Tmin", level = 0.95, method = c("asymptotic", "profile"))

#resampling bootstrap
Boot_fitbr1 <- nlsBoot(fitted_br1, niter = 999)
Boot_fitbr1$coefboot
Boot_fitbr1$bootCI
Boot_fitbr1$estiboot #interesante esto también
plot(Boot_fitbr1)
#intregions
nlsConfRegions (fitted_br1, length = 1000, exp = 1.5) #error :()
fitbr1Cont<-nlsContourRSS (fitted_br1, lseq = 100, exp = 2)
plot(fitbr1Cont)

#residuos: ver NORMALIDAD y HOMOCEDASTICIDAD
fitbr1_resid <- nlsResiduals(fitted_br1)
plot(fitbr1_resid,which=0) #todos
plot(fitbr1_resid,which=5) #histogram of residuals
plot(fitbr1_resid,which=6) #qqplot
test.nlsResiduals(fitbr1_resid) #test de normalidad Shapiro-Wilk

#preview GENIAL PARA VER initials...
preview (formula=r ~ briere1(a,temp,Tmin,Tmax),
         data = acari_test,
         start = list(a = 0.00014, #el que mejor me ajusta el RSS
                      Tmin =8.6,
                      Tmax= 39),variable = 1)

##¿funcionará para encontrar initials para coleoptera?
preview(formula= r~ briere1(a,temp,Tmin,Tmax),
        data = coleoptera_test,
        start = list(a = 0.000029,
                     Tmin =3.2,
                     Tmax= 37),variable=1)
#ir trasteando con valores que minimizan el RSS
## ver summaries
plotfit (fitted_br1) #error?
overview (fitted_br1)

#### 9. Assumption tests ####
# example acari
plot(fitted(fitted_br1_acari),residuals(fitted_br1_acari),
     xlab = "Fitted Values", ylab= "Residuals") #mala pinta
abline(a = 0, b = 0)
fitted_br1_acari_lm <- lm(r ~ temp, data = acari_test)
anova(fitted_br1_acari, fitted_br1_acari_lm) #diferencias significat.
# modelo no apropiado?

## ver homocedasticidad 
plot(fitted(fitted_br1_acari),abs(residuals(fitted_br1_acari)))
#claramente aumenta la varianza al desplazarse en la variable
library(car)
leveneTest(acari_test$r,acari_test$temp) #sale no homocedástico
# por tanto, esto no debería afectar a la estimación de parámetros,
# pero sí a la inferencia estadística
#normalidad
standardRes <- residuals(fitted_br1_acari)/summary(fitted_br1_acari)$sigma
qqnorm(standardRes, main= "")
abline(a=0,b=1)

#correlación
plot(residuals(fitted_br1_acari),c(residuals(fitted_br1_acari)[-1], NA),
     xlab = "Residuals",
     ylab = "Lagged residuals")
abline(a=0,b=1)     

## Dealing with assumptions violations
library(nlmer)
fitted_br1_var <- gnls(r ~ briere1(a,temp,Tmin,Tmax),
                       data = IR_data,
                       weights = varPower())
summary(fitted_br1_var)

#transformations box-cox
library(nlrwr)
BC_fitted_br1 <- boxcox.nls(fitted_br1)
summary(BC_fitted_br1)



#### 10. individual studies regression with Bootstrap ####
## first: assign a number to each unique article as ID
IR_data_ID <- IR_data %>%
  distinct(`Article Title`)%>%
  mutate(id=row_number())
IR_data_all <- inner_join(IR_data,IR_data_ID,by="Article Title")
startVals_list <- list() #para almacenar los starting
## one by one ##

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
sum_br1_1 <- summary(fitted_br1_1)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_1 <- c(sum_br1_1$coefficients[,1])
startVals_list[[1]] <- startVals_1 #llenar lista starting values

# probar con esos parámetros como iniciales:
fitted_br1_1 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID1,
                    start = startVals_1)


sum_br1_1 <- summary(fitted_br1_1)
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

#### _ _ Study 57 ####
IR_data_ID57 <- IR_data_all %>%
  filter(id==57)
par(mfrow=c(2,2))
plot(IR_data_ID57$temperature,IR_data_ID57$growth_rate)
plot(IR_data_ID57$temperature,briere1(a=0.0001,Tmin=16,Tmax=35,temp=IR_data_ID57$temperature))
grid_br1_57 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_57 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID57,
                      start = grid_br1_57,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_57 <- summary(fitted_br1_57)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_57 <- c(sum_br1_57$coefficients[,1])
startVals_list[[57]] <- startVals_57 #llenar lista starting values
fitted_br1_57 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID57,
                     start = startVals_57)
sum_br1_57 <- summary(fitted_br1_57)
sum_br1_57 
boot_br1_57 <- nlsBoot(fitted_br1_57, niter = 999)
coefs_57 <- data.frame(sum_br1_57$coefficients[,1:2],row.names = NULL)
boot_57 <- data.frame(boot_br1_57$estiboot)
Topt_est_57 <- Topt(Tmin=coefs_57[2,1],
                    Tmax=coefs_57[3,1],
                    m=2)
Topt_se_57 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_57), vcov(fitted_br1_57))

params_br1_57 <- data.frame(coefs_57[1,],
                            coefs_57[2,],
                            coefs_57[3,],
                            Topt_est_57,
                            Topt_se_57,
                            boot_57[1,],
                            boot_57[2,],
                            boot_57[3,],
                            startVals_57[1],
                            startVals_57[2],
                            startVals_57[3],
                            "all")
colnames(params_br1_57) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_57

#### _ _ Study 58 ####
IR_data_ID58 <- IR_data_all %>%
  filter(id==58)
par(mfrow=c(2,2))
plot(IR_data_ID58$temperature,IR_data_ID58$growth_rate)
plot(IR_data_ID58$temperature,briere1(a=0.0002,Tmin=17,Tmax=36,temp=IR_data_ID58$temperature))
plot(IR_data_ID58$temperature,briere1(a=0.00025,Tmin=17,Tmax=36,temp=IR_data_ID58$temperature))
grid_br1_58 <- expand.grid(list(a=seq(0.0002,0.0003,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_58 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID58,
                      start = grid_br1_58,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_58 <- summary(fitted_br1_58)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_58 <- c(sum_br1_58$coefficients[,1])
startVals_list[[58]] <- startVals_58 #llenar lista starting values
fitted_br1_58 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID58,
                     start = startVals_58)
sum_br1_58 <- summary(fitted_br1_58)
sum_br1_58 
boot_br1_58 <- nlsBoot(fitted_br1_58, niter = 999)
coefs_58 <- data.frame(sum_br1_58$coefficients[,1:2],row.names = NULL)
boot_58 <- data.frame(boot_br1_58$estiboot)
Topt_est_58 <- Topt(Tmin=coefs_58[2,1],
                    Tmax=coefs_58[3,1],
                    m=2)
Topt_se_58 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_58), vcov(fitted_br1_58))

params_br1_58 <- data.frame(coefs_58[1,],
                            coefs_58[2,],
                            coefs_58[3,],
                            Topt_est_58,
                            Topt_se_58,
                            boot_58[1,],
                            boot_58[2,],
                            boot_58[3,],
                            startVals_58[1],
                            startVals_58[2],
                            startVals_58[3],
                            "all")
colnames(params_br1_58) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_58

#### _ _ Study 59 ####
IR_data_ID59 <- IR_data_all %>%
  filter(id==59)
par(mfrow=c(2,2))
plot(IR_data_ID59$temperature,IR_data_ID59$growth_rate)
plot(IR_data_ID59$temperature,briere1(a=0.0001,Tmin=15,Tmax=31,temp=IR_data_ID59$temperature))
plot(IR_data_ID59$temperature,briere1(a=0.00009,Tmin=15,Tmax=31,temp=IR_data_ID59$temperature))
grid_br1_59 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_59 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID59,
                      start = grid_br1_59,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_59 <- summary(fitted_br1_59)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_59 <- c(sum_br1_59$coefficients[,1])
startVals_list[[59]] <- startVals_59 #llenar lista starting values
fitted_br1_59 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID59,
                     start = startVals_59)
sum_br1_59 <- summary(fitted_br1_59)
sum_br1_59 
boot_br1_59 <- nlsBoot(fitted_br1_59, niter = 999)
coefs_59 <- data.frame(sum_br1_59$coefficients[,1:2],row.names = NULL)
boot_59 <- data.frame(boot_br1_59$estiboot)
Topt_est_59 <- Topt(Tmin=coefs_59[2,1],
                    Tmax=coefs_59[3,1],
                    m=2)
Topt_se_59 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_59), vcov(fitted_br1_59))

params_br1_59 <- data.frame(coefs_59[1,],
                            coefs_59[2,],
                            coefs_59[3,],
                            Topt_est_59,
                            Topt_se_59,
                            boot_59[1,],
                            boot_59[2,],
                            boot_59[3,],
                            startVals_59[1],
                            startVals_59[2],
                            startVals_59[3],
                            "all")
colnames(params_br1_59) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_59

#### _ _ Study 60 ####
# NO, only 2 data

#### _ _ Study 61 ####
IR_data_ID61 <- IR_data_all %>%
  filter(id==61)
par(mfrow=c(2,2))
plot(IR_data_ID61$temperature,IR_data_ID61$growth_rate)
plot(IR_data_ID61$temperature,briere1(a=0.0001,Tmin=4,Tmax=26,temp=IR_data_ID61$temperature))
plot(IR_data_ID61$temperature,briere1(a=0.0002,Tmin=0,Tmax=26,temp=IR_data_ID61$temperature))
grid_br1_61 <- expand.grid(list(a=seq(0.00005,0.0003,by=0.00001),
                                Tmin=seq(0,10,by=0.5),
                                Tmax=seq(18,28,by=0.5)))
fitted_br1_61 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID61,
                      start = grid_br1_61,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_61 <- summary(fitted_br1_61)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_61 <- c(sum_br1_61$coefficients[,1])
startVals_list[[61]] <- startVals_61 #llenar lista starting values
fitted_br1_61 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID61,
                     start = startVals_61)
#NO sale

#### _ _ Study 62 ####
IR_data_ID62 <- IR_data_all %>%
  filter(id==62)
par(mfrow=c(2,2))
plot(IR_data_ID62$temperature,IR_data_ID62$growth_rate)
plot(IR_data_ID62$temperature,briere1(a=0.0001,Tmin=14,Tmax=34,temp=IR_data_ID62$temperature))
plot(IR_data_ID62$temperature,briere1(a=0.0002,Tmin=14,Tmax=34,temp=IR_data_ID62$temperature))
plot(IR_data_ID62$temperature,briere1(a=0.00024,Tmin=15,Tmax=35,temp=IR_data_ID62$temperature))
grid_br1_62 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_62 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID62,
                      start = grid_br1_62,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_62 <- summary(fitted_br1_62)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_62 <- c(sum_br1_62$coefficients[,1])
startVals_list[[62]] <- startVals_62 #llenar lista starting values
fitted_br1_62 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID62,
                     start = startVals_62)
sum_br1_62 <- summary(fitted_br1_62)
sum_br1_62 
boot_br1_62 <- nlsBoot(fitted_br1_62, niter = 999)
coefs_62 <- data.frame(sum_br1_62$coefficients[,1:2],row.names = NULL)
boot_62 <- data.frame(boot_br1_62$estiboot)
Topt_est_62 <- Topt(Tmin=coefs_62[2,1],
                    Tmax=coefs_62[3,1],
                    m=2)
Topt_se_62 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_62), vcov(fitted_br1_62))

params_br1_62 <- data.frame(coefs_62[1,],
                            coefs_62[2,],
                            coefs_62[3,],
                            Topt_est_62,
                            Topt_se_62,
                            boot_62[1,],
                            boot_62[2,],
                            boot_62[3,],
                            startVals_62[1],
                            startVals_62[2],
                            startVals_62[3],
                            "all")
colnames(params_br1_62) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_62

#### _ _ Study 63 ####
IR_data_ID63 <- IR_data_all %>%
  filter(id==63)
par(mfrow=c(2,2))
plot(IR_data_ID63$temperature,IR_data_ID63$growth_rate)
plot(IR_data_ID63$temperature,briere1(a=0.0001,Tmin=19,Tmax=38,temp=IR_data_ID63$temperature))
plot(IR_data_ID63$temperature,briere1(a=0.00009,Tmin=19.5,Tmax=38,temp=IR_data_ID63$temperature))
plot(IR_data_ID63$temperature,briere1(a=0.00005,Tmin=19,Tmax=38,temp=IR_data_ID63$temperature))
grid_br1_63 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(15,25,by=0.5),
                                Tmax=seq(32,46,by=0.5)))
fitted_br1_63 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID63,
                      start = grid_br1_63,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_63 <- summary(fitted_br1_63)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_63 <- c(sum_br1_63$coefficients[,1])
startVals_list[[63]] <- startVals_63 #llenar lista starting values
fitted_br1_63 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID63,
                     start = startVals_63)
#no sale

#### _ _ Study 64 ####
IR_data_ID64 <- IR_data_all %>%
  filter(id==64)
par(mfrow=c(2,2))
plot(IR_data_ID64$temperature,IR_data_ID64$growth_rate)
plot(IR_data_ID64$temperature,briere1(a=0.0001,Tmin=15,Tmax=37,temp=IR_data_ID64$temperature))
plot(IR_data_ID64$temperature,briere1(a=0.00018,Tmin=15,Tmax=36,temp=IR_data_ID64$temperature))
plot(IR_data_ID64$temperature,briere1(a=0.00015,Tmin=15,Tmax=35.5,temp=IR_data_ID64$temperature))
grid_br1_64 <- expand.grid(list(a=seq(0.0001,0.0002,by=0.00001),
                                Tmin=seq(6,20,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_64 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID64,
                      start = grid_br1_64,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_64 <- summary(fitted_br1_64)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_64 <- c(sum_br1_64$coefficients[,1])
startVals_list[[64]] <- startVals_64 #llenar lista starting values
fitted_br1_64 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID64,
                     start = startVals_64)
sum_br1_64 <- summary(fitted_br1_64)
sum_br1_64 
boot_br1_64 <- nlsBoot(fitted_br1_64, niter = 999)
coefs_64 <- data.frame(sum_br1_64$coefficients[,1:2],row.names = NULL)
boot_64 <- data.frame(boot_br1_64$estiboot)
Topt_est_64 <- Topt(Tmin=coefs_64[2,1],
                    Tmax=coefs_64[3,1],
                    m=2)
Topt_se_64 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_64), vcov(fitted_br1_64))

params_br1_64 <- data.frame(coefs_64[1,],
                            coefs_64[2,],
                            coefs_64[3,],
                            Topt_est_64,
                            Topt_se_64,
                            boot_64[1,],
                            boot_64[2,],
                            boot_64[3,],
                            startVals_64[1],
                            startVals_64[2],
                            startVals_64[3],
                            "a,Tmax")
colnames(params_br1_64) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_64

#### _ _ Study 65 ####
IR_data_ID65 <- IR_data_all %>%
  filter(id==65)
par(mfrow=c(2,2))
plot(IR_data_ID65$temperature,IR_data_ID65$growth_rate)
plot(IR_data_ID65$temperature,briere1(a=0.0001,Tmin=12,Tmax=41,temp=IR_data_ID65$temperature))
plot(IR_data_ID65$temperature,briere1(a=0.00009,Tmin=12,Tmax=41,temp=IR_data_ID65$temperature))
grid_br1_65 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(33,46,by=0.5)))
fitted_br1_65 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID65,
                      start = grid_br1_65,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_65 <- summary(fitted_br1_65)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_65 <- c(sum_br1_65$coefficients[,1])
startVals_list[[65]] <- startVals_65 #llenar lista starting values
fitted_br1_65 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID65,
                     start = startVals_65)
sum_br1_65 <- summary(fitted_br1_65)
sum_br1_65 
boot_br1_65 <- nlsBoot(fitted_br1_65, niter = 999)
coefs_65 <- data.frame(sum_br1_65$coefficients[,1:2],row.names = NULL)
boot_65 <- data.frame(boot_br1_65$estiboot)
Topt_est_65 <- Topt(Tmin=coefs_65[2,1],
                    Tmax=coefs_65[3,1],
                    m=2)
Topt_se_65 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_65), vcov(fitted_br1_65))

params_br1_65 <- data.frame(coefs_65[1,],
                            coefs_65[2,],
                            coefs_65[3,],
                            Topt_est_65,
                            Topt_se_65,
                            boot_65[1,],
                            boot_65[2,],
                            boot_65[3,],
                            startVals_65[1],
                            startVals_65[2],
                            startVals_65[3],
                            "all")
colnames(params_br1_65) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_65

#### _ _ Study 66 ####
IR_data_ID66 <- IR_data_all %>%
  filter(id==66)
par(mfrow=c(2,2))
plot(IR_data_ID66$temperature,IR_data_ID66$growth_rate)
plot(IR_data_ID66$temperature,briere1(a=0.0001,Tmin=14,Tmax=34,temp=IR_data_ID66$temperature))
plot(IR_data_ID66$temperature,briere1(a=0.00017,Tmin=14,Tmax=33,temp=IR_data_ID66$temperature))
plot(IR_data_ID66$temperature,briere1(a=0.00022,Tmin=16,Tmax=33.5,temp=IR_data_ID66$temperature))
grid_br1_66 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_66 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID66,
                      start = grid_br1_66,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_66 <- summary(fitted_br1_66)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_66 <- c(sum_br1_66$coefficients[,1])
startVals_list[[66]] <- startVals_66 #llenar lista starting values
fitted_br1_66 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID66,
                     start = startVals_66)
sum_br1_66 <- summary(fitted_br1_66)
sum_br1_66 
boot_br1_66 <- nlsBoot(fitted_br1_66, niter = 999)
coefs_66 <- data.frame(sum_br1_66$coefficients[,1:2],row.names = NULL)
boot_66 <- data.frame(boot_br1_66$estiboot)
Topt_est_66 <- Topt(Tmin=coefs_66[2,1],
                    Tmax=coefs_66[3,1],
                    m=2)
Topt_se_66 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_66), vcov(fitted_br1_66))

params_br1_66 <- data.frame(coefs_66[1,],
                            coefs_66[2,],
                            coefs_66[3,],
                            Topt_est_66,
                            Topt_se_66,
                            boot_66[1,],
                            boot_66[2,],
                            boot_66[3,],
                            startVals_66[1],
                            startVals_66[2],
                            startVals_66[3],
                            "all")
colnames(params_br1_66) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_66

#### _ _ Study 67 ####
IR_data_ID67 <- IR_data_all %>%
  filter(id==67)
par(mfrow=c(2,2))
#NO

#### _ _ Study 68 ####
IR_data_ID68 <- IR_data_all %>%
  filter(id==68)
par(mfrow=c(2,2))
plot(IR_data_ID68$temperature,IR_data_ID68$growth_rate)
plot(IR_data_ID68$temperature,briere1(a=0.0001,Tmin=14,Tmax=35.5,temp=IR_data_ID68$temperature))
plot(IR_data_ID68$temperature,briere1(a=0.00012,Tmin=14,Tmax=34,temp=IR_data_ID68$temperature))
grid_br1_68 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_68 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID68,
                      start = grid_br1_68,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_68 <- summary(fitted_br1_68)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_68 <- c(sum_br1_68$coefficients[,1])
startVals_list[[68]] <- startVals_68 #llenar lista starting values
fitted_br1_68 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID68,
                     start = startVals_68)
sum_br1_68 <- summary(fitted_br1_68)
sum_br1_68 
boot_br1_68 <- nlsBoot(fitted_br1_68, niter = 999)
coefs_68 <- data.frame(sum_br1_68$coefficients[,1:2],row.names = NULL)
boot_68 <- data.frame(boot_br1_68$estiboot)
Topt_est_68 <- Topt(Tmin=coefs_68[2,1],
                    Tmax=coefs_68[3,1],
                    m=2)
Topt_se_68 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_68), vcov(fitted_br1_68))

params_br1_68 <- data.frame(coefs_68[1,],
                            coefs_68[2,],
                            coefs_68[3,],
                            Topt_est_68,
                            Topt_se_68,
                            boot_68[1,],
                            boot_68[2,],
                            boot_68[3,],
                            startVals_68[1],
                            startVals_68[2],
                            startVals_68[3],
                            "all")
colnames(params_br1_68) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_68

#### _ _ Study 69 ####
IR_data_ID69 <- IR_data_all %>%
  filter(id==69)
par(mfrow=c(2,2))
plot(IR_data_ID69$temperature,IR_data_ID69$growth_rate)
plot(IR_data_ID69$temperature,briere1(a=0.0002,Tmin=14,Tmax=34,temp=IR_data_ID69$temperature))
plot(IR_data_ID69$temperature,briere1(a=0.0003,Tmin=14,Tmax=34,temp=IR_data_ID69$temperature))
plot(IR_data_ID69$temperature,briere1(a=0.00025,Tmin=8,Tmax=34,temp=IR_data_ID69$temperature))
grid_br1_69 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(0,15,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_69 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID69,
                      start = grid_br1_69,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_69 <- summary(fitted_br1_69)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_69 <- c(sum_br1_69$coefficients[,1])
startVals_list[[69]] <- startVals_69 #llenar lista starting values
fitted_br1_69 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID69,
                     start = startVals_69)
sum_br1_69 <- summary(fitted_br1_69)
sum_br1_69 
boot_br1_69 <- nlsBoot(fitted_br1_69, niter = 999)
coefs_69 <- data.frame(sum_br1_69$coefficients[,1:2],row.names = NULL)
boot_69 <- data.frame(boot_br1_69$estiboot)
Topt_est_69 <- Topt(Tmin=coefs_69[2,1],
                    Tmax=coefs_69[3,1],
                    m=2)
Topt_se_69 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_69), vcov(fitted_br1_69))

params_br1_69 <- data.frame(coefs_69[1,],
                            coefs_69[2,],
                            coefs_69[3,],
                            Topt_est_69,
                            Topt_se_69,
                            NA,NA,
                            NA,NA,
                            NA,NA,
                            startVals_69[1],
                            startVals_69[2],
                            startVals_69[3],
                            "all")
colnames(params_br1_69) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_69

#### _ _ Study 70 ####
IR_data_ID70 <- IR_data_all %>%
  filter(id==70)
par(mfrow=c(2,2))
plot(IR_data_ID70$temperature,IR_data_ID70$growth_rate)
plot(IR_data_ID70$temperature,briere1(a=0.00015,Tmin=8,Tmax=27,temp=IR_data_ID70$temperature))
plot(IR_data_ID70$temperature,briere1(a=0.0002,Tmin=8,Tmax=26.5,temp=IR_data_ID70$temperature))
plot(IR_data_ID70$temperature,briere1(a=0.0003,Tmin=8,Tmax=26.5,temp=IR_data_ID70$temperature))
grid_br1_70 <- expand.grid(list(a=seq(0.00010,0.00020,by=0.00001),
                                Tmin=seq(-5,15,by=0.5),
                                Tmax=seq(22,34,by=0.5)))
fitted_br1_70 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID70,
                      start = grid_br1_70,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_70 <- summary(fitted_br1_70)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_70 <- c(sum_br1_70$coefficients[,1])
startVals_list[[70]] <- startVals_70 #llenar lista starting values
fitted_br1_70 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID70,
                     start = startVals_70)
sum_br1_70 <- summary(fitted_br1_70)
sum_br1_70 
boot_br1_70 <- nlsBoot(fitted_br1_70, niter = 999)
coefs_70 <- data.frame(sum_br1_70$coefficients[,1:2],row.names = NULL)
boot_70 <- data.frame(boot_br1_70$estiboot)
Topt_est_70 <- Topt(Tmin=coefs_70[2,1],
                    Tmax=coefs_70[3,1],
                    m=2)
Topt_se_70 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_70), vcov(fitted_br1_70))

params_br1_70 <- data.frame(coefs_70[1,],
                            coefs_70[2,],
                            coefs_70[3,],
                            Topt_est_70,
                            Topt_se_70,
                            boot_70[1,],
                            boot_70[2,],
                            boot_70[3,],
                            startVals_70[1],
                            startVals_70[2],
                            startVals_70[3],
                            "Tmax")
colnames(params_br1_70) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_70

#### _ _ Study 71 ####
IR_data_ID71 <- IR_data_all %>%
  filter(id==71)
par(mfrow=c(2,2))
plot(IR_data_ID71$temperature,IR_data_ID71$growth_rate)
plot(IR_data_ID71$temperature,briere1(a=0.0001,Tmin=18,Tmax=37,temp=IR_data_ID71$temperature))
plot(IR_data_ID71$temperature,briere1(a=0.0002,Tmin=17,Tmax=38,temp=IR_data_ID71$temperature))
grid_br1_71 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(10,25,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_71 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID71,
                      start = grid_br1_71,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_71 <- summary(fitted_br1_71)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_71 <- c(sum_br1_71$coefficients[,1])
startVals_list[[71]] <- startVals_71 #llenar lista starting values
fitted_br1_71 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID71,
                     start = startVals_71)
sum_br1_71 <- summary(fitted_br1_71)
sum_br1_71 
boot_br1_71 <- nlsBoot(fitted_br1_71, niter = 999)
coefs_71 <- data.frame(sum_br1_71$coefficients[,1:2],row.names = NULL)
boot_71 <- data.frame(boot_br1_71$estiboot)
Topt_est_71 <- Topt(Tmin=coefs_71[2,1],
                    Tmax=coefs_71[3,1],
                    m=2)
Topt_se_71 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_71), vcov(fitted_br1_71))

params_br1_71 <- data.frame(coefs_71[1,],
                            coefs_71[2,],
                            coefs_71[3,],
                            Topt_est_71,
                            Topt_se_71,
                            boot_71[1,],
                            boot_71[2,],
                            boot_71[3,],
                            startVals_71[1],
                            startVals_71[2],
                            startVals_71[3],
                            "Tmax")
colnames(params_br1_71) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_71

#### _ _ Study 72 ####
IR_data_ID72 <- IR_data_all %>%
  filter(id==72)
par(mfrow=c(2,2))
plot(IR_data_ID72$temperature,IR_data_ID72$growth_rate)
plot(IR_data_ID72$temperature,briere1(a=0.0001,Tmin=16,Tmax=33,temp=IR_data_ID72$temperature))
plot(IR_data_ID72$temperature,briere1(a=0.00012,Tmin=15,Tmax=32.5,temp=IR_data_ID72$temperature))
plot(IR_data_ID72$temperature,briere1(a=0.0002,Tmin=15,Tmax=32.5,temp=IR_data_ID72$temperature))
grid_br1_72 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_72 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID72,
                      start = grid_br1_72,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_72 <- summary(fitted_br1_72)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_72 <- c(sum_br1_72$coefficients[,1])
startVals_list[[72]] <- startVals_72 #llenar lista starting values
fitted_br1_72 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID72,
                     start = startVals_72)
sum_br1_72 <- summary(fitted_br1_72)
sum_br1_72 
boot_br1_72 <- nlsBoot(fitted_br1_72, niter = 999)
coefs_72 <- data.frame(sum_br1_72$coefficients[,1:2],row.names = NULL)
boot_72 <- data.frame(boot_br1_72$estiboot)
Topt_est_72 <- Topt(Tmin=coefs_72[2,1],
                    Tmax=coefs_72[3,1],
                    m=2)
Topt_se_72 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_72), vcov(fitted_br1_72))

params_br1_72 <- data.frame(coefs_72[1,],
                            coefs_72[2,],
                            coefs_72[3,],
                            Topt_est_72,
                            Topt_se_72,
                            boot_72[1,],
                            boot_72[2,],
                            boot_72[3,],
                            startVals_72[1],
                            startVals_72[2],
                            startVals_72[3],
                            "all")
colnames(params_br1_72) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_72

#### _ _ Study 73 ####
IR_data_ID73 <- IR_data_all %>%
  filter(id==73)
par(mfrow=c(2,2))
plot(IR_data_ID73$temperature,IR_data_ID73$growth_rate)
plot(IR_data_ID73$temperature,briere1(a=0.0002,Tmin=14,Tmax=35,temp=IR_data_ID73$temperature))
grid_br1_73 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,42,by=0.5)))
fitted_br1_73 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID73,
                      start = grid_br1_73,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_73 <- summary(fitted_br1_73)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_73 <- c(sum_br1_73$coefficients[,1])
startVals_list[[73]] <- startVals_73 #llenar lista starting values
fitted_br1_73 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID73,
                     start = startVals_73)
sum_br1_73 <- summary(fitted_br1_73)
sum_br1_73 
boot_br1_73 <- nlsBoot(fitted_br1_73, niter = 999)
coefs_73 <- data.frame(sum_br1_73$coefficients[,1:2],row.names = NULL)
boot_73 <- data.frame(boot_br1_73$estiboot)
Topt_est_73 <- Topt(Tmin=coefs_73[2,1],
                    Tmax=coefs_73[3,1],
                    m=2)
Topt_se_73 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_73), vcov(fitted_br1_73))

params_br1_73 <- data.frame(coefs_73[1,],
                            coefs_73[2,],
                            coefs_73[3,],
                            Topt_est_73,
                            Topt_se_73,
                            boot_73[1,],
                            boot_73[2,],
                            boot_73[3,],
                            startVals_73[1],
                            startVals_73[2],
                            startVals_73[3],
                            "Tmin,Tmax")
colnames(params_br1_73) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_73

#### _ _ Study 74 ####
IR_data_ID74 <- IR_data_all %>%
  filter(id==74)
par(mfrow=c(2,2))
plot(IR_data_ID74$temperature,IR_data_ID74$growth_rate)
plot(IR_data_ID74$temperature,briere1(a=0.0002,Tmin=8,Tmax=32,temp=IR_data_ID74$temperature))
plot(IR_data_ID74$temperature,briere1(a=0.0002,Tmin=8,Tmax=33.5,temp=IR_data_ID74$temperature))
plot(IR_data_ID74$temperature,briere1(a=0.00018,Tmin=8,Tmax=33.5,temp=IR_data_ID74$temperature))
grid_br1_74 <- expand.grid(list(a=seq(0.00010,0.00025,by=0.00001),
                                Tmin=seq(0,15,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_74 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID74,
                      start = grid_br1_74,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_74 <- summary(fitted_br1_74)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_74 <- c(sum_br1_74$coefficients[,1])
startVals_list[[74]] <- startVals_74 #llenar lista starting values
fitted_br1_74 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID74,
                     start = startVals_74)
sum_br1_74 <- summary(fitted_br1_74)
sum_br1_74 
boot_br1_74 <- nlsBoot(fitted_br1_74, niter = 999)
coefs_74 <- data.frame(sum_br1_74$coefficients[,1:2],row.names = NULL)
boot_74 <- data.frame(boot_br1_74$estiboot)
Topt_est_74 <- Topt(Tmin=coefs_74[2,1],
                    Tmax=coefs_74[3,1],
                    m=2)
Topt_se_74 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_74), vcov(fitted_br1_74))

params_br1_74 <- data.frame(coefs_74[1,],
                            coefs_74[2,],
                            coefs_74[3,],
                            Topt_est_74,
                            Topt_se_74,
                            boot_74[1,],
                            boot_74[2,],
                            boot_74[3,],
                            startVals_74[1],
                            startVals_74[2],
                            startVals_74[3],
                            "a,Tmax")
colnames(params_br1_74) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_74

#### _ _ Study 75 ####
IR_data_ID75 <- IR_data_all %>%
  filter(id==75)
par(mfrow=c(2,2))
plot(IR_data_ID75$temperature,IR_data_ID75$growth_rate)
plot(IR_data_ID75$temperature,briere1(a=0.0002,Tmin=15,Tmax=36,temp=IR_data_ID75$temperature))
grid_br1_75 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,42,by=0.5)))
fitted_br1_75 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID75,
                      start = grid_br1_75,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_75 <- summary(fitted_br1_75)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_75 <- c(sum_br1_75$coefficients[,1])
startVals_list[[75]] <- startVals_75 #llenar lista starting values
fitted_br1_75 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID75,
                     start = startVals_75)
# NO

#### _ _ Study 76 ####
IR_data_ID76 <- IR_data_all %>%
  filter(id==76)
par(mfrow=c(2,2))
plot(IR_data_ID76$temperature,IR_data_ID76$growth_rate)
plot(IR_data_ID76$temperature,briere1(a=0.0002,Tmin=14,Tmax=36,temp=IR_data_ID76$temperature))
plot(IR_data_ID76$temperature,briere1(a=0.00023,Tmin=14,Tmax=36,temp=IR_data_ID76$temperature))
grid_br1_76 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_76 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID76,
                      start = grid_br1_76,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_76 <- summary(fitted_br1_76)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_76 <- c(sum_br1_76$coefficients[,1])
startVals_list[[76]] <- startVals_76 #llenar lista starting values
fitted_br1_76 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID76,
                     start = startVals_76)
sum_br1_76 <- summary(fitted_br1_76)
sum_br1_76 
boot_br1_76 <- nlsBoot(fitted_br1_76, niter = 999)
coefs_76 <- data.frame(sum_br1_76$coefficients[,1:2],row.names = NULL)
boot_76 <- data.frame(boot_br1_76$estiboot)
Topt_est_76 <- Topt(Tmin=coefs_76[2,1],
                    Tmax=coefs_76[3,1],
                    m=2)
Topt_se_76 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_76), vcov(fitted_br1_76))

params_br1_76 <- data.frame(coefs_76[1,],
                            coefs_76[2,],
                            coefs_76[3,],
                            Topt_est_76,
                            Topt_se_76,
                            boot_76[1,],
                            boot_76[2,],
                            boot_76[3,],
                            startVals_76[1],
                            startVals_76[2],
                            startVals_76[3],
                            "none")
colnames(params_br1_76) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_76

#### _ _ Study 77 ####
IR_data_ID77 <- IR_data_all %>%
  filter(id==77)
par(mfrow=c(2,2))
plot(IR_data_ID77$temperature,IR_data_ID77$growth_rate)
plot(IR_data_ID77$temperature,briere1(a=0.00015,Tmin=14,Tmax=32,temp=IR_data_ID77$temperature))
plot(IR_data_ID77$temperature,briere1(a=0.0002,Tmin=14,Tmax=33,temp=IR_data_ID77$temperature))
plot(IR_data_ID77$temperature,briere1(a=0.00024,Tmin=14,Tmax=33,temp=IR_data_ID77$temperature))
grid_br1_77 <- expand.grid(list(a=seq(0.00015,0.0003,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_77 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID77,
                      start = grid_br1_77,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_77 <- summary(fitted_br1_77)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_77 <- c(sum_br1_77$coefficients[,1])
startVals_list[[77]] <- startVals_77 #llenar lista starting values
fitted_br1_77 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID77,
                     start = startVals_77)
sum_br1_77 <- summary(fitted_br1_77)
sum_br1_77 
boot_br1_77 <- nlsBoot(fitted_br1_77, niter = 999)
coefs_77 <- data.frame(sum_br1_77$coefficients[,1:2],row.names = NULL)
boot_77 <- data.frame(boot_br1_77$estiboot)
Topt_est_77 <- Topt(Tmin=coefs_77[2,1],
                    Tmax=coefs_77[3,1],
                    m=2)
Topt_se_77 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_77), vcov(fitted_br1_77))

params_br1_77 <- data.frame(coefs_77[1,],
                            coefs_77[2,],
                            coefs_77[3,],
                            Topt_est_77,
                            Topt_se_77,
                            boot_77[1,],
                            boot_77[2,],
                            boot_77[3,],
                            startVals_77[1],
                            startVals_77[2],
                            startVals_77[3],
                            "all")
colnames(params_br1_77) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_77

#### _ _ Study 78 ####
IR_data_ID78 <- IR_data_all %>%
  filter(id==78)
par(mfrow=c(2,2))
plot(IR_data_ID78$temperature,IR_data_ID78$growth_rate)
#NO

#### _ _ Study 79 ####
IR_data_ID79 <- IR_data_all %>%
  filter(id==79)
par(mfrow=c(2,2))
plot(IR_data_ID79$temperature,IR_data_ID79$growth_rate)
plot(IR_data_ID79$temperature,briere1(a=0.0002,Tmin=12,Tmax=31,temp=IR_data_ID79$temperature))
plot(IR_data_ID79$temperature,briere1(a=0.00025,Tmin=12,Tmax=31,temp=IR_data_ID79$temperature))
plot(IR_data_ID79$temperature,briere1(a=0.00024,Tmin=15,Tmax=35,temp=IR_data_ID79$temperature))
grid_br1_79 <- expand.grid(list(a=seq(0.0002,0.0003,by=0.00001),
                                Tmin=seq(5,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_79 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID79,
                      start = grid_br1_79,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_79 <- summary(fitted_br1_79)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_79 <- c(sum_br1_79$coefficients[,1])
startVals_list[[79]] <- startVals_79 #llenar lista starting values
fitted_br1_79 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID79,
                     start = startVals_79)
sum_br1_79 <- summary(fitted_br1_79)
sum_br1_79 
boot_br1_79 <- nlsBoot(fitted_br1_79, niter = 999)
coefs_79 <- data.frame(sum_br1_79$coefficients[,1:2],row.names = NULL)
boot_79 <- data.frame(boot_br1_79$estiboot)
Topt_est_79 <- Topt(Tmin=coefs_79[2,1],
                    Tmax=coefs_79[3,1],
                    m=2)
Topt_se_79 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_79), vcov(fitted_br1_79))

params_br1_79 <- data.frame(coefs_79[1,],
                            coefs_79[2,],
                            coefs_79[3,],
                            Topt_est_79,
                            Topt_se_79,
                            boot_79[1,],
                            boot_79[2,],
                            boot_79[3,],
                            startVals_79[1],
                            startVals_79[2],
                            startVals_79[3],
                            "a,Tmax")
colnames(params_br1_79) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_79

#### _ _ Study 80 ####
IR_data_ID80 <- IR_data_all %>%
  filter(id==80)
par(mfrow=c(2,2))
plot(IR_data_ID80$temperature,IR_data_ID80$growth_rate)
plot(IR_data_ID80$temperature,briere1(a=0.0001,Tmin=15,Tmax=32,temp=IR_data_ID80$temperature))
plot(IR_data_ID80$temperature,briere1(a=0.0002,Tmin=14,Tmax=32,temp=IR_data_ID80$temperature))
plot(IR_data_ID80$temperature,briere1(a=0.00015,Tmin=14,Tmax=32,temp=IR_data_ID80$temperature))
grid_br1_80 <- expand.grid(list(a=seq(0.00005,0.0002,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_80 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID80,
                      start = grid_br1_80,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_80 <- summary(fitted_br1_80)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_80 <- c(sum_br1_80$coefficients[,1])
startVals_list[[80]] <- startVals_80 #llenar lista starting values
fitted_br1_80 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID80,
                     start = startVals_80)
sum_br1_80 <- summary(fitted_br1_80)
sum_br1_80 
boot_br1_80 <- nlsBoot(fitted_br1_80, niter = 999)
coefs_80 <- data.frame(sum_br1_80$coefficients[,1:2],row.names = NULL)
boot_80 <- data.frame(boot_br1_80$estiboot)
Topt_est_80 <- Topt(Tmin=coefs_80[2,1],
                    Tmax=coefs_80[3,1],
                    m=2)
Topt_se_80 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_80), vcov(fitted_br1_80))

params_br1_80 <- data.frame(coefs_80[1,],
                            coefs_80[2,],
                            coefs_80[3,],
                            Topt_est_80,
                            Topt_se_80,
                            boot_80[1,],
                            boot_80[2,],
                            boot_80[3,],
                            startVals_80[1],
                            startVals_80[2],
                            startVals_80[3],
                            "Tmax")
colnames(params_br1_80) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_80

#### _ _ Study 81 ####
IR_data_ID81 <- IR_data_all %>%
  filter(id==81)
par(mfrow=c(2,2))
plot(IR_data_ID81$temperature,IR_data_ID81$growth_rate)
plot(IR_data_ID81$temperature,briere1(a=0.0001,Tmin=15,Tmax=33,temp=IR_data_ID81$temperature))
plot(IR_data_ID81$temperature,briere1(a=0.0001,Tmin=13,Tmax=32,temp=IR_data_ID81$temperature))
plot(IR_data_ID81$temperature,briere1(a=0.00015,Tmin=13,Tmax=32,temp=IR_data_ID81$temperature))
grid_br1_81 <- expand.grid(list(a=seq(0.00005,0.0002,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_81 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID81,
                      start = grid_br1_81,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_81 <- summary(fitted_br1_81)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_81 <- c(sum_br1_81$coefficients[,1])
startVals_list[[81]] <- startVals_81 #llenar lista starting values
# no sale fitted_br1_81 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID81,
                     start = list(a=0.0001,Tmin=15,Tmax=34))
sum_br1_81 <- summary(fitted_br1_81)
sum_br1_81 
boot_br1_81 <- nlsBoot(fitted_br1_81, niter = 999)
coefs_81 <- data.frame(sum_br1_81$coefficients[,1:2],row.names = NULL)
boot_81 <- data.frame(boot_br1_81$estiboot)
Topt_est_81 <- Topt(Tmin=coefs_81[2,1],
                    Tmax=coefs_81[3,1],
                    m=2)
Topt_se_81 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_81), vcov(fitted_br1_81))

params_br1_81 <- data.frame(coefs_81[1,],
                            coefs_81[2,],
                            coefs_81[3,],
                            Topt_est_81,
                            Topt_se_81,
                            boot_81[1,],
                            boot_81[2,],
                            boot_81[3,],
                            startVals_81[1],
                            startVals_81[2],
                            startVals_81[3],
                            "a,Tmax")
colnames(params_br1_81) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_81

#### _ _ Study 82 #### 
#RETIRADO

#### _ _ Study 83 ####
IR_data_ID83 <- IR_data_all %>%
  filter(id==83)
par(mfrow=c(2,2))
plot(IR_data_ID83$temperature,IR_data_ID83$growth_rate)
plot(IR_data_ID83$temperature,briere1(a=0.0003,Tmin=14,Tmax=38,temp=IR_data_ID83$temperature))
plot(IR_data_ID83$temperature,briere1(a=0.0002,Tmin=14,Tmax=38,temp=IR_data_ID83$temperature))
plot(IR_data_ID83$temperature,briere1(a=0.00025,Tmin=14,Tmax=38,temp=IR_data_ID83$temperature))
grid_br1_83 <- expand.grid(list(a=seq(0.0002,0.0003,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_83 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID83,
                      start = grid_br1_83,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_83 <- summary(fitted_br1_83)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_83 <- c(sum_br1_83$coefficients[,1])
startVals_list[[83]] <- startVals_83 #llenar lista starting values
fitted_br1_83 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID83,
                     start = startVals_83)
sum_br1_83 <- summary(fitted_br1_83)
sum_br1_83 
boot_br1_83 <- nlsBoot(fitted_br1_83, niter = 999)
coefs_83 <- data.frame(sum_br1_83$coefficients[,1:2],row.names = NULL)
boot_83 <- data.frame(boot_br1_83$estiboot)
Topt_est_83 <- Topt(Tmin=coefs_83[2,1],
                    Tmax=coefs_83[3,1],
                    m=2)
Topt_se_83 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_83), vcov(fitted_br1_83))

params_br1_83 <- data.frame(coefs_83[1,],
                            coefs_83[2,],
                            coefs_83[3,],
                            Topt_est_83,
                            Topt_se_83,
                            boot_83[1,],
                            boot_83[2,],
                            boot_83[3,],
                            startVals_83[1],
                            startVals_83[2],
                            startVals_83[3],
                            "all")
colnames(params_br1_83) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_83

#### _ _ Study 84 ####
IR_data_ID84 <- IR_data_all %>%
  filter(id==84)
par(mfrow=c(2,2))
plot(IR_data_ID84$temperature,IR_data_ID84$growth_rate)
plot(IR_data_ID84$temperature,briere1(a=0.0003,Tmin=8,Tmax=28,temp=IR_data_ID84$temperature))
plot(IR_data_ID84$temperature,briere1(a=0.0006,Tmin=8,Tmax=28,temp=IR_data_ID84$temperature))
grid_br1_84 <- expand.grid(list(a=seq(0.0003,0.0009,by=0.0001),
                                Tmin=seq(0,15,by=0.5),
                                Tmax=seq(22,32,by=0.5)))
fitted_br1_84 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID84,
                      start = grid_br1_84,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_84 <- summary(fitted_br1_84)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_84 <- c(sum_br1_84$coefficients[,1])
startVals_list[[84]] <- startVals_84 #llenar lista starting values
fitted_br1_84 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID84,
                     start = startVals_84)
sum_br1_84 <- summary(fitted_br1_84)
sum_br1_84 
boot_br1_84 <- nlsBoot(fitted_br1_84, niter = 999)
coefs_84 <- data.frame(sum_br1_84$coefficients[,1:2],row.names = NULL)
boot_84 <- data.frame(boot_br1_84$estiboot)
Topt_est_84 <- Topt(Tmin=coefs_84[2,1],
                    Tmax=coefs_84[3,1],
                    m=2)
Topt_se_84 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_84), vcov(fitted_br1_84))

params_br1_84 <- data.frame(coefs_84[1,],
                            coefs_84[2,],
                            coefs_84[3,],
                            Topt_est_84,
                            Topt_se_84,
                            boot_84[1,],
                            boot_84[2,],
                            boot_84[3,],
                            startVals_84[1],
                            startVals_84[2],
                            startVals_84[3],
                            "Tmax")
colnames(params_br1_84) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_84

#### _ _ Study 85 ####
IR_data_ID85 <- IR_data_all %>%
  filter(id==85)
par(mfrow=c(2,2))
plot(IR_data_ID85$temperature,IR_data_ID85$growth_rate)
plot(IR_data_ID85$temperature,briere1(a=0.0001,Tmin=14,Tmax=32,temp=IR_data_ID85$temperature))
plot(IR_data_ID85$temperature,briere1(a=0.0001,Tmin=14,Tmax=31.5,temp=IR_data_ID85$temperature))
plot(IR_data_ID85$temperature,briere1(a=0.00015,Tmin=14,Tmax=31.5,temp=IR_data_ID85$temperature))
grid_br1_85 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(6,20,by=0.5),
                                Tmax=seq(28,38,by=0.5)))
fitted_br1_85 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID85,
                      start = grid_br1_85,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_85 <- summary(fitted_br1_85)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_85 <- c(sum_br1_85$coefficients[,1])
startVals_list[[85]] <- startVals_85 #llenar lista starting values
fitted_br1_85 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID85,
                     start = startVals_85)
sum_br1_85 <- summary(fitted_br1_85)
sum_br1_85 
boot_br1_85 <- nlsBoot(fitted_br1_85, niter = 999)
coefs_85 <- data.frame(sum_br1_85$coefficients[,1:2],row.names = NULL)
boot_85 <- data.frame(boot_br1_85$estiboot)
Topt_est_85 <- Topt(Tmin=coefs_85[2,1],
                    Tmax=coefs_85[3,1],
                    m=2)
Topt_se_85 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_85), vcov(fitted_br1_85))

params_br1_85 <- data.frame(coefs_85[1,],
                            coefs_85[2,],
                            coefs_85[3,],
                            Topt_est_85,
                            Topt_se_85,
                            boot_85[1,],
                            boot_85[2,],
                            boot_85[3,],
                            startVals_85[1],
                            startVals_85[2],
                            startVals_85[3],
                            "all")
colnames(params_br1_85) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_85

#### _ _ Study 86 ####
IR_data_ID86 <- IR_data_all %>%
  filter(id==86)
par(mfrow=c(2,2))
plot(IR_data_ID86$temperature,IR_data_ID86$growth_rate)
plot(IR_data_ID86$temperature,briere1(a=0.0001,Tmin=12,Tmax=32,temp=IR_data_ID86$temperature))
plot(IR_data_ID86$temperature,briere1(a=0.0001,Tmin=11,Tmax=32,temp=IR_data_ID86$temperature))
plot(IR_data_ID86$temperature,briere1(a=0.00012,Tmin=11,Tmax=33,temp=IR_data_ID86$temperature))
grid_br1_86 <- expand.grid(list(a=seq(0.00008,0.00018,by=0.00001),
                                Tmin=seq(5,15,by=0.5),
                                Tmax=seq(26,36,by=0.5)))
fitted_br1_86 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID86,
                      start = grid_br1_86,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_86 <- summary(fitted_br1_86)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_86 <- c(sum_br1_86$coefficients[,1])
startVals_list[[86]] <- startVals_86 #llenar lista starting values
fitted_br1_86 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID86,
                     start = startVals_86)
sum_br1_86 <- summary(fitted_br1_86)
sum_br1_86 
boot_br1_86 <- nlsBoot(fitted_br1_86, niter = 999)
coefs_86 <- data.frame(sum_br1_86$coefficients[,1:2],row.names = NULL)
boot_86 <- data.frame(boot_br1_86$estiboot)
Topt_est_86 <- Topt(Tmin=coefs_86[2,1],
                    Tmax=coefs_86[3,1],
                    m=2)
Topt_se_86 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_86), vcov(fitted_br1_86))

params_br1_86 <- data.frame(coefs_86[1,],
                            coefs_86[2,],
                            coefs_86[3,],
                            Topt_est_86,
                            Topt_se_86,
                            boot_86[1,],
                            boot_86[2,],
                            boot_86[3,],
                            startVals_86[1],
                            startVals_86[2],
                            startVals_86[3],
                            "all")
colnames(params_br1_86) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_86

#### _ _ Study 87 ####
IR_data_ID87 <- IR_data_all %>%
  filter(id==87)
par(mfrow=c(2,2))
plot(IR_data_ID87$temperature,IR_data_ID87$growth_rate)
plot(IR_data_ID87$temperature,briere1(a=0.0003,Tmin=5,Tmax=35,temp=IR_data_ID87$temperature))
grid_br1_87 <- expand.grid(list(a=seq(0.0002,0.0004,by=0.00001),
                                Tmin=seq(0,15,by=0.5),
                                Tmax=seq(35,45,by=0.5)))
fitted_br1_87 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID87,
                      start = grid_br1_87,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_87 <- summary(fitted_br1_87)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_87 <- c(sum_br1_87$coefficients[,1])
startVals_list[[87]] <- startVals_87 #llenar lista starting values
fitted_br1_87 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID87,
                     start = startVals_87)
#no sale

#### _ _ Study 88 ####
IR_data_ID88 <- IR_data_all %>%
  filter(id==88)
par(mfrow=c(2,2))
plot(IR_data_ID88$temperature,IR_data_ID88$growth_rate)
plot(IR_data_ID88$temperature,briere1(a=0.0003,Tmin=12,Tmax=30,temp=IR_data_ID88$temperature))
plot(IR_data_ID88$temperature,briere1(a=0.00035,Tmin=10,Tmax=29,temp=IR_data_ID88$temperature))
plot(IR_data_ID88$temperature,briere1(a=0.00035,Tmin=8,Tmax=29,temp=IR_data_ID88$temperature))
grid_br1_88 <- expand.grid(list(a=seq(0.0002,0.0004,by=0.00001),
                                Tmin=seq(0,20,by=0.5),
                                Tmax=seq(22,35,by=0.5)))
fitted_br1_88 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID88,
                      start = grid_br1_88,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_88 <- summary(fitted_br1_88)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_88 <- c(sum_br1_88$coefficients[,1])
startVals_list[[88]] <- startVals_88 #llenar lista starting values
fitted_br1_88 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID88,
                     start = startVals_88)
sum_br1_88 <- summary(fitted_br1_88)
sum_br1_88 
boot_br1_88 <- nlsBoot(fitted_br1_88, niter = 999)
coefs_88 <- data.frame(sum_br1_88$coefficients[,1:2],row.names = NULL)
boot_88 <- data.frame(boot_br1_88$estiboot)
Topt_est_88 <- Topt(Tmin=coefs_88[2,1],
                    Tmax=coefs_88[3,1],
                    m=2)
Topt_se_88 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_88), vcov(fitted_br1_88))

params_br1_88 <- data.frame(coefs_88[1,],
                            coefs_88[2,],
                            coefs_88[3,],
                            Topt_est_88,
                            Topt_se_88,
                            boot_88[1,],
                            boot_88[2,],
                            boot_88[3,],
                            startVals_88[1],
                            startVals_88[2],
                            startVals_88[3],
                            "a,Tmax")
colnames(params_br1_88) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_88

#### _ _ Study 89 ####
IR_data_ID89 <- IR_data_all %>%
  filter(id==89)
par(mfrow=c(2,2))
plot(IR_data_ID89$temperature,IR_data_ID89$growth_rate)
plot(IR_data_ID89$temperature,briere1(a=0.0001,Tmin=17,Tmax=27,temp=IR_data_ID89$temperature))
plot(IR_data_ID89$temperature,briere1(a=0.0002,Tmin=17,Tmax=26,temp=IR_data_ID89$temperature))
plot(IR_data_ID89$temperature,briere1(a=0.00025,Tmin=17,Tmax=26,temp=IR_data_ID89$temperature))
grid_br1_89 <- expand.grid(list(a=seq(0.00015,0.0003,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(22,32,by=0.5)))
fitted_br1_89 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID89,
                      start = grid_br1_89,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_89 <- summary(fitted_br1_89)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_89 <- c(sum_br1_89$coefficients[,1])
startVals_list[[89]] <- startVals_89 #llenar lista starting values
fitted_br1_89 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID89,
                     start = startVals_89)
#no

#### _ _ Study 90 ####
IR_data_ID90 <- IR_data_all %>%
  filter(id==90)
par(mfrow=c(2,2))
plot(IR_data_ID90$temperature,IR_data_ID90$growth_rate)
plot(IR_data_ID90$temperature,briere1(a=0.0001,Tmin=19,Tmax=34,temp=IR_data_ID90$temperature))
plot(IR_data_ID90$temperature,briere1(a=0.00022,Tmin=18,Tmax=35,temp=IR_data_ID90$temperature))
plot(IR_data_ID90$temperature,briere1(a=0.00017,Tmin=18,Tmax=34.5,temp=IR_data_ID90$temperature))
grid_br1_90 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(12,25,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_90 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID90,
                      start = grid_br1_90,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_90 <- summary(fitted_br1_90)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_90 <- c(sum_br1_90$coefficients[,1])
startVals_list[[90]] <- startVals_90 #llenar lista starting values
fitted_br1_90 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID90,
                     start = startVals_90)
sum_br1_90 <- summary(fitted_br1_90)
sum_br1_90 
boot_br1_90 <- nlsBoot(fitted_br1_90, niter = 999)
coefs_90 <- data.frame(sum_br1_90$coefficients[,1:2],row.names = NULL)
boot_90 <- data.frame(boot_br1_90$estiboot)
Topt_est_90 <- Topt(Tmin=coefs_90[2,1],
                    Tmax=coefs_90[3,1],
                    m=2)
Topt_se_90 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_90), vcov(fitted_br1_90))

params_br1_90 <- data.frame(coefs_90[1,],
                            coefs_90[2,],
                            coefs_90[3,],
                            Topt_est_90,
                            Topt_se_90,
                            boot_90[1,],
                            boot_90[2,],
                            boot_90[3,],
                            startVals_90[1],
                            startVals_90[2],
                            startVals_90[3],
                            "Tmax")
colnames(params_br1_90) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_90

#### _ _ Study 91 ####
IR_data_ID91 <- IR_data_all %>%
  filter(id==91)
par(mfrow=c(2,2))
plot(IR_data_ID91$temperature,IR_data_ID91$growth_rate)
plot(IR_data_ID91$temperature,briere1(a=0.0002,Tmin=10,Tmax=33,temp=IR_data_ID91$temperature))
plot(IR_data_ID91$temperature,briere1(a=0.0002,Tmin=10,Tmax=35,temp=IR_data_ID91$temperature))
plot(IR_data_ID91$temperature,briere1(a=0.0002,Tmin=10,Tmax=34,temp=IR_data_ID91$temperature))
grid_br1_91 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(5,15,by=0.5),
                                Tmax=seq(25,42,by=0.5)))
fitted_br1_91 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID91,
                      start = grid_br1_91,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_91 <- summary(fitted_br1_91)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_91 <- c(sum_br1_91$coefficients[,1])
startVals_list[[91]] <- startVals_91 #llenar lista starting values
fitted_br1_91 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID91,
                     start = startVals_91)
sum_br1_91 <- summary(fitted_br1_91)
sum_br1_91 
boot_br1_91 <- nlsBoot(fitted_br1_91, niter = 999)
coefs_91 <- data.frame(sum_br1_91$coefficients[,1:2],row.names = NULL)
boot_91 <- data.frame(boot_br1_91$estiboot)
Topt_est_91 <- Topt(Tmin=coefs_91[2,1],
                    Tmax=coefs_91[3,1],
                    m=2)
Topt_se_91 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_91), vcov(fitted_br1_91))

params_br1_91 <- data.frame(coefs_91[1,],
                            coefs_91[2,],
                            coefs_91[3,],
                            Topt_est_91,
                            Topt_se_91,
                            boot_91[1,],
                            boot_91[2,],
                            boot_91[3,],
                            startVals_91[1],
                            startVals_91[2],
                            startVals_91[3],
                            "all")
colnames(params_br1_91) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_91

#### _ _ Study 92 ####
IR_data_ID92 <- IR_data_all %>%
  filter(id==92)
par(mfrow=c(2,2))
plot(IR_data_ID92$temperature,IR_data_ID92$growth_rate)
plot(IR_data_ID92$temperature,briere1(a=0.0001,Tmin=13,Tmax=33,temp=IR_data_ID92$temperature))
plot(IR_data_ID92$temperature,briere1(a=0.00007,Tmin=13,Tmax=33,temp=IR_data_ID92$temperature))
grid_br1_92 <- expand.grid(list(a=seq(0.00004,0.0001,by=0.00001),
                                Tmin=seq(5,17,by=0.5),
                                Tmax=seq(25,37,by=0.5)))
fitted_br1_92 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID92,
                      start = grid_br1_92,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_92 <- summary(fitted_br1_92)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_92 <- c(sum_br1_92$coefficients[,1])
startVals_list[[92]] <- startVals_92 #llenar lista starting values
fitted_br1_92 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID92,
                     start = startVals_92)
sum_br1_92 <- summary(fitted_br1_92)
sum_br1_92 
boot_br1_92 <- nlsBoot(fitted_br1_92, niter = 999)
coefs_92 <- data.frame(sum_br1_92$coefficients[,1:2],row.names = NULL)
boot_92 <- data.frame(boot_br1_92$estiboot)
Topt_est_92 <- Topt(Tmin=coefs_92[2,1],
                    Tmax=coefs_92[3,1],
                    m=2)
Topt_se_92 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_92), vcov(fitted_br1_92))

params_br1_92 <- data.frame(coefs_92[1,],
                            coefs_92[2,],
                            coefs_92[3,],
                            Topt_est_92,
                            Topt_se_92,
                            boot_92[1,],
                            boot_92[2,],
                            boot_92[3,],
                            startVals_92[1],
                            startVals_92[2],
                            startVals_92[3],
                            "Tmax")
colnames(params_br1_92) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_92

#### _ _ Study 93 ####
IR_data_ID93 <- IR_data_all %>%
  filter(id==93)
par(mfrow=c(2,2))
plot(IR_data_ID93$temperature,IR_data_ID93$growth_rate)
plot(IR_data_ID93$temperature,briere1(a=0.0002,Tmin=15,Tmax=33,temp=IR_data_ID93$temperature))
plot(IR_data_ID93$temperature,briere1(a=0.00016,Tmin=15,Tmax=33.5,temp=IR_data_ID93$temperature))
plot(IR_data_ID93$temperature,briere1(a=0.0001,Tmin=10,Tmax=33.5,temp=IR_data_ID93$temperature))
grid_br1_93 <- expand.grid(list(a=seq(0.00005,0.0002,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(30,40,by=0.5)))
fitted_br1_93 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID93,
                      start = grid_br1_93,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_93 <- summary(fitted_br1_93)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_93 <- c(sum_br1_93$coefficients[,1])
startVals_list[[93]] <- startVals_93 #llenar lista starting values
fitted_br1_93 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID93,
                     start = startVals_93)
sum_br1_93 <- summary(fitted_br1_93)
sum_br1_93 
boot_br1_93 <- nlsBoot(fitted_br1_93, niter = 999)
coefs_93 <- data.frame(sum_br1_93$coefficients[,1:2],row.names = NULL)
boot_93 <- data.frame(boot_br1_93$estiboot)
Topt_est_93 <- Topt(Tmin=coefs_93[2,1],
                    Tmax=coefs_93[3,1],
                    m=2)
Topt_se_93 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_93), vcov(fitted_br1_93))

params_br1_93 <- data.frame(coefs_93[1,],
                            coefs_93[2,],
                            coefs_93[3,],
                            Topt_est_93,
                            Topt_se_93,
                            boot_93[1,],
                            boot_93[2,],
                            boot_93[3,],
                            startVals_93[1],
                            startVals_93[2],
                            startVals_93[3],
                            "all")
colnames(params_br1_93) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_93

#### _ _ Study 94 ####
IR_data_ID94 <- IR_data_all %>%
  filter(id==94)
par(mfrow=c(2,2))
plot(IR_data_ID94$temperature,IR_data_ID94$growth_rate)
plot(IR_data_ID94$temperature,briere1(a=0.0003,Tmin=13,Tmax=33,temp=IR_data_ID94$temperature))
plot(IR_data_ID94$temperature,briere1(a=0.0003,Tmin=13,Tmax=35,temp=IR_data_ID94$temperature))
grid_br1_94 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(8,18,by=0.5),
                                Tmax=seq(28,42,by=0.5)))
fitted_br1_94 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID94,
                      start = grid_br1_94,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_94 <- summary(fitted_br1_94)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_94 <- c(sum_br1_94$coefficients[,1])
startVals_list[[94]] <- startVals_94 #llenar lista starting values
fitted_br1_94 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID94,
                     start = startVals_94)
sum_br1_94 <- summary(fitted_br1_94)
sum_br1_94 
boot_br1_94 <- nlsBoot(fitted_br1_94, niter = 999)
coefs_94 <- data.frame(sum_br1_94$coefficients[,1:2],row.names = NULL)
boot_94 <- data.frame(boot_br1_94$estiboot)
Topt_est_94 <- Topt(Tmin=coefs_94[2,1],
                    Tmax=coefs_94[3,1],
                    m=2)
Topt_se_94 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_94), vcov(fitted_br1_94))

params_br1_94 <- data.frame(coefs_94[1,],
                            coefs_94[2,],
                            coefs_94[3,],
                            Topt_est_94,
                            Topt_se_94,
                            boot_94[1,],
                            boot_94[2,],
                            boot_94[3,],
                            startVals_94[1],
                            startVals_94[2],
                            startVals_94[3],
                            "none")
colnames(params_br1_94) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_94

#### _ _ Study 95 ####
IR_data_ID95 <- IR_data_all %>%
  filter(id==95)
par(mfrow=c(2,2))
plot(IR_data_ID95$temperature,IR_data_ID95$growth_rate)
plot(IR_data_ID95$temperature,briere1(a=0.0002,Tmin=20,Tmax=36,temp=IR_data_ID95$temperature))
plot(IR_data_ID95$temperature,briere1(a=0.0003,Tmin=20,Tmax=37,temp=IR_data_ID95$temperature))
plot(IR_data_ID95$temperature,briere1(a=0.00035,Tmin=19,Tmax=37,temp=IR_data_ID95$temperature))
grid_br1_95 <- expand.grid(list(a=seq(0.0001,0.0004,by=0.00001),
                                Tmin=seq(15,25,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_95 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID95,
                      start = grid_br1_95,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_95 <- summary(fitted_br1_95)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_95 <- c(sum_br1_95$coefficients[,1])
startVals_list[[95]] <- startVals_95 #llenar lista starting values
fitted_br1_95 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID95,
                     start = startVals_95)
#noo

#### _ _ Study 96 ####
IR_data_ID96 <- IR_data_all %>%
  filter(id==96)
par(mfrow=c(2,2))
plot(IR_data_ID96$temperature,IR_data_ID96$growth_rate)
plot(IR_data_ID96$temperature,briere1(a=0.00006,Tmin=15,Tmax=30.5,temp=IR_data_ID96$temperature))
grid_br1_96 <- expand.grid(list(a=seq(0.00001,0.0001,by=0.00001),
                                Tmin=seq(8,18,by=0.5),
                                Tmax=seq(25,35,by=0.5)))
fitted_br1_96 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID96,
                      start = grid_br1_96,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_96 <- summary(fitted_br1_96)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_96 <- c(sum_br1_96$coefficients[,1])
startVals_list[[96]] <- startVals_96 #llenar lista starting values
fitted_br1_96 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID96,
                     start = startVals_96)
sum_br1_96 <- summary(fitted_br1_96)
sum_br1_96 
boot_br1_96 <- nlsBoot(fitted_br1_96, niter = 999)
coefs_96 <- data.frame(sum_br1_96$coefficients[,1:2],row.names = NULL)
boot_96 <- data.frame(boot_br1_96$estiboot)
Topt_est_96 <- Topt(Tmin=coefs_96[2,1],
                    Tmax=coefs_96[3,1],
                    m=2)
Topt_se_96 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_96), vcov(fitted_br1_96))

params_br1_96 <- data.frame(coefs_96[1,],
                            coefs_96[2,],
                            coefs_96[3,],
                            Topt_est_96,
                            Topt_se_96,
                            NA,NA,
                            NA,NA,
                            NA,NA,
                            startVals_96[1],
                            startVals_96[2],
                            startVals_96[3],
                            "Tmax")
colnames(params_br1_96) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_96

#### _ _ Study 97 ####
IR_data_ID97 <- IR_data_all %>%
  filter(id==97)
par(mfrow=c(2,2))
#no

#### _ _ Study 98 ####
IR_data_ID98 <- IR_data_all %>%
  filter(id==98)
par(mfrow=c(2,2))
plot(IR_data_ID98$temperature,IR_data_ID98$growth_rate)
plot(IR_data_ID98$temperature,briere1(a=0.0002,Tmin=12,Tmax=37,temp=IR_data_ID98$temperature))
plot(IR_data_ID98$temperature,briere1(a=0.0003,Tmin=12,Tmax=37,temp=IR_data_ID98$temperature))
plot(IR_data_ID98$temperature,briere1(a=0.00025,Tmin=12,Tmax=37,temp=IR_data_ID98$temperature))
grid_br1_98 <- expand.grid(list(a=seq(0.0002,0.0003,by=0.00001),
                                Tmin=seq(5,15,by=0.5),
                                Tmax=seq(30,42,by=0.5)))
fitted_br1_98 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID98,
                      start = grid_br1_98,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_98 <- summary(fitted_br1_98)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_98 <- c(sum_br1_98$coefficients[,1])
startVals_list[[98]] <- startVals_98 #llenar lista starting values
fitted_br1_98 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID98,
                     start = startVals_98,
                     nls.control(warnOnly=TRUE))
sum_br1_98 <- summary(fitted_br1_98)
sum_br1_98 
boot_br1_98 <- nlsBoot(fitted_br1_98, niter = 999)
coefs_98 <- data.frame(sum_br1_98$coefficients[,1:2],row.names = NULL)
boot_98 <- data.frame(boot_br1_98$estiboot)
Topt_est_98 <- Topt(Tmin=coefs_98[2,1],
                    Tmax=coefs_98[3,1],
                    m=2)
Topt_se_98 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_98), vcov(fitted_br1_98))

params_br1_98 <- data.frame(coefs_98[1,],
                            coefs_98[2,],
                            coefs_98[3,],
                            Topt_est_98,
                            Topt_se_98,
                            NA,NA,
                            NA,NA,
                            NA,NA,
                            startVals_98[1],
                            startVals_98[2],
                            startVals_98[3],
                            "all")
colnames(params_br1_98) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_98

#### _ _ Study 99 ####
IR_data_ID99 <- IR_data_all %>%
  filter(id==99)
par(mfrow=c(2,2))
plot(IR_data_ID99$temperature,IR_data_ID99$growth_rate)
plot(IR_data_ID99$temperature,briere1(a=0.0002,Tmin=17,Tmax=35,temp=IR_data_ID99$temperature))
plot(IR_data_ID99$temperature,briere1(a=0.00025,Tmin=18,Tmax=34,temp=IR_data_ID99$temperature))
plot(IR_data_ID99$temperature,briere1(a=0.0003,Tmin=19,Tmax=34,temp=IR_data_ID99$temperature))
grid_br1_99 <- expand.grid(list(a=seq(0.0002,0.0004,by=0.00001),
                                Tmin=seq(12,24,by=0.5),
                                Tmax=seq(28,40,by=0.5)))
fitted_br1_99 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID99,
                      start = grid_br1_99,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_99 <- summary(fitted_br1_99)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_99 <- c(sum_br1_99$coefficients[,1])
startVals_list[[99]] <- startVals_99 #llenar lista starting values
fitted_br1_99 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID99,
                     start = startVals_99)
#no

#### _ _ Study 100 ####
IR_data_ID100 <- IR_data_all %>%
  filter(id==100)
par(mfrow=c(2,2))
plot(IR_data_ID100$temperature,IR_data_ID100$growth_rate)
plot(IR_data_ID100$temperature,briere1(a=0.0002,Tmin=12,Tmax=30.5,temp=IR_data_ID100$temperature))
plot(IR_data_ID100$temperature,briere1(a=0.0004,Tmin=12,Tmax=30.5,temp=IR_data_ID100$temperature))
plot(IR_data_ID100$temperature,briere1(a=0.0004,Tmin=14,Tmax=30.5,temp=IR_data_ID100$temperature))
grid_br1_100 <- expand.grid(list(a=seq(0.0002,0.0004,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(27,40,by=0.5)))
fitted_br1_100 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID100,
                      start = grid_br1_100,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_100 <- summary(fitted_br1_100)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_100 <- c(sum_br1_100$coefficients[,1])
startVals_list[[100]] <- startVals_100 #llenar lista starting values
fitted_br1_100 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID100,
                     start = startVals_100)
sum_br1_100 <- summary(fitted_br1_100)
sum_br1_100 
boot_br1_100 <- nlsBoot(fitted_br1_100, niter = 999)
coefs_100 <- data.frame(sum_br1_100$coefficients[,1:2],row.names = NULL)
boot_100 <- data.frame(boot_br1_100$estiboot)
Topt_est_100 <- Topt(Tmin=coefs_100[2,1],
                    Tmax=coefs_100[3,1],
                    m=2)
Topt_se_100 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_100), vcov(fitted_br1_100))

params_br1_100 <- data.frame(coefs_100[1,],
                            coefs_100[2,],
                            coefs_100[3,],
                            Topt_est_100,
                            Topt_se_100,
                            boot_100[1,],
                            boot_100[2,],
                            boot_100[3,],
                            startVals_100[1],
                            startVals_100[2],
                            startVals_100[3],
                            "Tmax")
colnames(params_br1_100) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_100

#### _ _ Study 101 ####
IR_data_ID101 <- IR_data_all %>%
  filter(id==101)
par(mfrow=c(2,2))
plot(IR_data_ID101$temperature,IR_data_ID101$growth_rate)
#no

#### _ _ Study 102 ####
IR_data_ID102 <- IR_data_all %>%
  filter(id==102)
par(mfrow=c(2,2))
plot(IR_data_ID102$temperature,IR_data_ID102$growth_rate)
#NO

#### _ _ Study 103 ####
IR_data_ID103 <- IR_data_all %>%
  filter(id==103)
par(mfrow=c(2,2))
plot(IR_data_ID103$temperature,IR_data_ID103$growth_rate)
plot(IR_data_ID103$temperature,briere1(a=0.0002,Tmin=12,Tmax=28,temp=IR_data_ID103$temperature))
plot(IR_data_ID103$temperature,briere1(a=0.0002,Tmin=11,Tmax=28,temp=IR_data_ID103$temperature))
plot(IR_data_ID103$temperature,briere1(a=0.0004,Tmin=8,Tmax=28.5,temp=IR_data_ID103$temperature))
grid_br1_103 <- expand.grid(list(a=seq(0.0002,0.0005,by=0.00001),
                                Tmin=seq(5,15,by=0.5),
                                Tmax=seq(24,34,by=0.5)))
fitted_br1_103 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID103,
                      start = grid_br1_103,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_103 <- summary(fitted_br1_103)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_103 <- c(sum_br1_103$coefficients[,1])
startVals_list[[103]] <- startVals_103 #llenar lista starting values
fitted_br1_103 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID103,
                     start = startVals_103)
plot(fitted_br1_103)
sum_br1_103 <- summary(fitted_br1_103)
sum_br1_103 
boot_br1_103 <- nlsBoot(fitted_br1_103, niter = 999)
coefs_103 <- data.frame(sum_br1_103$coefficients[,1:2],row.names = NULL)
boot_103 <- data.frame(boot_br1_103$estiboot)
Topt_est_103 <- Topt(Tmin=coefs_103[2,1],
                    Tmax=coefs_103[3,1],
                    m=2)
Topt_se_103 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_103), vcov(fitted_br1_103))

params_br1_103 <- data.frame(coefs_103[1,],
                            coefs_103[2,],
                            coefs_103[3,],
                            Topt_est_103,
                            Topt_se_103,
                            boot_103[1,],
                            boot_103[2,],
                            boot_103[3,],
                            startVals_103[1],
                            startVals_103[2],
                            startVals_103[3],
                            "all")
colnames(params_br1_103) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_103

#### _ _ Study 104 ####
IR_data_ID104 <- IR_data_all %>%
  filter(id==104)
par(mfrow=c(2,2))
plot(IR_data_ID104$temperature,IR_data_ID104$growth_rate)
plot(IR_data_ID104$temperature,briere1(a=0.0001,Tmin=14,Tmax=31,temp=IR_data_ID104$temperature))
plot(IR_data_ID104$temperature,briere1(a=0.0001,Tmin=14,Tmax=30.5,temp=IR_data_ID104$temperature))
grid_br1_104 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                 Tmin=seq(8,18,by=0.5),
                                 Tmax=seq(24,34,by=0.5)))
fitted_br1_104 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                       data = IR_data_ID104,
                       start = grid_br1_104,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_br1_104 <- summary(fitted_br1_104)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_104 <- c(sum_br1_104$coefficients[,1])
startVals_list[[104]] <- startVals_104 #llenar lista starting values
fitted_br1_104 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID104,
                      start = startVals_104)
sum_br1_104 <- summary(fitted_br1_104)
sum_br1_104 
boot_br1_104 <- nlsBoot(fitted_br1_104, niter = 999)
coefs_104 <- data.frame(sum_br1_104$coefficients[,1:2],row.names = NULL)
boot_104 <- data.frame(boot_br1_104$estiboot)
Topt_est_104 <- Topt(Tmin=coefs_104[2,1],
                     Tmax=coefs_104[3,1],
                     m=2)
Topt_se_104 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                           coef(fitted_br1_104), vcov(fitted_br1_104))

params_br1_104 <- data.frame(coefs_104[1,],
                             coefs_104[2,],
                             coefs_104[3,],
                             Topt_est_104,
                             Topt_se_104,
                             boot_104[1,],
                             boot_104[2,],
                             boot_104[3,],
                             startVals_104[1],
                             startVals_104[2],
                             startVals_104[3],
                             "all")
colnames(params_br1_104) <- c("a_est_nls","a_se_nls",
                              "Tmin_est_nls","Tmin_se_nls",
                              "Tmax_est_nls","Tmax_se_nls",
                              "Topt_est","Topt_se_delta",
                              "a_est_boot","a_se_boot",
                              "Tmin_est_boot","Tmin_se_boot",
                              "Tmax_est_boot","Tmax_se_boot",
                              "starting_a","starting_Tmin",
                              "starting_Tmax","which_signif")
params_br1_104

#### _ _ Study 105 ####
IR_data_ID105 <- IR_data_all %>%
  filter(id==105)
par(mfrow=c(2,2))
plot(IR_data_ID105$temperature,IR_data_ID105$growth_rate)
plot(IR_data_ID105$temperature,briere1(a=0.0001,Tmin=19,Tmax=39,temp=IR_data_ID105$temperature))
plot(IR_data_ID105$temperature,briere1(a=0.00013,Tmin=19,Tmax=39,temp=IR_data_ID105$temperature))
plot(IR_data_ID105$temperature,briere1(a=0.00011,Tmin=19,Tmax=39,temp=IR_data_ID105$temperature))
grid_br1_105 <- expand.grid(list(a=seq(0.00005,0.00019,by=0.00001),
                                 Tmin=seq(12,24,by=0.5),
                                 Tmax=seq(30,44,by=0.5)))
fitted_br1_105 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                       data = IR_data_ID105,
                       start = grid_br1_105,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_br1_105 <- summary(fitted_br1_105)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_105 <- c(sum_br1_105$coefficients[,1])
startVals_list[[105]] <- startVals_105 #llenar lista starting values
fitted_br1_105 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID105,
                      start = startVals_105)
sum_br1_105 <- summary(fitted_br1_105)
sum_br1_105 
boot_br1_105 <- nlsBoot(fitted_br1_105, niter = 999)
coefs_105 <- data.frame(sum_br1_105$coefficients[,1:2],row.names = NULL)
boot_105 <- data.frame(boot_br1_105$estiboot)
Topt_est_105 <- Topt(Tmin=coefs_105[2,1],
                     Tmax=coefs_105[3,1],
                     m=2)
Topt_se_105 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                           coef(fitted_br1_105), vcov(fitted_br1_105))

params_br1_105 <- data.frame(coefs_105[1,],
                             coefs_105[2,],
                             coefs_105[3,],
                             Topt_est_105,
                             Topt_se_105,
                             boot_105[1,],
                             boot_105[2,],
                             boot_105[3,],
                             startVals_105[1],
                             startVals_105[2],
                             startVals_105[3],
                             "Tmin,Tmax")
colnames(params_br1_105) <- c("a_est_nls","a_se_nls",
                              "Tmin_est_nls","Tmin_se_nls",
                              "Tmax_est_nls","Tmax_se_nls",
                              "Topt_est","Topt_se_delta",
                              "a_est_boot","a_se_boot",
                              "Tmin_est_boot","Tmin_se_boot",
                              "Tmax_est_boot","Tmax_se_boot",
                              "starting_a","starting_Tmin",
                              "starting_Tmax","which_signif")
params_br1_105

#### _ _ Study 106 ####
IR_data_ID106 <- IR_data_all %>%
  filter(id==106)
par(mfrow=c(2,2))
plot(IR_data_ID106$temperature,IR_data_ID106$growth_rate)
plot(IR_data_ID106$temperature,briere1(a=0.0001,Tmin=12,Tmax=32,temp=IR_data_ID106$temperature))
plot(IR_data_ID106$temperature,briere1(a=0.0008,Tmin=12,Tmax=30,temp=IR_data_ID106$temperature))
plot(IR_data_ID106$temperature,briere1(a=0.0009,Tmin=10,Tmax=29,temp=IR_data_ID106$temperature))
grid_br1_106 <- expand.grid(list(a=seq(0.0001,0.0012,by=0.0001),
                                 Tmin=seq(5,15,by=0.5),
                                 Tmax=seq(22,35,by=0.5)))
fitted_br1_106 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                       data = IR_data_ID106,
                       start = grid_br1_106,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_br1_106 <- summary(fitted_br1_106)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_106 <- c(sum_br1_106$coefficients[,1])
startVals_list[[106]] <- startVals_106 #llenar lista starting values
fitted_br1_106 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID106,
                      start = startVals_106)
sum_br1_106 <- summary(fitted_br1_106)
sum_br1_106 
boot_br1_106 <- nlsBoot(fitted_br1_106, niter = 999)
coefs_106 <- data.frame(sum_br1_106$coefficients[,1:2],row.names = NULL)
boot_106 <- data.frame(boot_br1_106$estiboot)
Topt_est_106 <- Topt(Tmin=coefs_106[2,1],
                     Tmax=coefs_106[3,1],
                     m=2)
Topt_se_106 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                           coef(fitted_br1_106), vcov(fitted_br1_106))

params_br1_106 <- data.frame(coefs_106[1,],
                             coefs_106[2,],
                             coefs_106[3,],
                             Topt_est_106,
                             Topt_se_106,
                             boot_106[1,],
                             boot_106[2,],
                             boot_106[3,],
                             startVals_106[1],
                             startVals_106[2],
                             startVals_106[3],
                             "Tmax")
colnames(params_br1_106) <- c("a_est_nls","a_se_nls",
                              "Tmin_est_nls","Tmin_se_nls",
                              "Tmax_est_nls","Tmax_se_nls",
                              "Topt_est","Topt_se_delta",
                              "a_est_boot","a_se_boot",
                              "Tmin_est_boot","Tmin_se_boot",
                              "Tmax_est_boot","Tmax_se_boot",
                              "starting_a","starting_Tmin",
                              "starting_Tmax","which_signif")
params_br1_106

#### _ _ Study 107 ####
IR_data_ID107 <- IR_data_all %>%
  filter(id==107)
par(mfrow=c(2,2))
plot(IR_data_ID107$temperature,IR_data_ID107$growth_rate)
plot(IR_data_ID107$temperature,briere1(a=0.0001,Tmin=14,Tmax=33,temp=IR_data_ID107$temperature))
plot(IR_data_ID107$temperature,briere1(a=0.00018,Tmin=14,Tmax=32,temp=IR_data_ID107$temperature))
grid_br1_107 <- expand.grid(list(a=seq(0.00005,0.0002,by=0.00001),
                                 Tmin=seq(10,20,by=0.5),
                                 Tmax=seq(28,38,by=0.5)))
fitted_br1_107 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                       data = IR_data_ID107,
                       start = grid_br1_107,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_br1_107 <- summary(fitted_br1_107)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_107 <- c(sum_br1_107$coefficients[,1])
startVals_list[[107]] <- startVals_107 #llenar lista starting values
fitted_br1_107 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID107,
                      start = startVals_107)
sum_br1_107 <- summary(fitted_br1_107)
sum_br1_107 
boot_br1_107 <- nlsBoot(fitted_br1_107, niter = 999)
coefs_107 <- data.frame(sum_br1_107$coefficients[,1:2],row.names = NULL)
boot_107 <- data.frame(boot_br1_107$estiboot)
Topt_est_107 <- Topt(Tmin=coefs_107[2,1],
                     Tmax=coefs_107[3,1],
                     m=2)
Topt_se_107 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                           coef(fitted_br1_107), vcov(fitted_br1_107))

params_br1_107 <- data.frame(coefs_107[1,],
                             coefs_107[2,],
                             coefs_107[3,],
                             Topt_est_107,
                             Topt_se_107,
                             boot_107[1,],
                             boot_107[2,],
                             boot_107[3,],
                             startVals_107[1],
                             startVals_107[2],
                             startVals_107[3],
                             "Tmax")
colnames(params_br1_107) <- c("a_est_nls","a_se_nls",
                              "Tmin_est_nls","Tmin_se_nls",
                              "Tmax_est_nls","Tmax_se_nls",
                              "Topt_est","Topt_se_delta",
                              "a_est_boot","a_se_boot",
                              "Tmin_est_boot","Tmin_se_boot",
                              "Tmax_est_boot","Tmax_se_boot",
                              "starting_a","starting_Tmin",
                              "starting_Tmax","which_signif")
params_br1_107

#### _ _ Study 108 ####
IR_data_ID108 <- IR_data_all %>%
  filter(id==108)
par(mfrow=c(2,2))
plot(IR_data_ID108$temperature,IR_data_ID108$growth_rate)
#no

#### _ _ Study 109 ####
IR_data_ID109 <- IR_data_all %>%
  filter(id==109)
par(mfrow=c(2,2))
plot(IR_data_ID109$temperature,IR_data_ID109$growth_rate)
plot(IR_data_ID109$temperature,briere1(a=0.0001,Tmin=19,Tmax=35.5,temp=IR_data_ID109$temperature))
plot(IR_data_ID109$temperature,briere1(a=0.00012,Tmin=19,Tmax=35.5,temp=IR_data_ID109$temperature))
plot(IR_data_ID109$temperature,briere1(a=0.00015,Tmin=19,Tmax=35.5,temp=IR_data_ID109$temperature))
grid_br1_109 <- expand.grid(list(a=seq(0.0001,0.0002,by=0.00001),
                                 Tmin=seq(12,24,by=0.5),
                                 Tmax=seq(30,44,by=0.5)))
fitted_br1_109 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                       data = IR_data_ID109,
                       start = grid_br1_109,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_br1_109 <- summary(fitted_br1_109)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_109 <- c(sum_br1_109$coefficients[,1])
startVals_list[[109]] <- startVals_109 #llenar lista starting values
fitted_br1_109 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID109,
                      start = startVals_109)
sum_br1_109 <- summary(fitted_br1_109)
sum_br1_109 
boot_br1_109 <- nlsBoot(fitted_br1_109, niter = 999)
coefs_109 <- data.frame(sum_br1_109$coefficients[,1:2],row.names = NULL)
boot_109 <- data.frame(boot_br1_109$estiboot)
Topt_est_109 <- Topt(Tmin=coefs_109[2,1],
                     Tmax=coefs_109[3,1],
                     m=2)
Topt_se_109 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                           coef(fitted_br1_109), vcov(fitted_br1_109))

params_br1_109 <- data.frame(coefs_109[1,],
                             coefs_109[2,],
                             coefs_109[3,],
                             Topt_est_109,
                             Topt_se_109,
                             boot_109[1,],
                             boot_109[2,],
                             boot_109[3,],
                             startVals_109[1],
                             startVals_109[2],
                             startVals_109[3],
                             "all")
colnames(params_br1_109) <- c("a_est_nls","a_se_nls",
                              "Tmin_est_nls","Tmin_se_nls",
                              "Tmax_est_nls","Tmax_se_nls",
                              "Topt_est","Topt_se_delta",
                              "a_est_boot","a_se_boot",
                              "Tmin_est_boot","Tmin_se_boot",
                              "Tmax_est_boot","Tmax_se_boot",
                              "starting_a","starting_Tmin",
                              "starting_Tmax","which_signif")
params_br1_109


#### _ _ Study 110 ####
IR_data_ID110 <- IR_data_all %>%
  filter(id==110)
par(mfrow=c(2,2))
plot(IR_data_ID110$temperature,IR_data_ID110$growth_rate)
plot(IR_data_ID110$temperature,briere1(a=0.0003,Tmin=12,Tmax=40,temp=IR_data_ID110$temperature))
plot(IR_data_ID110$temperature,briere1(a=0.0002,Tmin=12,Tmax=40,temp=IR_data_ID110$temperature))
grid_br1_110 <- expand.grid(list(a=seq(0.00015,0.00025,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(32,42,by=0.5)))
fitted_br1_110 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID110,
                      start = grid_br1_110,
                      algorithm = "brute-force",
                      trace = TRUE)
sum_br1_110 <- summary(fitted_br1_110)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_110 <- c(sum_br1_110$coefficients[,1])
startVals_list[[110]] <- startVals_110 #llenar lista starting values
fitted_br1_110 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID110,
                     start = startVals_110)
sum_br1_110 <- summary(fitted_br1_110)
sum_br1_110 
boot_br1_110 <- nlsBoot(fitted_br1_110, niter = 999)
coefs_110 <- data.frame(sum_br1_110$coefficients[,1:2],row.names = NULL)
boot_110 <- data.frame(boot_br1_110$estiboot)
Topt_est_110 <- Topt(Tmin=coefs_110[2,1],
                    Tmax=coefs_110[3,1],
                    m=2)
Topt_se_110 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                          coef(fitted_br1_110), vcov(fitted_br1_110))

params_br1_110 <- data.frame(coefs_110[1,],
                            coefs_110[2,],
                            coefs_110[3,],
                            Topt_est_110,
                            Topt_se_110,
                            boot_110[1,],
                            boot_110[2,],
                            boot_110[3,],
                            startVals_110[1],
                            startVals_110[2],
                            startVals_110[3],
                            "all")
colnames(params_br1_110) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "starting_a","starting_Tmin",
                             "starting_Tmax","which_signif")
params_br1_110

#### _ _ Study 111 ####
IR_data_ID111 <- IR_data_all %>%
  filter(id==111)
par(mfrow=c(2,2))
plot(IR_data_ID111$temperature,IR_data_ID111$growth_rate)
plot(IR_data_ID111$temperature,briere1(a=0.00025,Tmin=12,Tmax=39,temp=IR_data_ID111$temperature))
plot(IR_data_ID111$temperature,briere1(a=0.00023,Tmin=11,Tmax=38.5,temp=IR_data_ID111$temperature))
grid_br1_111 <- expand.grid(list(a=seq(0.00017,0.00027,by=0.00001),
                                 Tmin=seq(8,18,by=0.5),
                                 Tmax=seq(32,42,by=0.5)))
fitted_br1_111 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                       data = IR_data_ID111,
                       start = grid_br1_111,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_br1_111 <- summary(fitted_br1_111)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_111 <- c(sum_br1_111$coefficients[,1])
startVals_list[[111]] <- startVals_111 #llenar lista starting values
fitted_br1_111 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID111,
                      start = startVals_111)
sum_br1_111 <- summary(fitted_br1_111)
sum_br1_111 
boot_br1_111 <- nlsBoot(fitted_br1_111, niter = 999)
coefs_111 <- data.frame(sum_br1_111$coefficients[,1:2],row.names = NULL)
boot_111 <- data.frame(boot_br1_111$estiboot)
Topt_est_111 <- Topt(Tmin=coefs_111[2,1],
                     Tmax=coefs_111[3,1],
                     m=2)
Topt_se_111 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                           coef(fitted_br1_111), vcov(fitted_br1_111))

params_br1_111 <- data.frame(coefs_111[1,],
                             coefs_111[2,],
                             coefs_111[3,],
                             Topt_est_111,
                             Topt_se_111,
                             boot_111[1,],
                             boot_111[2,],
                             boot_111[3,],
                             startVals_111[1],
                             startVals_111[2],
                             startVals_111[3],
                             "all")
colnames(params_br1_111) <- c("a_est_nls","a_se_nls",
                              "Tmin_est_nls","Tmin_se_nls",
                              "Tmax_est_nls","Tmax_se_nls",
                              "Topt_est","Topt_se_delta",
                              "a_est_boot","a_se_boot",
                              "Tmin_est_boot","Tmin_se_boot",
                              "Tmax_est_boot","Tmax_se_boot",
                              "starting_a","starting_Tmin",
                              "starting_Tmax","which_signif")
params_br1_111

#### _ _ Study 112 ####
IR_data_ID112 <- IR_data_all %>%
  filter(id==112)
par(mfrow=c(2,2))
plot(IR_data_ID112$temperature,IR_data_ID112$growth_rate)
plot(IR_data_ID112$temperature,briere1(a=0.00025,Tmin=12,Tmax=36,temp=IR_data_ID112$temperature))
plot(IR_data_ID112$temperature,briere1(a=0.0002,Tmin=11,Tmax=36.5,temp=IR_data_ID112$temperature))
plot(IR_data_ID112$temperature,briere1(a=0.00025,Tmin=11,Tmax=36.5,temp=IR_data_ID112$temperature))
grid_br1_112 <- expand.grid(list(a=seq(0.0002,0.0003,by=0.00001),
                                 Tmin=seq(8,18,by=0.5),
                                 Tmax=seq(32,42,by=0.5)))
fitted_br1_112 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                       data = IR_data_ID112,
                       start = grid_br1_112,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_br1_112 <- summary(fitted_br1_112)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_112 <- c(sum_br1_112$coefficients[,1])
startVals_list[[112]] <- startVals_112 #llenar lista starting values
fitted_br1_112 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID112,
                      start = startVals_112)
sum_br1_112 <- summary(fitted_br1_112)
sum_br1_112 
boot_br1_112 <- nlsBoot(fitted_br1_112, niter = 999)
coefs_112 <- data.frame(sum_br1_112$coefficients[,1:2],row.names = NULL)
boot_112 <- data.frame(boot_br1_112$estiboot)
Topt_est_112 <- Topt(Tmin=coefs_112[2,1],
                     Tmax=coefs_112[3,1],
                     m=2)
Topt_se_112 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                           coef(fitted_br1_112), vcov(fitted_br1_112))

params_br1_112 <- data.frame(coefs_112[1,],
                             coefs_112[2,],
                             coefs_112[3,],
                             Topt_est_112,
                             Topt_se_112,
                             boot_112[1,],
                             boot_112[2,],
                             boot_112[3,],
                             startVals_112[1],
                             startVals_112[2],
                             startVals_112[3],
                             "all")
colnames(params_br1_112) <- c("a_est_nls","a_se_nls",
                              "Tmin_est_nls","Tmin_se_nls",
                              "Tmax_est_nls","Tmax_se_nls",
                              "Topt_est","Topt_se_delta",
                              "a_est_boot","a_se_boot",
                              "Tmin_est_boot","Tmin_se_boot",
                              "Tmax_est_boot","Tmax_se_boot",
                              "starting_a","starting_Tmin",
                              "starting_Tmax","which_signif")
params_br1_112

#### _ _ Study 113 ####
IR_data_ID113 <- IR_data_all %>%
  filter(id==113)
par(mfrow=c(2,2))
plot(IR_data_ID113$temperature,IR_data_ID113$growth_rate)
plot(IR_data_ID113$temperature,briere1(a=0.0003,Tmin=13,Tmax=39,temp=IR_data_ID113$temperature))
plot(IR_data_ID113$temperature,briere1(a=0.00027,Tmin=13,Tmax=41,temp=IR_data_ID113$temperature))
plot(IR_data_ID113$temperature,briere1(a=0.00022,Tmin=13,Tmax=40.5,temp=IR_data_ID113$temperature))
grid_br1_113 <- expand.grid(list(a=seq(0.0002,0.0003,by=0.00001),
                                 Tmin=seq(8,18,by=0.5),
                                 Tmax=seq(34,44,by=0.5)))
fitted_br1_113 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                       data = IR_data_ID113,
                       start = grid_br1_113,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_br1_113 <- summary(fitted_br1_113)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_113 <- c(sum_br1_113$coefficients[,1])
startVals_list[[113]] <- startVals_113 #llenar lista starting values
fitted_br1_113 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID113,
                      start = startVals_113)
sum_br1_113 <- summary(fitted_br1_113)
sum_br1_113 
boot_br1_113 <- nlsBoot(fitted_br1_113, niter = 999)
coefs_113 <- data.frame(sum_br1_113$coefficients[,1:2],row.names = NULL)
boot_113 <- data.frame(boot_br1_113$estiboot)
Topt_est_113 <- Topt(Tmin=coefs_113[2,1],
                     Tmax=coefs_113[3,1],
                     m=2)
Topt_se_113 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                           coef(fitted_br1_113), vcov(fitted_br1_113))

params_br1_113 <- data.frame(coefs_113[1,],
                             coefs_113[2,],
                             coefs_113[3,],
                             Topt_est_113,
                             Topt_se_113,
                             boot_113[1,],
                             boot_113[2,],
                             boot_113[3,],
                             startVals_113[1],
                             startVals_113[2],
                             startVals_113[3],
                             "all")
colnames(params_br1_113) <- c("a_est_nls","a_se_nls",
                              "Tmin_est_nls","Tmin_se_nls",
                              "Tmax_est_nls","Tmax_se_nls",
                              "Topt_est","Topt_se_delta",
                              "a_est_boot","a_se_boot",
                              "Tmin_est_boot","Tmin_se_boot",
                              "Tmax_est_boot","Tmax_se_boot",
                              "starting_a","starting_Tmin",
                              "starting_Tmax","which_signif")
params_br1_113

#### _ _ Study 114 ####
IR_data_ID114 <- IR_data_all %>%
  filter(id==114)
par(mfrow=c(2,2))
plot(IR_data_ID114$temperature,IR_data_ID114$growth_rate)
plot(IR_data_ID114$temperature,briere1(a=0.00025,Tmin=13,Tmax=39,temp=IR_data_ID114$temperature))
plot(IR_data_ID114$temperature,briere1(a=0.00027,Tmin=13,Tmax=41,temp=IR_data_ID114$temperature))
plot(IR_data_ID114$temperature,briere1(a=0.00022,Tmin=13,Tmax=40.5,temp=IR_data_ID114$temperature))
grid_br1_114 <- expand.grid(list(a=seq(0.0002,0.0003,by=0.00001),
                                 Tmin=seq(8,18,by=0.5),
                                 Tmax=seq(34,44,by=0.5)))
fitted_br1_114 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                       data = IR_data_ID114,
                       start = grid_br1_114,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_br1_114 <- summary(fitted_br1_114)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_114 <- c(sum_br1_114$coefficients[,1])
startVals_list[[114]] <- startVals_114 #llenar lista starting values
fitted_br1_114 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID114,
                      start = startVals_114)
sum_br1_114 <- summary(fitted_br1_114)
sum_br1_114 
boot_br1_114 <- nlsBoot(fitted_br1_114, niter = 999)
coefs_114 <- data.frame(sum_br1_114$coefficients[,1:2],row.names = NULL)
boot_114 <- data.frame(boot_br1_114$estiboot)
Topt_est_114 <- Topt(Tmin=coefs_114[2,1],
                     Tmax=coefs_114[3,1],
                     m=2)
Topt_se_114 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                           coef(fitted_br1_114), vcov(fitted_br1_114))

params_br1_114 <- data.frame(coefs_114[1,],
                             coefs_114[2,],
                             coefs_114[3,],
                             Topt_est_114,
                             Topt_se_114,
                             boot_114[1,],
                             boot_114[2,],
                             boot_114[3,],
                             startVals_114[1],
                             startVals_114[2],
                             startVals_114[3],
                             "all")
colnames(params_br1_114) <- c("a_est_nls","a_se_nls",
                              "Tmin_est_nls","Tmin_se_nls",
                              "Tmax_est_nls","Tmax_se_nls",
                              "Topt_est","Topt_se_delta",
                              "a_est_boot","a_se_boot",
                              "Tmin_est_boot","Tmin_se_boot",
                              "Tmax_est_boot","Tmax_se_boot",
                              "starting_a","starting_Tmin",
                              "starting_Tmax","which_signif")
params_br1_114

#### _ _ Study 115 ####
IR_data_ID115 <- IR_data_all %>%
  filter(id==115)
par(mfrow=c(2,2))
plot(IR_data_ID115$temperature,IR_data_ID115$growth_rate)
plot(IR_data_ID115$temperature,briere1(a=0.0003,Tmin=13,Tmax=37,temp=IR_data_ID115$temperature))
plot(IR_data_ID115$temperature,briere1(a=0.00027,Tmin=13,Tmax=38.5,temp=IR_data_ID115$temperature))
plot(IR_data_ID115$temperature,briere1(a=0.00022,Tmin=13,Tmax=39,temp=IR_data_ID115$temperature))
grid_br1_115 <- expand.grid(list(a=seq(0.00015,0.00025,by=0.00001),
                                 Tmin=seq(8,15,by=0.5),
                                 Tmax=seq(34,44,by=0.5)))
fitted_br1_115 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                       data = IR_data_ID115,
                       start = grid_br1_115,
                       algorithm = "brute-force",
                       trace = TRUE)
sum_br1_115 <- summary(fitted_br1_115)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
startVals_115 <- c(sum_br1_115$coefficients[,1])
startVals_list[[115]] <- startVals_115 #llenar lista starting values
fitted_br1_115 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID115,
                      start = startVals_115)
sum_br1_115 <- summary(fitted_br1_115)
sum_br1_115 
boot_br1_115 <- nlsBoot(fitted_br1_115, niter = 999)
coefs_115 <- data.frame(sum_br1_115$coefficients[,1:2],row.names = NULL)
boot_115 <- data.frame(boot_br1_115$estiboot)
Topt_est_115 <- Topt(Tmin=coefs_115[2,1],
                     Tmax=coefs_115[3,1],
                     m=2)
Topt_se_115 <- deltamethod(~ ((2*2*x3+(2+1)*x2)+sqrt(4*(2^2)*(x3^2)+((2+1)^2)*(x2^2)-4*(2^2)*x2*x3))/(4*2+2), 
                           coef(fitted_br1_115), vcov(fitted_br1_115))

params_br1_115 <- data.frame(coefs_115[1,],
                             coefs_115[2,],
                             coefs_115[3,],
                             Topt_est_115,
                             Topt_se_115,
                             boot_115[1,],
                             boot_115[2,],
                             boot_115[3,],
                             startVals_115[1],
                             startVals_115[2],
                             startVals_115[3],
                             "all")
colnames(params_br1_115) <- c("a_est_nls","a_se_nls",
                              "Tmin_est_nls","Tmin_se_nls",
                              "Tmax_est_nls","Tmax_se_nls",
                              "Topt_est","Topt_se_delta",
                              "a_est_boot","a_se_boot",
                              "Tmin_est_boot","Tmin_se_boot",
                              "Tmax_est_boot","Tmax_se_boot",
                              "starting_a","starting_Tmin",
                              "starting_Tmax","which_signif")
params_br1_115


#### 11. List of starting values ####
startVals_list
startVals_df <- startVals_list %>%
  as_tibble(.name_repair = "minimal")
colnames(startVals_df) <- seq(1,109,1)
rownames(startVals_df) <- c("a","Tmin","Tmax")
startVals_df_long <- startVals_df %>%
 t()%>%
  as_tibble()%>%
  write_csv(file = "startVals.csv")

#### 12. Ensemble data frame with parameter estimates####
list_pattern <- ls(pattern = "params_br1_")#list objects with that word
df_list <- mget(list_pattern) #combines it into a df

ids_raw <- str_sub(list_pattern,12,15) #extract ids
ids <-as.numeric(ids_raw) #ordenar
parameters <- df_list %>%
  bind_rows() %>%
  bind_cols(ids) %>%
  mutate(id=...19) %>%
  select(-...19)%>%
  group_by_all()%>%
  summarise(id=sort(id))

  write_csv(file="parameters_briere1_by_study.csv")

#### 13. Merge parameters df with info####
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

