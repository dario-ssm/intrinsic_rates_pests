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
IR_data_all
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
summary(fitted_br1_1)
# probar con esos parámetros como iniciales:
fitted_br1_1 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID1,
                    start = list(a = 0.00008,
                                 Tmin =4.5,
                                 Tmax= 33.5))


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
                           "a,Tmax")
colnames(params_br1_1) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "which_signif")

params_br1_1

#### _ _ Study 2 ####
IR_data_ID2 <- IR_data_all %>%
  filter(id==2)
par(mfrow=c(2,2))
plot(IR_data_ID2$temperature,IR_data_ID2$growth_rate)
plot(IR_data_ID2$temperature,briere1(a=0.00003,Tmin=8,Tmax=35,temp=IR_data_ID2$temperature))# <-
plot(IR_data_ID2$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID2$temperature))
plot(IR_data_ID2$temperature,briere1(a=0.00001,Tmin=8,Tmax=35,temp=IR_data_ID2$temperature))

fitted_br1_2 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID2,
                    start = list(a = 0.00003,
                                 Tmin =10,
                                 Tmax= 40))
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
                           "all")
colnames(params_br1_2) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "which_signif")

params_br1_2

#### _ _ Study 3 ####
IR_data_ID3 <- IR_data_all %>%
  filter(id==3)
par(mfrow=c(2,2))
plot(IR_data_ID1$temperature,IR_data_ID1$growth_rate)
plot(IR_data_ID1$temperature,briere1(a=0.0002,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))
plot(IR_data_ID1$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))# <-
plot(IR_data_ID1$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))
grid_br1_3 <- expand.grid(list(a=seq(0.00001,0.00015,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(28,50,by=0.5)))
fitted_br1_3 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID3,
                     start = grid_br1_3,
                     algorithm = "brute-force",
                     trace = TRUE)
summary(fitted_br1_3)
fitted_br1_3 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID3,
                    start = list(a=0.00015,Tmin=6.5,Tmax=29.5))

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
                           "Tmax")
colnames(params_br1_3) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "which_signif")

params_br1_3

#### _ _ Study 4 ####
IR_data_ID4 <- IR_data_all %>%
  filter(id==4)
par(mfrow=c(2,2))
plot(IR_data_ID4$temperature,IR_data_ID4$growth_rate)
plot(IR_data_ID4$temperature,briere1(a=0.00003,Tmin=8,Tmax=35,temp=IR_data_ID4$temperature))
plot(IR_data_ID4$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID4$temperature))# <-
plot(IR_data_ID4$temperature,briere1(a=0.00001,Tmin=8,Tmax=35,temp=IR_data_ID4$temperature))

fitted_br1_4 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID4,
                    start = list(a = 0.00009,
                                 Tmin =8,
                                 Tmax= 35))

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
                           "none")
colnames(params_br1_4) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "which_signif")
params_br1_4

#### _ _ Study 5 ####
IR_data_ID5 <- IR_data_all %>%
  filter(id==5)
par(mfrow=c(2,2))
plot(IR_data_ID5$temperature,IR_data_ID5$growth_rate)
plot(IR_data_ID5$temperature,briere1(a=0.00023,Tmin=8,Tmax=35,temp=IR_data_ID5$temperature))# <-
plot(IR_data_ID5$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID5$temperature))
plot(IR_data_ID5$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID5$temperature))
grid_br1_5 <- expand.grid(list(a=seq(0.0002,0.00025,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(28,50,by=0.5)))
fitted_br1_5 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID5,
                     start = grid_br1_5,
                     algorithm = "brute-force",
                     trace = TRUE)
summary(fitted_br1_5)
#con esos valores introducimos nuevos iniciales
#fitted_br1_5 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
 #                   data = IR_data_ID5,
  #                  start = list(a=0.00024,Tmin=5.5,Tmax=33))

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
                           "all")
colnames(params_br1_5) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "Topt_est","Topt_se_delta",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot",
                            "which_signif")
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
summary(fitted_br1_6)
#con esos valores introducimos nuevos iniciales
fitted_br1_6 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID6,
                    start = list(a=0.00013,Tmin=13,Tmax=34))


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
                            "Tmin,Tmax")
colnames(params_br1_6) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_6
#### _ _ Study 7 ####
IR_data_ID7 <- IR_data_all %>%
  filter(id==7)
par(mfrow=c(2,2))
plot(IR_data_ID7$temperature,IR_data_ID7$growth_rate)
plot(IR_data_ID7$temperature,briere1(a=0.00023,Tmin=8,Tmax=35,temp=IR_data_ID7$temperature))
plot(IR_data_ID7$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID7$temperature))
plot(IR_data_ID7$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID7$temperature))# <-
grid_br1_7 <- expand.grid(list(a=seq(0.00007,0.00012,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(28,45,by=0.5)))
fitted_br1_7 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID7,
                     start = list(a=0.00009,Tmin=5,Tmax=35))

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
                            "Tmin,Tmax")
colnames(params_br1_7) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_7
#### _ _ Study 8 ####
IR_data_ID8 <- IR_data_all %>%
  filter(id==8)
par(mfrow=c(2,2))
plot(IR_data_ID8$temperature,IR_data_ID8$growth_rate)
plot(IR_data_ID8$temperature,briere1(a=0.00023,Tmin=8,Tmax=35,temp=IR_data_ID8$temperature))
plot(IR_data_ID8$temperature,briere1(a=0.00025,Tmin=8,Tmax=35,temp=IR_data_ID8$temperature))
plot(IR_data_ID8$temperature,briere1(a=0.00027,Tmin=8,Tmax=35,temp=IR_data_ID8$temperature))# <-
grid_br1_8 <- expand.grid(list(a=seq(0.00021,0.00029,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(28,45,by=0.5)))

fitted_br1_8 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID8,
                     start = grid_br1_8,
                     algorithm = "brute-force",
                     trace = TRUE)
summary(fitted_br1_8) #usamos esos como starting
#fitted_br1_8 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
 #                   data = IR_data_ID8,
  #                  start = list(a=0.00021,Tmin=1.5,Tmax=32.5))
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
                            boot_8[1,],
                            boot_8[2,],
                            boot_8[3,],
                            "Tmax")
colnames(params_br1_8) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(28,45,by=0.5)))
fitted_br1_9 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID9,
                     start = grid_br1_9,
                     algorithm = "brute-force",
                     trace = TRUE)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
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
                            "all")
colnames(params_br1_9) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_9

#### _ _ Study 10 ####
IR_data_ID10 <- IR_data_all %>%
  filter(id==10)
par(mfrow=c(2,2))
plot(IR_data_ID10$temperature,IR_data_ID10$growth_rate)
plot(IR_data_ID10$temperature,briere1(a=0.00023,Tmin=8,Tmax=35,temp=IR_data_ID10$temperature))
plot(IR_data_ID10$temperature,briere1(a=0.00015,Tmin=8,Tmax=35,temp=IR_data_ID10$temperature))
plot(IR_data_ID10$temperature,briere1(a=0.00015,Tmin=16,Tmax=37,temp=IR_data_ID10$temperature))# <-
grid_br1_10 <- expand.grid(list(a=seq(0.00010,0.00020,by=0.00001),
                               Tmin=seq(8,20,by=0.5),
                               Tmax=seq(28,40,by=0.5)))
fitted_br1_10 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID10,
                     start = grid_br1_10,
                     algorithm = "brute-force",
                     trace = TRUE)
summary(fitted_br1_10)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_10 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID10,
                    start = list(a=0.00020,Tmin=19.5,Tmax=35.5))

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
                            "Tmax")
colnames(params_br1_10) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
                               Tmin=seq(10,16,by=0.5),
                               Tmax=seq(28,35,by=0.5)))
fitted_br1_11 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID11,
                     start = grid_br1_11,
                     algorithm = "brute-force",
                     trace = TRUE)
summary(fitted_br1_11)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_11 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID11,
                    start = list(a=0.00019,Tmin=10,Tmax=31))

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
                            "all")
colnames(params_br1_11) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_12)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_12 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID12,
                     start = list(a=0.00018,Tmin=13,Tmax=38))

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
                            "all")
colnames(params_br1_12) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
                                Tmin=seq(5,15,by=0.5),
                                Tmax=seq(30,39,by=0.5)))
fitted_br1_13 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID13,
                      start = grid_br1_13,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_13)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_13 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID13,
                     start = list(a=0.00014,Tmin=5,Tmax=35.5))

sum_br1_13 <- summary(fitted_br1_13)
sum_br1_13 #Tmin no significativo
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
                            "a,Tmax")
colnames(params_br1_13) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_14)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_14 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID14,
                     start = list(a=0.00015,Tmin=12.5,Tmax=33))

sum_br1_14 <- summary(fitted_br1_14)
sum_br1_14 #Tmin no significativo, a no significativo
boot_br1_14 <- nlsBoot(fitted_br1_14, niter = 999) #sale mal...
coefs_14 <- data.frame(sum_br1_14$coefficients[,1:2],row.names = NULL)
boot_14 <- data.frame(boot_br1_14$estiboot)
params_br1_14 <- data.frame(coefs_14[1,],
                            NA,NA,
                            coefs_14[3,],
                            boot_14[1,],
                            boot_14[2,],
                            boot_14[3,])
colnames(params_br1_13) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot")
params_br1_13

#### _ _ Study 15 ####
IR_data_ID15 <- IR_data_all %>%
  filter(id==15)
par(mfrow=c(2,2))
plot(IR_data_ID15$temperature,IR_data_ID15$growth_rate)
plot(IR_data_ID15$temperature,briere1(a=0.00023,Tmin=13,Tmax=33,temp=IR_data_ID15$temperature))
plot(IR_data_ID15$temperature,briere1(a=0.00022,Tmin=13,Tmax=30,temp=IR_data_ID15$temperature))
plot(IR_data_ID15$temperature,briere1(a=0.00035,Tmin=13,Tmax=30,temp=IR_data_ID15$temperature))# <-
grid_br1_15 <- expand.grid(list(a=seq(0.00022,0.00037,by=0.00001),
                                Tmin=seq(8,18,by=0.5),
                                Tmax=seq(28,35,by=0.5)))
fitted_br1_15 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID15,
                      start = grid_br1_15,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_15)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_15 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID15,
                     start = list(a=0.00022,Tmin=12,Tmax=30.5))

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
                            "Tmax")
colnames(params_br1_15) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_16)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_16 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID16,
                     start = list(a=0.00023,Tmin=15,Tmax=35))

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
                            "all")
colnames(params_br1_16) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_16

#### _ _ Study 17 ####
IR_data_ID17 <- IR_data_all %>%
  filter(id==17)
par(mfrow=c(2,2))
plot(IR_data_ID17$temperature,IR_data_ID17$growth_rate)
plot(IR_data_ID17$temperature,briere1(a=0.0001,Tmin=8,Tmax=33,temp=IR_data_ID17$temperature))
plot(IR_data_ID17$temperature,briere1(a=0.0001,Tmin=15,Tmax=30,temp=IR_data_ID17$temperature))
plot(IR_data_ID17$temperature,briere1(a=0.00023,Tmin=17,Tmax=30,temp=IR_data_ID17$temperature))
grid_br1_17 <- expand.grid(list(a=seq(0.00010,0.00028,by=0.00001),
                                Tmin=seq(12,20,by=0.5),
                                Tmax=seq(28,35,by=0.5)))
fitted_br1_17 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID17,
                      start = grid_br1_17,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_17)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_17 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID17,
                     start = list(a=0.00028,Tmin=18.5,Tmax=30)) # NO SALE

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
                            "Tmax")
colnames(params_br1_17) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_18)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_18 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID18,
                     start = list(a=0.00004,Tmin=8.5,Tmax=30))#no hay forma

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
                            "Tmax")
colnames(params_br1_18) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_18

#### _ _ Study 19 ####
IR_data_ID19 <- IR_data_all %>%
  filter(id==19)
par(mfrow=c(2,2))
plot(IR_data_ID19$temperature,IR_data_ID19$growth_rate)
plot(IR_data_ID19$temperature,briere1(a=0.00015,Tmin=15,Tmax=35,temp=IR_data_ID19$temperature))
plot(IR_data_ID19$temperature,briere1(a=0.00017,Tmin=15,Tmax=35,temp=IR_data_ID19$temperature))
plot(IR_data_ID19$temperature,briere1(a=0.0002,Tmin=17,Tmax=35,temp=IR_data_ID19$temperature))# <-
grid_br1_19 <- expand.grid(list(a=seq(0.00012,0.00025,by=0.00001),
                                Tmin=seq(11,19,by=0.5),
                                Tmax=seq(30,38,by=0.5)))
fitted_br1_19 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID19,
                      start = grid_br1_19,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_19)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_19 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID19,
                     start = list(a=0.00025,Tmin=17.5,Tmax=33))

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
                            "all")
colnames(params_br1_19) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_20)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_20 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID20,
                     start = list(a=0.00008,Tmin=13.5,Tmax=31))#no hay forma

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
                            "Tmax")
colnames(params_br1_20) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_21)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_21 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID21,
                     start = list(a=1.100e-04,Tmin=1.600e+01,Tmax=3.550e+01))#no hay forma

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
                            "all")
colnames(params_br1_21) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_22)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_22 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID22,
                     start = list(a=0.00003,Tmin=6.5,Tmax=43)) #NO sale
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
                            "Tmax")
colnames(params_br1_22) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_23)
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
                            "Tmax")
colnames(params_br1_23) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_24)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_24 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID24,
                     start = list(a=0.00015,Tmin=3.5,Tmax=30)) #NO sale
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
                            "Tmax")
colnames(params_br1_24) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_25)
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
grid_br1_26 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(5,18,by=0.5),
                                Tmax=seq(30,40,by=0.5)))
fitted_br1_26 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID26,
                      start = grid_br1_26,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_26)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
#fitted_br1_26 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
 #                    data = IR_data_ID26,
  #                   start = list(a=0.0001,Tmin=13.5,Tmax=37.5))
sum_br1_26 <- summary(fitted_br1_26)
sum_br1_26 
#Sale mucho mejor la del grid, así que pasamos los del grid
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
                            "Tmax")
colnames(params_br1_26) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_26

#### _ _ Study 27 ####
IR_data_ID27 <- IR_data_all %>%
  filter(id==27)
par(mfrow=c(2,2))
plot(IR_data_ID27$temperature,IR_data_ID27$growth_rate)
plot(IR_data_ID27$temperature,briere1(a=0.0001,Tmin=20,Tmax=34,temp=IR_data_ID27$temperature))
plot(IR_data_ID27$temperature,briere1(a=0.0002,Tmin=22,Tmax=36,temp=IR_data_ID27$temperature))
plot(IR_data_ID27$temperature,briere1(a=0.0002,Tmin=21,Tmax=35,temp=IR_data_ID27$temperature))# <-
grid_br1_27 <- expand.grid(list(a=seq(0.0001,0.0003,by=0.00001),
                                Tmin=seq(15,25,by=0.5),
                                Tmax=seq(30,40,by=0.5)))
fitted_br1_27 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID27,
                      start = grid_br1_27,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_27)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_27 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID27,
                   start = list(a=0.00024,Tmin=23,Tmax=34.5),
                   nls.control(warnOnly=TRUE))
sum_br1_27 <- summary(fitted_br1_27)
sum_br1_27 
## NO SALE

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
summary(fitted_br1_28)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_28 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID28,
                     start = list(a=0.00008,Tmin=12.5,Tmax=35.5)) #NO sale
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
                            "Tmin,Tmax")
colnames(params_br1_28) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_28

#### _ _ Study 29 ####
IR_data_ID29 <- IR_data_all %>%
  filter(id==29)
par(mfrow=c(2,2))
plot(IR_data_ID29$temperature,IR_data_ID29$growth_rate)
plot(IR_data_ID29$temperature,briere1(a=0.0001,Tmin=15,Tmax=35,temp=IR_data_ID29$temperature))
plot(IR_data_ID29$temperature,briere1(a=0.00045,Tmin=15,Tmax=35,temp=IR_data_ID29$temperature))
plot(IR_data_ID29$temperature,briere1(a=0.00045,Tmin=16,Tmax=30.5,temp=IR_data_ID29$temperature))


grid_br1_29 <- expand.grid(list(a=seq(0.0001,0.0005,by=0.0001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(28,35,by=0.5)))
fitted_br1_29 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID29,
                      start = grid_br1_29,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_29)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_29 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID29,
                     start = list(a=0.0005,Tmin=17,Tmax=30.5)) #NO sale
# no sale

#### _ _ Study 30 #### 
excluido por solo 2 datos
IR_data_ID30 <- IR_data_all %>%
  filter(id==30)
par(mfrow=c(2,2))
plot(IR_data_ID30$temperature,IR_data_ID30$growth_rate)
plot(IR_data_ID30$temperature,briere1(a=0.0001,Tmin=15,Tmax=35,temp=IR_data_ID30$temperature))
plot(IR_data_ID30$temperature,briere1(a=0.0001,Tmin=14,Tmax=35,temp=IR_data_ID30$temperature))
grid_br1_30 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(10,20,by=0.5),
                                Tmax=seq(30,40,by=0.5)))
fitted_br1_28 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID28,
                      start = grid_br1_28,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_28)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_28 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID28,
                     start = list(a=0.00008,Tmin=12.5,Tmax=35.5)) #NO sale
sum_br1_28 <- summary(fitted_br1_28)
sum_br1_28 
boot_br1_28 <- nlsBoot(fitted_br1_28, niter = 999)
coefs_28 <- data.frame(sum_br1_28$coefficients[,1:2],row.names = NULL)
boot_28 <- data.frame(boot_br1_28$estiboot)
params_br1_28 <- data.frame(coefs_28[1,],
                            coefs_28[2,],
                            coefs_28[3,],
                            boot_28[1,],
                            boot_28[2,],
                            boot_28[3,],
                            "Tmin,Tmax")
colnames(params_br1_28) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_28

#### _ _ Study 31 ####
IR_data_ID31 <- IR_data_all %>%
  filter(id==31)
par(mfrow=c(2,2))
plot(IR_data_ID31$temperature,IR_data_ID31$growth_rate)
plot(IR_data_ID31$temperature,briere1(a=0.0001,Tmin=15,Tmax=35,temp=IR_data_ID31$temperature))
plot(IR_data_ID31$temperature,briere1(a=0.0002,Tmin=15,Tmax=32,temp=IR_data_ID31$temperature))
plot(IR_data_ID31$temperature,briere1(a=0.00023,Tmin=12,Tmax=32,temp=IR_data_ID31$temperature)) # < --
grid_br1_31 <- expand.grid(list(a=seq(0.00015,0.00025,by=0.00001),
                                Tmin=seq(5,15,by=0.5),
                                Tmax=seq(28,36,by=0.5)))
fitted_br1_31 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID31,
                      start = grid_br1_31,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_31)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_31 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID31,
                     start = list(a=0.00015,Tmin=6.5,Tmax=33)) #NO sale
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
                            "Tmax")
colnames(params_br1_31) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_deltha",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_31

#### _ _ Study 32 ####
IR_data_ID32 <- IR_data_all %>%
  filter(id==32)
par(mfrow=c(2,2))
plot(IR_data_ID32$temperature,IR_data_ID32$growth_rate)
plot(IR_data_ID32$temperature,briere1(a=0.0001,Tmin=15,Tmax=35,temp=IR_data_ID32$temperature))
plot(IR_data_ID32$temperature,briere1(a=0.00016,Tmin=14,Tmax=37,temp=IR_data_ID32$temperature))
plot(IR_data_ID32$temperature,briere1(a=0.00013,Tmin=15,Tmax=38,temp=IR_data_ID32$temperature)) # < --
grid_br1_32 <- expand.grid(list(a=seq(0.0001,0.00025,by=0.00001),
                                Tmin=seq(8,20,by=0.5),
                                Tmax=seq(30,42,by=0.5)))
fitted_br1_32 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID32,
                      start = grid_br1_32,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_32)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_32 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID32,
                     start = list(a=0.00013,Tmin=16,Tmax=39)) #NO sale
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
                            "all")
colnames(params_br1_32) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_deltha",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_33)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_33 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID33,
                     start = list(a=0.00016,Tmin=18,Tmax=25.5)) #NO sale
sum_br1_33 <- summary(fitted_br1_33)
sum_br1_33 
boot_br1_33 <- nlsBoot(fitted_br1_33, niter = 999)
coefs_33 <- data.frame(sum_br1_33$coefficients[,1:2],row.names = NULL)
boot_33 <- data.frame(boot_br1_33$estiboot)
# what about Topt?
## first we need to check if parameters are normally distributed:
hist(boot_br1_33$coefboot[,"a"]) #buena pinta
hist(boot_br1_33$coefboot[,"Tmin"]) #meeh
hist(boot_br1_33$coefboot[,"Tmax"]) #buena pinta
shapiro.test(boot_br1_33$coefboot[,"a"]) #no salen
shapiro.test(boot_br1_33$coefboot[,"Tmin"]) #no salen
shapiro.test(boot_br1_33$coefboot[,"Tmax"]) #no salen
# pongamos que asumimos normalidad...
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
                            "all")
colnames(params_br1_33) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")


#### _ _ Study 34 ####
IR_data_ID34 <- IR_data_all %>%
  filter(id==34)
par(mfrow=c(2,2))
plot(IR_data_ID34$temperature,IR_data_ID34$growth_rate)
plot(IR_data_ID34$temperature,briere1(a=0.0001,Tmin=19,Tmax=31,temp=IR_data_ID34$temperature))
plot(IR_data_ID33$temperature,briere1(a=0.0002,Tmin=19,Tmax=30.5,temp=IR_data_ID34$temperature))
plot(IR_data_ID33$temperature,briere1(a=0.0002,Tmin=19,Tmax=30,temp=IR_data_ID34$temperature)) # < --
grid_br1_34 <- expand.grid(list(a=seq(0.0001,0.00025,by=0.00001),
                                Tmin=seq(13,23,by=0.5),
                                Tmax=seq(25,35,by=0.5)))
fitted_br1_34 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID34,
                      start = grid_br1_34,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_34)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_34 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID34,
                     start = list(a=0.00025,Tmin=18.5,Tmax=30)) #NO sale
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
                            "Tmin,Tmax")
colnames(params_br1_34) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_35)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_35 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID35,
                     start = list(a=0.00006,Tmin=13.5,Tmax=37.5))
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
                            "Tmax")
colnames(params_br1_35) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_36)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_36 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID36,
                     start = list(a=0.00006,Tmin=5.5,Tmax=30))
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
                            "a,Tmax")
colnames(params_br1_36) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_37)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_37 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID37,
                     start = list(a=0.00006,Tmin=9,Tmax=46))
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
                            "none")
colnames(params_br1_37) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_38)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_38 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID38,
                     start = list(a=0.0001,Tmin=13,Tmax=25.5))
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
                            "all")
colnames(params_br1_38) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_38

#### _ _ Study 39 ####
IR_data_ID39 <- IR_data_all %>%
  filter(id==39)
par(mfrow=c(2,2))
plot(IR_data_ID39$temperature,IR_data_ID39$growth_rate)
plot(IR_data_ID39$temperature,briere1(a=0.0001,Tmin=8,Tmax=27,temp=IR_data_ID39$temperature))
plot(IR_data_ID39$temperature,briere1(a=0.0001,Tmin=8,Tmax=26,temp=IR_data_ID39$temperature))
plot(IR_data_ID39$temperature,briere1(a=0.0001,Tmin=8,Tmax=25.5,temp=IR_data_ID39$temperature))
grid_br1_39 <- expand.grid(list(a=seq(0.00005,0.00015,by=0.00001),
                                Tmin=seq(0,15,by=0.5),
                                Tmax=seq(20,30,by=0.5)))
fitted_br1_39 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                      data = IR_data_ID39,
                      start = grid_br1_39,
                      algorithm = "brute-force",
                      trace = TRUE)
summary(fitted_br1_39)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_39 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID39,
                     start = list(a=0.00005,Tmin=3,Tmax=28))
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
                            "a,Tmax")
colnames(params_br1_39) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
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
summary(fitted_br1_40)
# viendo esos parametros, buscamos starters cercanos para facilitar el ajuste
fitted_br1_40 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID40,
                     start = list(a=0.00017,Tmin=9,Tmax=33.5))
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
                            "all")
colnames(params_br1_40) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "Topt_est","Topt_se_delta",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot",
                             "which_signif")
params_br1_40
