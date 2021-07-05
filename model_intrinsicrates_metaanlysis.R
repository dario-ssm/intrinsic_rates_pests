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
x11() 
preview(formula=r~briere1(a,temp,Tmin,Tmax),
        data=IR_test,
        start=list(a=0.0001,Tmin=8,Tmax=40)

fitt_br1_pooled <- nls(r~briere1(a,temp,Tmin,Tmax),
                       data= IR_test,
                       start=list(a=0.0001,Tmin=8,Tmax=40))

summary(fitt_br1_pooled)
coef(fitt_br1_pooled) #no realista el -6.78ºC de LDT
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
Tmin = -6.77ºC
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
  

#### 9. individual studies regression with Bootstrap ####
## first: assign a number to each unique article as ID
IR_data_ID <- IR_data %>%
  distinct(`Article Title`)%>%
  mutate(id=row_number())
IR_data_all <- inner_join(IR_data,IR_data_ID,by="Article Title")
View(IR_data_all)
## one by one ##

#### _ _ Study 1 ####
IR_data_ID1 <- IR_data_all %>%
  filter(id==1)
x11()
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
                    algorithm = "brute-force")
sum_br1_1 <- summary(fitted_br1_1)
boot_br1_1 <- nlsBoot(fitted_br1_1, niter = 999)
coefs_1 <- data.frame(sum_br1_1$coefficients[,1:2],row.names = NULL)
boot_1 <- data.frame(boot_br1_1$estiboot)
params_br1_1 <- data.frame(coefs_1[1,],
                           coefs_1[2,],
                           coefs_1[3,],
                           boot_1[1,],
                           boot_1[2,],
                           boot_1[3,])
colnames(params_br1_1) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot")

#### _ _ Study 2 ####
IR_data_ID2 <- IR_data_all %>%
  filter(id==2)
x11()
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
params_br1_2 <- data.frame(coefs_2[1,],
                           coefs_2[2,],
                           coefs_2[3,],
                           boot_2[1,],
                           boot_2[2,],
                           boot_2[3,])
colnames(params_br1_2) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot")
#### _ _ Study 3 ####
IR_data_ID3 <- IR_data_all %>%
  filter(id==3)
x11()
par(mfrow=c(2,2))
plot(IR_data_ID1$temperature,IR_data_ID1$growth_rate)
plot(IR_data_ID1$temperature,briere1(a=0.0002,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))
plot(IR_data_ID1$temperature,briere1(a=0.0001,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))# <-
plot(IR_data_ID1$temperature,briere1(a=0.00009,Tmin=8,Tmax=35,temp=IR_data_ID1$temperature))
grid_br1_3 <- expand.grid(list(a=seq(0.00001,0.00015,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(28,50,by=0.5)))
fitted_br1_3 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID3,
                     start = list(a=0.00015,Tmin=6.5,Tmax=30))

fitted_br1_3 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID3,
                     start = grid_br1_3,
                     algorithm = "brute-force",
                     trace = TRUE)
sum_br1_3 <- summary(fitted_br1_3)
boot_br1_3 <- nlsBoot(fitted_br1_3, niter = 999)
coefs_3 <- data.frame(sum_br1_3$coefficients[,1:2],row.names = NULL)
boot_3 <- data.frame(boot_br1_3$estiboot)
params_br1_3 <- data.frame(NA,NA,
                           NA,NA,
                           coefs_3[3,],
                           boot_3[1,],
                           boot_3[2,],
                           boot_3[3,])
colnames(params_br1_3) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot")
#### _ _ Study 4 ####
IR_data_ID4 <- IR_data_all %>%
  filter(id==4)
x11()
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
grid_br1_4 <- expand.grid(list(a=seq(0.00001,0.00015,by=0.00001),
                               Tmin=seq(0,15,by=0.5),
                               Tmax=seq(25,45,by=0.5)))
fitted_br1_4 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID4,
                     start = grid_br1_4,
                     algorithm = "brute-force",
                     trace = FALSE)

sum_br1_4 <- summary(fitted_br1_4)
# NO SIGNIFICATIVO
#### _ _ Study 5 ####
IR_data_ID5 <- IR_data_all %>%
  filter(id==5)
x11()
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
sum_br1_5 <- summary(fitted_br1_5)
boot_br1_5 <- nlsBoot(fitted_br1_5, niter = 999)
coefs_5 <- data.frame(sum_br1_5$coefficients[,1:2],row.names = NULL)
boot_5 <- data.frame(boot_br1_5$estiboot)
params_br1_5 <- data.frame(coefs_5[1,],
                           coefs_5[2,],
                           coefs_5[3,],
                           boot_5[1,],
                           boot_5[2,],
                           boot_5[3,])
colnames(params_br1_5) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot")
#### _ _ Study 6 ####
IR_data_ID6 <- IR_data_all %>%
  filter(id==6)
x11()
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
boot_br1_6 <- nlsBoot(fitted_br1_6, niter = 999)
coefs_6 <- data.frame(sum_br1_6$coefficients[,1:2],row.names = NULL)
boot_6 <- data.frame(boot_br1_6$estiboot)
params_br1_6 <- data.frame(NA,NA,
                           coefs_6[2,],
                           coefs_6[3,],
                           boot_6[1,],
                           boot_6[2,],
                           boot_6[3,])
colnames(params_br1_6) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot")
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
params_br1_7 <- data.frame(NA,NA,
                           coefs_7[2,],
                           coefs_7[3,],
                           boot_7[1,],
                           boot_7[2,],
                           boot_7[3,])
colnames(params_br1_7) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot")
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
fitted_br1_8 <- nls(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                    data = IR_data_ID8,
                    start = list(a=0.00027,Tmin=5,Tmax=30))
fitted_br1_8 <- nls2(formula= growth_rate ~ briere1(a,temp = temperature,Tmin,Tmax),
                     data = IR_data_ID8,
                     start = grid_br1_8,
                     algorithm = "brute-force",
                     trace = TRUE)

sum_br1_8 <- summary(fitted_br1_8)
boot_br1_8 <- nlsBoot(fitted_br1_8, niter = 999)
coefs_8 <- data.frame(sum_br1_8$coefficients[,1:2],row.names = NULL)
boot_8 <- data.frame(boot_br1_8$estiboot)
params_br1_8 <- data.frame(NA,NA,
                           NA,NA,
                           coefs_8[3,],
                           boot_8[1,],
                           boot_8[2,],
                           boot_8[3,])
colnames(params_br1_8) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot")
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
params_br1_9 <- data.frame(coefs_9[1,],
                           coefs_9[2,],
                           coefs_9[3,],
                           boot_9[1,],
                           boot_9[2,],
                           boot_9[3,])
colnames(params_br1_9) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot")

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
params_br1_10 <- data.frame(NA,NA,
                           NA,NA,
                           coefs_10[3,],
                           boot_10[1,],
                           boot_10[2,],
                           boot_10[3,])
colnames(params_br1_10) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot")

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
params_br1_11 <- data.frame(coefs_11[1,],
                           coefs_11[2,],
                           coefs_11[3,],
                           boot_11[1,],
                           boot_11[2,],
                           boot_11[3,])
colnames(params_br1_11) <- c("a_est_nls","a_se_nls",
                            "Tmin_est_nls","Tmin_se_nls",
                            "Tmax_est_nls","Tmax_se_nls",
                            "a_est_boot","a_se_boot",
                            "Tmin_est_boot","Tmin_se_boot",
                            "Tmax_est_boot","Tmax_se_boot")
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
params_br1_12 <- data.frame(coefs_12[1,],
                            coefs_12[2,],
                            coefs_12[3,],
                            boot_12[1,],
                            boot_12[2,],
                            boot_12[3,])
colnames(params_br1_12) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot")
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
params_br1_13 <- data.frame(coefs_13[1,],
                            NA,NA,
                            coefs_13[3,],
                            boot_13[1,],
                            boot_13[2,],
                            boot_13[3,])
colnames(params_br1_13) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot")
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
boot_13 <- data.frame(boot_br1_13$estiboot)
params_br1_13 <- data.frame(coefs_13[1,],
                            NA,NA,
                            coefs_13[3,],
                            boot_13[1,],
                            boot_13[2,],
                            boot_13[3,])
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
params_br1_15 <- data.frame(NA,NA,
                            NA,NA,
                            coefs_15[3,],
                            boot_15[1,],
                            boot_15[2,],
                            boot_15[3,])
colnames(params_br1_15) <- c("a_est_nls","a_se_nls",
                             "Tmin_est_nls","Tmin_se_nls",
                             "Tmax_est_nls","Tmax_se_nls",
                             "a_est_boot","a_se_boot",
                             "Tmin_est_boot","Tmin_se_boot",
                             "Tmax_est_boot","Tmax_se_boot")
params_br1_15
