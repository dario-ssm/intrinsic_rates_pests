rm(list=ls())

## exploring filters

#load tidyverse
library(tidyverse)

## read document and filter by "yes" in column "Filter" 
pooled_yes <- read_csv("/home/dario/Documentos/salvar pen/UAH/Beca investigacion/Resultados_busqueda/busqueda_definitiva/pooled.csv",
                       col_names = TRUE) %>% 
  filter(Filter=="yes") %>%
  glimpse()


## plot by years
pooled_yes %>% ggplot(aes(Year))+
  geom_bar(aes(colour=Year))

## select samples of 50 to show titles in order to assure good review
pooled_subsample_titles <- pooled_yes %>% 
  slice(sample(1:1890, 95)) %>%
  select("Article Title")%>%
  view()

#### second filter ####
# After first filter, now we have different categories
filtered2 <- read_csv("/home/dario/Documentos/salvar pen/UAH/Beca investigacion/Resultados_busqueda/busqueda_definitiva/filtered_review.csv", 
                      col_names=TRUE) %>%
  filter(Filter == "yes")%>%
  glimpse()


#See counts of combinations of categories:
filtered2 %>%
  count(Approach, Subapproach,Lab_response)%>%
  view()
# See grouping categories:
filtered2 %>%
  count(Approach)%>%
  view()
filtered2 %>%
  count(Approach, Subapproach)%>%
  view()

# Extract  a csv with each filtered dataset
intrinsic_growth <- filtered2 %>%
  filter(Approach == "experimental" &
           Subapproach == "laboratory" &
           Lab_response == "rm" )%>%
  write_csv(path = "/home/dario/Documentos/salvar pen/UAH/Beca investigacion/Resultados_busqueda/busqueda_definitiva/intrinsic_growth.csv")

#### aggregate and exploratory of where, when and how for rm ####

intrinsic_classified <- read_csv("/home/dario/Documentos/salvar pen/UAH/Beca investigacion/Resultados_busqueda/busqueda_definitiva/intrinsic_growth.csv",
                                 col_names = TRUE) %>%
  filter(Filter_3 == "yes" & n==1) # pass filter
# 1s are for first row of each study)

#by order
intrinsic_classified %>%
  count(order)%>%
  view()
#by family
intrinsic_classified %>%
  count(family)%>%
  view()
#by feeding_guilds
intrinsic_classified %>%
  count(feeding_guild)%>%
  view()
#by host plant family
intrinsic_classified %>%
  count(h_p_family)%>%
  view()
#by diet class
intrinsic_classified %>%
  count(diet_family)%>%
  view()
#by both order and family
intrinsic_classified %>%
  count(order,family)%>%
  view()

#by both family and feeding guild
intrinsic_classified %>%
  count(family,feeding_guild)%>%
  view()
intrinsic_classified %>%
  count(order,feeding_guild)%>%
  view()

#### mapping where primary studies have been carried out ####
library(rworldmap)
library(osmdata)
library(sf)
library(ggmap)
#create map along coordinates
sbbox <- intrinsic_classified %>%
  make_bbox(lon = lon, lat =lat, f = .1)

sq_map<- get_stamenmap(sbbox, maptype = "terrain-background", zoom = 2) 
sq_map<-get_map(location = sbbox, maptype = "toner-background", source = "google")

# plot map with ggplot2
ggmap(sq_map) +
  geom_point(data = intrinsic_classified, #base de datos
             mapping = aes(x = lon, y = lat),#coordenadas x e y del mapa
             size=2,#tamaÃ±o de los puntos (mÃ¡s pequeÃ±o el nÂº, mÃ¡s pequeÃ±o el punto)
             color = "darkred", #color de los puntos, ver ColorChart 
             alpha=0.3)#transparencia (entre 0 y 1)
# alternative using Rbase
newmap <- getMap(resolution="high") ## obteniendo mapa mundial a alta resoluciÃ³n
plot(newmap, col="lightgrey", border="white")
points(intrinsic_classified$lon,intrinsic_classified$lat,
       pch=19,cex=1,col=adjustcolor("darkred",0.27))
## gganimate (not working well)
library(maps)
library(ggplot2)
map <- ggplot(data = intrinsic_classified, aes(x = lon, y = lat)) +
  borders("world", colour = "gray90", fill = "gray90") +
  geom_point(aes(colour = order), alpha = 0.35, size = 3)+theme_classic()+
  theme(legend.position = "bottom")


library(gganimate)
map_states <- map +
  transition_states(states = Year,   # variable for movement
                    transition_length = 1, # relative time of transition
                    state_length = 2,      # relative time of pause
                    wrap = TRUE)   

#### FIRSTS EXAMINATIONS AND VARIABLE ANALYSIS ####
#load fixed dataset
## with r_m values paired with temperature values.
## third_filtered by whether r_m are extractable or not

library(readr)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(svglite)
IR_data <-read_csv("/Users/Ecologia/Documents/USUARIOS/DARÍO/intrinsic_growth2.csv") %>%
  filter(Filter_3 == "yes") %>% #select only those which have passed the filter 3 in the source csv
  mutate(growth_rate = replace(growth_rate, growth_rate <= 0,0))%>%
  filter(as.numeric(temperature)<50)%>%
  glimpse() # str()
#since some numerical variables are read as characters, we change them to dbl
IR_data$temperature <- as.numeric(IR_data$temperature)
IR_data$growth_rate <- as.numeric(IR_data$growth_rate)
IR_data$error <- as.numeric(IR_data$n_1)
IR_data$RH <- as.numeric(IR_data$RH)
IR_data$lat <- as.numeric(IR_data$lat)

IR_data %>% 
  count(`Article Title`)%>%
  View()
#### a) scatterplots grouped para ver relaciones ####

scatter_IR_grouped <- ggplot(data=IR_data, aes(x=temperature, y=growth_rate))+
  geom_point(alpha=0.2,aes(color=order))+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette="Dark2")
ggsave(filename = "scatter_rm_all.png",units = "cm" ,height = 20, width = 15)

#### b) linearity,normality,regressions ####
#histograms
hist(IR_data$temperature)
hist(IR_data$growth_rate) #poisson with lambda = 1 ?
shapiro.test(IR_data$temperature) #temperature not-normal? by Saphiro test
shapiro.test(IR_data$growth_rate) #r not-normal? by Saphiro test
ks.test(IR_data$temperature,"pnorm")#temperature not-normal? by Kolmogorov-Smirnov
ks.test(IR_data$growth_rate,"pnorm") #r not normal by Kolmogorov-Smirnov

#skewness and kurtosis
library(moments)
skewness(IR_data$temperature,na.rm=TRUE) # just slightly left-skewed
skewness (IR_data$growth_rate,na.rm=TRUE) #highly right-skewed
kurtosis(IR_data$temperature,na.rm=TRUE) #a bit flattened
kurtosis(IR_data$growth_rate,na.rm=TRUE) #strongly peaked

#linear regression requirements
par(mfrow=c(2,2))
lm1 <- lm(growth_rate ~ temperature, data=IR_data)
summary(lm1)
plot(lm1) #not normal (deviation at edges), not sure if homocedastic(more variance in the centre)
#other function
library(performance)
library(qqplotr)
check_model(lm1) 

# log-transform
log_rm_data <- data.frame(log(IR_data$growth_rate),IR_data$temperature)
colnames(log_rm_data) <- c("rm","temp")
hist(log_rm_data$rm)
shapiro.test(as.numeric(log_rm_data$rm))
ks.test(log_rm_data$rm,"pnorm") 
#subset: checking covariables
ir_subdata <- IR_data %>%
  select(temperature,growth_rate,RH,lat,daylength)
pairs(ir_subdata) #the only one that might occur is growth_rate ~ lat

#let's examine lat
lat <- abs(IR_data$lat)
temp <- IR_data$temperature
rm <- IR_data$growth_rate
hist(lat) #poco normales
skewness(lat,na.rm=TRUE)#left-skewed
kurtosis(lat,na.rm=TRUE)#well
ks.test(lat,"pnorm") #not-normal

covars <- tibble(temp,rm,lat) #build a tibble to facilitate correlation analysis
pairs(covars) #lat ~ temp has no interest
GGally::ggpairs()
llm2 <- lm(rm~lat,data=covars)
summary(lm2) #a correlation seems reasonable (significant)
check_model(lm2) #at extremes: nonlinear, non-homocedastic, ab-normal

#homocedasticity
fligner.test(growth_rate~temperature,data = IR_data) #no homocedástico
#### EXPLORING GROUPING VARIABLES ####
#### a) Order ####
scatter_IR_grouped_se <- ggplot(data=IR_data, aes(x=temperature, y=growth_rate,color=order,fill = order))+
  geom_point(alpha=0.4)+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette="Dark2")
ggsave(filename = "scatter_rm_all.png",units = "cm" ,height = 20, width = 15)

# now grided
scatter_IR_grided <- ggplot(data=IR_data, aes(x=temperature, y=growth_rate,color=order,fill = order))+
  geom_point(alpha=0.2)+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  facet_wrap(.~order,nrow = 2)+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave(filename = "bonf_order.svg",units = "cm" ,
       height = 20, width = 15)

## a.1) See what happens into different groups
# All acari as one
acari <- IR_data %>%
  filter(order == "Acari>Prostigmata" | 
         order == "Acari>Trombidiformes")%>%
  mutate(order = "Acari")%>%
  glimpse()

acari_curve <- ggplot(acari, aes(x=temperature,y=growth_rate))+
  geom_point(color="lightcoral",alpha=0.5)+
  labs(title = "Intrinsic rate of increase values for Acari",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth(color="lightcoral",fill="lightcoral")+
  theme_bw()
ggsave(filename = "acari_all.png",units = "cm" ,height = 20, width = 15)

#acari by family:
acari_curve_fam <- ggplot(acari, aes(x=temperature,y=growth_rate,color=family,fill=family))+
  geom_point(alpha=0.5)+
  labs(title = "Intrinsic rate of increase values for Acari",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()
ggsave(filename = "acari_by_family.png",units = "cm" ,height = 20, width = 15)

## hemiptera by family
hemiptera <- IR_data %>%
  filter(order == "Hemiptera")%>%
  glimpse()

hemiptera_curve_fam <- ggplot(hemiptera, aes(x=temperature,y=growth_rate,color=family,fill=family))+
  geom_point(alpha=0.5)+
  labs(title = "Intrinsic rate of increase values for Hemiptera",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()
ggsave(filename = "hemiptera_family_all.png",units = "cm" ,height = 20, width = 15)

hemiptera_curve_grid<- ggplot(hemiptera, aes(x=temperature,y=growth_rate,color=family,fill=family))+
  geom_point(alpha=0.25)+
  labs(title = "Intrinsic rate of increase values for Hemiptera",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  facet_wrap(.~family,nrow=2)+
  theme_bw()+
  theme(legend.position = "none")
ggsave(filename = "hemiptera_family_grided.png",units = "cm" ,height = 20, width = 15)

#como siguen siendo muy variables para Aphididae, probemos a dividir por especies
aphids <- hemiptera %>%
  filter(family == "Aphididae")
#se parte en dos el código para que funcione
 aphids <- aphids %>% mutate(Species_complete= paste(aphids$genus,aphids$species))%>%
  glimpse()
  
aphids_curve_spp <- ggplot(aphids, aes(x=temperature,y=growth_rate,color=Species_complete,fill=Species_complete))+
  geom_point(alpha=0.5)+
  labs(title = "Intrinsic rate of increase values for aphids",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  facet_wrap(.~Species_complete)+
  geom_smooth()+
  theme_bw()+
  theme(legend.position="none")
ggsave(filename = "aphids_by_species.png",units = "cm" ,height = 20, width = 15)

#vamos a ver por géneros de lepidópteras
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
ggsave(filename = "lepidoptera_by_fam.png",units = "cm" ,height = 20, width = 15)

# vamos a ver las Noctuidae
noctuids <- lepidoptera %>%
  filter(family== "Noctuidae")
noctuids <- noctuids %>%
  mutate(Species_complete= paste(noctuids$genus,noctuids$species))%>%
  glimpse()

noctuids_curve_spp <- ggplot(noctuids, aes(x=temperature,y=growth_rate,color=Species_complete,fill=Species_complete))+
  geom_point(alpha=0.5)+
  labs(title = "Intrinsic rate of increase values for noctuids",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()
#meh, muy poco


#### b) Feeding guilds ####
rm_feeding <- IR_data %>%
  select(temperature,growth_rate,lat,feeding_guild,order)

scatter_rmfeed_grided_by_order <- ggplot(data=rm_feeding, aes(x=temperature, y=growth_rate))+
  geom_point(aes(color=feeding_guild),alpha=0.4)+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  facet_wrap(.~order,nrow = 2)+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave(filename = "scatter_feeding_ord2.png",units = "cm" ,height = 20, width = 15)

scatter_rmfeed_grided_by_feedinguild <- ggplot(data=rm_feeding, aes(x=temperature, y=growth_rate))+
  geom_point(aes(color=feeding_guild),alpha=0.4)+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth(aes(color=feeding_guild))+
  facet_wrap(.~feeding_guild)+
  theme_bw()
ggsave(filename = "scatter_feeding2.png",units = "cm" ,height = 20, width = 15)

scatter_rmorder_grided_by_feedinguild <- ggplot(data=rm_feeding, aes(x=temperature, y=growth_rate))+
  geom_point(aes(color=order),alpha=0.4)+
  labs(subtitle = "Intrinsic rate of increase values for the entire dataset",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth(aes(color=order,fill=order))+
  facet_wrap(.~feeding_guild)+
  theme_bw()
ggsave(filename = "scatter_feeding_order.png",units = "cm" ,height = 20, width = 15)

#### c) latitude ####
scatter_IR_grouped_lat <- ggplot(data=IR_data, aes(x=temperature, y=growth_rate))+
  geom_point(alpha=0.5,size=3,aes(color=abs(lat)),position = position_jitter(w = 2, h = 0))+
     labs(subtitle = "Intrinsic rate of increase values and latitude",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()+
  scale_colour_gradient(low = "lightcoral",high = "cyan4")
ggsave(filename = "scatter_latitude.png",units = "cm" ,height = 20, width = 15)

#¿efectos de latitud según orden?
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
ggsave(filename = "scatter_latitude_gridord",units = "cm" ,height = 20, width = 15)

# vamos a discretizar la variable para tirar las smooths
IR_data_disLat <- IR_data %>%
  mutate(lat=abs(lat)) %>%
  mutate(lat= cut_number(lat,n=4,
                         labels=c("Tropical","warm","mild","cold temperate")))
  #drop_na(lat)
write_csv(IR_data_disLat,file="intrinsic_dislat.csv")
# grid por orden, smooths para categorías de latitud
scatter_grid_disLat <- ggplot(data=IR_data_disLat, aes(x=temperature, y=growth_rate,color=lat,fill=lat))+
  geom_point(alpha=0.3,size=3,aes(color=lat),position = position_jitter(w = 2, h = 0))+
  labs(subtitle = "Intrinsic rate of increase values and latitude",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()+
  facet_wrap(.~order,nrow=2)+
  theme(legend.position = "bottom")
ggsave(filename = "smooths_lat_gridord.png",units = "cm" ,height = 20, width = 15)


scatter_grid_disLat2 <- ggplot(data=IR_data_disLat,aes(x=temperature, y=growth_rate))+
  geom_point(alpha=0.3,size=3,aes(color=order),position=position_jitter(w=2,h=0))+
  labs(subtitle = "Intrinsic rate of increase values and latitude",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()+
  facet_wrap(.~lat)
ggsave(filename = "smooths_lat_gridbylat2.png",units = "cm" ,height = 20, width = 15)

#parece que hay tendencia, vamos a verlo para hemiptera y lepidoptera
lepi_dislat <- IR_data_disLat %>%
  filter(order=="Lepidoptera")

lepidoptera_rm_to_lat<- ggplot(lepi_dislat, aes(x=temperature,y=growth_rate,color=lat,fill=lat))+
  geom_point(alpha=0.5,size=3,position=position_jitter(w=2,h=0))+
  labs(title = "Intrinsic rate of increase values for Lepidoptera",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()

hemi_dislat <- IR_data_disLat %>%
  filter(order=="Hemiptera")
hemiptera_rm_to_lat<- ggplot(hemi_dislat, aes(x=temperature,y=growth_rate,color=lat,fill=lat))+
  geom_point(alpha=0.5,size=3,position=position_jitter(w=2,h=0))+
  labs(title = "Intrinsic rate of increase values for Hemiptera",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  facet_wrap(.~lat)+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave(filename = "hemiptera_lat2.png",units = "cm" ,height = 20, width = 15)

#en aphids
aphids_dislat <- IR_data_disLat %>%
  filter(family=="Aphididae")

aphids_rm_to_lat<- ggplot(aphids_dislat, aes(x=temperature,y=growth_rate,color=lat,fill=lat))+
  geom_point(alpha=0.5,size=3,position=position_jitter(w=2,h=0))+
  labs(title = "Intrinsic rate of increase values for Aphids",
       x= "Temperature (ºC)",
       y= "intrinsic rate of increase (r)")+
  geom_smooth()+
  theme_bw()

#### REGRESSIONS ####
# since data are not normally distributed neither for growth rate 
# (Poisson with lambda =1, probably) nor temperature (predictor),
# we'll apply generalized linear models to explore correlations 
# between variables that we have scanned in those scatterplots
#### a) GLMs ####
glm1 <- glm(growth_rate~temperature,
            data=IR_data,
            family = "poisson") #significativo

summary(glm1)
check_model(glm1)
par(mfrow=c(2,2))
plot(glm1)
E1 <- resid(glm1, type="pearson")
N <- nrow(IR_data)
p <- length(coef(glm1))
OvDisp_glm1 <- sum(E1^2)/(N-p)
# muy pequeño (no hay sobredispersión)
# AER::dispersiontest()
# AER: sim <- simulateResiduals(fmnb,refit=T,n=99)
#deviance explained 
D_2_glm1 <-(59.633-40.151)/59.633 # divianza explicada (null - resid)/null
# D^2=0.326698, not bad?
# para aplicar un equivalente a ANOVA, ver condiciones

#media
mean(IR_data$growth_rate,na.rm = TRUE)

#AIC
AIC(glm1)
# da un AIC infinito...¿no es la mejor distribución?
#### b) LMs transformations ####
# probamos a convertir en normal con sqrt()
lm1.2 <- lm(sqrt(growth_rate)~temperature, data=IR_data)#modelo
summary(lm1.2)#ver
AIC(lm1.2)
BIC(lm1.2)
#añadimos orden
lm1.3<- lm(sqrt(growth_rate)~temperature+order,data=IR_data)
summary(lm1.3)
anova(lm1.3)
AIC(lm1.3) #tiene un AIC más alto que el del 1.2
BIC(lm1.3)

lm1.4 <- lm(sqrt(growth_rate)~temperature+order*feeding_guild,data=IR_data)
summary(lm1.4)
AIC(lm1.4)

lm1.5 <- lm(sqrt(growth_rate)~feeding_guild, data=IR_data)
summary(lm1.5)
AIC(lm1.5)
anova(glm1.5)
bonf_feed <- pairwise.t.test(IR_data$growth_rate,IR_data$feeding_guild,data=IR_data, 
                              p.adjust.method = "bonferroni")
print(bonf_feed)

## vamos a transformar con sqrt() las dos variables y ver linealidad, normalidad y homocedasticidad
lm_transf<-lm(sqrt(growth_rate)~sqrt(temperature),IR_data)
performance::check_model(lm_transf)
lm_transf2<-lm(sqrt(growth_rate)~temperature,IR_data)
performance::check_model(lm_transf2)
cbrt <- function(x){
  y= x^(1/3)
y}
df <- data.frame(cbrt(IR_data$growth_rate),IR_data$temperature)
colnames(df) <- c("rm_cbrt","temp")
lm_transf3<-lm(rm_cbrt~temp,data=df)
performance::check_model(lm_transf3)

df2 <- data.frame(df$rm_cbrt,cbrt(df$temp))
colnames(df2) <- c("rm_cbrt","temp_cbrt")
lm_transf4 <- lm(rm_cbrt~temp_cbrt,data=df2)
performance::check_model(lm_transf4)

#ANOVAs
lm_order <- lm(sqrt(growth_rate)~order,data=IR_data)
anova(lm_order)
summary(lm_order)
bonf_order <- pairwise.t.test(IR_data$growth_rate,IR_data$order,data=IR_data, 
                              p.adjust.method = "bonferroni")
bonf_plot_order <- ggplot(IR_data, aes(y=growth_rate))+
  geom_boxplot(aes(fill=order))+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(title= "Variación de r en distintos taxones",
       x= "Orden", y= "Tasa de crecimiento intrínseco (r)")
ggsave("bonf_plot_order2.svg",units = "cm" ,height = 18, width = 20)

lm_feed <- lm(sqrt(growth_rate)~feeding_guild,data=IR_data)
anova(lm_feed)
summary(lm_feed)
bonf_feed <- pairwise.t.test(IR_data$growth_rate,IR_data$feeding_guild,data=IR_data, 
                              p.adjust.method = "bonferroni")
bonf_plot_feed <- ggplot(IR_data, aes(y=growth_rate))+
  geom_boxplot(aes(fill=feeding_guild))+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(title= "Variación de r en distintos feeding guilds",
       x= "Orden", y= "Tasa de crecimiento intrínseco (r)")
ggsave("bonf_plot_feed.svg",units = "cm" ,height = 18, width = 20)

####_ _ b.1. compare slopes####
# 1) Entre ordenes
#lepidoptera
lm_lepi <- lm(sqrt(growth_rate)~temperature, lepidoptera)
summary(lm_lepi)
slope_lepi <- (lm_lepi$coefficients[1])^(2)
#hemiptera
lm_hemip <- lm(sqrt(growth_rate)~temperature, hemiptera)
summary(lm_hemip)
slope_hemip <- (lm_hemip$coefficients[1])^(2)
#acari
lm_acari <- lm(sqrt(growth_rate)~temperature, acari)
summary(lm_acari)
slope_acari <- (lm_acari$coefficients[1])^(2)
slopes_order <- data.frame(slope_hemip,slope_acari,slope_lepi)

## por feeding_guild
#borers
borer <- IR_data%>%
  filter(feeding_guild=="borer")
lm_borer <- lm(sqrt(growth_rate)~temperature, borer)
summary(lm_borer)
slope_borer <- (lm_borer$coefficients[1])^(2)
#chewer
chewer <- IR_data%>%
  filter(feeding_guild=="chewer")
lm_chewer <- lm(sqrt(growth_rate)~temperature, chewer)
summary(lm_chewer)
slope_chewer <- (lm_chewer$coefficients[1])^(2)
#miner
miner <- IR_data%>%
  filter(feeding_guild=="miner"|feeding_guild=="leafminer")
lm_miner <- lm(sqrt(growth_rate)~temperature, miner)
summary(lm_miner)
slope_miner <- (lm_miner$coefficients[1])^(2)
#sucker
sucker <- IR_data%>%
  filter(feeding_guild=="sucker")
lm_sucker <- lm(sqrt(growth_rate)~temperature, sucker)
summary(lm_sucker)
slope_sucker <- (lm_sucker$coefficients[1])^(2)
slopes_feed <- data.frame(slope_sucker,slope_chewer,slope_borer,
                          slope_miner)

#lat?
lm_lat <- lm(sqrt(growth_rate)~abs(lat),IR_data)
summary(lm_lat)
#assuming poisson while mean < 5, and assuming less than 3
# random effects: Laplace or GHQ
# functions glmer or glmmML, LR test (random) or Wald Z (fixed)


## incluir la variable orden
glm2 <-  glm(growth_rate~temperature+order,
             data=IR_data,
             family = "poisson") #significativo
summary(glm2)

AIC(glm1)
#### c) GLMM ####
####  c.1 glmmML #### 
library(glmmML)
#### c.2 glmer ####
# para hacer un mixed effects incluyendo variable aleatoria
library(lme4)
#primero crear variable identificadora de artículo
##NO SOY CAPAZ
IR_data <- IR_data %>%
  rename(title = `Article Title`) #cambiar el título
  
glmm1 <- glmer(growth_rate~temperature+(1|title), data=IR_data,
      poisson)


print(glmm1) #??????
#se supone que deberíamos crear un modelo para growth~temp y otro para growth~temp+random y comparar con chi-cuadrado
glm2 <- glm(growth_rate~Authors,data=IR_data,
            family="poisson")



glm2 <- glm(growth_rate~lat,
            data=IR_data,
            family="poisson")
summary(glm2)
E1 <- resid(glm2, type="pearson")
N <- nrow(IR_data)
p <- length(coef(glm2))
OvDisp_glm2 <- sum(E1^2)/(N-p)
# muy pequeño (no hay sobredispersión)
#deviance explained 
D_2_glm2 <-(54.418-52.800)/54.418 # divianza explicada (null - resid)/null
# D^2=0.0297, poco poder explicativo
#¿los AIC son infinitamente altos?

#ver ahora metiendo familias etc
glm3 <- glm(growth_rate~temperature+order,
            data=IR_data,
            family = "poisson")

summary(glm3)

#### ELIMINAR OUTLIERS ####
# Se eliminan los outliers de áfidos
IR_data <- IR_data%>%
filter(`UT (Unique WOS ID)` != "WOS:000071638600063"&
         `UT (Unique WOS ID)`!= "WOS:A1995RM24300002" &
         growth_rate <= 0.6)
#Repetir análisis anteriores con el dataset sobreescrito
#### CONCLUSIONS ####
# 1) we conclude that it might be useful to include lat as a covariable.
# 2) data are not NORMAL, so the regression should follow generalized methods
# 3) data are not LINEAR, so we will have to use a non-linear approach (ex. Brière-1)
# 4) data have random variables such as Study (var. "Article Title"). Latitude as a covariate?
# 5) MODEL: generalized non-linear mixed-effects (hierarchichal) regression (GNLMER).

