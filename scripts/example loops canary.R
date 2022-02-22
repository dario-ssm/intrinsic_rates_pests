library(dplyr)
islas <- c("La Palma","Tenerife","El Hierro","La Gomera","Gran Canaria", "Lanzarote", "Fuerteventura")
fincas_tabaibal <- c(10,45,4,6,32,12,15)
mean_masa_euphorbias <- c(20,30,26,34,45,41,23) 
sd_masa <- c(6,2,9,5,3,5,4)
euphorbia_canariensis <- data.frame(islas,
                                numero = fincas_tabaibal,
                                mean_masa_euphorbias,
                                sd_masa) %>% 
  glimpse()

#para cada isla, genera una simulación rnorm con n= número de fincas, mean = masa_euphorbia, sd = sd_masa
islas_new <- data.frame(isla = rep(islas, fincas_tabaibal),
                    peso = 0)
for (i in islas){
  n_isla <- euphorbia_canariensis[euphorbia_canariensis$islas ==i,"numero"]
  mean_isla <- euphorbia_canariensis[euphorbia_canariensis$islas ==i,"mean_masa_euphorbias"]
  sd_isla <- euphorbia_canariensis[euphorbia_canariensis$islas ==i,"sd_masa"]
  simulated_isla <- rnorm(n_isla,mean_isla,sd_isla)
  islas_new[islas_new$isla == i,2] <- simulated_isla
  }
  