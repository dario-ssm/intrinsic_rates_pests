#### Script for thermal suitability simulation ####
# Authors: Dario San Segundo Molina & Ignacio Morales Castilla
# Date: 06/06/2022

library(raster)
library(stringr)
library(ggplot2)
library(dplyr)
library(lubridate)
library(rworldmap)
library(purrr)
library(tidyr)
library(readr)
library(cowplot)
library(performance)
library(envirem)
library(abind)
library(rasterVis)
library(scales)
# 1. generic example -------------------------------------------------------
set.seed(110722)
daily_min <- c(rnorm(30,5,3), rnorm(30,10,3),rnorm(30,15,3))
daily_max <- c(rnorm(30,15,3), rnorm(30,25,3),rnorm(30,35,3))
daily_min2 <- c(rnorm(30,7,3), rnorm(30,11,3),rnorm(30,16,3))
daily_max2 <- c(rnorm(30,14,4), rnorm(30,24,5),rnorm(30,33,4))
tbhigh <- 21
tblow <- 18

sim_day_temps <- function(tmin, tmax, timestep){
  temp <- ((tmax-tmin)/2) * cos((pi * (timestep - 0.5)) / 48) + ((tmax + tmin)/2)
}
sum_a <- c(NULL)
a <- for(i in 1:48){
  aaa <- sim_day_temps(daily_min, daily_max, i)
  sum_a <- c(sum_a, aaa)
}

tmin <- c(daily_min, daily_min2)
tmin_rep <- rep(tmin, each = 2)
tmax <- c(daily_max, daily_max2)
tmax_rep <- rep(tmax, each = 2)
day <- rep(seq(1,90),2)
cell <- rep(1:2, each = 90)

season_length <- 90
tmin_desph <- c(NULL)
for(pixel in 1:2){
  tmin_cell <- tmin_rep[(season_length*(pixel-1)*2+1):(season_length*pixel*2)]
  index_length <- length(tmin_cell)
    for(i in 1:index_length){
    tmin_cell[i] <- if_else(i == index_length,
                            true = tmin_cell[i],
                            if_else(condition = i %% 2 == 0,
                                    true = tmin_cell[i+1],
                                    false = tmin_cell[i]))
  } 
  tmin_desph <- c(tmin_desph, tmin_cell)
}

n_cells = 2
ph_temps <- tibble(tmin = tmin_desph,
                   tmax = tmax_rep,
                   day = rep(day, each = n_cells),
                   cell = rep(c(1,2), each = 180),
                   mark = rep(c("left","right"), 180))

suitability <- ph_temps %>%
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, n_cells*season_length*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh &
                                                             sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
  

# 2. applied example for Spain -------------------------------------------------------
## we use linear mixed effects results of thermal breadth
## from meta-analysis (thermal_traits_analysis_sideB.R)
tblow <- 20.46
tbhigh <- 34.47

## climatic data
tmax_brick <- brick(x = "/Users/dario-ssm/Downloads/Spain02_v5.0_010reg_aa3d/Spain02_v5.0_010reg_aa3d/Spain02_v5.0_DD_010reg_aa3d_tasmax.nc")
tmin_brick <- brick(x = "/Users/dario-ssm/Downloads/Spain02_v5.0_010reg_aa3d/Spain02_v5.0_010reg_aa3d/Spain02_v5.0_DD_010reg_aa3d_tasmin.nc")
daily_tmin <- as_tibble(tmin_brick[,]) %>%
  tibble::rownames_to_column(var="cell")%>%
  mutate(cell=as.numeric(cell)) %>%
  pivot_longer(cols = starts_with("X"),
               names_to = "date",
               values_to = "tmin") %>% 
  mutate(date = str_sub(date, 2, 11),
         date = ymd(date),
         year = year(date))

daily_tmax<- as_tibble(tmax_brick[,]) %>%
  tibble::rownames_to_column(var="cell")%>%
  mutate(cell=as.numeric(cell)) %>%
  pivot_longer(cols = starts_with("X"),
               names_to = "date",
               values_to = "tmax") %>% 
  mutate(date = str_sub(date, 2, 11),
         date = ymd(date),
         year = year(date))
daily_temps <- inner_join(daily_tmin, daily_tmax)

#take just one year (2015)
daily_temps_2015_rep <- daily_temps %>% 
  filter(year == 2015) %>% 
  slice(rep(1:n(), each=2))
  
n_cells = 10902

# compute separately
### take a template for 1 day raster:
template_one_year <-  subset(tmax_brick,
                            which(str_sub(names(tmax_brick),2,5) ==  2015))
template_last_doy <-  template_one_year[[365]]
  
### prepare lagged data set
season_length <- 365
tmin_lagged <- c(NULL)

for(pixel in 1:n_cells){
  tmin_cell <- daily_temps_2015_rep$tmin[(season_length*(pixel-1)*2+1):(season_length*pixel*2)]
  index_length <- length(tmin_cell)
  for(i in 1:index_length){
    tmin_cell[i] <- if_else(i == index_length,
                            true = tmin_cell[i],
                            if_else(condition = i %% 2 == 0,
                                    true = tmin_cell[i+1],
                                    false = tmin_cell[i]))
  } 
  tmin_lagged <- c(tmin_lagged, tmin_cell)
  print(paste("cell", pixel, "/", n_cells))
}

temps_lagged <- tibble(tmin = tmin_lagged,
                       tmax = daily_temps_2015_rep$tmax,
                       day = day(daily_temps_2015_rep$date),
                       cell = daily_temps_2015_rep$cell,
                       mark = rep(c("left","right"), 365*10902))

## Eleven partitions are made based on cell index to avoid allocating huge vectors
### group 1
suitability_1 <- temps_lagged %>%
  filter(cell %in% c(1:1000)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*1000*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
### group 2
suitability_2 <- temps_lagged %>%
  filter(cell %in% c(1001:2000)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*1000*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
### group 3
suitability_3 <- temps_lagged %>%
  filter(cell %in% c(2001:3000)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*1000*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
### group 4
suitability_4 <- temps_lagged %>%
  filter(cell %in% c(3001:4000)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*1000*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
### group 5
suitability_5 <- temps_lagged %>%
  filter(cell %in% c(4001:5000)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*1000*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
### group 6
suitability_6 <- temps_lagged %>%
  filter(cell %in% c(5001:6000)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*1000*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
### group 7
suitability_7 <- temps_lagged %>%
  filter(cell %in% c(6001:7000)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*1000*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
### group 8
suitability_8 <- temps_lagged %>%
  filter(cell %in% c(7001:8000)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*1000*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
### group 9
suitability_9 <- temps_lagged %>%
  filter(cell %in% c(8001:9000)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*1000*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
### group 10
suitability_10 <- temps_lagged %>%
  filter(cell %in% c(9001:10000)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*1000*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))
### group 11
suitability_11 <- temps_lagged %>%
  filter(cell %in% c(10001:10902)) %>% 
  slice(rep(1:n(), each=48)) %>%
  mutate(minutes = rep(1:48, 365*902*2)) %>% 
  group_by(day, cell, mark) %>% 
  mutate(sim = sim_day_temps(tmin = tmin,
                             tmax = tmax,
                             timestep = minutes)) %>% 
  group_by(cell, day) %>% 
  summarise(daily_suitable = sum(series_logical <- if_else(sim <= tbhigh & sim >= tblow,
                                                           1,
                                                           0))/4) %>% 
  group_by(cell) %>% 
  summarise(suitability = sum(daily_suitable, na.rm = TRUE))

### join all subgroups:
suitability_pooled <- bind_rows(suitability_1,
                                suitability_2,
                                suitability_3,
                                suitability_4,
                                suitability_5,
                                suitability_6,
                                suitability_7,
                                suitability_8,
                                suitability_9,
                                suitability_10,
                                suitability_11) %>% 
  pivot_wider(names_from = cell, #return into raster-like matrix format
              values_from = suitability) %>% 
  t() %>% 
  as_tibble()
dim(template_last_doy)
dim(suitability_pooled)
values(template_last_doy) <- as.matrix(suitability_pooled)
map_suitability <- template_last_doy
plot(map_suitability)
writeRaster(map_suitability,
            "map_suitability.tif", 
            format = "GTiff",
            overwrite = TRUE)
template_last_doy