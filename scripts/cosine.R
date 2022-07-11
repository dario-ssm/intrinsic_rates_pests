#### Script for thermal suitability simulation ####
# Authors: Dario San Segundo Molina & Ignacio Morales Castilla
# Date: 06/06/2022


# 1. generic example -------------------------------------------------------
set.seed(110722)
daily_min <- c(rnorm(30,5,3), rnorm(30,10,3),rnorm(30,15,3))
daily_max <- c(rnorm(30,15,3), rnorm(30,25,3),rnorm(30,35,3))
tbhigh <- 21
tblow <- 18
sim_day_temps <- function(tmin, tmax, timestep){
  temp <- ((tmax-tmin)/2) * cos((pi * (timestep - 0.5)) / 48) + ((tmax + tmin)/2)
}

suitable_period <- c(NULL)
for(day in 60:70){
  series_temps_day_left <- rep(NA, 48)
  series_temps_day_right <- rep(NA, 48)
  for(i in 1:48){
    series_temps_day_left[i] <- sim_day_temps(daily_min[day],
                                              daily_max[day],
                                              timestep = i)
    series_temps_day_right[i] <- if_else(day < length(daily_min),
                                         sim_day_temps(daily_min[day+1],
                                                       daily_max[day],
                                                       timestep = i),
                                         sim_day_temps(daily_min[day],
                                                       daily_max[day],
                                                       timestep = i))
  }
  series_temps_day <- c(rev(series_temps_day_left), series_temps_day_right)
  series_temps <- c(series_temps, series_temps_day)
  daily_suitable <- sum(series_logical <- if_else(series_temps_day <= tbhigh &
                                                    series_temps_day >= tblow,
                                                  1,
                                                  0))/4
  suitable_period <- c(suitable_period, daily_suitable)
}
sum(suitable_period)

# 2. applied example for Spain -------------------------------------------------------
## we use linear mixed effects results of thermal breadth
