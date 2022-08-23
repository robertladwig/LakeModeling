#' Created on Thu Aug 19 11:29:34 2021
#' 
#' @author: Robert Ladwig
#' @email: rladwig2@wisc.edu

## remove everything from workspace
rm(list = ls())

# set wd to current dir of script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## colors for plotting
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(LakeMetabolizer)
library(reshape2)
library(lubridate)

source('1D_HeatMixing_functions.R')

# Sky Pond
# https://github.com/bellaoleksy/FABs-DO-concept

## lake configurations
zmax = 5 # maximum lake depth
nx = 15 # number of layers we will have
dt = 3600  # 24 hours times 60 min/hour times 60 seconds/min
dx = zmax/nx # spatial step

## area and depth values of our lake 
hyps_all <- get_hypsography(hypsofile = 'bc/LakeEnsemblR_bathymetry_standard.csv',
                            dx = dx, nx = nx)
# hyps_all[[1]][15] = hyps_all[[1]][15] /1e2
# hyps_all[[3]][15] = hyps_all[[1]][15] * dx

df.obs <- read_csv('FABs-DO-concept/data/loch_buoy_long_temp_DO_2016-07-19_to_2018-09-11.csv')
df_obs = df.obs %>%
  dplyr::filter(parameter == 'temp') %>%
  arrange(dateTime) %>%
  select(datetime = dateTime, Depth_meter = depth, Water_Temperature_celsius = value)

df_obs_do = df.obs %>%
  dplyr::filter(parameter == 'DO') %>%
  arrange(dateTime) %>%
  select(datetime = dateTime, Depth_meter = depth, Dissolved_oxygen_gram_per_cubicMeter = value)


# colnames(df.obs) <- c('datetime', '0.5', '2.0', '4.0', '6.0', '7.0')
# df_obs <- melt(df.obs, id.vars = c("datetime")) %>%
#   arrange(datetime) %>%
#   rename(Depth_meter = variable, Water_Temperature_celsius = value)

write_csv(x = df_obs, file = 'bc/obs.txt')

# wnd <- read_csv('FABs-DO-concept/data/sky_2016_windSpeed.txt')
# par <- read_csv('FABs-DO-concept/data/sky_2016_PAR.txt')
# range_meteo <- range(wnd$dateTime)
dat <- read_csv('bc/WY1992to2019_LochVale_Climate.csv')

meteo <- read_csv('FABs-DO-concept/data/WY2016to2022_LochVale_WX.csv')
idx <- na.contiguous(meteo$SWin)
meteo <- meteo[match( as.POSIXct(as.Date(min(df_obs$datetime))), (as.POSIXct((meteo$dateTime),
                                                                             format = '%m/%d/%y %H:%M'))):
                 match( as.POSIXct(as.Date(max(df_obs$datetime))), (as.POSIXct((meteo$dateTime),
                                                                               format = '%m/%d/%y %H:%M'))) ,]
range_meteo <- range(as.POSIXct((meteo$dateTime),  format = '%m/%d/%y %H:%M'))

interpolated_sw <- approx(x = as.numeric(as.POSIXct((dat$timestamp), format = '%Y/%m/%d %H:%M:%S')), 
                          y = dat$SWin,
                               xout = as.numeric(as.POSIXct((meteo$dateTime),  format = '%m/%d/%y %H:%M')), rule = 2)$y
interpolated_rh <- approx(x = as.numeric(as.POSIXct((dat$timestamp),format = '%Y/%m/%d %H:%M:%S')), 
                          y = dat$RH,
                          xout = as.numeric(as.POSIXct((meteo$dateTime),format = '%m/%d/%y %H:%M')), rule = 2)$y
interpolated_lw <- approx(x = as.numeric(as.POSIXct((dat$timestamp),format = '%Y/%m/%d %H:%M:%S')), 
                          y = dat$LWin,
                          xout = as.numeric(as.POSIXct((meteo$dateTime),format = '%m/%d/%y %H:%M')), rule = 2)$y
interpolated_airt <- approx(x = as.numeric(as.POSIXct((meteo$dateTime),format = '%m/%d/%y %H:%M')), 
                          y = meteo$T_air,
                          xout = as.numeric(as.POSIXct((meteo$dateTime),format = '%m/%d/%y %H:%M')), rule = 2)$y


# airT <- read_csv('bc/mainwx_airT_2m_6m_daily_19821001_20180118.csv') 
# airT$date <- as.Date(airT$date, format = '%m/%d/%y')
# airT = airT %>%
#   dplyr::filter(date >= range_meteo[1] & date <= range_meteo[2])
# interpolated_2m_airT <- approx(x = as.numeric(as.POSIXct(paste(airT$date, '00:00:00'))), y = airT$Tave2M_C, 
#                                xout = as.numeric(par$dateTime[which(!is.na(par$PAR))]), rule = 2)$y

# df <- read_csv('bc/NLDAS2_Mendota_1979_2016_cell_5_GLMReady_cut_timezonecorr.csv') %>%
#   dplyr::filter(Date >= range_meteo[1] & Date <= range_meteo[2])
df = data.frame(matrix(ncol = 9, nrow = nrow(meteo)))
colnames(df) = c("datetime","Shortwave_Radiation_Downwelling_wattPerMeterSquared",
                 "Longwave_Radiation_Downwelling_wattPerMeterSquared",
                "Air_Temperature_celsius",
                "Relative_Humidity_percent",
                "Ten_Meter_Elevation_Wind_Speed_meterPerSecond",
                "Precipitation_millimeterPerDay", "Snowfall_millimeterPerDay",
                "Surface_Level_Barometric_Pressure_pascal")
df$Shortwave_Radiation_Downwelling_wattPerMeterSquared <-interpolated_sw #/ 2.16
df$Ten_Meter_Elevation_Wind_Speed_meterPerSecond <- (meteo$WSpd)
df$Air_Temperature_celsius <- interpolated_airt * 0.5
df$Relative_Humidity_percent <- interpolated_rh
df$Surface_Level_Barometric_Pressure_pascal = 98393
df$Longwave_Radiation_Downwelling_wattPerMeterSquared = interpolated_lw
df$Precipitation_millimeterPerDay = -999
df$Snowfall_millimeterPerDay =-999
df$datetime = as.POSIXct((meteo$dateTime), format = '%m/%d/%y %H:%M')

write_csv(x = df, file = 'bc/LakeEnsemblR_meteo_standard.csv')

## atmospheric boundary conditions
meteo_data <- provide_meteorology(meteofile = 'bc/LakeEnsemblR_meteo_standard.csv',
                    secchifile = NULL, 
                    windfactor = 1.)

meteo_all = list()

meteo_all[[1]] = meteo_data[[1]] %>%
  dplyr::filter(datetime >= range_meteo[1]  & datetime <= range_meteo[2] )

meteo_all[[1]]$dt <- as.POSIXct(meteo_all[[1]]$datetime) - (as.POSIXct(meteo_all[[1]]$datetime)[1]) + 1
# meteo_all[[1]]$Surface_Level_Barometric_Pressure_pascal = meteo_all[[1]]$Surface_Level_Barometric_Pressure_pascal + 3000



# df.obs1 <- read_csv('FABs-DO-concept/data/sky_2016_DO_0.5m.txt') %>%
#   dplyr::filter(dateTime >= range_meteo[1] & dateTime <= range_meteo[2])
# df.obs2 <- read_csv('FABs-DO-concept/data/sky_2016_DO_6.5m.txt') %>%
#   dplyr::filter(dateTime >= range_meteo[1] & dateTime <= range_meteo[2])
# df.obs_do <- data.frame('datetime' = df.obs1$dateTime, '0.5' = df.obs1$DO, '6.5' = df.obs2$DO)
# colnames(df.obs_do) <- c('datetime', '0.5', '6.5')
# df_obs_do <- melt(df.obs_do, id.vars = c("datetime")) %>%
#   arrange(datetime) %>%
#   rename(Depth_meter = variable, Dissolved_oxygen_gram_per_cubicMeter = value)

# df_obs_do_tofile <- df_obs_do
# colnames(df_obs_do_tofile) <- c("datetime","Depth_meter","Water_Temperature_celsius")
df_obs_do_tofile = df_obs_do
write_csv(x = df_obs_do_tofile, file = 'bc/obs_do.txt')
  

## here we define our initial profile
u_ini <- initial_profile(initfile = 'bc/obs.txt', nx = nx, dx = dx,
                     depth = hyps_all[[2]],
                     processed_meteo = meteo_all[[1]])

meteo_all[[2]] = meteo_data[[2]]

### EXAMPLE RUNS
# do the meteo interpolation all at the start
hydrodynamic_timestep = 24 * 3600 #24/4 * dt
total_runtime <- as.numeric(floor(range_meteo[2]-range_meteo[1]))
startingDate <- meteo_all[[1]]$datetime[1]

meteo = get_interp_drivers(meteo_all=meteo_all, total_runtime=total_runtime, 
                           hydrodynamic_timestep=hydrodynamic_timestep, dt=dt, method="integrate",
                           secchi =F)



temp <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
              nrow = nx)
avgtemp <- matrix(NA, ncol = 6,
                nrow = total_runtime * hydrodynamic_timestep/ dt)
temp_diff <- matrix(NA, ncol = (total_runtime * hydrodynamic_timestep/ dt) ,
               nrow = nx)
temp_mix <- matrix(NA, ncol = (total_runtime * hydrodynamic_timestep/ dt) ,
                    nrow = nx)
temp_conv <- matrix(NA, ncol = (total_runtime * hydrodynamic_timestep/ dt) ,
                    nrow = nx)
temp_ice <- matrix(NA, ncol = (total_runtime * hydrodynamic_timestep/ dt) ,
                    nrow = nx)
diff <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
               nrow = nx)
# meteo <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
#                nrow = 9)
buoyancy <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
               nrow = nx)
dissoxygen <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
               nrow = nx)

dissoxygen <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
                     nrow = nx)
icethickness <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
                         nrow = 1)
icelog <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
                       nrow = 1)


peri_kd <- rep(0, length( hyps_all[[2]]))
peri_kd[length(peri_kd)] = 0.8
km = peri_kd * 0# 0 * peri_kd

do_ini = (14.7 - 0.0017 * 4800) * exp(-0.0225*u_ini)
do_ini = initial_profile(initfile = 'bc/obs_do.txt', nx = nx, dx = dx,
                         depth = hyps_all[[2]],
                         processed_meteo = meteo_all[[1]])

if (exists('res')) {remove('res')}

for (i in 1:total_runtime){
  if (exists('res')){
    u = res$temp[, ncol(res$temp)]
    do = res$dissoxygen[, ncol(res$temp)]
    startTime = res$endtime
    endTime =  res$endtime + hydrodynamic_timestep -1
    ice = res$icefla
    Hi = res$icethickness 
    iceT = res$icemovAvg
    supercooled = res$supercooled
    kd_light = 0.3
    matrix_range = max(1, (startTime/dt)):((endTime/dt)) # this returns floats, not ints, after first round through? seems to cause issue down below in setting up avgtime
    matrix_range_start = max(1, round(startTime/dt) + 1)
    matrix_range_end = round(endTime/dt)
  } else {
    u = u_ini
    do = do_ini
    startTime = 1
    endTime = hydrodynamic_timestep - 1
    ice = FALSE
    Hi = 0
    iceT = 20#6
    supercooled = 0
    kd_light = 0.3
    matrix_range = max(1, (startTime/dt)):((endTime/dt)+1)
    matrix_range_start = max(1, round(startTime/dt))
    matrix_range_end = round(endTime/dt)
  }
  res <-  run_thermalmodel(u = u, 
                            startTime = startTime, 
                            endTime =  endTime, 
                            ice = ice, 
                            Hi = Hi, 
                            iceT = iceT,
                            supercooled = supercooled,
                           dt_iceon_avg = 0.1, 
                           Hgeo = 0.1,
                            kd_light = kd_light,
                            zmax = zmax,
                            nx = nx,
                            dt = dt,
                            dx = dx,
                            area = hyps_all[[1]], # area
                            depth = hyps_all[[2]], # depth
                            volume = hyps_all[[3]], # volume
                            daily_meteo = meteo[,matrix_range_start:matrix_range_end],
                            # secview = meteo_all[[2]],
                            Cd = 0.0013,
                            pgdl_mode = 'off',
                           scheme = 'implicit',
                           km = km,
                           do = do,
                           Fvol = 0.2, #0.01
                           Fred = 1.5, #0.005,##0.36,
                           Do2 = 1.08 * 10^(-4),
                           delta_DBL = 1/1000,
                           eff_area = seq(from = 1e-10,to = 1e-4,length.out = length(hyps_all[[1]]))
                           #  seq(from = 1e-20,to = 0.0002,length.out = length(hyps_all[[1]]))
                           )

  temp[, matrix_range_start:matrix_range_end] =  res$temp
  diff[, matrix_range_start:matrix_range_end] =  res$diff
  avgtemp[matrix_range_start:matrix_range_end,] <- as.matrix(res$average)
  buoyancy[, matrix_range_start:matrix_range_end] =  res$temp
  dissoxygen[, matrix_range_start:matrix_range_end] =  res$dissoxygen
  icethickness[, matrix_range] =  res$icethickness_matrix
  icelog[, matrix_range] = res$iceflag
  
  average <- res$average %>%
    mutate(datetime = as.POSIXct(startingDate + time),
           Date = as.Date(datetime, format = "%m/%d/%Y")) %>%
    summarise_all(mean)
}

time =  startingDate + seq(1, ncol(temp), 1) * dt


surf_z = which.min( abs(seq(0, zmax, length = nx)  - min(df_obs$Depth_meter) ))
bottom_z = which.min( abs(seq(0, zmax, length = nx)  - max(df_obs$Depth_meter) ))
surf_temp = df_obs %>% filter(Depth_meter == min(Depth_meter))
bottom_temp = df_obs %>% filter(Depth_meter == max(Depth_meter))

plot(time, temp[surf_z,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(-10,20), lwd = 2)
points(surf_temp$datetime, surf_temp$Water_Temperature_celsius)

plot(time, temp[bottom_z,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(-10,20), lwd = 2)
points(bottom_temp$datetime, bottom_temp$Water_Temperature_celsius)

surf_z = which.min( abs(seq(0, zmax, length = nx)  - min( df_obs_do$Depth_meter) ))
bottom_z = which.min( abs(seq(0, zmax, length = nx)  - max( df_obs_do$Depth_meter) ))
surf_do =  df_obs_do %>% filter(Depth_meter == min(Depth_meter))
bottom_do =  df_obs_do %>% filter(Depth_meter == max(Depth_meter))

plot(time, dissoxygen[surf_z,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='DO (g/m3)', ylim=c(0,20), lwd = 2)
points(surf_do$datetime, surf_do$Dissolved_oxygen_gram_per_cubicMeter)

plot(time, dissoxygen[bottom_z,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='DO (g/m3)', ylim=c(0,20), lwd = 2)
points(bottom_do$datetime, bottom_do$Dissolved_oxygen_gram_per_cubicMeter)

idx = match(surf_temp$datetime, bottom_temp$datetime)
plot(bottom_temp$datetime, 100 *(calc_dens(bottom_temp$Water_Temperature_celsius) -  calc_dens( surf_temp$Water_Temperature_celsius[!is.na(idx)])), col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temp)', ylim=c(0,20), lwd = 2) 
points(bottom_do$datetime, bottom_do$Dissolved_oxygen_gram_per_cubicMeter)

icelog_ts = ifelse(icelog == T,1,0)
plot(seq(1, ncol(temp))*dt/24/3600,icelog_ts)



plot(seq(1, ncol(temp))*dt/24/3600, icethickness)

# plotting for checking model output and performance
plot(seq(1, ncol(temp))*dt/24/3600, temp[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(-1,30), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(temp))*dt/24/3600, temp[i,], 
        lty = 'dashed',lwd =2)
}

plot(seq(1, ncol(dissoxygen))*dt/24/3600, dissoxygen[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Diss. Oxygen (g/m3)', ylim=c(0,25), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(dissoxygen))*dt/24/3600, dissoxygen[i,], 
        lty = 'dashed',lwd =2)
}

plot(seq(1, ncol(temp))*dt/24/3600, (14.7 - 0.0017 * 4800) * exp(-0.0225*temp[1,]), col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Saturated Diss. Oxygen (g/m3)', ylim=c(3,10), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(temp))*dt/24/3600, (14.7 - 0.0017 * 4800) * exp(-0.0225*temp[i,]), 
        lty = 'dashed',lwd =2)
}




avgtemp = as.data.frame(avgtemp)
colnames(avgtemp) = c('time', 'epi', 'hyp', 'tot', 'stratFlag', 'thermoclineDep')

ggplot(avgtemp) +
  geom_line(aes(time, epi, col = 'epilimnion')) +
  geom_line(aes(time, hyp, col = 'hypolimnion')) +
  geom_line(aes(time, tot, col = 'total')) +
  theme_minimal() 
ggplot(avgtemp) +
  geom_line(aes(time, stratFlag)) +
  theme_minimal() 
ggplot(avgtemp) +
  geom_line(aes(time, thermoclineDep)) +
  geom_point(aes(time, stratFlag, col = as.factor(stratFlag))) +
  scale_y_reverse() +
  theme_minimal() 

df <- data.frame(cbind(time, t(temp)) )
colnames(df) <- c("time", as.character(paste0(seq(1,nrow(temp)))))
m.df <- reshape2::melt(df, "time")
m.df$time <- time

g1 <- ggplot(m.df, aes((time), as.numeric(variable)*dx)) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(-2,15),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Temp [degC]')+
  scale_y_reverse() ; g1

df.do <- data.frame(cbind(time, t(dissoxygen)) )
colnames(df.do) <- c("time", as.character(paste0(seq(1,nrow(dissoxygen)))))
m.df.do <- reshape2::melt(df.do, "time")
m.df.do$time <- time

g2 <- ggplot(m.df.do, aes((time), as.numeric(variable)*dx)) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(0,15),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Oxygen [g/m3]')+
  scale_y_reverse() ;g2

g1 / g2


## vertical temperature profiles
for (i in seq(1,ncol(temp), length.out = 300)){
  n = i
  i = floor(i)
  
  sim = m.df %>% 
    dplyr::filter(time == time[i]) %>%
    mutate(dosat = (14.7 - 0.0017 * 4800) * exp(-0.0225*value))
  sim.do = m.df.do %>% 
    dplyr::filter(time == time[i]) 
  
  obs = df_obs %>% 
    dplyr::filter(datetime == time[i]) %>%
    mutate(dosat = (14.7 - 0.0017 * 4800) * exp(-0.0225*Water_Temperature_celsius))
  if (length(obs$Water_Temperature_celsius) == 0){
    obs = data.frame('datetime' = time[i], 'Depth_meter' = 0, 'Water_Temperature_celsius' = NA,
                     'dosat' = NA)
  }
  obs.do = df_obs_do %>% 
    dplyr::filter(datetime == time[i]) 
  if (length(obs.do$Dissolved_oxygen_gram_per_cubicMeter) == 0){
    obs.do = data.frame('datetime' = time[i], 'Depth_meter' = 0, 'Dissolved_oxygen_gram_per_cubicMeter' = NA)
  }

  ggplot() +
    geom_path(data = sim, aes(value, 
                              as.numeric(variable) * dx, col = 'T (degC)'), size = 1.2) +
    geom_point(data = obs, aes(as.numeric(Water_Temperature_celsius), 
                               as.numeric(as.character(Depth_meter)), col = 'obs T (degC)'), size = 1.2) +
    geom_path(data = sim, aes(as.numeric(dosat) * 2, 
                              as.numeric(variable) * dx, col = 'saturated DO (g/m3)'), size = 1.2,
              linetype  = 'dashed') +
    geom_point(data = obs, aes(as.numeric(dosat), 
                               as.numeric(as.character(Depth_meter)), col = 'obs saturated DO (g/m3)'), size = 1.2) +
    geom_path(data = sim.do, aes(as.numeric(value) * 2, 
                              as.numeric(variable) * dx, col = 'DO (g/m3)'), size = 1.2) +
    geom_point(data = obs.do, aes(as.numeric(Dissolved_oxygen_gram_per_cubicMeter) * 2, 
                                  as.numeric(as.character(Depth_meter)) , col = 'obs DO (g/m3)'), size = 1.2) +
    xlab('temp. (deg C)') + ylab('depth (m)')+
    scale_y_reverse() +
    scale_color_manual(values = c("#56B4E9", "#56B4E9", "lightblue",'red', 'lightblue', 'red')) +
    scale_x_continuous(sec.axis = sec_axis(~ . /2, name = "diss. oxygen (g/m3)"), limits = c(-2,30)) +
    ggtitle( time[i]) + 
    labs(col='') +
    # xlim(3, 30) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(paste0('animation/pic_',match(n, seq(1,ncol(temp), length.out = 300)),'.png'),
         width = 4, height = 5, units = 'in')
  
}
