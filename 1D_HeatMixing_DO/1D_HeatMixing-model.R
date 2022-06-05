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

source('1D_HeatMixing_functions.R')

## lake configurations
zmax = 3 # maximum lake depth
nx = 10 # number of layers we will have
dt = 3600  # 24 hours times 60 min/hour times 60 seconds/min
dx = zmax/nx # spatial step

## area and depth values of our lake 
hyps_all <- get_hypsography(hypsofile = 'bc/LakeEnsemblR_bathymetry_standard.csv',
                            dx = dx, nx = nx)

## atmospheric boundary conditions
meteo_data <- provide_meteorology(meteofile = 'bc/LakeEnsemblR_meteo_standard.csv',
                    secchifile = 'bc/light.csv', 
                    windfactor = 0.8)

meteo_all = list()

meteo_all[[1]] = meteo_data[[1]] %>%
  dplyr::filter(datetime >= '2009-07-01 00:00:00' & datetime <= '2009-07-07 23:30:00')

meteo_all[[1]]$dt <- as.POSIXct(meteo_all[[1]]$datetime) - (as.POSIXct(meteo_all[[1]]$datetime)[1]) + 1
meteo_all[[1]]$Surface_Level_Barometric_Pressure_pascal = meteo_all[[1]]$Surface_Level_Barometric_Pressure_pascal + 3000

## here we define our initial profile
u_ini <- initial_profile(initfile = 'bc/obs.txt', nx = nx, dx = dx,
                     depth = hyps_all[[2]],
                     processed_meteo = meteo_all[[1]])

meteo_all[[2]] = meteo_data[[2]]

### EXAMPLE RUNS
hydrodynamic_timestep = 24 * 3600 #24/4 * dt
total_runtime <- 7
startingDate <- meteo_all[[1]]$datetime[1]

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

# do the meteo interpolation all at the start
meteo = get_interp_drivers(meteo_all=meteo_all, total_runtime=total_runtime, 
                           hydrodynamic_timestep=hydrodynamic_timestep, dt=dt, method="integrate")

peri_kd <- rep(0, length( hyps_all[[2]]))
peri_kd[length(peri_kd)] = 1
km = 2.0 * peri_kd

do_ini = rep(10, length(u_ini))

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
    kd_light = 0.1
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
    iceT = 6
    supercooled = 0
    kd_light = 0.1
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
                            Cd = 0.0008,
                            pgdl_mode = 'off',
                           scheme = 'implicit',
                           km = km,
                           do = do)

  temp[, matrix_range_start:matrix_range_end] =  res$temp
  diff[, matrix_range_start:matrix_range_end] =  res$diff
  avgtemp[matrix_range_start:matrix_range_end,] <- as.matrix(res$average)
  buoyancy[, matrix_range_start:matrix_range_end] =  res$temp
  dissoxygen[, matrix_range_start:matrix_range_end] =  res$dissoxygen
  
  average <- res$average %>%
    mutate(datetime = as.POSIXct(startingDate + time),
           Date = as.Date(datetime, format = "%m/%d/%Y")) %>%
    summarise_all(mean)
  
}


# plotting for checking model output and performance
plot(seq(1, ncol(temp))*dt/24/3600, temp[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(17,30), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(temp))*dt/24/3600, temp[i,], 
        lty = 'dashed',lwd =2)
}

plot(seq(1, ncol(dissoxygen))*dt/24/3600, dissoxygen[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Diss. Oxygen (g/m3)', ylim=c(0,10), lwd = 2)
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


time =  startingDate + seq(1, ncol(temp), 1) * dt
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

ggplot(m.df, aes((time), as.numeric(variable)*dx)) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(-2,35),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Temp [degC]')+
  scale_y_reverse() 

df.do <- data.frame(cbind(time, t(dissoxygen)) )
colnames(df.do) <- c("time", as.character(paste0(seq(1,nrow(dissoxygen)))))
m.df.do <- reshape2::melt(df.do, "time")
m.df.do$time <- time

ggplot(m.df.do, aes((time), as.numeric(variable)*dx)) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(0,10),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Oxygen [g/m3]')+
  scale_y_reverse() 


## vertical temperature profiles
for (i in seq(1,ncol(temp), length.out = 50)){
  n = i
  i = floor(i)
  
  sim = m.df %>% 
    filter(time == time[i]) %>%
    mutate(dosat = (14.7 - 0.0017 * 4800) * exp(-0.0225*value))
  sim.do = m.df.do %>% 
    filter(time == time[i]) 

  ggplot() +
    geom_path(data = sim, aes(value, 
                              as.numeric(variable) * dx, col = 'T (degC)'), size = 1.2) +
    geom_path(data = sim, aes(dosat * 2, 
                              as.numeric(variable) * dx, col = 'saturated DO (g/m3)'), size = 1.2,
              linetype  = 'dashed') +
    geom_path(data = sim.do, aes(value * 2, 
                              as.numeric(variable) * dx, col = 'DO (g/m3)'), size = 1.2) +
    xlab('temp. (deg C)') + ylab('depth (m)')+
    scale_y_reverse() +
    scale_color_manual(values = c("#56B4E9", "lightblue", "red")) +
    scale_x_continuous(sec.axis = sec_axis(~ . /2, name = "diss. oxygen (g/m3)"), limits = c(0,30)) +
    ggtitle( time[i]) + 
    labs(col='') +
    # xlim(3, 30) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(paste0('animation/pic_',match(n, seq(1,ncol(temp), length.out = 50)),'.png'),
         width = 4, height = 5, units = 'in')
  
}
