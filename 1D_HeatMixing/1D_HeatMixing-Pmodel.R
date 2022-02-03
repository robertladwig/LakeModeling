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

source('1D_HeatMixing_functions.R')

## lake configurations
zmax = 25 # maximum lake depth
nx = 25 # number of layers we will have
dt = 3600 # 24 hours times 60 min/hour times 60 seconds/min
dx = zmax/nx # spatial step

## area and depth values of our lake 
hyps_all <- get_hypsography(hypsofile = 'bc/LakeEnsemblR_bathymetry_standard.csv',
                            dx = dx, nx = nx)

## here we define our initial profile
u <- initial_profile(initfile = 'bc/obs.txt', nx = nx, dx = dx,
                     depth = hyps_all[[3]])

## atmo8,spheric boundary conditions
meteo_all <- provide_meteorology(meteofile = 'bc/LakeEnsemblR_meteo_standard.csv',
                    secchifile = 'bc/light.csv', windfactor = 0.8)

### EXAMPLE RUNS
# 1 day
temp <- c()
avgtemp <- c()
res <- run_thermalmodel(u = u, 
                        startTime = 1, 
                        endTime = 24*3600, 
                        kd_light = NULL,
                        zmax = zmax,
                        nx = nx,
                        dt = dt,
                        dx = dx,
                        area = hyps_all[[1]], # area
                        depth = hyps_all[[2]], # depth
                        volume = hyps_all[[3]], # volume
                        daily_meteo = meteo_all[[1]],
                        secview = meteo_all[[2]])
temp <-cbind(temp, res$temp)
avgtemp <- rbind(avgtemp, res$average)
# doing another day
for (i in 1:365){
  res <-  run_thermalmodel(u = res$temp[, ncol(res$temp)], 
                            startTime = res$endtime, 
                            endTime =  res$endtime + 24*3600, 
                            ice = res$iceflag, 
                            Hi = res$icethickness, 
                            iceT = res$icemovAvg,
                            supercooled = res$supercooled,
                            kd_light = NULL,
                            zmax = zmax,
                            nx = nx,
                            dt = dt,
                            dx = dx,
                            area = hyps_all[[1]], # area
                            depth = hyps_all[[2]], # depth
                            volume = hyps_all[[3]], # volume
                            daily_meteo = meteo_all[[1]],
                            secview = meteo_all[[2]])
  temp <-cbind(temp, res$temp[,-1])
  avgtemp <- rbind(avgtemp, res$average[-1,])
}


plot(seq(1, ncol(temp))*dt/24/3600, temp[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(-1,35), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(temp))*dt/24/3600, temp[i,], 
        lty = 'dashed',lwd =2)
}
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
  theme_minimal() 


