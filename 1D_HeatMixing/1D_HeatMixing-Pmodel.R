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

## atmospheric boundary conditions
meteo_all <- provide_meteorology(meteofile = 'bc/LakeEnsemblR_meteo_standard.csv',
                    secchifile = 'bc/light.csv', windfactor = 0.8)


### EXAMPLE RUNS
# 1 day
temp <- c()
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
temp <- rbind(temp, res$average)
# doing another day
res2 <-  run_thermalmodel(u = res$temp[, ncol(res$temp)], 
                          startTime = res$endtime, 
                          endTime =  res$endtime + 24*3600, 
                          ice = res$iceflag, 
                          Hi = res$icethickness[length(res$icethickness)], 
                          iceT = res$icemovAvg,
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
temp <- rbind(temp, res2$average[-1,])
# doing another, another day
res3 <-  run_thermalmodel(u = res2$temp[, ncol(res2$temp)], 
                          startTime = res2$endtime, 
                          endTime =  res2$endtime + 24*3600, 
                          ice = res2$iceflag, 
                          Hi = res2$icethickness[length(res2$icethickness)], 
                          iceT = res2$icemovAvg,
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
temp <- rbind(temp, res3$average[-1,])

ggplot(temp) +
  geom_line(aes(time, epi, col = 'epilimnion')) +
  geom_line(aes(time, hyp, col = 'hypolimnion')) +
  geom_line(aes(time, tot, col = 'total')) +
  theme_minimal() /
ggplot(temp) +
  geom_line(aes(time, stratFlag)) +
  theme_minimal() 


