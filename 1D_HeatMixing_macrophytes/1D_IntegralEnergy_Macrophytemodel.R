#' Created on Thu Aug 19 11:29:34 2021
#' 
#' @author: Robert Ladwig
#' @email: rladwig2@wisc.edu
#' 
#' Temperature transport equation as Tt = 1/A K Tzz
#' To run the model you will need (a) one initial water temperature profile,  
#' (b) time series data of meteorological drivers: air temperature, wind speed 
#' and short-wave radiation, and
#' (c) your lake's hypsography (area over depth), or at least maximum deoth and
#' surface area
#' 
#' Diffusion code is based on 12 steps to Navier-Stokes by (c) Lorena A. Barba, 
#' Gilbert F. Forsyth 2017.
#' https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/ 
#' 
#' Eddy diffusivity is estimated from buoyancy frequency according to Hondzo and 
#' Stefan (1993) Lake Water Temperature Simulation Model. ASCE
#' https://doi.org/10.1061/(ASCE)0733-9429(1993)119:11(1251) 
#' 
#' Mixing dynamics code is taken from Herb & Stefan (2004) Temperature 
#' stratification and Mixing Dynamics in a Shallow Lake with Submersed 
#' Macrophytes. Lake & Reservoir Management
#' https://www.tandfonline.com/doi/pdf/10.1080/07438140409354159
#' 
#' Convective overturn algorithm is taken from Saloranta & Andersen (2007) 
#' MyLakeâA multi-year lake simulation model code suitable for uncertainty 
#' and sensitivity analysis simulations. Ecological Modeling
#' https://doi.org/10.1016/j.ecolmodel.2007.03.018 
#' 
#' Ice formation and growth/decay code is taken from Saloranta & Andersen (2007) 
#' MyLakeâA multi-year lake simulation model code suitable for uncertainty 
#' and sensitivity analysis simulations. Ecological Modeling
#' https://doi.org/10.1016/j.ecolmodel.2007.03.018 

## remove everything from workspace
#rm(list = ls())

# set wd to current dir of script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source('1D_macrophyte-functions.R')

#devtools::install_github("aemon-j/gotmtools")

## colors for plotting
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(rLakeAnalyzer)

## lake configurations
zmax = 2 # maximum depth
nx = 20 # number of layers we will have
dt = 3600 # time step, 60 min times 60 seconds/min
dx = zmax/nx # spatial step


# area and depth values of our lake 
hyps_all <- get_hypsography(hypsofile = 'bc/hyps.csv',
                            dx = dx, nx = nx)


## atmospheric boundary conditions
meteo_all <- provide_meteorology(meteofile = 'bc/ZOO799_GLM_MetData.csv',
                                 secchifile = 'bc/secchi.csv',
                                 pressurefile = 'bc/pressure.csv')

# initial water temperature profile
u <- initial_profile(initfile = 'bc/temp_profiles_2020.csv',
                     nx = nx,
                     dx = dx,
                     processed_meteo = meteo_all[[1]])




# macrophyte
macro_all <- get_macrophyte(canopyfile = 'bc/canopy_F19.csv',
                            biomassfile = 'bc/biomass_F19.csv',
                            processed_meteo = meteo_all[[1]])

nt = as.double(max(meteo_all[[1]]$dt)) # maximum time duration


# u = u
# startTime = 1
# endTime = nt
# kd_light = NULL
# zmax = zmax
# nx = nx
# dt = dt
# dx = dx
# area = hyps_all[[1]] # area
# depth = hyps_all[[2]] # depth
# volume = hyps_all[[3]] # volume
# meteo = meteo_all[[1]] # meteorology
# light = meteo_all[[2]] # light
# pressure = meteo_all[[3]] # pressure
# canpy = macro_all[[1]] # canopy height
# biomass = macro_all[[2]] # macrophyte density
# Cd = 0.0013 # wind momentum drag
# km = 0 #0.04 # macrophyte light extinction
# Cdplant = 1e3 # macrophyte momentum drag
# ahat = 0.02 # macrophyte area to volume
# reflect = 0.6
# infra = 0.4
# windfactor = 1.2 # wind multiplier
# shortwavefactor = 1 # shortwave radiation multiplier
# diffusionfactor = 1 # diffusion multiplier
# diffmethod = 2
# ice = FALSE
# Hi = 0
# iceT = 6
# supercooled = 0
# scheme = 'explicit'
# kd_light = NULL
# densThresh = 1e-3
# reflect = 0.3
# infra = 0.7
# eps = 0.97
# emissivity = 0.97
# sigma = 5.67 * 10^(-8)
# p2 = 1
# B = 0.61
# g = 9.81
# meltP = 1
# dt_iceon_avg = 0.8
# Hgeo = 0.1 # geothermal heat
# KEice = 1/1000
# Ice_min = 0.1

### RUNS
temp <- c()
diff <- c()
buoy <- c()
macroz <- c()
mixing <- c()
res <- run_thermalmodel(u = u, 
                        startTime = 1, 
                        endTime = nt, 
                        kd_light = NULL,
                        zmax = zmax,
                        nx = nx,
                        dt = dt,
                        dx = dx,
                        area = hyps_all[[1]], # area
                        depth = hyps_all[[2]], # depth
                        volume = hyps_all[[3]], # volume
                        meteo = meteo_all[[1]], # meteorology
                        light = meteo_all[[2]], # light
                        pressure = meteo_all[[3]], # pressure
                        canpy = macro_all[[1]], # canopy height
                        biomass = macro_all[[2]], # macrophyte density
                        Cd = 0.0013, # wind momentum drag
                        km = 0.2, #0.008,#0.04, # macrophyte light extinction
                        Cdplant = 1,#1e3, # macrophyte momentum drag
                        ahat = 0.5, # macrophyte area to volume
                        reflect = 0.6, 
                        infra = 0.4,
                        windfactor = 1, # wind multiplier 
                        shortwavefactor = 0.8, # shortwave radiation multiplier
                        diffusionfactor = 1, # diffusion multiplier
                        diffmethod = 2,
                        densThresh = 1e-2,
                        Hgeo = 1,
                        rho_plant = 30
                        )
temp <-cbind(temp, res$temp)
diff <-cbind(diff, res$diff)
buoy <-cbind(buoy, res$buoyancy)
macroz <-cbind(macroz, res$macroheight)
mixing <-cbind(mixing, res$mixing)

time =  meteo_all[[1]]$Date[1] + seq(1, ncol(temp))*dt
df <- data.frame(cbind(time, t(temp)) )
colnames(df) <- c("time", as.character(paste0(seq(1,nrow(temp)))))
m.df <- reshape2::melt(df, "time")
m.df$time <- time

df.kz <- data.frame(cbind(time, t(diff)) )
colnames(df.kz) <- c("time", as.character(paste0(seq(1,nrow(diff))*dx)))
m.df.kz <- reshape2::melt(df.kz, "time")
m.df.kz$time <- time

df.n2 <- data.frame(cbind(time, t(buoy)) )
colnames(df.n2) <- c("time", as.character(paste0(seq(1,nrow(buoy))*dx)))
m.df.n2 <- reshape2::melt(df.n2, "time")
m.df.n2$time <- time

g1 <- ggplot(m.df, aes((time), as.numeric(as.character(variable)))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(10,45),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Temp [degC]')+
  scale_y_reverse() 
g2 <- ggplot(m.df.kz, aes((time), as.numeric(as.character(variable)))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(#limits = c(-20,35),
    colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Diffusion [m2/s]')+
  scale_y_reverse() 
g3 <- ggplot(m.df.n2, aes((time), as.numeric(as.character(variable)))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(#limits = c(-20,35),
    colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'N2 [s-2]')+
  scale_y_reverse() 
g <- g1 / g2 / g3 ; g
ggsave(filename = paste0('heatmap_temp_diff_buoy','.png'),plot = g, width = 15, height = 8, units = 'in')



## model goodness
obs <- read.csv('bc/temp_profiles_2020.csv')
obs <- obs %>%
  filter(pond == 'F', site_id == '19') %>%
  select(datetime, temp_depth_m, temp_c)

df.sim <- df
colnames(df.sim) <- c("datetime", as.character(paste0('temp_c.',seq(1,nrow(temp))*dx)))
df.sim$datetime <-  time

idx <- (match(as.POSIXct(obs$datetime), as.POSIXct(df.sim$datetime) ))


obs <- obs[which(!is.na(idx)), ]

deps <-  hyps_all[[2]]
if (min(unique(obs$temp_depth_m)) < min(deps)){
  deps[which.min(deps)] <- min(unique(obs$temp_depth_m)) 
}
if (max(unique(obs$temp_depth_m)) > max(deps)){
  deps[which.max(deps)] <- max(unique(obs$temp_depth_m)) 
}

df.sim.interp <- NULL
for (i in 1:nrow(df.sim)){
  df.sim.interp <- rbind(df.sim.interp,
                         approx(deps, df.sim[i, -1], unique(obs$temp_depth_m))$y)
}
df.sim.interp <- as.data.frame(df.sim.interp)
df.sim.interp$datetime <-   time
colnames(df.sim.interp) <- c(as.character(paste0('temp_c.',unique(obs$temp_depth_m))), 'datetime')

wide.obs <- reshape(obs, idvar = "datetime", timevar = "temp_depth_m", direction = "wide")
m.obs <- reshape2::melt(wide.obs, id = 'datetime')
m.obs$datetime <- as.POSIXct(m.obs$datetime)
m.obs$group <- 'obs'
m.df.sim.interp <- reshape2::melt(df.sim.interp, id = 'datetime')
m.df.sim.interp$group <- 'sim'


rmse <- data.frame('variable' = NULL, 'fit' = NULL)
for (i in unique(as.character(m.obs$variable))){
  o <- m.obs
  o$variable <- as.character(o$variable)
  o = o %>%
    filter(variable == i) 
  s <- m.df.sim.interp
  s$variable <- as.character(s$variable)
  s = s %>%
    filter(variable == i)
  rmse <- rbind(rmse, data.frame('variable' = i,
                                 'fit' = sqrt((sum((o$value-s$value)**2))/nrow(o))))
}

ggplot() +
  geom_line(data = m.obs,aes(datetime, value, col = group)) +
  geom_line(data = m.df.sim.interp, aes(datetime, value, col = group)) +
  geom_text(data=rmse, aes( as.POSIXct("2020-05-26 10:30:00 CDT"), y=17, label=round(fit,2)),                 
            color="black", size =3) +
  facet_wrap(~ variable) +
  xlab('') + ylab('Temp. (deg C)')+
  theme_bw()
ggsave(filename = paste0('fieldcomparison.png'), width = 15, height = 8, units = 'in')



## vertical temperature profiles
# for (i in seq(1,ncol(temp), length.out = 200)){
#   n = i
#   i = floor(i)
# 
#   as.numeric(gsub(".*temp_c.","",as.character(factor(m.df.sim.interp$variable))))
# 
#   mc <- data.frame('datetime' = rep(time[1],length(seq(2,(2-macroz[i]),-0.1))),
#                    'canopy' = seq(2,(2-macroz[i]),-0.1),
#                    'value' =rep(30,length(seq(2,(2-macroz[i]),-0.1))))
#   m.mc <- reshape2::melt(mc, id = c('datetime','value'))
# 
#   sim = m.df.sim.interp %>%
#     filter(datetime == time[i]) %>%
#     mutate(group = 'modeled') %>%
#     mutate(depth =  unique(  as.numeric(gsub(".*temp_c.","",as.character(factor(variable))))))
#   obs = m.obs %>%
#     filter(datetime == time[i]) %>%
#     mutate(group = 'observed') %>%
#     mutate(depth =  unique(  as.numeric(gsub(".*temp_c.","",as.character(factor(variable))))))
# 
#   ggplot() +
#     geom_path(data = sim, aes(value,
#                               depth, col = group), size = 1.2) +
#     # facet_wrap(~ factor(variable, level = c(paste0('wtemp.',seq(0,24,1)))), scales = 'free') +
#     geom_line(data = mc, aes(value, canopy, col = 'macro'), size = 1.2, linetype=2) +
#     geom_point(data = obs ,aes(value, depth, col = group), size =1.2) +
#     xlab('temp. (deg C)') + ylab('depth (m)')+
#     scale_y_reverse() +
#     scale_color_manual(values = c('#69b3a2',"#E69F00", "#56B4E9")) +
#     ggtitle( time[i]) +
#     labs(col='') +
#     xlim(15, 40) +
#     theme_bw()
# 
#   ggsave(paste0('../../animation_macrophyte/pic_',match(n, seq(1,ncol(temp),length.out=200)),'.png'),
#          width = 4, height = 5, units = 'in')
# 
# }

