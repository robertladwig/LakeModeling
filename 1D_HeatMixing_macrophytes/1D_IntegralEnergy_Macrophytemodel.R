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
rm(list = ls())

# set wd to current dir of script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source('1D_macrophyte-functions.R')

#devtools::install_github("aemon-j/gotmtools")

## colors for plotting
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(rLakeAnalyzer)
library(adagio)

## lake configurations
zmax = 2 # maximum depth
nx = 8 # number of layers we will have
dt = 3600 # time step, 60 min times 60 seconds/min
dx=zmax/nx


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

# dx = 0.5 # reduce spatial step
# try reducing number of layers (dx=0.1, 0.2, to 0.5)


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


nt = 2045*  3600

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
                        km = 0.025, #0.008,#0.04, # macrophyte light extinction
                        Cdplant = 1e3,#1e3, # macrophyte momentum drag
                        ahat = 0.2, # macrophyte area to volume
                        reflect = 0.6, 
                        infra = 0.4,
                        windfactor = 1., # wind multiplier 
                        shortwavefactor = 0.7, # shortwave radiation multiplier
                        diffusionfactor = 1, # diffusion multiplier
                        tempfactor = 1.1, 
                        diffmethod = 1,
                        densThresh = 1e-2,
                        Hgeo = 0,
                        rho_plant = 300,
                        scheme = 'implicit'
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

m.df2<-m.df %>% 
  mutate(variable=as.numeric(variable)) %>% 
  mutate(depth=variable*dx)

df.kz <- data.frame(cbind(time, t(diff)) )
colnames(df.kz) <- c("time", as.character(paste0(seq(1,nrow(diff))*dx)))
m.df.kz <- reshape2::melt(df.kz, "time")
m.df.kz$time <- time

df.n2 <- data.frame(cbind(time, t(buoy)) )
colnames(df.n2) <- c("time", as.character(paste0(seq(1,nrow(buoy))*dx)))
m.df.n2 <- reshape2::melt(df.n2, "time")
m.df.n2$time <- time

g1 <- ggplot(m.df, aes((time), dx*as.numeric(as.character(variable)))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(10,30),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Temp [degC]')+
  scale_y_reverse() 
g2 <- ggplot(m.df.kz, aes((time), dx * as.numeric(as.character(variable)))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(#limits = c(-20,35),
    colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Diffusion [m2/s]')+
  scale_y_reverse() 
g3 <- ggplot(m.df.n2, aes((time), dx * as.numeric(as.character(variable)))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(#limits = c(-20,35),
    colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'N2 [s-2]')+
  scale_y_reverse() 
g <- g1 / g2 / g3 ; g
ggsave(filename = paste0('heatmap_temp_diff_buoy','.png'),plot = g, width = 15, height = 8, units = 'in')


### MANUSCRIPT FIGURE ---------------------------------------------------------------------------------------
# MAIN HEAT MAP FIGURE PANEL
# data frame = m.df2
# convert time to day of year for consistency 
library(lubridate)
m.df2$doy<-yday(m.df2$time)
m.df2$doy_frac<-hour(m.df2$time)
m.df2$min<-minute(m.df2$time)
m.df2$minfrac[m.df2$min=="30"] <- 0.5
m.df2$minfrac[m.df2$min=="0"] <- 0
m.df2$hourfrac<-(m.df2$doy_frac + m.df2$minfrac)/24
m.df2$doy_frac<-m.df2$doy+m.df2$hourfrac

m.df2<-m.df2 %>% 
  select(time,doy,doy_frac,value, depth) %>%
  rename(temp_c=value)

output_heatmap<-
ggplot(m.df2, aes((doy_frac), depth )) +
  geom_raster(aes(fill = as.numeric(temp_c)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(10,36),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_linedraw()  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  xlab('Day of Year') + 
  ylab('Depth (m)') + theme(legend.position="none")+
  labs(fill = 'Temperature (C)')+  scale_y_reverse() 

ggsave("no_mac_heatmap.png", output_heatmap, width=6, height=3.5, units="in")  
### -----------------------------------------------------------------------------------------------------------

## model goodness
obs <- read.csv('bc/temp_profiles_2020.csv')
obs <- obs %>%
  filter(pond == 'B', site_id == '20') %>%
  select(datetime, temp_depth_m, temp_c)

df.sim <- df
colnames(df.sim) <- c("datetime", as.character(paste0('temp_c.',seq(1,nrow(temp))*dx)))
df.sim$datetime <-  time

idx <- (match(as.POSIXct(obs$datetime), as.POSIXct(df.sim$datetime) ))


obs <- obs[which(!is.na(idx)), ]

deps <- seq(1,nrow(temp))*dx
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

# ORIGINAL PLOT
ggplot() +
  geom_line(data = m.obs,aes(datetime, value, col = group)) +
  geom_line(data = m.df.sim.interp, aes(datetime, value, col = group)) +
  geom_text(data=rmse, aes( as.POSIXct("2020-05-26 10:30:00 CDT"), y=17, label=round(fit,2)),                 
            color="black", size =3) +
  facet_wrap(~ variable) +
  xlab('') + ylab('Temp. (deg C)')+
  theme_bw()
ggsave(filename = paste0('fieldcomparison.png'), width = 15, height = 8, units = 'in')

### MANUSCRIPT MODEL FIT ------------------------------------------------------------------------------------
# Data frames: observed values: m.obs; model = m.df.sim.interp
# convert time to day of year for consistency 
#library(lubridate)
m.obs$doy<-yday(m.obs$datetime)
m.obs$doy_frac<-hour(m.obs$datetime)
m.obs$min<-minute(m.obs$datetime)
m.obs$minfrac[m.obs$min=="30"] <- 0.5
m.obs$minfrac[m.obs$min=="0"] <- 0
m.obs$hourfrac<-(m.obs$doy_frac + m.obs$minfrac)/24
m.obs$doy_frac<-m.obs$doy+m.obs$hourfrac

m.obs.fig<-m.obs %>% 
  select(datetime,doy,doy_frac,value, variable,group) %>%
  rename(temp_c=value) %>% 
  mutate(variable=(recode(variable, "temp_c.0" = "0 m", "temp_c.0.25"="0.25 m", "temp_c.0.5"="0.5 m",
                "temp_c.0.75"="0.75 m", "temp_c.1.25"="1.25 m", "temp_c.1.5"= "1.5 m", "temp_c.2"="2 m")))

m.df.sim.interp$doy<-yday(m.df.sim.interp$datetime)
m.df.sim.interp$doy_frac<-hour(m.df.sim.interp$datetime)
m.df.sim.interp$min<-minute(m.df.sim.interp$datetime)
m.df.sim.interp$minfrac[m.df.sim.interp$min=="30"] <- 0.5
m.df.sim.interp$minfrac[m.df.sim.interp$min=="0"] <- 0
m.df.sim.interp$hourfrac<-(m.df.sim.interp$doy_frac + m.df.sim.interp$minfrac)/24
m.df.sim.interp$doy_frac<-m.df.sim.interp$doy+m.df.sim.interp$hourfrac

m.df.sim.interp.fig<-m.df.sim.interp %>% 
  select(datetime,doy,doy_frac,value, variable,group) %>%
  rename(temp_c=value) %>% 
  mutate(variable=(recode(variable, "temp_c.0" = "0 m", "temp_c.0.25"="0.25 m", "temp_c.0.5"="0.5 m",
                "temp_c.0.75"="0.75 m", "temp_c.1.25"="1.25 m", "temp_c.1.5"= "1.5 m", "temp_c.2"="2 m")))
write.csv(m.df.sim.interp.fig, "model_out_calibrated.csv",row.names = F)

write.csv(m.obs.fig, "model_out_observed.csv",row.names = F)


# determination of change in temp over heatwave
m.df.sim_daily<-m.df.sim.interp.fig %>% 
  filter(variable=="2 m") %>% 
  group_by(doy) %>% 
  summarize(across(c(temp_c), ~mean(.,na.rm=T))) %>% 
  ungroup()

rmse.fig<-rmse %>% 
  mutate(variable=(recode(variable, "temp_c.0" = "0 m", "temp_c.0.25"="0.25 m", "temp_c.0.5"="0.5 m",
                          "temp_c.0.75"="0.75 m", "temp_c.1.25"="1.25 m", "temp_c.1.5"= "1.5 m", "temp_c.2"="2 m")))

# PLOT
modelgoodnessplot<-
  ggplot() +
  geom_line(data = m.obs.fig,aes(doy_frac, temp_c), col="grey40") +
  geom_line(data = m.df.sim.interp.fig, aes(doy_frac, temp_c), col="#FF0066") +
  geom_text(data=rmse.fig, aes( x=150, y=35, label=round(fit,2)),                 
            color="black", size =3.5) +
  facet_wrap(~ variable, nrow=4) +
  xlab('Day of Year') + ylab('Water Temperature (C)')+
  theme_linedraw()  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("Model_Fit.png",modelgoodnessplot, width = 7, height = 11, units = 'in')

### -----------------------------------------------------------------------------------------------------------------




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

#### OPTIMIZATION ROUTINE --------------------------------------------------------------------------------------

optim_macro <- function(p, scaling = TRUE, lb, ub){

  if(!is.numeric(p)){
    p = values.optim
  }
  if (scaling == TRUE){
    p <- wrapper_scales(p, lb, ub)
  }

  temp <- c()
  

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
                          Cd = p[1], # wind momentum drag
                          km = p[2], #0.008,#0.04, # macrophyte light extinction
                          Cdplant = p[3],#1e3, # macrophyte momentum drag
                          ahat = p[4], # macrophyte area to volume
                          reflect = 0.6, 
                          infra = 0.4,
                          windfactor = p[5], # wind multiplier 
                          shortwavefactor = p[6], # shortwave radiation multiplier
                          diffusionfactor = p[7], # diffusion multiplier
                          diffmethod = 2,
                          densThresh = 1e-3,
                          Hgeo = p[8],
                          rho_plant = p[9]
  )
  temp <-cbind(temp, res$temp)
  
  time =  meteo_all[[1]]$Date[1] + seq(1, ncol(temp))*dt
  df <- data.frame(cbind(time, t(temp)) )
  colnames(df) <- c("time", as.character(paste0(seq(1,nrow(temp)))))
  m.df <- reshape2::melt(df, "time")
  m.df$time <- time
  
  
  ## model goodness
  obs <- read.csv('bc/temp_profiles_2020.csv')
  obs <- obs %>%
    filter(pond == 'B', site_id == '20') %>%
    select(datetime, temp_depth_m, temp_c)
  
  df.sim <- df
  colnames(df.sim) <- c("datetime", as.character(paste0('temp_c.',seq(1,nrow(temp))*dx)))
  df.sim$datetime <-  time
  
  idx <- (match(as.POSIXct(obs$datetime), as.POSIXct(df.sim$datetime) ))
  
  
  obs <- obs[which(!is.na(idx)), ]
  
  deps <- seq(1,nrow(temp))*dx
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
  return(sum(rmse$fit, na.rm = T))
}

#can play around with sum versus median


scaling = TRUE
calib_setup <- data.frame('x0' = c(1.3e-3, 0.005, 1, 0.5, 1, 1, 1, 0.1, 30),
                          'ub' = c(4.3e-3, 0.9, 1e3, 0.9, 1.3, 1.3, 1.2, 0.2, 3e3),
                          'lb' = c(1.3e-4, 1e-2, 1e-2, 1e-2, 0.7, 0.4, 0.8, 1e-2, 3e-2))

if (scaling){
  init.val <- (calib_setup$x0 - calib_setup$lb) *10 /(calib_setup$ub-calib_setup$lb) 
}

opt <- pureCMAES(par = init.val, fun = optim_macro, lower = rep(0,length(init.val)), 
          upper = rep(10,length(init.val)), 
          sigma = 0.5, 
          stopfitness = 1, 
          stopeval = 30,
          scaling = scaling, lb = calib_setup$lb, ub = calib_setup$ub)
  
val <- wrapper_scales(opt$xmin, lb = calib_setup$lb, ub = calib_setup$ub)
write.csv(val, 'calib_results.txt', quote = F, row.names = F)

### OPTIMISED RUN
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
                        Cd = val[1], # wind momentum drag
                        km = val[2], #0.008,#0.04, # macrophyte light extinction
                        Cdplant = val[3],#1e3, # macrophyte momentum drag
                        ahat = val[4], # macrophyte area to volume
                        reflect = 0.6, 
                        infra = 0.4,
                        windfactor = val[5], # wind multiplier 
                        shortwavefactor = val[6], # shortwave radiation multiplier
                        diffusionfactor = val[7], # diffusion multiplier
                        diffmethod = 2,
                        densThresh = 1e-2,
                        Hgeo = val[8],
                        rho_plant = val[9]
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
  filter(pond == 'B', site_id == '20') %>%
  select(datetime, temp_depth_m, temp_c)

df.sim <- df
colnames(df.sim) <- c("datetime", as.character(paste0('temp_c.',seq(1,nrow(temp))*dx)))
df.sim$datetime <-  time

idx <- (match(as.POSIXct(obs$datetime), as.POSIXct(df.sim$datetime) ))


obs <- obs[which(!is.na(idx)), ]

deps <- seq(1,nrow(temp))*dx #deps <-  hyps_all[[2]]
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
                                 'fit' = sqrt((sum((o$value-s$value)**2, na.rm = TRUE))/nrow(o))))
}

ggplot() +
  geom_line(data = m.obs,aes(datetime, value, col = group)) +
  geom_line(data = m.df.sim.interp, aes(datetime, value, col = group)) +
  geom_text(data=rmse, aes( as.POSIXct("2020-05-26 10:30:00 CDT"), y=17, label=round(fit,2)),                 
            color="black", size =3) +
  facet_wrap(~ variable, scales = 'free') +
  xlab('') + ylab('Temp. (deg C)')+
  theme_bw()
ggsave(filename = paste0('fieldcomparison.png'), width = 15, height = 8, units = 'in')

