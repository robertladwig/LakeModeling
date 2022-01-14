#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#' Created on Thu Aug 19 11:29:34 2021
#' 
#' @author: robert
#' @email: rladwig2@wisc.edu
#' 
#' Temperature equation as Tt = 1/A K Tzz
#' Diffusion code is based on 12 steps to Navier-Stokes by (c) Lorena A. Barba, Gilbert F. Forsyth 2017.
#' https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/ 
#' 
#' Eddy diffusivity is estimated from buoyancy frequency according to Hondzo and Stefan (1993)
#' Lake Water Temperature Simulation Model. ASCE
#' https://doi.org/10.1061/(ASCE)0733-9429(1993)119:11(1251) 
#' 
#' Mixing dynamics code is taken from Herb & Stefan (2004) Temperature stratification and Mixing 
#' Dynamics in a Shallow Lake with Submersed Macrophytes. Lake & Reservoir Management
#' https://www.tandfonline.com/doi/pdf/10.1080/07438140409354159
#' 
#' Convective overturn algorithm is taken from Salorante & Andersen (2007) MyLake—A multi-year 
#' lake simulation model code suitable for uncertainty and sensitivity analysis simulations. 
#' Ecological Modeling
#' https://doi.org/10.1016/j.ecolmodel.2007.03.018 

# remove everything from workspace
rm(list = ls())

# set wd to current dir of script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(tidyverse)
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

zmax = 1.7
nx = 10 # number of layers we will have
dt = 30 * 60 # time step, 30 min times 60 seconds/min
dx = zmax/nx # spatial step

# area and depth values of our lake 
area = seq(-350,-1e-1, length.out = nx) * (-1)
depth = seq(0,nx, length.out = nx)

# function to calculate density from temperature
calc_dens <-function(wtemp){
  dens = 999.842594 + (6.793952 * 1e-2 * wtemp) - (9.095290 * 1e-3 *wtemp**2) + 
    (1.001685 * 1e-4 * wtemp**3) - (1.120083 * 1e-6* wtemp**4) + 
    (6.536336 * 1e-9 * wtemp**5)
  return(dens)
}

# initial water temperature profile
init.df <- read.csv('bc/initialprofile.txt') %>%
  arrange(depth_m)
u = approx(init.df$depth_m, init.df$temp_c,
           seq(0, nx * dx, length.out= nx))$y
# u[(0):(10)] = 7

rho = calc_dens(u)

# this is our attempt for turbulence closure, estimating eddy diffusivity
eddy_diffusivity <-function(rho, depth, g, rho_0){
  buoy = rep(1, (nx)) * 7e-5
  for (i in seq(1, nx-1)){#range(0, nx - 1):
    buoy[i] = sqrt( abs(rho[i+1] - rho[i]) / (depth[i+1] - depth[i]) * g/rho_0 )
  }
  buoy[nx] = sqrt( abs(rho[nx-1] - rho[nx]) / abs(depth[nx-1] - depth[nx]) * g/rho_0 )
  
  low_values_flags = buoy < 7e-5  # Where values are low
  buoy[low_values_flags] = 7e-5
  
  kz = 0.00706 *( 3.8 * 1e1)**(0.56) * (buoy)**(-0.43)
  return(kz)
}

kz = eddy_diffusivity(rho, depth, 9.81, 998.2) / 86400# 1e4

# atmospheric boundary conditions
air.df <- read.csv('bc/Ames_weather_station.csv')
airT = data.frame('time' = as.POSIXct( air.df$valid[which(!is.na(as.double(air.df$tmpc)))], format = '%m/%d/%Y %H:%M'),
  'atr' = as.double(na.omit(as.double(air.df$tmpc))) ,
  'dtr' = as.double(na.omit(as.double(air.df$tmpc))) - ((100 - as.double(air.df$relh[which(!is.na(as.double(air.df$tmpc)))]))/5.))
met.df <- read.csv('bc/weather_station_ponds.csv')
wQ = data.frame('time' = as.POSIXct(met.df$date_time, format = '%m/%d/%Y %H:%M'),
                'Jsw' = met.df$par/4.6 , #/4.6, PAR is 50% of total short-wave
                'vW' = met.df$wind_speed)
wQ$Uw <- 1225*0.0013* wQ$vW^2 # 19.0 + 0.95 * (wQ$vW)^2 # wind shear stress
wQ$Uw <- 19.0 + 0.95 * (wQ$vW)^2 # wind shear stress

minDates <- which(!is.na(match(wQ$time, airT$time)))
maxDates <- match(wQ$time, airT$time)[which(!is.na(match(wQ$time, airT$time)))]

airT = airT %>%
  filter(time >= airT$time[maxDates[1]] & time <= airT$time[maxDates[length(maxDates)]])
wQ = wQ %>%
  filter(time >= wQ$time[minDates[1]] & time <= wQ$time[minDates[length(minDates)]])

airT$dt = airT$time - (airT$time[1]) + 1
wQ$dt = wQ$time - (wQ$time[1]) +1

nt = as.double(max(airT$dt)) # maximum simulation length

# linearization of driver data, so model can have dynamic step
Jsw <- approxfun(x = wQ$dt, y = wQ$Jsw, method = "linear", rule = 2)
Tair <- approxfun(x = airT$dt, y = airT$atr, method = "linear", rule = 2)
Dew <- approxfun(x = airT$dt, y = airT$dtr, method = "linear", rule = 2)
Uw <- approxfun(x = wQ$dt, y = wQ$Uw, method = "linear", rule = 2)
vW <- approxfun(x = wQ$dt, y = wQ$vW, method = "linear", rule = 2)

## additional parameters to run the model
# meteorology
Rl <- 0.03 # 0.3
Acoeff <- 0.6 # coefficient between 0.5 - 0.7
sigma <-  4.9 * 10^(-3) / (24 * 3600) # 11.7 * 10^(-8) # cal / (cm2 d K4) or: 4.9 * 10^(-3) # Stefan-Boltzmann constant in [J (m2 d K4)-1]
eps <- 0.97 # emissivity of water
cp <- 4184 # specific heat (J/kg/C)
c1 <- 0.47 # Bowen's coefficient
a <- 7 # constant
c <- 9e4 # empirical constant
g <- 9.81  # gravity (m/s2)


macro <- read_csv('bc/ModelTest_MacrophyteSummary.csv')
macroheight <- mean(macro$stand_height_m)
macrobiomass <- mean(macro$biomass_gDWperM3, na.rm = T)
# vertical heating
reflect <- 0.6 # fraction of reflected solar radiation
infra = 0.3 # fraction infrared radiation
kd = 0.1# 1.0 #0.2 # light attenuation coefficient
km = 0.06 #0.4 # specific light attenuation coefficient for macrophytes
P = macrobiomass # macrophyte biomass per unit volume in gDW m-3, e.g. 100

# dissipative turbulent energy by macrophytes
Cd = 1.0 # plant form drag coefficient
ahat = 0.4 # plant surface area per unit volume
Hmacrophytes <- seq(dx, zmax, length.out = nx) 
Hmacrophytes <- ifelse(Hmacrophytes < ((max(Hmacrophytes) - macroheight)), 0, 1) #c(rep(0,2),rep(1,8)) # height of macrophytes (abundance)
rho_mp = 70 # biomass density

longwave <- function(sigma, Tair, Acoeff, eair, Rl){  # longwave radiation into the lake
  lw = (sigma * (Tair(n * dt) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) 
  return(lw)
}
backscattering <- function(eps, sigma, Twater){ # backscattering longwave radiation from the lake
  back = (eps * sigma * (Twater + 273)^4) 
  return((-1) * back)
}
latent <- function(c1, wind, Twater, Tair){ # convection / latent heat
  latent <- (c1 * wind * (Twater - Tair))
  return((-1) * latent)
}
sensible <- function(wind, esat, eair){ # evaporation / sensible heat 
  sensible = (wind * ((esat) - (eair))) 
  return((-1) * sensible)
}

# plot initial profile
plot( u, seq(0, nx * dx, length.out=(nx)),  
      ylim = rev(range(seq(0, dx * nx, length.out=(nx)))), xlim = c(0,35), type ='l');

um <- c()
Hts <- c()
# modeling code for vertical 1D mixing and heat transport
for (n in 1:floor(nt/dt)){  #iterate through time
  un = u # prior temperature values
  kz = eddy_diffusivity(calc_dens(un), depth, 9.81, 998.2) / 86400
  kzn = kz   
  
  eair <- (4.596 * exp((17.27 * Dew(n * dt)) / (237.3 + Dew(n * dt)))) # air vapor pressure
  esat <- 4.596 * exp((17.27 * Tair(n * dt)) / (237.3 + Tair(n * dt))) # saturation vapor pressure
  RH <- eair/esat *100 # relative humidity
  es <- 4.596 * exp((17.27 * u[1])/ (273.3+u[1]))
  Q <- (Jsw(n * dt) + 
          longwave(sigma = sigma, Tair = Tair(n * dt), Acoeff = Acoeff, eair = eair, Rl = Rl) +
          backscattering(eps = eps, sigma = sigma, Twater = un[1]) +
          latent(c1 = c1, wind = Uw(n * dt), Twater = un[1], Tair = Tair(n *dt)) + 
          sensible(wind = Uw(n * dt), esat = esat, eair = eair))
  
  # heat addition over depth
  H = (1- reflect) * (1- infra) * (Jsw(n * dt)+
                                     (sigma * (Tair(n * dt) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)))  * 
    exp(-(kd + km * P) *seq(dx,nx*dx,length.out=nx)) 
  
  ## (1) DIFFUSION
  # surface layer
  u[1] = un[1] + 
     Q * area[1]/(area[1]*dx)*1/(4181 * calc_dens(un[1]) ) * dt#+
  
 Hts <- append(Hts, Q *  area[1]/(area[1]*dx)*1/(4181 * calc_dens(un[1]) ))
  # all other layers in between
  for (i in 2:(nx-1)){
    u[i] = un[i] + 
      kzn[i] * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1]) +
      H[i] * area[i]/(area[i]*dx) * 1/(4181 * calc_dens(un[i]) ) *dt
  }
  
  # bottom layer
  u[nx] = un[nx] +
    H[nx] * area[nx]/(area[nx]*dx) * 1/(4181 * calc_dens(un[nx]) )* dt
  
  ## (2) TURBULENT MIXING OF MIXED LAYER
  # the mixed layer depth is determined for each time step by comparing kinetic energy available
  # from wind and the potential energy required to completely mix the water column to a given depth
  Zcv <- seq(1, nx) %*% area / sum(area) # center of volume
  tau = 1225 * 0.0013 * vW(n * dt)^2 # wind shear is air density times shear coefficient times wind velocity squared
  KE = vW(n * dt) *  tau * dt # kinetic energy as function of wind
  maxdep = 1
  for (dep in 1:(nx)){
    if (dep == 1){
      DKE = Hmacrophytes[dep]*(rho_mp* ahat * Cd) *Uw(n * dt)^3 * dt  *dx
      PE = abs(g *  ( seq(1,nx)[dep] - Zcv)  * calc_dens(u[dep]) * dx)
      KE = KE - DKE
    } else {
      PEprior = PE
      DKEprior = DKE
      DKE = Hmacrophytes[dep]*(rho_mp * ahat * Cd) *Uw(n * dt)^3 * dt  *dx + DKEprior
      PE = abs(g *  ( seq(1,nx)[dep] - Zcv)  * calc_dens(u[dep]) * dx +
                 PEprior) 
      KE = KE - DKE
    }
      if (PE > KE){
        maxdep = dep
        break
      }
    maxdep = dep
  }
  u[1:maxdep] = max(u[1:maxdep])
  
  ## (3) DENSITY INSTABILITIES
  # convective overturn: Convective mixing is induced by an unstable density profile. 
  # All groups of water layers where the vertical density profile is unstable are mixed with the 
  # first stable layer below the unstable layer(s) (i.e., a layer volume weighed means of 
  # temperature and other variables are calculated for the mixed water column). 
  # This procedure is continued until the vertical density profile in the whole water column becomes neutral or stable.
  dens_u = calc_dens(u) 
  diff_dens_u <- (diff(dens_u)) 
  diff_dens_u[abs(diff(dens_u)) < 1e-4] = 0
  while (any(diff_dens_u < 0)){
    dens_u = calc_dens(u) 
    for (dep in 1:(nx-1)){
      if (dens_u[dep+1] < dens_u[dep] & abs(dens_u[dep+1] - dens_u[dep]) > 1e-4){
        u[dep:(dep+1)] = mean(u[dep:(dep+1)])
        break
      }
    }
    dens_u = calc_dens(u) 
    diff_dens_u <- (diff(dens_u)) 
    diff_dens_u[abs(diff(dens_u)) < 1e-4] = 0
  }

  um <- cbind(um, u)
  
  lines( u, seq(0, dx * nx, length.out=(nx)),
          ylim = rev(range(seq(0, dx * nx, length.out=(nx)))), lty = 'dashed');
}

str(um)
plot(seq(1, ncol(um))*dt/24/3600, um[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(20,40), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(um))*dt/24/3600, um[i,], col = sample(col_vector,1), lty = 'dashed', lwd =2)
}

df = data.frame('1' =NULL)
name = NULL
for (i in seq(1,ncol(um), length.out = 20)){
  i = floor(i)
  name = append(name, round((i*dt)/24/3600,1))
  df = append(df, data.frame('1' = um[,i] + round((i*dt)/24/3600,1)))
}
df.m = matrix(unlist(df), ncol = 20)
colnames(df.m) = name
depth = seq(dx, zmax, length.out= nx)

df.m.m = reshape2::melt(df.m)
ggplot2::ggplot(df.m.m) +
  geom_point(aes(x = value/40, y=Var1*dx,  col = value-Var2, group = Var2)) +
  geom_line(aes(x = value/40, y=Var1*dx,  col = value-Var2, group = Var2)) +
  scale_y_reverse() + xlab('Time') + ylab('Depth')+labs(col='Temp')+
  scale_color_gradient(low = "lightblue", high = "red") +
  theme_minimal()
