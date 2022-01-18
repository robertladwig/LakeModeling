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

## colors for plotting
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(rLakeAnalyzer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals)))

## lake configurations
zmax = 2#1.7
nx = 10 # number of layers we will have
dt = 3600 # time step, 30 min times 60 seconds/min
dx = zmax/nx # spatial step

# nt = 10  * 24 * 60 * 60

# area and depth values of our lake 
area = seq(-350,-1e-1, length.out = nx) * (-1)
depth = seq(0,nx, length.out = nx)

## function to calculate density from temperature
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

rho = calc_dens(u)

## this is our attempt for turbulence closure, estimating eddy diffusivity
eddy_diffusivity <-function(rho, depth, g, rho_0){
  buoy = rep(1, (nx)) * 7e-5
  for (i in seq(1, nx-1)){#range(0, nx - 1):
    buoy[i] = sqrt( abs(rho[i+1] - rho[i]) / (depth[i+1] - depth[i]) * g/rho_0 )
  }
  buoy[nx] = sqrt( abs(rho[nx-1] - rho[nx]) / abs(depth[nx-1] - depth[nx]) * 
                     g/rho_0 )
  
  low_values_flags = buoy < 7e-5  # Where values are low
  buoy[low_values_flags] = 7e-5
  
  kz = 0.00706 *( 3.8 * 1e1)**(0.56) * (buoy)**(-0.43)
  return(kz)
}

kz = eddy_diffusivity(rho, depth, 9.81, 998.2) / 86400# 1e4

## atmospheric boundary conditions
## create daily meteorological variables
air.df <- read.csv('bc/Ames_weather_station.csv')
airT = data.frame('time' = as.POSIXct( air.df$valid[which(!is.na(as.double(air.df$tmpc)))], format = '%m/%d/%Y %H:%M'),
                  'atr' = as.double(na.omit(as.double(air.df$tmpc))) ,
                  'dtr' = as.double(na.omit(as.double(air.df$tmpc))) - ((100 - as.double(air.df$relh[which(!is.na(as.double(air.df$tmpc)))]))/5.),
                  'relhum' = as.double(na.omit(as.double(air.df$relh))))
airT$ea <- (101.325 * exp(13.3185 * (1 - (373.15 / (airT$atr + 273.15))) -
                            1.976 * (1 - (373.15 / (airT$atr  + 273.15)))**2 -
                            0.6445 * (1 - (373.15 / (airT$atr  + 273.15)))**3 -
                            0.1229 * (1 - (373.15 / (airT$atr  + 273.15)))**4)) *airT$relhum /100
airT$dt <- as.POSIXct(airT$time) - (as.POSIXct(airT$time)[1]) + 1

met.df <- read.csv('bc/weather_station_ponds.csv')
wQ = data.frame('time' = as.POSIXct(met.df$date_time, format = '%m/%d/%Y %H:%M'),
                'Jsw' = met.df$par/2 , #/4.6, PAR is 50% of total short-wave
                'Uw' = met.df$wind_speed)

minDates <- which(!is.na(match(wQ$time, airT$time)))
maxDates <- match(wQ$time, airT$time)[which(!is.na(match(wQ$time, airT$time)))]

startDate <- airT$time[maxDates[1]]

airT = airT %>%
  filter(time >= airT$time[maxDates[1]] & time <= airT$time[maxDates[length(maxDates)]])
wQ = wQ %>%
  filter(time >= wQ$time[minDates[1]] & time <= wQ$time[minDates[length(minDates)]])

airT$Cloud_Cover <- gotmtools::calc_cc(date = as.POSIXct(airT$time),
                     airt = airT$atr,
                     relh = airT$relhum,
                     swr = approx(wQ$time, wQ$Jsw, airT$time)$y,
                     lat = 43, lon = -89.41,
                     elev =  mean(air.df$elevation))
airT$lw <- gotmtools::calc_in_lwr(cc = airT$Cloud_Cover,
                                  airt = airT$atr,
                                  relh = airT$relhum)
airT$Pa = 100000

airT$dt <- as.POSIXct(airT$time) - (as.POSIXct(airT$time)[1]) + 1
wQ$dt <- as.POSIXct(wQ$time) - (as.POSIXct(wQ$time)[1]) + 1


## linearization of driver data, so model can have dynamic step
Jsw <- approxfun(x = wQ$dt, y = wQ$Jsw, method = "linear", rule = 2)
Jlw <- approxfun(x = airT$dt, y = airT$lw, method = "linear", rule = 2)
Tair <- approxfun(x = airT$dt, y = airT$atr, method = "linear", rule = 2)
ea <- approxfun(x = airT$dt, y = airT$ea, method = "linear", rule = 2)
Uw <- approxfun(x = wQ$dt, y = wQ$Uw, method = "linear", rule = 2)
CC <- approxfun(x = airT$dt, y = airT$Cloud_Cover, method = "linear", rule = 2)
Pa <- approxfun(x = airT$dt, y = airT$Pa, method = "linear", rule = 2)

## additional parameters to run the model
# meteorology
Rl <- 0.03 # 0.3
Acoeff <- 0.6 # coefficient between 0.5 - 0.7
sigma <-  4.9 * 10^(-3) / (24 * 3600) # 11.7 * 10^(-8) # cal / (cm2 d K4) or: 
# 4.9 * 10^(-3) # Stefan-Boltzmann constant in [J (m2 d K4)-1]
eps <- 0.97 # emissivity of water
cp <- 4184 # specific heat (J/kg/C)
c1 <- 0.47 # Bowen's coefficient
a <- 7 # constant
c <- 9e4 # empirical constant
g <- 9.81  # gravity (m/s2)

# vertical heating
reflect <- 0.6 # fraction of reflected solar radiation
infra = 0.3 # fraction infrared radiation
kd = 1# 0.2# 1.0 #0.2 # light attenuation coefficient

emissivity = 0.97
sigma = 5.67 * 10^(-8)
p2 = 1
B = 0.61

# macrophyte
macro <- read_csv('bc/ModelTest_MacrophyteSummary.csv')
macro$dt <- ((macro$doy) - ((macro$doy)[1]) + 1)* 86400
macro$dt[1] = 1
# macroheight <- mean(macro$stand_height_m)
# macrobiomass <- mean(macro$biomass_gDWperM3, na.rm = T)

# macro$stand_height_m <- mean(macro$stand_height_m)
# macro$biomass_gDWperM3 <- mean(macro$biomass_gDWperM3, na.rm = T)

macroheight <- approxfun(x = macro$dt, y = macro$stand_height_m, method = "linear", rule = 2)
macrobiomss <- approxfun(x = macro$dt, y = macro$biomass_gDWperM3, method = "linear", rule = 2)

# vertical heating
reflect <- 0.6 # fraction of reflected solar radiation
infra = 0.3 # fraction infrared radiation
kd = 0.1# 1.0 #0.2 # light attenuation coefficient
km = 0.06#0.06 #0.4 # specific light attenuation coefficient for macrophytes
# P = macrobiomass # macrophyte biomass per unit volume in gDW m-3, e.g. 100

# dissipative turbulent energy by macrophytes
Cd = 1.0 #1.0 # plant form drag coefficient
ahat = 0.04#0.4 # 0.4 plant surface area per unit volume
# Hmacrophytes <- seq(dx, zmax, length.out = nx) 
# Hmacrophytes <- ifelse(Hmacrophytes < ((max(Hmacrophytes) - macroheight)), 0, 1) #c(rep(0,2),rep(1,8)) # height of macrophytes (abundance)
rho_mp = 1.0 # 70 #70 # biomass density

# longwave <- function(cc, sigma, Tair, ea, emissivity, Jlw){  # longwave radiation into
#   lw = emissivity * Jlw 
#   return(lw)
# }
longwave <- function(cc, sigma, Tair, ea, emissivity, Jlw){  # longwave radiation into
  Tair = Tair + 273.15
  Ea <- 1.24 * (1 + 0.17 * cc**2) * (ea/Tair)^(1/7)
  lw <- emissivity * Ea *sigma * Tair**4
  return(lw)
}
backscattering <- function(emissivity, sigma, Twater){ # backscattering longwave 
  # radiation from the lake
  Twater = Twater + 273.15
  back = (eps * sigma * (Twater )^4) 
  return((-1) * back)
}
sensible <- function(p2, B, Tair, Twater, Uw){ # convection / sensible heat
  Twater = Twater + 273.15
  Tair = Tair + 273.15
  fu = 4.4 + 1.82 * Uw + 0.26 *(Twater - Tair)
  sensible <- ( p2 * B * fu * (Twater - Tair)) 
  return((-1) * sensible)
}
latent <- function(Tair, Twater, Uw, p2, pa, ea){ # evaporation / latent heat 
  Twater = Twater + 273.15
  Tair = Tair + 273.15
  Pressure = pa / 100
  fu = 4.4 + 1.82 * Uw + 0.26 *(Twater - Tair)
  fw = 0.61 * (1 + 10^(-6) * Pressure * (4.5 + 6 * 10^(-5) * Twater**2))
  ew = fw * 10 * ((0.7859+0.03477* Twater)/(1+0.00412* Twater))
  latent = fu * p2 * (ew - ea) 
  return((-1) * latent)
}

## plot initial profile
plot( u, seq(0, nx * dx, length.out=(nx)),  
      ylim = rev(range(seq(0, dx * nx, length.out=(nx)))), xlim = c(0,35), 
      type ='l');

ice = FALSE
Hi= 0
nt = as.double(max(airT$dt)) # maximum simulation length

um <- matrix(NA, ncol = floor(nt/dt), nrow = nx)
kzm <- matrix(NA, ncol = floor(nt/dt), nrow = nx)
n2m <- matrix(NA, ncol = floor(nt/dt), nrow = nx)
Hts <- rep(NA, length = floor(nt/dt))
Swf <- rep(NA, length = floor(nt/dt))
Lwf <- rep(NA, length = floor(nt/dt))
BLwf <- rep(NA, length = floor(nt/dt))
Lf <- rep(NA, length = floor(nt/dt))
Sf <- rep(NA, length = floor(nt/dt))
mix <- rep(NA, length = floor(nt/dt))
therm.z <- rep(NA, length = floor(nt/dt))
mix.z <- rep(NA, length = floor(nt/dt))
Him <- rep(NA, length = floor(nt/dt))
macroz <- rep(NA, length = floor(nt/dt))

## modeling code for vertical 1D mixing and heat transport
for (n in 1:floor(nt/dt)){  #iterate through time
  print(paste0(round(n*100/floor(nt/dt),3),' %'))
  un = u # prior temperature values
  kz = eddy_diffusivity(calc_dens(un), depth, 9.81, 998.2) / 86400
  
  if (ice){
    kzn = rep(1E-6, length = length(kz))
    absorp = 0.85
  } else {
    kzn = kz   
    absorp = 1-reflect# 0.3
  }
  kzm[, n] <- kzn
  
  # surface heat flux
  Q <- (absorp * Jsw(n * dt) + longwave(cc = CC(n * dt), sigma = sigma, Tair = Tair(n * dt), ea = ea(n * dt), emissivity = emissivity, Jlw = Jlw(n * dt)) + # longwave(emissivity = emissivity, Jlw = Jlw(n * dt)) +
          backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) +
          latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt)) + 
          sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n * dt))) 
  
  # heat addition over depth
  P = macrobiomss(n * dt)
  H = (1- reflect) * (1- infra) * (Jsw(n * dt))  * #
    exp(-(kd + km * P) *seq(dx,nx*dx,length.out=nx)) 

  ## (1) DIFFUSION
  # surface layer
  u[1] = un[1] +
     Q * area[1]/(area[1]*dx)*1/(4184 * calc_dens(un[1]) ) *dt +
    H[1] * area[1]/(area[1]*dx) * 1/(4184 * calc_dens(un[1]) )* dt
  
 Hts[n] <-  Q *  area[1]/(area[1]*dx)*1/(4181 * calc_dens(un[1]) ) # append(Hts, Q *  area[1]/(area[1]*dx)*1/(4181 * calc_dens(un[1]) ))
 Swf[n] <-  absorp * Jsw(n * dt) # append(Swf, 0.3 * Jsw(n * dt))
 Lwf[n] <-  longwave(cc = CC(n * dt), sigma = sigma, Tair = Tair(n * dt), ea = ea(n * dt), emissivity = emissivity, Jlw = Jlw(n * dt))# longwave(emissivity = emissivity, Jlw = Jlw(n * dt))  #append(Lwf, longwave(emissivity = emissivity, Jlw = Jlw(n * dt)) )
 BLwf[n] <-  backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) #append(BLwf, backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1])) 
 Lf[n] <- latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt)) #append(Lf, latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt)) )
 Sf[n] <- sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n * dt)) #append(Sf, sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n * dt)))
 # 
  # all other layers in between
  for (i in 2:(nx-1)){
    u[i] = un[i] + 
      kzn[i] * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1]) +
      H[i] * area[i]/(area[i]*dx) * 1/(4184 * calc_dens(un[i]) )* dt
  }
  
  # bottom layer
  u[nx] = un[nx] + 
    H[nx] * area[nx]/(area[nx]*dx) * 1/(4181 * calc_dens(un[nx]) ) * dt
  
  ## (2) TURBULENT MIXING OF MIXED LAYER
  # the mixed layer depth is determined for each time step by comparing kinetic 
  # energy available from wind and the potential energy required to completely 
  # mix the water column to a given depth
  Zcv <- seq(1, nx) %*% area / sum(area) # center of volume
  tau = 1.225 * 0.0013 * Uw(n * dt)^2 # wind shear is air density times shear 
  if (Uw(n * dt) <= 15) {
    c10 = 0.0005 * sqrt(Uw(n * dt))
  } else {
    c10 = 0.0026
  }
  shear = sqrt((c10 * calc_dens(un[1]))/1.225) *  Uw(n * dt)# shear velocity
  # coefficient times wind velocity squared
  KE = shear *  tau * dt # kinetic energy as function of wind
  
  Hmacrophytes <- seq(dx, zmax, length.out = nx) 
  Hmacrophytes <- ifelse(Hmacrophytes < ((max(Hmacrophytes) - macroheight(n*dt) )), 0, 1) #c(rep(0,2),rep(1,8)) # height of macrophytes (abundance)
  macroz[n] = macroheight(n *dt)
  
  maxdep = 1
  for (dep in 1:(nx-1)){
    if (dep == 1){
      PE = abs(g *  ( seq(1,nx)[dep] - Zcv)  * calc_dens(u[dep]) * dx)
      # PE = abs(g *   seq(1,nx)[dep] *( seq(1,nx)[dep+1] - Zcv)  * 
                 # abs(calc_dens(u[dep+1])- calc_dens(u[dep])))
      DKE = Hmacrophytes[dep]*(rho_mp* ahat * Cd) *Uw(n * dt)^3 * dt  *dx
    } else {
      PEprior = PE
      PE = abs(g *  ( seq(1,nx)[dep] - Zcv)  * calc_dens(u[dep]) * dx +
                 PEprior)
      # PE = abs(g *   seq(1,nx)[dep] *( seq(1,nx)[dep+1] - Zcv)  * 
                 # abs(calc_dens(u[dep+1])- calc_dens(u[dep]))) + PEprior
      DKEprior = DKE
      DKE = Hmacrophytes[dep]*(rho_mp * ahat * Cd) *Uw(n * dt)^3 * dt  *dx + DKEprior
      KE = KE - DKE
    }
      if (PE > KE){
        maxdep = dep
        break
      }
    maxdep = dep
  }
  u[1:maxdep] = mean(u[1:maxdep])
  mix[n] <- KE/PE #append(mix, KE/PE)
  therm.z[n] <- maxdep #append(therm.z, maxdep)
  
  ## (3) DENSITY INSTABILITIES
  # convective overturn: Convective mixing is induced by an unstable density 
  # profile. All groups of water layers where the vertical density profile is 
  # unstable are mixed with the first stable layer below the unstable layer(s) 
  # (i.e., a layer volume weighed means of temperature and other variables are 
  # calculated for the mixed water column). This procedure is continued until 
  # the vertical density profile in the whole water column becomes neutral or stable.
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
  
  dens_u_n2 = calc_dens(u) 
  n2 <- 9.81/mean(calc_dens(u)) * (lead(dens_u_n2) - lag(dens_u_n2))/dx
  max.n2 <- ifelse(max(n2, na.rm = T) > 1E-4, which.max(n2) * dx, dx * nx)
  mix.z[n] <- max.n2
  
  n2m[, n] <- n2
  um[, n] <- u
  
  lines( u, seq(0, dx * nx, length.out=(nx)),
          ylim = rev(range(seq(0, dx * nx, length.out=(nx)))), lty = 'dashed');
  
  ## (4) ICE FORMATION
  # according to Hostetler & Bartlein (1990): 
  # (1) ice forms when surface water temp <= 1 deg C and melts when > 1 deg
  # (2) rate of ice formation/melting is exponential function of ice thickness
  # (the thicker the ice, the slower the formation rate, and vice versa)
  # (3) heat of fusion is added/subtracted from surface energy balance
  # (4) diffusion below ice only happens on molecular level
  # (5) with ice, surface absorption of incoming solar radiation increases to 85 %
  if (any(u <= 0) == TRUE){
    supercooled <- which(u < 0)
    initEnergy <- sum((0-u[supercooled])*hyps$Area_meterSquared[supercooled] * dx * 4.18E6)
    
    if (ice != TRUE) {
      Hi <- 0.01+(initEnergy/(910*333500))/max(hyps$Area_meterSquared)
    } else {
      if (Tair(n*dt) > 0){
        Tice <- 0
        Hi = Hi -max(c(0, dt*(((1-0.4)*Jsw(n * dt))+(backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) +
                                                      latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt)) + 
                                                      sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n * dt))) )/(1000*333500)))
      } else {
        Tice <-  ((1/(10 * Hi)) * 0 +  Tair(n*dt)) / (1 + (1/(10 * Hi))) 
        Hi <- sqrt(Hi**2 + 2 * 2.1/(910 * 333500)* (0 - Tice) * dt)
      }
    }
    ice = TRUE
    if (Hi > 0){
      u[supercooled] = 0
      u[1] = 0
    }
    Him[n] <- Hi
  } else if (ice == TRUE & Hi > 0) {
        if (Tair(n*dt) > 0){
          Tice <- 0
          Hi = Hi -max(c(0, dt*(((1-0.2)*Jsw(n * dt))+(longwave(cc = CC(n * dt), sigma = sigma, Tair = Tair(n * dt), ea = ea(n * dt), emissivity = emissivity, Jlw = Jlw(n * dt)) + backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) +
                                                    latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt)) + 
                                                    sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n * dt))) )/(1000*333500)))
        } else {
          Tice <-  ((1/(10 * Hi)) * 0 +  Tair(n*dt)) / (1 + (1/(10 * Hi))) 
          Hi <- sqrt(Hi**2 + 2 * 2.1/(910 * 333500)* (0 - Tice) * dt)
        }
        u[supercooled] = 0
        u[1] = 0
        Him[n] <- Hi
  } else if (ice == TRUE & Hi <= 0){
    ice = FALSE 
  }
}

str(um)
## water temperature time series at different depths
plot(seq(1, ncol(um))*dt/24/3600, um[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(17,40), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(um))*dt/24/3600, um[i,], col = sample(col_vector,1), 
        lty = 'dashed',lwd =2)
}

# schmidt stability (processed model output)
time =  startDate + seq(1, ncol(um))*dt#/24/3600
df <- data.frame(cbind(time, t(um)) )
colnames(df) <- c("datetime", as.character(paste0('wtr_',seq(1,nrow(um))*dx)))
df$datetime <- time
df.h <- data.frame('depths' = depth, 'areas' = area)
SI <- rLakeAnalyzer::ts.schmidt.stability(wtr = df, bathy = df.h)

ggplot(SI, aes(datetime, schmidt.stability)) +
  geom_line() +
  ylab('Schmidt Stability (J/m2)') + xlab('')+
  theme_minimal()

# lake number (processed model output)
df.u <- data.frame('datetime' = df$datetime,
                   'wind' = approx(wQ$time, wQ$Uw, df$datetime)$y)
LN <- rLakeAnalyzer::ts.lake.number(wtr = df, wnd = df.u,
                                    wnd.height = 10,bathy = df.h)

ggplot(LN, aes(datetime, (lake.number))) +
  geom_line() +
  scale_y_continuous(trans='log10') +
  ylab('log10 Lake Number ()') + xlab('')+
  theme_minimal()

g.stab <- ggplot(LN, aes(datetime, (lake.number))) +
  geom_line() +
  ylab('log10 Lake Number ()') + xlab('')+
  geom_line(data = SI, aes(y = schmidt.stability*1000 ),
            linetype = 2, col = 'blue') +
  scale_y_continuous(trans='log10',
    sec.axis = sec_axis(trans = ~ . /1000,
                        name = "Schmidt Stability (J/m2)")) +
  theme_minimal()

# ice thickness (direct model output)
plot(seq(1, ncol(um))*dt/24/3600, Him, type = 'l', 
     col = 'red', xlab = 'Time', 
     ylab= 'Ice thickness (m)', lwd= 3)

# macrophyte height (direct model output)
plot(seq(1, ncol(um))*dt/24/3600, macroz, type = 'l', 
     col = 'red', xlab = 'Time', 
     ylab= 'Macrophyte height (m)', lwd= 3)
g.macro <- ggplot(data.frame('time' = time,
                             'macro' =macroz))+
  geom_line(aes(time, macro)) +
  xlab('') + ylab('Macrophyte height (m)')+
  theme_minimal(); g.macro
 
## surface mixed layer depth (direct model output)
therm.z.roll = zoo::rollmean(therm.z, 14)
plot(seq(1, ncol(um))*dt/24/3600, therm.z, type = 'l',ylim = rev(range( seq(0,zmax, length.out=nx))), 
     col = 'red', xlab = 'Time', 
        ylab= 'Mixed Layer Depth (m)', lwd= 3)
lines(therm.z.roll, lty ='dashed', lwd =3)

# decision if lake is stratified or not: 1 deg C criterium
strat.state <- um[1,] - um[nx,]
strat.state <- ifelse(strat.state > 1, 1, 0)
plot(seq(1, ncol(um))*dt/24/3600, strat.state, col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Stratified conditions', ylim=c(0,1), lwd = 2)

## Max. buoyancy frequency layer depth (direct model output)
mix.z.roll = zoo::rollmean(mix.z, 14)
plot(mix.z, type = 'l',ylim = rev(range( seq(0,zmax, length.out=nx))), 
     col = 'red', xlab = 'Time', 
     ylab= 'Max. Buoyancy Freq. Depth (m)', lwd= 3)
lines(mix.z.roll, lty ='dashed', lwd =3)


# meteorological heat fluxes
fluxes.df <- data.frame('time' = seq(1, ncol(um))*dt/24/3600,
                        'shortwave' = Swf,
                        'longwavein' = Lwf,
                        'longwaveout' = BLwf,
                        'latent' = Lf,
                        'sensible' = Sf,
                        'total' = Swf + Lwf + BLwf + Lf + Sf)
m.fluxes.df <- reshape2::melt(fluxes.df, id = 'time')
ggplot(m.fluxes.df) +
  geom_line(aes(time, value)) +
  facet_wrap(~ variable, scales = 'free')

## vertical water temperature profiles over time
df = data.frame('1' =NULL)
name = NULL
for (i in seq(1,ncol(um), length.out = 20)){
  i = floor(i)
  name = append(name, round((i*dt)/24/3600,1))
  df = append(df, data.frame('1' = um[,i] + round((i*dt)/24/3600,1)))
}
df.m = matrix(unlist(df), ncol = 20)
colnames(df.m) = name

df.m.m = reshape2::melt(df.m)
ggplot2::ggplot(df.m.m) +
  geom_point(aes(x = value/40, y=Var1,  col = value-Var2, group = Var2)) +
  geom_line(aes(x = value/40, y=Var1,  col = value-Var2, group = Var2)) +
  scale_y_reverse() + xlab('Time') + ylab('Depth')+labs(col='Temp')+
  scale_color_gradient(low = "lightblue", high = "red") +
  theme_minimal()

## vertical temperature profiles
# for (i in seq(1,ncol(um), length.out = 200)){
#   n = i
#   i = floor(i)
#   png(paste0('../../animation_macrophyte/pic_',match(n, seq(1,ncol(um), 
# length.out=200)),'.png'))
#   plot(um[,i], seq(dx,zmax, length.out=nx),
#        ylim = rev(range( seq(dx,zmax, length.out=nx))), col = 'red', 
# type = 'l', xlab = 'Temperature (degC)', 
#        ylab='Depth (m)',
#        xlim = c(0,35), main = paste0('time (d): ',round((i*dt)/24/3600,1)),
#        lwd = 3)
#   dev.off()
# }
# 

## contour plot of water temperature
time =  seq(1, ncol(um))*dt/24/3600
df <- data.frame(cbind(time, t(um)) )
colnames(df) <- c("time", as.character(paste0(seq(1,nrow(um))*dx)))
m.df <- reshape2::melt(df, "time")
m.df$time <- time

df.kz <- data.frame(cbind(time, t(kzm)) )
colnames(df.kz) <- c("time", as.character(paste0(seq(1,nrow(kzm))*dx)))
m.df.kz <- reshape2::melt(df.kz, "time")
m.df.kz$time <- time

df.n2 <- data.frame(cbind(time, t(n2m)) )
colnames(df.n2) <- c("time", as.character(paste0(seq(1,nrow(n2m))*dx)))
m.df.n2 <- reshape2::melt(df.n2, "time")
m.df.n2$time <- time

g1 <- ggplot(m.df, aes(as.numeric(time), as.numeric(as.character(variable)))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(15,35),
                         colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Temp [degC]')+
  scale_y_reverse() 
g2 <- ggplot(m.df.kz, aes(as.numeric(time), as.numeric(as.character(variable)))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(#limits = c(-20,35),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Diffusion [m2/s]')+
  scale_y_reverse() 
g3 <- ggplot(m.df.n2, aes(as.numeric(time), as.numeric(as.character(variable)))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(#limits = c(-20,35),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'N2 [s-2]')+
  scale_y_reverse() 
g <- g1 / g2 / g3 / g.stab / g.macro; g
ggsave(filename = 'heatmap.png',plot = g, width = 15, height = 8, units = 'in')


## model goodness
obs <- read.csv('bc/temp_profiles_2020.csv')
obs <- obs %>%
  filter(pond == 'F') %>%
  select(datetime, temp_depth_m, temp_c)

df.sim <- df
colnames(df.sim) <- c("datetime", as.character(paste0('temp_c.',seq(1,nrow(um))*dx)))
df.sim$datetime <-   startDate + seq(1, ncol(um))*dt#/24/3600

# idx <- na.omit(match(as.POSIXct(df.sim$datetime), as.POSIXct(obs$datetime) ))
idx <- (match(as.POSIXct(obs$datetime), as.POSIXct(df.sim$datetime) ))

# df.sim <- df.sim[idx, ]
obs <- obs[which(!is.na(idx)), ]

deps <- seq(1,nrow(um))*dx
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
# df.sim.interp <- apply(df.sim[,-1], 1, function(x) approx(deps, x, unique(obs$temp_depth_m))$y)
df.sim.interp$datetime <-   startDate + seq(1, ncol(um))*dt#/24/3600
colnames(df.sim.interp) <- c(as.character(paste0('temp_c.',unique(obs$temp_depth_m))), 'datetime')

wide.obs <- reshape(obs, idvar = "datetime", timevar = "temp_depth_m", direction = "wide")
m.obs <- reshape2::melt(wide.obs, id = 'datetime')
m.obs$datetime <- as.POSIXct(m.obs$datetime)
m.obs$group <- 'obs'
m.df.sim.interp <- reshape2::melt(df.sim.interp, id = 'datetime')
m.df.sim.interp$group <- 'sim'

# df.rmse <- merge(m.df.sim.interp, m.obs, by = c('datetime', 'variable'))
# df.rmse %>%
#   group_by(variable, datetime) %>%
#   mutate(rmse = sqrt((sum(value.x-value.y)**2)/length(value.y))) %>%
#   group_by(variable) %>%
#   summarise(sum(rmse))

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
  # annotate(data=rmse, geom="text", x=as.POSIXct("2020-05-25 10:30:00 CDT"), y=16, label=(rmse$fit),
           # color="red") +
  geom_text(data=rmse, aes( as.POSIXct("2020-05-26 10:30:00 CDT"), y=17, label=round(fit,2)),                 
            color="black", size =3) +
  facet_wrap(~ variable) +
  xlab('') + ylab('Temp. (deg C)')+
  theme_bw()
ggsave(filename = 'fieldcomp.png', width = 15, height = 8, units = 'in')







               