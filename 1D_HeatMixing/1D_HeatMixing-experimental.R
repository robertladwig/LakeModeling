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
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals)))

## lake configurations
zmax = 25 # maximum lake depth
nx = 25 # number of layers we will have
dt = 3600 # 24 hours times 60 min/hour times 60 seconds/min
dx = zmax/nx # spatial step

nyear = 3
nt = nyear * 365* 24 * 60 * 60 # as.double(max(wQ$dt)) # maximum simulation length

## area and depth values of our lake 
hyps <- read_csv('bc/LakeEnsemblR_bathymetry_standard.csv')
area = approx(hyps$Depth_meter,hyps$Area_meterSquared,seq(1,nx*dx, 
                                                          length.out= nx))$y
area[which.min(area)] <- 1e-2
depth = depth= seq(1,nx*dx, length.out = nx)
volume <- c(rev(diff(pracma::cumtrapz(area, depth))*(-1)),0)
volume[which(volume == 0)] = min(volume[-which(volume == 0)])
volume <- rep(0, (length(depth)-1))
for (p in 1:length(volume)){
  volume[p] <- pracma::trapz(depth[p:(p+1)],area[p:(p+1)])
}
volume <- c(volume, 1000)

## function to calculate density from temperature
calc_dens <-function(wtemp){
  dens = 999.842594 + (6.793952 * 1e-2 * wtemp) - (9.095290 * 1e-3 *wtemp**2) + 
    (1.001685 * 1e-4 * wtemp**3) - (1.120083 * 1e-6* wtemp**4) + 
    (6.536336 * 1e-9 * wtemp**5)
  return(dens)
}

## here we define our initial profile
obs <- read_csv('bc/obs.txt')
init.df <- obs %>% 
  filter(datetime == min(datetime)) %>%
  arrange(Depth_meter)
if (max(depth) > max(init.df$Depth_meter)){
  init.df <- rbind(init.df, init.df[nrow(init.df),])
  init.df$Depth_meter[nrow(init.df)] <- max(depth)
}
u = approx(init.df$Depth_meter, init.df$Water_Temperature_celsius,
           seq(0, nx * dx, length.out= nx))$y

rho = calc_dens(u)

## this is our attempt for turbulence closure, estimating eddy diffusivity
eddy_diffusivity <-function(rho, depth, g, rho_0, ice){
  buoy = rep(1, (nx)) * 7e-5
  for (i in seq(1, nx-1)){#range(0, nx - 1):
    buoy[i] = sqrt( abs(rho[i+1] - rho[i]) / (depth[i+1] - depth[i]) * g/rho_0 )
  }
  buoy[nx] = sqrt( abs(rho[nx-1] - rho[nx]) / abs(depth[nx-1] - depth[nx]) * 
                     g/rho_0 )
  
  low_values_flags = buoy < 7e-5  # Where values are low
  buoy[low_values_flags] = 7e-5
  
  if (ice){
    ak <- 0.000898
  } else{
    ak <- 0.00706 *( max(area)/1E6)**(0.56)
  }
  
  kz = ak * (buoy**2)**(-0.43)
  return(kz)
}

kz = eddy_diffusivity(rho, depth, g = 9.81, rho_0 = 009.2, ice = FALSE) / 86400# 1e4

## atmospheric boundary conditions
## create daily meteorological variables
meteo <- read_csv('bc/LakeEnsemblR_meteo_standard.csv')
# daily_meteo <-  meteo %>%
#   mutate(date = as.Date(datetime)) %>%
#   group_by(date) %>%
#   summarise_all(mean)
daily_meteo <- meteo
daily_meteo$date = daily_meteo$datetime
daily_meteo$Cloud_Cover <- gotmtools::calc_cc(date = as.POSIXct(daily_meteo$date),
                     airt = daily_meteo$Air_Temperature_celsius,
                     relh = daily_meteo$Relative_Humidity_percent,
                     swr = daily_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared,
                     lat = 43, lon = -89.41,
                     elev = 258)
daily_meteo$dt <- as.POSIXct(daily_meteo$date) - (as.POSIXct(daily_meteo$date)[1]) + 1
daily_meteo$ea <- (daily_meteo$Relative_Humidity_percent * (4.596 * exp((17.27*(daily_meteo$Air_Temperature_celsius))/
                                                                          (237.3 + (daily_meteo$Air_Temperature_celsius) )))/100)
daily_meteo$ea <- (101.325 * exp(13.3185 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15))) -
                1.976 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15)))**2 -
                0.6445 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15)))**3 -
                0.1229 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15)))**4)) *daily_meteo$Relative_Humidity_percent/100
daily_meteo$ea <- (daily_meteo$Relative_Humidity_percent/100) * 10^(9.28603523 - 2322.37885/(daily_meteo$Air_Temperature_celsius + 273.15))
startDate <- daily_meteo$datetime[1]

## calibration parameters
daily_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared <-
  daily_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared 
daily_meteo$Ten_Meter_Elevation_Wind_Speed_meterPerSecond <-
  daily_meteo$Ten_Meter_Elevation_Wind_Speed_meterPerSecond * 0.8# wind speed multiplier
Cd <- 0.0009 #0.0006 # wind shear drag coefficient, usually set at 0.0013 because 'appropriate for most engineering solutions' (Fischer 1979)
meltP <- 5 # melt energy multiplier
dt_iceon_avg =  0.8 # moving average modifier for ice onset
Hgeo <- 0.1 # 0.1 W/m2 geothermal heat flux
KEice <- 1/1000
Ice_min <- 0.1
# kd = 0.4 # 0.2# 1.0 #0.2 # light attenuation coefficient

## light
# Package ID: knb-lter-ntl.31.30 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Secchi Disk Depth; Other Auxiliary Base Crew Sample Data 1981 - current.
secview <- read_csv('bc/light.csv') %>%
  filter(sampledate >= startDate)
if (secview$sampledate[1] >= startDate){
  secview <- rbind(data.frame('sampledate' = startDate,
                              'secnview' = secview$secnview[1]),
                   secview)
}
secview$dt <- as.POSIXct(secview$sampledate) - (as.POSIXct(secview$sampledate)[1]) + 1
secview$kd <- 1.7 / secview$secnview
secview$kd  <- zoo::na.approx(secview$kd)

## linearization of driver data, so model can have dynamic step
Jsw <- approxfun(x = daily_meteo$dt, y = daily_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared, method = "linear", rule = 2)
Jlw <- approxfun(x = daily_meteo$dt, y = daily_meteo$Longwave_Radiation_Downwelling_wattPerMeterSquared, method = "linear", rule = 2)
Tair <- approxfun(x = daily_meteo$dt, y = daily_meteo$Air_Temperature_celsius, method = "linear", rule = 2)
ea <- approxfun(x = daily_meteo$dt, y = daily_meteo$ea, method = "linear", rule = 2)
Uw <- approxfun(x = daily_meteo$dt, y = daily_meteo$Ten_Meter_Elevation_Wind_Speed_meterPerSecond, method = "linear", rule = 2)
CC <- approxfun(x = daily_meteo$dt, y = daily_meteo$Cloud_Cover, method = "linear", rule = 2)
Pa <- approxfun(x = daily_meteo$dt, y = daily_meteo$Surface_Level_Barometric_Pressure_pascal, method = "linear", rule = 2)
kd <- approxfun(x = secview$dt, y = secview$kd, method = "constant", rule = 2)
RH <- approxfun(x = daily_meteo$dt, y = daily_meteo$Relative_Humidity_percent, method = "constant", rule = 2)

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
reflect <- 0.3 # fraction of reflected solar radiation
infra = 0.7 # fraction infrared radiation

emissivity = 0.97
sigma = 5.67 * 10^(-8)
p2 = 1
B = 0.61

# longwave <- function(cc, sigma, Tair, ea, emissivity, Jlw){  # longwave radiation into
#   lw = emissivity * Jlw
#   return(lw)
# }
longwave <- function(cc, sigma, Tair, ea, emissivity, Jlw){  # longwave radiation into
  Tair = Tair + 273.15
  p <- (1.33 * ea/Tair)
  Ea <- 1.24 * (1 + 0.17 * cc**2) * p**(1/7)
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
latent <- function(Tair, Twater, Uw, p2, pa, ea, RH){ # evaporation / latent heat
  Twater = Twater + 273.15
  Tair = Tair + 273.15
  Pressure = pa / 100
  fu = 4.4 + 1.82 * Uw + 0.26 *(Twater - Tair)
  fw = 0.61 * (1 + 10^(-6) * Pressure * (4.5 + 6 * 10^(-5) * Twater**2))
  ew = fw * 10 * ((0.7859+0.03477* Twater)/(1+0.00412* Twater))
  latent = fu * p2 * (ew - ea)# * 1.33) #* 1/6
  return((-1) * latent)
}
# latent <- function(Tair, Twater, Uw, p2, pa, ea, RH){ # evaporation / latent heat
#   Twater = Twater
#   Tair = Tair
#   Pressure = pa / 100
#   Lv = 2.501 * 10^6 - 2370 * Twater
#   es = 6.11 * exp((17.27 * Tair)/(237.3 + Tair))
#   ez = (RH * es)/100
#   qz = (0.622 * ez) / Pressure
#   esat = 6.11 * exp((17.27 * Twater)/(237.3 + Twater))
#   q0 = (0.622 * esat)/ Pressure
#   Ra = 2.87 * (1 + 0.608 * qz)
#   rho_z = (100 * Pressure)/(Ra * (Tair + 273.16))
#   Ce = 0.0013
#   latent = rho_z * Lv * Ce * Uw * (q0-qz)
#   return((-1) * latent)
# }

## plot initial profile
plot( u, seq(0, nx * dx, length.out=(nx)),  
      ylim = rev(range(seq(0, dx * nx, length.out=(nx)))), xlim = c(0,35), 
      type ='l');

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


densThresh <- 1e-3
ice = FALSE
Hi= 0
iceT <- 6
scheme = 'explicit' # options are 'explicit' (FTCS, Forward Time Centered Space) or 'implicit' (Crank-Nicholson scheme)

start.time <- Sys.time()
## modeling code for vertical 1D mixing and heat transport
for (n in 1:floor(nt/dt)){  #iterate through time
  print(paste0(round(n*100/floor(nt/dt),3),' %'))
  un = u # prior temperature values
  kz = eddy_diffusivity(calc_dens(un), depth, 9.81, 998.2, ice) / 86400
  
  if (ice & Tair(n*dt) <= 0){
    kzn = kz
    absorp = 1 - 0.7
    infra = 1 - absorp
  } else if (ice & Tair(n*dt) >= 0){
    kzn = kz
    absorp = 1 - 0.3
    infra = 1 - absorp
    } else if (!ice) {
    kzn = kz   
    absorp = 1 - reflect# 0.3
    infra = 1 - absorp
  }
  kzm[, n] <- kzn
  
  ## (1) Heat addition
  # surface heat flux
  Q <- (absorp * Jsw(n * dt) + longwave(cc = CC(n * dt), sigma = sigma, Tair = Tair(n * dt), ea = ea(n * dt), emissivity = emissivity, Jlw = Jlw(n * dt)) + #longwave(emissivity = emissivity, Jlw = Jlw(n * dt)) +
          backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) +
          latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt), RH = RH(n * dt)) + 
          sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n * dt)))  
  
  # heat addition over depth
  # H = (1- reflect) * (1- infra) * (Jsw(n * dt))  * #
  #   exp(-(kd ) *seq(dx,nx*dx,length.out=nx)) 
  H =  (1- infra) * (Jsw(n * dt))  * #
    exp(-(kd(n * dt) ) *seq(dx,nx*dx,length.out=nx)) 
  
  Hg <- (area-lead(area))/dx * Hgeo/(4181 * calc_dens(un[1])) 
  Hg[which(is.na(Hg))] <- min(Hg, na.rm = TRUE)

  # add heat to all layers
  if (scheme == 'implicit'){
    ## (2) DIFFUSION
    # surface layer
    un[1] = un[1] +    Q * area[1]/(dx)*1/(4184 * calc_dens(un[1]) ) * dt/area[1]
    # # bottom layer
    un <- un  + (H * area/(dx) * 1/(4184 * calc_dens(un) ))* dt/area
    
    j <- length(kzn)
    y <- array(0, c(j,j))
    
    # all other layers in between
    # Linearized heat conservation equation matrix (diffusion only)
    az <- (dt/dx**2) * kzn                                         #coefficient for i-1
    cz <- (dt/dx**2) * kzn                #coefficient for i+1
    bz <- 1 + 2 * (dt/dx**2) * kzn                                                         #coefficient for i+1
    #Boundary conditions, surface
    az[1] <- 0
    #cz(1) remains unchanged
    bz[1]<- 1 #+ az[1] + (dt/dx**2) * kzn    
    #Boundary conditions, bottom
    #az(end) remains unchanged
    cz[length(cz)] <- 0
    bz[length(bz)] <- 1 #+ (dt/dx**2) * kzn    + cz[length(cz)]
    y[0 + 1:(j - 1) * (j + 1)] <- -cz[-length(bz)]	# superdiagonal
    y[1 + 0:(j - 1) * (j + 1)] <- bz	# diagonal
    y[2 + 0:(j - 2) * (j + 1)] <- -az[-1] 	# subdiagonal
    
    y[1,2] <- 0#- 2 * (dt/dx**2) * kzn[1]           
    y[nrow(y), (ncol(y)-1)] = 0#-2 * (dt/dx**2) * kzn[ncol(y)]           
    
    u <- solve(y, un)
  }

  
  # save variables for model diagnostics
  Hts[n] <-  Q *  area[1]/(area[1]*dx)*1/(4181 * calc_dens(un[1]) ) # append(Hts, Q *  area[1]/(area[1]*dx)*1/(4181 * calc_dens(un[1]) ))
  Swf[n] <-  absorp * Jsw(n * dt) # append(Swf, 0.3 * Jsw(n * dt))
  Lwf[n] <-  longwave(cc = CC(n * dt), sigma = sigma, Tair = Tair(n * dt), ea = ea(n * dt), emissivity = emissivity, Jlw = Jlw(n * dt))#longwave(emissivity = emissivity, Jlw = Jlw(n * dt))  #append(Lwf, longwave(emissivity = emissivity, Jlw = Jlw(n * dt)) )
  BLwf[n] <-  backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) #append(BLwf, backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1])) 
  Lf[n] <- latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt), RH = RH(n * dt)) #append(Lf, latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt)) )
  Sf[n] <- sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n * dt)) #append(Sf, sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n * dt)))
  
  
  ## (2) DIFFUSION
  # surface layer
  if (scheme == 'explicit'){
    u[1] = un[1] +
      (Q * area[1]/(dx)*1/(4184 * calc_dens(un[1]) ) +
         abs(H[1+1]-H[1]) * area[1]/(dx) * 1/(4184 * calc_dens(un[1]) ) +
         Hg[1]) * dt/area[1]
    
    # all other layers in between
    for (i in 2:(nx-1)){
      u[i] = un[i] +
        (area[i] * kzn[i] * 1 / dx**2 * (un[i+1] - 2 * un[i] + un[i-1]) +
           abs(H[i+1]-H[i]) * area[i]/(dx) * 1/(4184 * calc_dens(un[i]) ) +
           Hg[i])* dt/area[i]
    }
    # bottom layer
    u[nx] = un[nx] +
      abs(H[nx]-H[nx-1]) * area[nx]/(area[nx]*dx) * 1/(4181 * calc_dens(un[nx]) +
                                                         Hg[nx]/area[nx]) * dt
  }

  ## (3) TURBULENT MIXING OF MIXED LAYER
  # the mixed layer depth is determined for each time step by comparing kinetic 
  # energy available from wind and the potential energy required to completely 
  # mix the water column to a given depth
  Zcv <- seq(1, nx) %*% area / sum(area) # center of volume
  tau = 1.225 * Cd * Uw(n * dt)^2 # wind shear is air density times wind velocity 
  if (Uw(n * dt) <= 15) {
    c10 = 0.0005 * sqrt(Uw(n * dt))
  } else {
    c10 = 0.0026
  }
  shear = sqrt((c10 * calc_dens(un[1]))/1.225) *  Uw(n * dt) # shear velocity
  # coefficient times wind velocity squared
  KE = shear *  tau * dt # kinetic energy as function of wind
  
  if (ice){
    KE = KE * KEice
  }
  maxdep = 1
  for (dep in 1:(nx-1)){
    if (dep == 1){
      # PE = abs(g *  ( seq(1,nx)[dep] - Zcv)  * calc_dens(u[dep]) * dx)
      PE = abs(g *   seq(1,nx)[dep] *( seq(1,nx)[dep+1] - Zcv)  *
                 # abs(calc_dens(u[dep+1])- calc_dens(u[dep])))
                 abs(calc_dens(u[dep+1])- mean(calc_dens(u[1:dep]))))
    } else {
      PEprior = PE
      # PE = abs(g *  ( seq(1,nx)[dep] - Zcv)  * calc_dens(u[dep]) * dx +
      #            PEprior)
      PE = abs(g *   seq(1,nx)[dep] *( seq(1,nx)[dep+1] - Zcv)  *
      # abs(calc_dens(u[dep+1])- calc_dens(u[dep]))) + PEprior
        abs(calc_dens(u[dep+1])- mean(calc_dens(u[1:dep])))) + PEprior
      
    }
      if (PE > KE){
        maxdep = dep-1
        break
      } else if (dep>1 & PE < KE ){
        u[(dep-1):dep] = (u[(dep-1):dep] %*% volume[(dep-1):dep])/sum(volume[(dep-1):dep])
      }
    maxdep = dep
  }
  # u[1:maxdep] = (u[1:(maxdep)] %*% volume[1:(maxdep)])/sum(volume[1:(maxdep)]) #mean(u[1:maxdep])
  mix[n] <- KE/PE #append(mix, KE/PE)
  therm.z[n] <- maxdep #append(therm.z, maxdep)
  
  ## (4) DENSITY INSTABILITIES
  # convective overturn: Convective mixing is induced by an unstable density 
  # profile. All groups of water layers where the vertical density profile is 
  # unstable are mixed with the first stable layer below the unstable layer(s) 
  # (i.e., a layer volume weighed means of temperature and other variables are 
  # calculated for the mixed water column). This procedure is continued until 
  # the vertical density profile in the whole water column becomes neutral or stable.
  dens_u = calc_dens(u) 
  diff_dens_u <- (diff(dens_u)) 
  diff_dens_u[abs(diff(dens_u)) <= densThresh] = 0
  while (any(diff_dens_u < 0)){
    dens_u = calc_dens(u) 
    for (dep in 1:(nx-1)){
      if (dens_u[dep+1] < dens_u[dep] & abs(dens_u[dep+1] - dens_u[dep]) >= densThresh){
        u[dep:(dep+1)] = (u[dep:(dep+1)] %*% volume[dep:(dep+1)])/sum(volume[dep:(dep+1)]) #mean(u[dep:(dep+1)])
        break
      }
    }
    dens_u = calc_dens(u) 
    diff_dens_u <- (diff(dens_u)) 
    diff_dens_u[abs(diff(dens_u)) <= densThresh] = 0
  }
  
  dens_u_n2 = calc_dens(u) 
  n2 <- 9.81/mean(calc_dens(u)) * (lead(dens_u_n2) - lag(dens_u_n2))/dx
  max.n2 <- ifelse(max(n2, na.rm = T) > 1E-4, which.max(n2) * dx, dx * nx)
  mix.z[n] <- max.n2
  

  
  ## (5) ICE FORMATION
  # according to Hostetler & Bartlein (1990): 
  # (1) ice forms when surface water temp <= 1 deg C and melts when > 1 deg
  # (2) rate of ice formation/melting is exponential function of ice thickness
  # (the thicker the ice, the slower the formation rate, and vice versa)
  # (3) heat of fusion is added/subtracted from surface energy balance
  # (4) diffusion below ice only happens on molecular level
  # (5) with ice, surface absorption of incoming solar radiation increases to 85 %
  icep  = max(dt_iceon_avg,  (dt/86400))
  x = (dt/86400) / icep
  iceT = iceT * (1 - x) + u[1] * x
  if ((iceT <= 0) == TRUE){
  # if (any(u <= 0) == TRUE){
    supercooled <- which(u < 0)
    initEnergy <- sum((0-u[supercooled])*hyps$Area_meterSquared[supercooled] * dx * 4.18E6)
    
    if (ice != TRUE) {
      Hi <- Ice_min+(initEnergy/(910*333500))/max(hyps$Area_meterSquared)
    } else {
      if (Tair(n*dt) > 0){
        Tice <- 0
        Hi = Hi -max(c(0, meltP * dt*((absorp*Jsw(n * dt))+(longwave(cc = CC(n * dt), sigma = sigma, Tair = Tair(n * dt), ea = ea(n * dt), emissivity = emissivity, Jlw = Jlw(n * dt)) +
                                                      backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) +
                                                      latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt),  RH = RH(n * dt)) + 
                                                      sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n * dt))) )/(1000*333500)))
      } else {
        Tice <-  ((1/(10 * Hi)) * 0 +  Tair(n*dt)) / (1 + (1/(10 * Hi))) 
        Hi <- min(Ice_min, sqrt(Hi**2 + 2 * 2.1/(910 * 333500)* (0 - Tice) * dt))
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
          Hi = Hi -max(c(0, meltP * dt*((absorp*Jsw(n * dt))+(backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) +
                                                    latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt),  RH = RH(n * dt)) + 
                                                    sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n * dt))) )/(1000*333500))) 
        } else {
          Tice <-  ((1/(10 * Hi)) * 0 +  Tair(n*dt)) / (1 + (1/(10 * Hi))) 
          Hi <- min(Ice_min, sqrt(Hi**2 + 2 * 2.1/(910 * 333500)* (0 - Tice) * dt))
        }
        u[supercooled] = 0
        u[1] = 0
        Him[n] <- Hi
  } else if (ice == TRUE & Hi <= 0){
    ice = FALSE 
  }
  
  n2m[, n] <- n2
  um[, n] <- u
  
  lines( u, seq(0, dx * nx, length.out=(nx)),
         ylim = rev(range(seq(0, dx * nx, length.out=(nx)))), lty = 'dashed');
  
  # print(un)
  # print(u)
  # cat ("Press [enter] to continue")
  # line <- readline()
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

str(um)
## water temperature time series at different depths
plot(seq(1, ncol(um))*dt/24/3600, um[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(-1,35), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(um))*dt/24/3600, um[i,], col = sample(col_vector,1), 
        lty = 'dashed',lwd =2)
}

# ice thickness (direct model output)
plot(seq(1, ncol(um))*dt/24/3600, Him, type = 'l', 
     col = 'red', xlab = 'Time', 
     ylab= 'Ice thickness (m)', lwd= 3)
 
## surface mixed layer depth (direct model output)
therm.z.roll = zoo::rollmean(therm.z, 14)
plot(seq(1, ncol(um))*dt/24/3600, therm.z, type = 'l',ylim = rev(range( seq(0,zmax, length.out=nx))), 
     col = 'red', xlab = 'Time', 
        ylab= 'Mixed Layer Depth (m)', lwd= 3)
lines(therm.z.roll, lty ='dashed', lwd =3)

## decision if lake is stratified or not: 1 deg C criterium
strat.state <- um[1,] - um[nx,]
strat.state <- ifelse(strat.state > 1, 1, 0)
plot(seq(1, ncol(um))*dt/24/3600, strat.state, col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Stratified conditions', ylim=c(0,1), lwd = 2)

## max. buoyancy frequency layer depth (direct model output)
mix.z.roll = zoo::rollmean(mix.z, 14)
plot(mix.z, type = 'l',ylim = rev(range( seq(0,zmax, length.out=nx))), 
     col = 'red', xlab = 'Time', 
     ylab= 'Max. Buoyancy Freq. Depth (m)', lwd= 3)
lines(mix.z.roll, lty ='dashed', lwd =3)

## check stratification durations
df.ice <- data.frame('time' = startDate + seq(1, ncol(um))*dt,
                     'ice' = Him,
                     'stratified' = strat.state)
g.ice <- ggplot(df.ice) +
  geom_line(aes(time, ice)) +
  ylab('ice thickness (m)') + xlab('')+
  geom_line(aes(time, stratified),
            linetype = 2, col = 'blue') +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . & 1,
                                         name = "stratified (yes/no)")) +
  theme_minimal();g.ice


# meteorological heat fluxes
fluxes.df <- data.frame('time' =  startDate + seq(1, ncol(um))*dt,#/24/3600seq(1, ncol(um))*dt/24/3600,
                        'shortwave' = Swf,
                        'longwave' = Lwf + BLwf,
                        'latent' = Lf,
                        'sensible' = Sf,
                        'total' = Swf + Lwf + BLwf + Lf + Sf)
m.fluxes.df <- reshape2::melt(fluxes.df, id = 'time')
m.fluxes.df$type = 'model'
ggplot(m.fluxes.df) +
  geom_line(aes(time, value)) +
  facet_wrap(~ variable, scales = 'free')

glm.meteo <- read.csv('bc/lake.csv') %>%
  rename(shortwave = Daily.Qsw, latent = Daily.Qe,
         sensible = Daily.Qh, longwave = Daily.Qlw) %>%
  mutate(total = shortwave + longwave + sensible + latent) %>%
  select(time, shortwave, latent, sensible, longwave, total)
m.glm.meteo <- reshape2::melt(glm.meteo, id = 'time')
m.glm.meteo$type = 'GLM'
m.glm.meteo$time <- as.POSIXct(m.glm.meteo$time)

time =  startDate + seq(1, ncol(um))*dt
g.heat <- ggplot(m.fluxes.df, aes(time, value)) +
  geom_line(aes(time, value, col = type)) +
  geom_line(aes(y=zoo::rollmean(value, 24, na.pad=TRUE))) +
  geom_line(data =  m.glm.meteo, aes(time, value, col = type)) +
  xlim(min(time), max(time))+
  facet_wrap(~ variable, scales = 'free'); g.heat
ggsave(filename = paste0('boundaryconditions','_',scheme,'_',nyear,'.png'),g.heat,  width = 15, height = 8, units = 'in')


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



## contour plot of water temperature
# time =  seq(1, ncol(um))*dt/24/3600
time =  startDate + seq(1, ncol(um))*dt
df <- data.frame(cbind(time, t(um)) )
colnames(df) <- c("time", as.character(paste0(seq(1,nrow(um)))))
m.df <- reshape2::melt(df, "time")
m.df$time <- as.POSIXct(time)

df.kz <- data.frame(cbind(time, t(kzm)) )
colnames(df.kz) <- c("time", as.character(paste0(seq(1,nrow(kzm)))))
m.df.kz <- reshape2::melt(df.kz, "time")
m.df.kz$time <- as.POSIXct(time)

df.n2 <- data.frame(cbind(time, t(n2m)) )
colnames(df.n2) <- c("time", as.character(paste0(seq(1,nrow(n2m)))))
m.df.n2 <- reshape2::melt(df.n2, "time")
m.df.n2$time <- as.POSIXct(time)

g1 <- ggplot(m.df, aes((time), as.numeric(variable))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(-2,40),
                         colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Temp [degC]')+
  scale_y_reverse() 
g2 <- ggplot(m.df.kz, aes((time), as.numeric(variable))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(#limits = c(-20,35),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Diffusion [m2/s]')+
  scale_y_reverse() 
g3 <- ggplot(m.df.n2, aes((time), as.numeric(variable))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(#limits = c(-20,35),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'N2 [s-2]')+
  scale_y_reverse() 
g <- g1 / g2 / g3 / g.ice; g
ggsave(filename = paste0('heatmaps','_',scheme,'_',nyear,'.png'),plot = g, width = 15, height = 8, units = 'in')



# Package ID: knb-lter-ntl.130.29 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: High Frequency Water Temperature Data - Lake  Mendota Buoy 2006 - current.
inUrl2  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/130/29/63d0587cf326e83f57b054bf2ad0f7fe" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")

dt2 <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "sampledate",     
                 "year4",     
                 "month",     
                 "daynum",     
                 "hour",     
                 "depth",     
                 "wtemp",     
                 "flag_wtemp"    ), check.names=TRUE)

unlink(infile2)

# attempting to convert dt2$sampledate dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp2sampledate<-as.Date(dt2$sampledate,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp2sampledate) == length(tmp2sampledate[!is.na(tmp2sampledate)])){dt2$sampledate <- tmp2sampledate } else {print("Date conversion failed for dt2$sampledate. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp2sampledate) 
if (class(dt2$year4)=="factor") dt2$year4 <-as.numeric(levels(dt2$year4))[as.integer(dt2$year4) ]               
if (class(dt2$year4)=="character") dt2$year4 <-as.numeric(dt2$year4)
if (class(dt2$month)=="factor") dt2$month <-as.numeric(levels(dt2$month))[as.integer(dt2$month) ]               
if (class(dt2$month)=="character") dt2$month <-as.numeric(dt2$month)
if (class(dt2$daynum)=="factor") dt2$daynum <-as.numeric(levels(dt2$daynum))[as.integer(dt2$daynum) ]               
if (class(dt2$daynum)=="character") dt2$daynum <-as.numeric(dt2$daynum)
if (class(dt2$depth)=="factor") dt2$depth <-as.numeric(levels(dt2$depth))[as.integer(dt2$depth) ]               
if (class(dt2$depth)=="character") dt2$depth <-as.numeric(dt2$depth)
if (class(dt2$wtemp)=="factor") dt2$wtemp <-as.numeric(levels(dt2$wtemp))[as.integer(dt2$wtemp) ]               
if (class(dt2$wtemp)=="character") dt2$wtemp <-as.numeric(dt2$wtemp)
if (class(dt2$flag_wtemp)!="factor") dt2$flag_wtemp<- as.factor(dt2$flag_wtemp)


dt2$bhour <- ifelse(dt2$hour %/% 100 >= 1, dt2$hour/100, dt2$hour)
dt2$datetime <- as.POSIXct(paste0(dt2$sampledate,' ',dt2$bhour,':00:00'), format = "%Y-%m-%d %H:%M:%S")

# Package ID: knb-lter-ntl.29.29 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Physical Limnology of Primary Study Lakes 1981 - current.

# inUrl3 <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/29/29/03e232a1b362900e0f059859abe8eb97"
# infile3 <- tempfile()
# download.file(inUrl3, infile3, method = "curl")
# dt1 <- read_csv(infile3, skip = 1, quote = "\"", guess_max = 1e+05, 
#                 col_names = c("lakeid", "year4", "daynum", "sampledate", 
#                               "depth", "rep", "sta", "event", "wtemp", "o2", "o2sat", 
#                               "deck", "light", "frlight", "flagdepth", "flagwtemp", 
#                               "flago2", "flago2sat", "flagdeck", "flaglight", "flagfrlight"))
# dt1

# time =  startDate + seq(1, ncol(um))*dt#/24/3600
# obs <- dt1 %>%
#   filter(lakeid == 'ME') %>%
#   rename(datetime = sampledate) %>%
#   select(datetime, depth, wtemp)
time =  startDate + seq(1, ncol(um))*dt#/24/3600
obs <- dt2 %>%
  select(datetime, depth, wtemp)

df.sim <- df
colnames(df.sim) <- c("datetime", as.character(paste0('wtemp.',seq(1,nrow(um))*dx)))
df.sim$datetime <-   startDate + seq(1, ncol(um))*dt#/24/3600

# idx <- na.omit(match(as.POSIXct(df.sim$datetime), as.POSIXct(obs$datetime) ))
idx <- (match(as.POSIXct(obs$datetime), as.POSIXct(df.sim$datetime) ))

# df.sim <- df.sim[idx, ]
obs <- obs[which(!is.na(idx)), ]

idz <- which(obs$depth %in% seq(0,24,1))
obs = obs[idz,]

deps <- seq(1,nrow(um))*dx
if (min(unique(obs$depth)) < min(deps)){
  deps[which.min(deps)] <- min(unique(obs$depth)) 
}
if (max(unique(obs$depth)) > max(deps)){
  deps[which.max(deps)] <- max(unique(obs$depth)) 
}

df.sim.interp <- NULL
for (i in 1:nrow(df.sim)){
  df.sim.interp <- rbind(df.sim.interp,
                         approx(deps, df.sim[i, -1], unique(obs$depth))$y)
}
df.sim.interp <- as.data.frame(df.sim.interp)
# df.sim.interp <- apply(df.sim[,-1], 1, function(x) approx(deps, x, unique(obs$temp_depth_m))$y)
df.sim.interp$datetime <-   startDate + seq(1, ncol(um))*dt#/24/3600
colnames(df.sim.interp) <- c(as.character(paste0('wtemp.',unique(obs$depth))), 'datetime')

obs <- data.frame(obs)
obs$depth <- factor(obs$depth)

wide.obs <- reshape(obs, idvar = "datetime", timevar = "depth", direction = "wide")
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
  id.r <-  (match(as.POSIXct(s$datetime), as.POSIXct(o$datetime) ))
  s <- s[which(!is.na(id.r)),]
  rmse <- rbind(rmse, data.frame('variable' = i,
                                 'fit' = sqrt((sum((o$value-s$value)**2, na.rm = T))/nrow(o))))
}

m.obs$variable <-  factor(m.obs$variable, levels=paste0('wtemp.',seq(0,24,1)))
m.df.sim.interp$variable <-  factor(m.df.sim.interp$variable, levels=paste0('wtemp.',seq(0,24,1)))

ggplot() +
  geom_point(data = m.obs,aes(datetime, value, col = group), size =0.3) +
  geom_line(data = m.df.sim.interp, aes(datetime, value, col = group)) +
  geom_text(data=rmse, aes( as.POSIXct("2010-01-01 10:30:00 CDT"), y=17, label=round(fit,2)),                 
            color="black", size =3) +
  # ylim(-5,40)+
  facet_wrap(~ factor(variable, level = c(paste0('wtemp.',seq(0,24,1)))), scales = 'free') +
  xlab('') + ylab('Temp. (deg C)')+
  theme_bw()
ggsave(filename = paste0('fielcomparison','_',scheme,'_',nyear,'.png'), width = 15, height = 8, units = 'in')



## averaged responses
bf.obs <- apply(wide.obs[,-1], 1, function(x) rLakeAnalyzer::buoyancy.freq(wtr = x, depths = as.numeric(unique(obs$depth))))
bf.sim <- apply(df.sim.interp[,-22], 1, function(x) rLakeAnalyzer::buoyancy.freq(wtr = x, depths = as.numeric(unique(obs$depth))))

z.bf.obs <- apply(bf.obs,2, function(x) which.max(x))
z.bf.sim <- apply(bf.sim,2, function(x) which.max(x))
df.z.df.obs <- data.frame('time' = wide.obs$datetime, 'z' = as.numeric(z.bf.obs))
df.z.df.sim <- data.frame('time' = df.sim.interp$datetime, 'z' = z.bf.sim)

g.therm <- ggplot() +
  geom_line(data = df.z.df.obs,
            aes(time, z, col = 'observed'), alpha = 0.7) +
  geom_line(data = df.z.df.sim,
            aes(time, z, col = 'sim'), alpha = 0.7) +
    scale_y_reverse() + xlab('Time') + ylab('Thermocline depth') +
  theme_minimal()

avg.epi.obs <- NULL
avg.hyp.obs <- NULL
for (j in 1:nrow(df.z.df.obs)){
  d = wide.obs[,-1]
  if (is.na(df.z.df.obs$z[j])){
    df.z.df.obs$z[j] = 1
  }
  avg.epi.obs <- append(avg.epi.obs,mean(as.numeric(d[j,1:df.z.df.obs$z[j]], na.rm = T)))
  avg.hyp.obs <- append(avg.hyp.obs,mean(as.numeric(d[j,df.z.df.obs$z[j]:ncol(d)], na.rm = T)))
}

avg.epi.sim <- NULL
avg.hyp.sim <- NULL
for (j in 1:nrow(df.z.df.sim)){
  d = df.sim.interp[,-22]
  if (is.na(df.z.df.sim$z[j])){
    df.z.df.sim$z[j] = 1
  }
  avg.epi.sim <- append(avg.epi.sim,mean(as.numeric(d[j,1:df.z.df.sim$z[j]], na.rm = T)))
  avg.hyp.sim <- append(avg.hyp.sim,mean(as.numeric(d[j,df.z.df.sim$z[j]:ncol(d)], na.rm = T)))
}

df.avg.obs <- data.frame('time' = wide.obs$datetime,
                         'epi' = avg.epi.obs,
                         'hyp' = avg.hyp.obs,
                         'type' = 'obs')
df.avg.sim <- data.frame('time' = df.sim.interp$datetime,
                         'epi' = avg.epi.sim,
                         'hyp' = avg.hyp.sim,
                         'type' = 'sim')

g.avg <- ggplot() +
  geom_point(data = df.avg.obs,
            aes(time, epi, col = 'observed epi'), alpha = 0.7) +
  geom_point(data = df.avg.obs,
            aes(time, hyp, col = 'observed hyp'), alpha = 0.7) +
  geom_line(data = df.avg.sim,
            aes(time, epi, col = 'simulated epi'), alpha = 0.7) +
  geom_line(data = df.avg.sim,
            aes(time, hyp, col = 'simulated hyp'), alpha = 0.7) +
  xlab('Time') + ylab('Average temp.') +
  theme_minimal()

g.average <- g.therm / g.avg; g.average
ggsave(filename = paste0('averaged','_',scheme,'_',nyear,'.png'),plot = g.average, width = 15, height = 8, units = 'in')


# stratification dates
obs_dens <- abs(calc_dens(wide.obs$wtemp.1) - calc_dens(wide.obs$wtemp.20))
sim_dens <- abs(calc_dens(df.sim.interp$wtemp.1) - calc_dens(df.sim.interp$wtemp.20))

bindt <- ifelse(obs_dens >= 0.1, 1, 0)
bindt[which(is.na(bindt))] <- 0
bindy <- ifelse(sim_dens >= 0.1, 1, 0)

df.obs.dens <- data.frame('time' =  wide.obs$datetime,
                          'grad' =  bindt, 
                          'year' = lubridate::year(wide.obs$datetime))
df.sim.dens <- data.frame('time' =  df.sim.interp$datetime,
                          'grad' =  bindy, 
                          'year' = lubridate::year(df.sim.interp$datetime))

ggplot() +
  geom_line(data = df.obs.dens, aes(time, grad, col ='obs'))+
  geom_line(data = df.sim.dens, aes(time, grad,col='sim')) +
  facet_wrap(~ year, scales = 'free')
for (t in unique(lubridate::year(df.obs.dens$time))){
  d = df.obs.dens %>%
    filter(year == t)
  for (i in (1+24*15):(nrow(d) -24*15)){
    if (all(bindt[(i-1-24*15):(i-1)] != 1) & all(bindt[i:(i+(24*15))] == 1))
      print(paste0('start: ',d$datetime[i]))
    if (all(bindt[(i-1-24*15):(i-1)] == 1) & all(bindt[i:(i+(24*15))] != 1))
      print(paste0('end: ',d$datetime[i]))
  }

}


for (i in (1+24*30):(length(sim_dens)-(60*24))){
  if (all(bindy[(i-1-24*30):(i-1)] != 1) & all(bindy[i:(i+(24*30))] == 1))
    print(paste0('start: ',df.sim.interp$datetime[i]))
  if (all(bindy[(i-1-24*30):(i-1)] == 1) & all(bindy[i:(i+(24*30))] != 1))
    print(paste0('end: ',df.sim.interp$datetime[i]))
}



## vertical temperature profiles
for (i in seq(1,ncol(um), length.out = 200)){
  n = i
  i = floor(i)

  sim = m.df.sim.interp %>% 
    filter(datetime == time[i]) %>%
    mutate(group = 'modeled') %>%
    mutate(depth =  as.numeric(gsub(".*wtemp.","",as.character(factor(variable, level = c(paste0('wtemp.',seq(0,24,1))))))))
  obs = m.obs %>%
    filter(datetime == time[i]) %>%
    mutate(group = 'observed') %>%
    mutate(depth =  as.numeric(gsub(".*wtemp.","",as.character(factor(variable, level = c(paste0('wtemp.',seq(0,24,1))))))))
  
  ggplot() +
    geom_path(data = sim, aes(value, 
                              depth, col = group), size = 1.2) +
    # facet_wrap(~ factor(variable, level = c(paste0('wtemp.',seq(0,24,1)))), scales = 'free') +
    geom_point(data = obs ,aes(value, depth, col = group), size =1.2) +
     xlab('temp. (deg C)') + ylab('depth (m)')+
    scale_y_reverse() +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    ggtitle( time[i]) + 
    labs(col='') +
    xlim(-5, 35) + 
    theme_bw()
  
  ggsave(paste0('../../animation_mendota/pic_',match(n, seq(1,ncol(um),length.out=200)),'.png'),
         width = 4, height = 5, units = 'in')
  
}
