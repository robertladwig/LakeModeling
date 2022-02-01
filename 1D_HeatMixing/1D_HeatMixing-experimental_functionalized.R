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
  
  kz = ak * (buoy)**(-0.43)
  return(kz)
}

kz = eddy_diffusivity(rho, depth, g = 9.81, rho_0 = 009.2, ice = FALSE) / 86400# 1e4

## atmospheric boundary conditions
## create daily meteorological variables
meteo <- read_csv('bc/LakeEnsemblR_meteo_standard.csv')

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

densThresh <- 1e-3
ice = FALSE
Hi= 0
iceT <- 6
scheme = 'explicit' # options are 'explicit' (FTCS, Forward Time Centered Space) or 'implicit' (Crank-Nicholson scheme)


run_thermalmodel <- function(u, startTime, endTime, ice, Hi, iceT, scheme, kd_light){
  
  um <- matrix(NA, ncol =length( seq(startTime, endTime, dt)/dt), nrow = nx)
  kzm <- matrix(NA, ncol = length( seq(startTime, endTime, dt)/dt), nrow = nx)
  n2m <- matrix(NA, ncol = length( seq(startTime, endTime, dt)/dt), nrow = nx)
  mix <- rep(NA, length = length( seq(startTime, endTime, dt)/dt))#(floor(endTime/dt - startTime/dt)))
  therm.z <- rep(NA, length =length( seq(startTime, endTime, dt)/dt))
  mix.z <- rep(NA, length = length( seq(startTime, endTime, dt)/dt))
  Him <- rep(NA, length = length( seq(startTime, endTime, dt)/dt))
  
  ## plot initial profile
  plot( u, seq(0, nx * dx, length.out=(nx)),  
        ylim = rev(range(seq(0, dx * nx, length.out=(nx)))), xlim = c(0,35), 
        type ='l');
  
  if (!is.null(kd_light)){
    kd <- approxfun(x = seq(1, endTime, 1), y = rep(kd_light, (endTime)), method = "constant", rule = 2)
  } else {
    kd <- approxfun(x = secview$dt, y = secview$kd, method = "constant", rule = 2)
  }
  
  start.time <- Sys.time()
  ## modeling code for vertical 1D mixing and heat transport
  for (n in seq(startTime, endTime, dt)){#1:(floor(endTime/dt - startTime/dt))){  #iterate through time 1:floor(nt/dt)

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
    kzm[, match(n, seq(startTime, endTime, dt))] <- kzn
    
    ## (1) Heat addition
    # surface heat flux
    Q <- (absorp * Jsw(n) + longwave(cc = CC(n), sigma = sigma, Tair = Tair(n), ea = ea(n), emissivity = emissivity, Jlw = Jlw(n)) + #longwave(emissivity = emissivity, Jlw = Jlw(n)) +
            backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) +
            latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt), RH = RH(n)) + 
            sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n)))  
    
    # heat addition over depth
    H =  (1- infra) * (Jsw(n))  * #
      exp(-(kd(n) ) *seq(dx,nx*dx,length.out=nx)) 
    
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
    tau = 1.225 * Cd * Uw(n)^2 # wind shear is air density times wind velocity 
    if (Uw(n) <= 15) {
      c10 = 0.0005 * sqrt(Uw(n))
    } else {
      c10 = 0.0026
    }
    shear = sqrt((c10 * calc_dens(un[1]))/1.225) *  Uw(n) # shear velocity
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
    mix[match(n, seq(startTime, endTime, dt))] <- KE/PE #append(mix, KE/PE)
    therm.z[match(n, seq(startTime, endTime, dt))] <- maxdep #append(therm.z, maxdep)
    
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
    mix.z[match(n, seq(startTime, endTime, dt))] <- max.n2
    
    
    
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
          Hi = Hi -max(c(0, meltP * dt*((absorp*Jsw(n))+(longwave(cc = CC(n), sigma = sigma, Tair = Tair(n), ea = ea(n), emissivity = emissivity, Jlw = Jlw(n)) +
                                                                backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) +
                                                                latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt),  RH = RH(n)) + 
                                                                sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n))) )/(1000*333500)))
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
        Hi = Hi -max(c(0, meltP * dt*((absorp*Jsw(n))+(backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1]) +
                                                              latent(Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n *dt), p2 = p2, pa = Pa(n*dt), ea=ea(n*dt),  RH = RH(n)) + 
                                                              sensible(p2 = p2, B = B, Tair = Tair(n*dt), Twater = un[1], Uw = Uw(n))) )/(1000*333500))) 
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
    
    n2m[, match(n, seq(startTime, endTime, dt))] <- n2
    um[, match(n, seq(startTime, endTime, dt))] <- u
    
    lines( u, seq(0, dx * nx, length.out=(nx)),
           ylim = rev(range(seq(0, dx * nx, length.out=(nx)))), lty = 'dashed');
    
    # print(un)
    # print(u)
    # cat ("Press [enter] to continue")
    # line <- readline()
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  return(list('temp'  = um,
              'diff' = kzm,
              'mixing' = mix,
              'buoyancy' = n2m,
              'icethickness' = Him,
              'iceflag' = ice,
              'icemovAvg' = iceT,
              'mixingdepth' = mix.z,
              'thermoclinedepth' = therm.z,
              'endtime' = endTime))
}

### EXAMPLE RUNS
# 1 day
temp <- c()
res <- run_thermalmodel(u = u, 
                        startTime = 1, 
                        endTime = 24*3600, 
                        ice = ice, 
                        Hi = Hi, 
                        iceT = iceT, 
                        scheme = scheme,
                        kd_light = NULL)
temp <- cbind(temp, res$temp)
# doing another day
res2 <-  run_thermalmodel(u = res$temp[, ncol(res$temp)], 
                          startTime = res$endtime, 
                          endTime =  res$endtime + 24*3600, 
                          ice = res$iceflag, 
                          Hi = res$icethickness[length(res$icethickness)], 
                          iceT = res$icemovAvg, scheme,
                          kd_light = NULL)
temp <- cbind(temp, res2$temp[,-1])
# doing another, another day
res3 <-  run_thermalmodel(u = res2$temp[, ncol(res2$temp)], 
                          startTime = res2$endtime, 
                          endTime =  res2$endtime + 24*3600, 
                          ice = res2$iceflag, 
                          Hi = res2$icethickness[length(res2$icethickness)], 
                          iceT = res2$icemovAvg, scheme,
                          kd_light = NULL)
temp <- cbind(temp, res3$temp[,-1])

um = temp
str(um)
## water temperature time series at different depths
plot(seq(1, ncol(um))*dt/24/3600, um[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(-1,35), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(um))*dt/24/3600, um[i,], col = sample(col_vector,1), 
        lty = 'dashed',lwd =2)
}

time =  startDate + seq(1, ncol(um))*dt#/24/3600

df <- data.frame(cbind(time, t(um)) )
df.sim <- df
colnames(df.sim) <- c("datetime", as.character(paste0('wtemp.',seq(1,nrow(um))*dx)))
df.sim$datetime <-   startDate + seq(1, ncol(um))*dt#/24/3600


## averaged responses
bf.sim <- apply(df.sim[,-1], 1, function(x) rLakeAnalyzer::buoyancy.freq(wtr = x, depths = as.numeric(unique(obs$depth))))

z.bf.sim <- apply(bf.sim,2, function(x) which.max(x))

df.z.df.sim <- data.frame('time' = df.sim.interp$datetime, 'z' = z.bf.sim)

g.therm <- ggplot() +
  geom_line(data = df.z.df.sim,
            aes(time, z, col = 'sim'), alpha = 0.7) +
    scale_y_reverse() + xlab('Time') + ylab('Thermocline depth') +
  theme_minimal()


avg.epi.sim <- NULL
avg.hyp.sim <- NULL
for (j in 1:nrow(df.z.df.sim)){
  d = df.sim[,-1]
  if (is.na(df.z.df.sim$z[j])){
    df.z.df.sim$z[j] = 1
  }
  avg.epi.sim <- append(avg.epi.sim,((as.numeric(d[j,1:df.z.df.sim$z[j]], na.rm = T) %*% 
                                       area[1:df.z.df.sim$z[j]] )/ 
                                       sum(area[1:df.z.df.sim$z[j]])))
  avg.hyp.sim <- append(avg.hyp.sim,((as.numeric(d[j,df.z.df.sim$z[j]:ncol(d)], na.rm = T)%*% 
                                        area[df.z.df.sim$z[j]:ncol(d)] )/ 
                                       sum(area[df.z.df.sim$z[j]:ncol(d)])))
}

df.avg.sim <- data.frame('time' = df.sim.interp$datetime,
                         'epi' = avg.epi.sim,
                         'hyp' = avg.hyp.sim,
                         'type' = 'sim')

g.avg <- ggplot() +
  geom_line(data = df.avg.sim,
            aes(time, epi, col = 'simulated epi'), alpha = 0.7) +
  geom_line(data = df.avg.sim,
            aes(time, hyp, col = 'simulated hyp'), alpha = 0.7) +
  xlab('Time') + ylab('Average temp.') +
  theme_minimal()

g.average <- g.therm / g.avg; g.average
ggsave(filename = paste0('averaged','_',scheme,'.png'),plot = g.average, width = 15, height = 8, units = 'in')

