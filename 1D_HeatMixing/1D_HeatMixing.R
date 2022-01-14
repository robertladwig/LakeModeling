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
#' MyLake—A multi-year lake simulation model code suitable for uncertainty 
#' and sensitivity analysis simulations. Ecological Modeling
#' https://doi.org/10.1016/j.ecolmodel.2007.03.018 

## remove everything from workspace
rm(list = ls())

# set wd to current dir of script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## colors for plotting
library(tidyverse)
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals)))

## lake configurations
zmax = 25 # maximum lake depth
nx = 25 # number of layers we will have
dt = 24* 3600 # 24 hours times 60 min/hour times 60 seconds/min
dx = zmax/nx # spatial step

## area and depth values of our lake 
hyps <- read_csv('bc/LakeEnsemblR_bathymetry_standard.csv')
area = approx(hyps$Depth_meter,hyps$Area_meterSquared,seq(1,nx*dx, 
                                                          length.out= nx))$y
area[which.min(area)] <- 1e-2
depth = depth= seq(1,nx*dx, length.out = nx)

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
meteo <-  read_delim('bc/meteo.txt', delim ='\t')

wQ = data.frame('time' = as.POSIXct(meteo$date, format = '%m/%d/%Y'),
                'Jsw' = meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared, 
                'vW' = meteo$Ten_Meter_Elevation_Wind_Speed_meterPerSecond,
                'airT' = meteo$Air_Temperature_celsius,
                'dewT' = meteo$Dewpoint_Air_Temperature_Celsius)
wQ$vW <- wQ$vW * 0.3 #* 0.3
wQ$Uw <- 1225*0.0013* wQ$vW^2 # 19.0 + 0.95 * (wQ$vW)^2 # wind shear stress
wQ$Uwalt <- 19.0 + 0.95 * (wQ$vW)^2

wQ$dt = wQ$time - (wQ$time[1]) +1

nyear = 5
nt = nyear * 365 * 86400 # as.double(max(wQ$dt)) # maximum simulation length

## linearization of driver data, so model can have dynamic step
Jsw <- approxfun(x = wQ$dt, y = wQ$Jsw, method = "linear", rule = 2)
Tair <- approxfun(x = wQ$dt, y = wQ$airT, method = "linear", rule = 2)
Dew <- approxfun(x = wQ$dt, y = wQ$dewT, method = "linear", rule = 2)
Uw <- approxfun(x = wQ$dt, y = wQ$Uw, method = "linear", rule = 2)
vW <- approxfun(x = wQ$dt, y = wQ$vW, method = "linear", rule = 2)
Uwalt <- approxfun(x = wQ$dt, y = wQ$Uwalt, method = "linear", rule = 2)

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

longwave <- function(sigma, Tair, Acoeff, eair, Rl){  # longwave radiation into 
  # the lake
  lw = (sigma * (Tair(n * dt) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) 
  return(lw)
}
backscattering <- function(eps, sigma, Twater){ # backscattering longwave 
  # radiation from the lake
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

## plot initial profile
plot( u, seq(0, nx * dx, length.out=(nx)),  
      ylim = rev(range(seq(0, dx * nx, length.out=(nx)))), xlim = c(0,35), 
      type ='l');

um <- c()
Hts <- c()
Swf <- c()
Lwf <- c()
BLwf <- c()
Lf <- c()
Sf <- c()
mix <- c()
therm.z <- c()
mix.z <- c()

## modeling code for vertical 1D mixing and heat transport
for (n in 1:floor(nt/dt)){  #iterate through time
  un = u # prior temperature values
  kz = eddy_diffusivity(calc_dens(un), depth, 9.81, 998.2) / 86400
  kzn = kz   
  
  eair <- (4.596 * exp((17.27 * Dew(n * dt)) / (237.3 + Dew(n * dt)))) # air 
  # vapor pressure
  esat <- 4.596 * exp((17.27 * Tair(n * dt)) / (237.3 + Tair(n * dt))) # saturation 
  # vapor pressure
  RH <- eair/esat *100 # relative humidity
  es <- 4.596 * exp((17.27 * u[1])/ (273.3+u[1]))
  # surface heat flux
  Q <- (Jsw(n * dt) + 
          longwave(sigma = sigma, Tair = Tair(n * dt), Acoeff = Acoeff, 
                   eair = eair, Rl = Rl) +
          backscattering(eps = eps, sigma = sigma, Twater = un[1]) +
          latent(c1 = c1, wind = Uw(n * dt), Twater = un[1], Tair = 
                   Tair(n *dt)) + 
          sensible(wind = Uw(n * dt), esat = esat, eair = eair)) #Uwalt
  
  # heat addition over depth
  H = (1- reflect) * (1- infra) * (Jsw(n * dt))  * #
    exp(-(kd ) *seq(dx,nx*dx,length.out=nx)) 

  ## (1) DIFFUSION
  # surface layer
  u[1] = un[1] +
     Q * area[1]/(area[1]*dx)*1/(4184 * calc_dens(un[1]) ) *dt
  
 Hts <- append(Hts, Q *  area[1]/(area[1]*dx)*1/(4181 * calc_dens(un[1]) ))
 Swf <- append(Swf, Jsw(n * dt))
 Lwf <- append(Lwf,  (sigma * (Tair(n * dt) + 273)^4 * (Acoeff + 0.031 * 
                                                          sqrt(eair)) * (1 - Rl)))
 BLwf <- append(BLwf, (-1)*  (eps * sigma * (un[1] + 273)^4))
 Lf <- append(Lf, (-1) * (c1 * Uw(n * dt) * (un[1] - Tair(n * dt))) )
 Sf <- append(Sf, (-1)* (Uw(n * dt) * ((esat) - (eair))) )
 
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
  tau = 1225 * 0.0013 * vW(n * dt)^2 # wind shear is air density times shear 
  # coefficient times wind velocity squared
  KE = vW(n * dt) *  tau * dt # kinetic energy as function of wind
  maxdep = 1
  for (dep in 1:(nx)){
    if (dep == 1){
      PE = abs(g *  ( seq(1,nx)[dep] - Zcv)  * calc_dens(u[dep]) * dx)
    } else {
      PEprior = PE
      PE = abs(g *  ( seq(1,nx)[dep] - Zcv)  * calc_dens(u[dep]) * dx +
                 PEprior) 
    }
      if (PE > KE){
        maxdep = dep
        break
      }
    maxdep = dep
  }
  u[1:maxdep] = mean(u[1:maxdep])
  mix <- append(mix, KE/PE)
  therm.z <- append(therm.z, maxdep)
  
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
  mix.z <- append(mix.z, max.n2)
  
  um <- cbind(um, u)
  
  lines( u, seq(0, dx * nx, length.out=(nx)),
          ylim = rev(range(seq(0, dx * nx, length.out=(nx)))), lty = 'dashed');
}

str(um)
## water temperature time series at different depths
plot(seq(1, ncol(um))*dt/24/3600, um[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(0,40), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(um))*dt/24/3600, um[i,], col = sample(col_vector,1), 
        lty = 'dashed',lwd =2)
}
 
## surface mixed layer depth (direct model output)
therm.z.roll = zoo::rollmean(therm.z, 14)
plot(therm.z, type = 'l',ylim = rev(range( seq(0,zmax, length.out=nx))), 
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


## meteorological heat fluxes
# plot(Hts )
# plot(Swf, col = 'red', type = 'l', ylim = c(-1000, 1000))
# lines(Lwf, col = 'blue')
# lines(BLwf, col = 'cyan')
# lines(Lf, col = 'yellow')
# lines(Sf, col = 'green')
# lines(Swf+Lwf+BLwf+Lf+Sf, col = 'black', lty =2)
# 

## vertical water temperature profiles over time
# df = data.frame('1' =NULL)
# name = NULL
# for (i in seq(1,ncol(um), length.out = 20)){
#   i = floor(i)
#   name = append(name, round((i*dt)/24/3600,1))
#   df = append(df, data.frame('1' = um[,i] + round((i*dt)/24/3600,1)))
# }
# df.m = matrix(unlist(df), ncol = 20)
# colnames(df.m) = name
# 
# df.m.m = reshape2::melt(df.m)
# ggplot2::ggplot(df.m.m) +
#   geom_point(aes(x = value/40, y=Var1,  col = value-Var2, group = Var2)) +
#   geom_line(aes(x = value/40, y=Var1,  col = value-Var2, group = Var2)) +
#   scale_y_reverse() + xlab('Time') + ylab('Depth')+labs(col='Temp')+
#   scale_color_gradient(low = "lightblue", high = "red") +
#   theme_minimal()

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
# filled.contour(x = seq(1, ncol(um))*dt/24/3600,
#         y = seq(1, nx),
#         z = t(um))

