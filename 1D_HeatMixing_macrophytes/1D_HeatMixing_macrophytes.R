#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#' Created on Thu Aug 19 11:29:34 2021
#' 
#' @author: robert
#' @email: rladwig2@wisc.edu
#' 
#' temperature equation as Tt = 1/A K Tzz
#' Diffusion code is based on 12 steps to Navier-Stokes by (c) Lorena A. Barba, Gilbert F. Forsyth 2017.
#' 
#' Eddy diffusivity is estimated from buoyancy frequency according to Hondzo and Stefan (1993)
#' 
#' Mixing dynamics code is taken from Herb & Stefan (2004) Temperature stratification and Mixing 
#' Dynamics in a Shallow Lake with Submersed Macrophytes. Lake & Reservoir Management
#' https://www.tandfonline.com/doi/pdf/10.1080/07438140409354159
#' 
#' Convective overturn algorithm is taken from Salorante & Andersen (2007) MyLake—A multi-year 
#' lake simulation model code suitable for uncertainty and sensitivity analysis simulations. 
#' Ecological Modeling
#' https://doi.org/10.1016/j.ecolmodel.2007.03.018 


nx = 30
dx = 2 / (nx - 1)
days = 365*1
nt = 86400 * days    #the number of timesteps we want to calculate
nu = 1e-6 #the value of viscosity
sigma = .2 #sigma is a parameter, we'll learn more about it later
dt = sigma * dx**2 / max(nu) #dt is defined using sigma ... more later!
# dt = 1 
dx = 1.03

# area and depth values of our lake 
area = seq(-1e2,0, length.out = nx) * (-1)
depth = seq(0,nx, length.out = nx)

calc_dens <-function(wtemp){
  dens = 999.842594 + (6.793952 * 1e-2 * wtemp) - (9.095290 * 1e-3 *wtemp**2) + 
    (1.001685 * 1e-4 * wtemp**3) - (1.120083 * 1e-6* wtemp**4) + 
    (6.536336 * 1e-9 * wtemp**5)
  return(dens)
}

# here we define our initial profile
u = rep(1,max(nx))  * 5#10    
# u[(0):(10)] = 7

rho = calc_dens(u)

# this is our attempt for turbulence closure, estimating eddy diffusivity
eddy_diffusivity <-function(rho, depth, g, rho_0){
  buoy = rep(1, length(nx)) * 7e-5
  for (i in seq(1, nx-1)){#range(0, nx - 1):
    buoy[i] = sqrt( abs(rho[i+1] - rho[i]) / (depth[i+1] - depth[i]) * g/rho_0 )
  }
  
  low_values_flags = buoy < 7e-5  # Where values are low
  buoy[low_values_flags] = 7e-5
  
  kz = 0.00706 *( 3.8 * 1e1)**(0.56) * (buoy)**(-0.43)
  return(kz)
}

kz = eddy_diffusivity(rho, depth, 9.81, 998.2) / 86400# 1e4

# plot initial profile
plot( u, seq(0, 30, length.out=(nx)),  
      ylim = rev(range(seq(0, 30, length.out=(nx)))), xlim = c(0,35), type ='l');
# lines( u, seq(0, 30, length.out=(nx)),  
       # ylim = rev(range(seq(0, 30, length.out=(nx)))) );

# atmospheric boundary conditions
bc =c(seq(0,350, length.out = 43200) ,seq(350,0, length.out=43200)) 
bc.approx = approxfun(x = seq(1,86400), y = bc, method = "linear", rule = 2)
# linearize bc to get bigger time steps

bound <- matrix(c(seq(1,12,1),
                  169, 274, 414, 552, 651, 684, 642, 537, 397, 259, 160, 127,
                  8.3, 9., 13.5,13.9,21.8,24.7,29.4,26.6,24.9,15.,9.7,6.6,
                  2.8,3.3,4.9,4.,5.3,7.8,11.8,11.5,7.7,6.8,6.5,2.4,
                  11.6,11.7,16.4,15.6,16.6,16.7,12.7,11.7,14.,12.9,14.8,11.6), nrow = 12, byrow = FALSE)
bound <- as.data.frame(bound)
colnames(bound) <- c('Month','Jsw','Tair','Dew','vW')
bound$Uw <- 19.0 + 0.95 * (bound$vW * 1000/3600)^2 # function to calculate wind shear stress (and transforming wind speed from km/h to m/s)
bound$vW <- bound$vW * 1000/3600
bound$Day <- cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30))

Jsw <- approxfun(x = bound$Day * 24 * 3600, y = bound$Jsw, method = "linear", rule = 2)
Tair <- approxfun(x = bound$Day* 24 * 3600, y = bound$Tair, method = "linear", rule = 2)
Dew <- approxfun(x = bound$Day* 24 * 3600, y = bound$Dew, method = "linear", rule = 2)
Uw <- approxfun(x = bound$Day* 24 * 3600, y = bound$Uw, method = "linear", rule = 2)
vW <- approxfun(x = bound$Day* 24 * 3600, y = bound$vW, method = "linear", rule = 2)

Rl <- 0.3
Acoeff <- 0.6 # coefficient between 0.5 - 0.7
sigma <- 11.7 * 10^(-8) # cal / (cm2 d K4) or: 4.9 * 10^(-3) # Stefan-Boltzmann constant in [J (m2 d K4)-1]
eps <- 0.97 # emissivity of water
rho <- 0.9982 # density (g per cm3)
cp <- 0.99 # specific heat (cal per gram per deg C)
c1 <- 0.47 # Bowen's coefficient
a <- 7 # constant
c <- 9e4 # empirical constant
g <- 9.81  # gravity (m/s2)

reflect <- 0.6 # fraction of reflected solar radiation
infra = 0.3 # fraction infrared radiation
kd = 0.2 # light attenuation coefficient
km = 0.4 # specific light attenuation coefficient for macrophytes
P = 0 # macrophyte biomass per unit volume in gDW m-3

um <- c()
# modeling code for vertical 1D diffusion
for (n in 1:floor(nt/dt)){  #iterate through time
  un = u ##copy the existing values of u into un
  kz = eddy_diffusivity(calc_dens(un), depth, 9.81, 998.2) / 86400#1e4
  kzn = kz     # u[0] = un[0] + 1/area[0] * kzn[0] * dt / dx**3 * (2 * un[0] - 5 * un[0+1] + 4 * un[0+2] - un[0+3]) + bc[n]/(depth[0+1]-depth[0]) * 1/(4181 * calc_dens(un[0]))
  # u[0] = un[0] + bc[n]/(depth[0+1]-depth[0]) * 1/(4181 * calc_dens(un[0]) * area[0])
  eair <- (4.596 * exp((17.27 * Dew(n * dt)) / (237.3 + Dew(n * dt)))) # air vapor pressure
  esat <- 4.596 * exp((17.27 * Tair(n * dt)) / (237.3 + Tair(n * dt))) # saturation vapor pressure
  RH <- eair/esat *100 # relative humidity
  es <- 4.596 * exp((17.27 * u[1])/ (273.3+u[1]))
  Q <- (Jsw(n * dt) + 
          (sigma * (Tair(n * dt) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
    (eps * sigma * (u[1] + 273)^4)  - # backscattering longwave radiation from the lake
    (c1 * Uw(n * dt) * (u[1] - Tair(n * dt))) - # convection
    (Uw(n * dt) * ((es) - (eair))) ) # evaporation
  
  H = (1- reflect) * (1- infra) * Jsw(n * dt) * exp(-(kd + km * P) *seq(1,nx)) 
  
  u[1] = un[1] +  kzn[1] * dt / dx**2 *  (un[2] - un[1]) + #1/area[1] *
     Q * area[1]/(4181 * calc_dens(un[1]) ) +
     H[1] * 1/(4181 * calc_dens(un[1]) ) #* area[0]) bc.approx(n*dt)/(depth[1+1]-depth[1])
  # u[nx] = un[nx] + 1/area[nx] * kzn[nx] * dt / dx**3 * (2 * un[dx] - 5 * un[dx-1] + 4 * un[dx-2] - un[dx-3])
  for (i in 2:(nx-1)){
    u[i] = un[i] + #1/area[i] *# (area[i]-area[i+1])/(depth[i+1]-depth[i]) * 
      kzn[i] * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1]) +
      H[i] * 1/(4181 * calc_dens(un[i]) )
  }
  
  # the mixed layer depth is determined for each time step by comparing kinetic enery available
  # from wind and the potential energy required to completely mix the water column to a given depth
  Zcv <- seq(1, nx) %*% area / sum(area)
  KE = Uw(n * dt) *  vW(n * dt) * dt
  maxdep = 1
  for (dep in 1:(nx)){
    if (dep == 1){
      # PE = seq(1,nx)[dep] * g * ( seq(1,nx)[dep+1] - Zcv) * (
        # calc_dens(un[dep+1]) - calc_dens(un[dep]))
      PE = abs(g/area[1] *  ( seq(1,nx)[dep] - Zcv) * area[dep] * calc_dens(u[dep]) * 1 )
    } else {
      PEprior = PE
      # PE = seq(1,nx)[dep] * g * ( seq(1,nx)[dep+1] - Zcv) * (
        # calc_dens(un[dep+1]) - calc_dens(un[dep])) + PEprior
      PE = abs(g/area[1] *  ( seq(1,nx)[dep] - Zcv) * area[dep] * calc_dens(u[dep]) * 1 +
        PEprior)
      # PE = abs(g/area[1] *  ( seq(1,nx)[1:dep] - Zcv) * area[1:dep] * calc_dens(un[1:dep]) + seq(1,nx)[1:dep])
    }
    # PE = abs(g/area[1] *  ( seq(1,nx)[1:dep] - Zcv) * area[1:dep] * calc_dens(un[1:dep]) + seq(1,nx)[1:dep]) 
      if (PE > KE){
        maxdep = dep
        break
      }
    maxdep = dep
  }
  # if (maxdep != 1){print('mixing!')}
  u[1:maxdep] = rep(u[1],maxdep)#max(u[1:maxdep])
  
  # convective overturn: Convective mixing is induced by an unstable density profile. 
  # All groups of water layers where the vertical density profile is unstable are mixed with the 
  # first stable layer below the unstable layer(s) (i.e., a layer volume weighed means of 
  # temperature and other variables are calculated for the mixed water column). 
  # This procedure is continued until the vertical density profile in the whole water column becomes neutral or stable.
  dens_u = calc_dens(u) 
  diff_dens_u <- (diff(dens_u)) 
  diff_dens_u[abs(diff(dens_u)) < 1e-1] = 0
  while (any(diff_dens_u < 0)){
    dens_u = calc_dens(u) 
    for (dep in 1:(nx-1)){
      if (dens_u[dep+1] < dens_u[dep] & abs(dens_u[dep+1] - dens_u[dep]) > 1e-1){
        u[dep:(dep+1)] = mean(u[dep:(dep+1)])#max(u[1:maxdep])
        print(dep)
        break
      }
    }
    dens_u = calc_dens(u) 
    diff_dens_u <- (diff(dens_u)) 
    diff_dens_u[abs(diff(dens_u)) < 1e-1] = 0
  }

  
  um <- cbind(um, u)
  
  
  lines( u, seq(0, 30, length.out=(nx)),
          ylim = rev(range(seq(0, 30, length.out=(nx)))), lty = 'dashed');
}

str(um)
plot(seq(1, ncol(um))*dt/24/3600, um[1,], col = 'red', type = 'l', xlab = 'Time (d)', ylab='Temeprature (degC)')
lines(seq(1, ncol(um))*dt/24/3600, um[15,], col = 'orange', lty = 'dashed')
lines(seq(1, ncol(um))*dt/24/3600, um[18,], col = 'green', lty = 'dashed')
lines(seq(1, ncol(um))*dt/24/3600, um[20,], col = 'magenta', lty = 'dashed')
lines(seq(1, ncol(um))*dt/24/3600, um[25,], col = 'blue', lty = 'dashed')

for (i in seq(1,ncol(um), length.out = 100)){
  i = floor(i)
  png(paste0('Projects/animation_macrophyte/',i,'.png'))
  plot(um[,i], seq(0,30, length.out=nx),
       ylim = rev(range(seq(0, 30, length.out=(nx)))), col = 'red', type = 'l', xlab = 'Time (d)', ylab='Temeprature (degC)',
       xlim = c(-1,30), main = paste0('time (d): ',round((i*dt)/24/3600),1))
  dev.off()
}

filled.contour(x = seq(1, ncol(um))*dt/24/3600,
        y = seq(1, nx),
        z = t(um))

