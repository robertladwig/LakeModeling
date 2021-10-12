#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#' Created on Thu Aug 19 11:29:34 2021
#' 
#' @author: robert
#' @email: rladwig2@wisc.edu
#' 
#' temperature equation as Tt = 1/A K Tzz
#' code is based on 12 steps to Navier-Stokes by (c) Lorena A. Barba, Gilbert F. Forsyth 2017.
#' eddy diffusivity is estimated from buoyancy frequency according to Hondzo and Stefan (1993)
#' 
#' Mixing dynamics code is taken from Herb & Stefan (2004) Temperature stratification and Mixing 
#' Dynamics in a Shallow Lake with Submersed Macrophytes. Lake & Reservoir Management
#' https://www.tandfonline.com/doi/pdf/10.1080/07438140409354159


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
u = rep(1,max(nx))  * 10    
u[(0):(10)] = 25 

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
      ylim = rev(range(seq(0, 30, length.out=(nx)))), xlim = c(0,40), type ='l');
# lines( u, seq(0, 30, length.out=(nx)),  
       # ylim = rev(range(seq(0, 30, length.out=(nx)))) );

# atmospheric boundary conditions
bc =c(seq(0,350, length.out = 43200) ,seq(350,0, length.out=43200)) 
bc.approx = approxfun(x = seq(1,86400), y = bc, method = "linear", rule = 2)
# linearize bc to get bigger time steps

# modeling code
for (n in 1:floor(nt/dt)){  #iterate through time
  un = u ##copy the existing values of u into un
  kz = eddy_diffusivity(calc_dens(un), depth, 9.81, 998.2) / 86400#1e4
  kzn = kz     # u[0] = un[0] + 1/area[0] * kzn[0] * dt / dx**3 * (2 * un[0] - 5 * un[0+1] + 4 * un[0+2] - un[0+3]) + bc[n]/(depth[0+1]-depth[0]) * 1/(4181 * calc_dens(un[0]))
  # u[0] = un[0] + bc[n]/(depth[0+1]-depth[0]) * 1/(4181 * calc_dens(un[0]) * area[0])
  u[1] = un[1] + 1/area[1] * kzn[1] * dt / dx**2 *  (un[2] - un[1]) + 
    bc.approx(n*dt)/(depth[1+1]-depth[1]) * 1/(4181 * calc_dens(un[1]) )#* area[0])
  # u[nx] = un[nx] + 1/area[nx] * kzn[nx] * dt / dx**3 * (2 * un[dx] - 5 * un[dx-1] + 4 * un[dx-2] - un[dx-3])
  for (i in 2:(nx-1)){
    u[i] = un[i] + 1/area[i] * (area[i]-area[i+1])/(depth[i+1]-depth[i]) * 
      kzn[i] * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])
  }
  lines( u, seq(0, 30, length.out=(nx)),
          ylim = rev(range(seq(0, 30, length.out=(nx)))), lty = 'dashed');
}