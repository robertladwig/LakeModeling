#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 11:29:34 2021

@author: robert

temperature equation as Tt = 1/A K Tzz
code is based on 12 steps to Navier-Stokes by (c) Lorena A. Barba, Gilbert F. Forsyth 2017.
eddy diffusivity is estimated from buoyancy frequency according to Hondzo and Stefan (1993)
"""
import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
%matplotlib inline

nx = 30
dx = 2 / (nx - 1)
nt = 90    #the number of timesteps we want to calculate
nu = 1e-5 #the value of viscosity
sigma = .2 #sigma is a parameter, we'll learn more about it later
dt = sigma * dx**2 / numpy.max(nu) #dt is defined using sigma ... more later!

area = numpy.linspace(1e2,0,nx)
depth = numpy.linspace(0,nx,nx)

def calc_dens(wtemp):
    dens = 999.842594 + (6.793952 * 1e-2 * wtemp) - (9.095290 * 1e-3 *wtemp**2) + (1.001685 * 1e-4 * wtemp**3) - (1.120083 * 1e-6* wtemp**4) + (6.536336 * 1e-9 * wtemp**5)
    return dens

u = numpy.ones(nx)  * 10    #a numpy array with nx elements all equal to 1.
u[int(0):int(10)] = 25 #setting u = 2 between 0.5 and 1 as per our I.C.s

rho = calc_dens(u)

def eddy_diffusivity(rho, depth, g, rho_0):
    buoy = numpy.ones(nx) * 7e-5
    
    for i in range(0, nx - 1):
        buoy[i] = numpy.sqrt( numpy.abs(rho[i+1] - rho[i]) / (depth[i+1] - depth[i]) * g/rho_0 )
        
    low_values_flags = buoy < 7e-5  # Where values are low
    buoy[low_values_flags] = 7e-5
    
    kz = 0.00706 *( 3.8 * 1e1)**(0.56) * (buoy)**(-0.43)
    return(kz)

kz = eddy_diffusivity(rho, depth, 9.81, 998.2) / 1e4

pyplot.plot(numpy.linspace(0, 30, nx), u);

un = numpy.ones(nx) #our placeholder array, un, to advance the solution in time

for n in range(nt):  #iterate through time
#    print(u)
#    print(n)
    un = u.copy() ##copy the existing values of u into un
    kz = eddy_diffusivity(calc_dens(un), depth, 9.81, 998.2) / 1e4
    kzn = kz.copy() 
    for i in range(1, nx - 1):
        u[i] = un[i] + 1/area[i] * kzn[i] * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])
        
pyplot.plot(numpy.linspace(0, 30, nx), u);