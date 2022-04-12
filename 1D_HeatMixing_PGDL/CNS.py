# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 20:26:13 2022

@author: rladwig
"""

import numpy as np

# INPUT DATA FROM PREVIOUS MODULE
t = np.arange(25, 5, -1) # temperature profile from previous module output
dt = 3600 # model time step - fixed
dx = 1 # model space step - fixed

# OUTPUT FROM MLP
d = np.array([1e-5] * len(t)) # estimated diffusivity values

# IMPLEMENTATION OF CRANK-NICHOLSON SCHEME
j = len(t)
y = np.zeros((len(t), len(t)))

alpha = (dt/dx**2) * d    

az = alpha # subdiagonal
bz = 2 * (1 + alpha) # diagonal
cz = -alpha # superdiagonal

bz[0] = 1
az[len(az)-2] = 0
bz[len(bz)-1] = 1
cz[0] = 0

# tridiagonal matrix
for k in range(j-1):
    y[k][k] = bz[k]
    y[k][k+1] = cz[k]
    y[k+1][k] = az[k]
    
y[j-1, j-1] = 1

mn = t * 0.0    
mn[0] = t[0]
mn[len(mn)-1] = t[len(t)-1]

for k in range(1,j-1):
    mn[k] = alpha[k] * t[k-1] + 2 * (1 - alpha[k]) * t[k] + alpha[k] * t[k]

# DERIVED TEMPERATURE OUTPUT FOR NEXT MODULE
output = np.linalg.solve(y, mn)