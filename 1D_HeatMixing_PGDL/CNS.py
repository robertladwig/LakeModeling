# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 20:26:13 2022

@author: rladwig
"""

import numpy as np

# INPUT DATA FROM PREVIOUS MODULE
t = np.array([11.4673, 11.6500, 11.6500, 11.3945, 11.1238, 10.9389, 10.8103, 10.7024,
        10.5719, 10.3305,  9.8574,  9.3268,  8.8159,  8.3246,  7.8447,  7.3750,
         6.9170,  6.4724,  6.0415,  5.6232,  5.2154,  4.8159,  4.4228,  4.0340,
         3.6476]) # temperature profile from previous module output
dt = 3600 # model time step - fixed
dx = 1 # model space step - fixed

# OUTPUT FROM MLP
d = np.array([3.7190e-05,
        3.7190e-05,
        2.0770e-05,
        2.0541e-05,
        2.4498e-05,
        2.8890e-05,
        3.1367e-05,
        2.9106e-05,
        2.2591e-05,
        1.7295e-05,
        1.7028e-05,
        1.7991e-05,
        1.9071e-05,
        2.0154e-05,
        2.1378e-05,
        2.2856e-05,
        2.4662e-05,
        2.6884e-05,
        2.9676e-05,
        3.3355e-05,
        3.7190e-05,
        3.7190e-05,
        3.7190e-05,
        3.7190e-05,
        3.7190e-05]) # estimated diffusivity values

# IMPLEMENTATION OF CRANK-NICHOLSON SCHEME
j = len(t)
y = np.zeros((len(t), len(t)))

alpha = (dt/dx**2) * d    

az = - alpha # subdiagonal
bz = 2 * (1 + alpha) # diagonal
cz = - alpha # superdiagonal

az[0] = 0
bz[0] = 1
cz[len(az)-1] = 0
bz[len(bz)-1] = 1

az =  np.delete(az,0)
cz =  np.delete(cz,len(cz)-1)

# tridiagonal matrix
for k in range(j-1):
    y[k][k] = bz[k]
    y[k][k+1] = cz[k]
    y[k+1][k] = az[k]

y[0,1] = 0    
y[j-1, j-1] = 1
y[j-1, j-2] = 0

mn = t * 0.0    
mn[0] = t[0]
mn[len(mn)-1] = t[len(t)-1]

for k in range(1,j-1):
    mn[k] = alpha[k] * t[k-1] + 2 * (1 - alpha[k]) * t[k] + alpha[k] * t[k]

# DERIVED TEMPERATURE OUTPUT FOR NEXT MODULE
output = np.linalg.solve(y, mn)

print(output)