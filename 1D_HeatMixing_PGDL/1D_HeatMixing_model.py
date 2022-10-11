import numpy as np
import pandas as pd
import os
from math import pi, exp, sqrt
from scipy.interpolate import interp1d
from copy import deepcopy
import datetime
import matplotlib.pyplot as plt
import seaborn as sns


# Robert: can add your path(s) here explicitly; otherwise if it's not me I'm assuming it's the one from your scripts
if os.environ.get("USERNAME") == 'cal':
     os.chdir("/home/cal/Documents/PostDoc/Projects/LakeModeling/1D_HeatMixing_PGDL")
else:
     os.chdir("C:/Users/ladwi/Documents/Projects/R/LakeModeling/1D_HeatMixing_PGDL")
from oneD_HeatMixing_Functions import get_hypsography, provide_meteorology, initial_profile, run_thermalmodel

## lake configurations
zmax = 25 # maximum lake depth
nx = 25 # number of layers we will have
dt = 3600 # 24 hours times 60 min/hour times 60 seconds/min
dx = zmax/nx # spatial step

## area and depth values of our lake 
hyps_all = get_hypsography(hypsofile = 'bc/LakeEnsemblR_bathymetry_standard.csv',
                            dx = dx, nx = nx)
                            
## atmospheric boundary conditions
meteo_all = provide_meteorology(meteofile = 'bc/LakeEnsemblR_meteo_standard.csv',
                    secchifile = 'bc/light.csv', 
                    windfactor = 0.8)

## here we define our initial profile
u_ini = initial_profile(initfile = 'bc/obs.txt', nx = nx, dx = dx,
                     depth = hyps_all[1],
                     processed_meteo = meteo_all[0])
                     
hydrodynamic_timestep = 24 * dt
total_runtime = 365 * 1

startingDate = meteo_all[0]['date'][0]

nTotalSteps = int(total_runtime * hydrodynamic_timestep/ dt)
temp = np.full([nx, nTotalSteps], np.nan)
avgtemp = np.full([nTotalSteps, 6], np.nan)
temp_diff = np.full([nx, nTotalSteps], np.nan)
temp_mix = np.full([nx, nTotalSteps], np.nan)
temp_conv = np.full([nx, nTotalSteps], np.nan)
temp_ice = np.full([nx, nTotalSteps], np.nan)
diff = np.full([nx, nTotalSteps], np.nan)
meteo = np.full([9, nTotalSteps], np.nan)
buoyancy = np.full([nx, nTotalSteps], np.nan)

Start = datetime.datetime.now()
if 'res' in locals() or 'res' in globals():
  del res
  
for i in range(total_runtime):
  if 'res' in locals() or 'res' in globals():
    u = res['temp'][:,-1]
    startTime = res['endtime']
    endTime = res['endtime'] + hydrodynamic_timestep - 1
    ice = res['iceflag']
    Hi = res['icethickness']
    iceT = res['icemovAvg']
    supercooled = res['supercooled']
    kd_light = None
    matrix_range_start = deepcopy(matrix_range_end)# max(0, round(startTime/dt))
    matrix_range_end = matrix_range_start + 24# round(endTime/dt) 
  else:
    u = deepcopy(u_ini)
    startTime = 1
    endTime = hydrodynamic_timestep - 1
    ice = False
    Hi = 0
    iceT = 6
    supercooled = 0
    kd_light = None
    matrix_range_start = 0 #max(0, round(startTime/dt))
    matrix_range_end = 24 #round(endTime/dt)
    
  res = run_thermalmodel(
    u = u,
    startTime = startTime, 
    endTime =  endTime,
    area = hyps_all[0],
    volume = hyps_all[2],
    depth = hyps_all[1],
    zmax = zmax,
    nx = nx,
    dt = dt,
    dx = dx,
    daily_meteo = meteo_all[0],
    secview = meteo_all[1],
    ice = ice,
    Hi = Hi,
    iceT = iceT,
    supercooled = supercooled,
    scheme='implicit',
    kd_light = kd_light,
    denThresh=1e-3,
    reflect=0.3,
    infra=0.7,
    eps=0.97,
    emissivity=0.97,
    sigma=5.67e-8,
    p2=1,
    B=0.61,
    g=9.81,
    Cd = 0.0008, # momentum coeff (wind)
    meltP=5,
    dt_iceon_avg=0.8,
    Hgeo=0.1, # geothermal heat
    KEice=1/1000,
    Ice_min=0.1,
    pgdl_mode = 'on')
  
  temp[:, matrix_range_start:(matrix_range_end)] =  res['temp']
  diff[:, matrix_range_start:matrix_range_end] =  res['diff']
  avgtemp[matrix_range_start:matrix_range_end,:] = res['average'].values
  temp_diff[:, matrix_range_start:matrix_range_end] =  res['temp_diff']
  temp_mix[:, matrix_range_start:matrix_range_end] =  res['temp_mix']
  temp_conv[:, matrix_range_start:matrix_range_end] =  res['temp_conv']
  temp_ice[:, matrix_range_start:matrix_range_end] =  res['temp_ice']
  buoyancy[:, matrix_range_start:matrix_range_end] =  res['temp']
  meteo[:, matrix_range_start:matrix_range_end] =  res['meteo_input']
  buoyancy[:, matrix_range_start:matrix_range_end] = res['buoyancy_pgdl']
  
# convert averages from array to data frame
avgtemp_df = pd.DataFrame(avgtemp, columns=["time", "thermoclineDep", "epi", "hypo", "tot", "stratFlag"])

End = datetime.datetime.now()
print(End - Start)

# epi/hypo/total
colors = ['#F8766D', '#00BA38', '#619CFF']
avgtemp_df.plot(x='time', y=['epi', 'hypo', 'tot'], color=colors, kind='line')
plt.show()

# stratflag
avgtemp_df.plot(x='time', y=['stratFlag'], kind='line', color="black")
plt.show()

# thermocline depth
avgtemp_df.plot(x='time', y=['thermoclineDep'], color="black")
plt.gca().invert_yaxis()
plt.scatter(avgtemp_df.time, avgtemp_df.stratFlag, c=avgtemp_df.stratFlag)
plt.show()

# heatmap of temps  
plt.subplots(figsize=(40,40))
sns.heatmap(temp, cmap=plt.cm.get_cmap('Spectral_r'), xticklabels=1000, yticklabels=2)
plt.show()


