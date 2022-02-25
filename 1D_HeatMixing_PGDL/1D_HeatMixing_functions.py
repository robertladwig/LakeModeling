import numpy as np
import pandas as pd
from math import pi, exp
from scipy.interpolate import interp1d

## function to calculate density from temperature
def calc_dens(wtemp):
    dens = (999.842594 + (6.793952 * 1e-2 * wtemp) - (9.095290 * 1e-3 *wtemp**2) +
      (1.001685 * 1e-4 * wtemp**3) - (1.120083 * 1e-6* wtemp**4) + 
      (6.536336 * 1e-9 * wtemp**5))
    return dens

## this is our attempt for turbulence closure, estimating eddy diffusivity
def eddy_diffusivity(rho, depth, g, rho_0, ice, area):
    buoy = np.ones(nx) * 7e-5
    buoy[:-1] = np.abs(rho[1:] - rho[:-1]) / (depth[1:] - depth[:-1]) * g / rho_0
    buoy[-1] = buoy[-2]
        
    low_values_flags = buoy < 7e-5  # Where values are low
    buoy[low_values_flags] = 7e-5
    
    if ice == True:
      ak = 0.000898
    else:
      ak = 0.00706 *( max(area)/1E6)**(0.56)
    
    kz = ak * (buoy)**(-0.43)
    return(kz)
  
def provide_meteorology(meteofile, secchifile, windfactor):
    meteo = pd.read_csv(meteofile)
    daily_meteo = meteo
    daily_meteo['date'] = pd.to_datetime(daily_meteo['datetime'])
    # TODO !!!!!!!!!!!!!!!!!! Implement actual cloud cover function
    daily_meteo['Cloud_Cover'] = 0
    daily_meteo['dt'] = (daily_meteo['date'] - daily_meteo['date'][0]).astype('timedelta64[s]') + 1
    daily_meteo['ea'] = (daily_meteo['Relative_Humidity_percent'] * 
      (4.596 * np.exp((17.27*(daily_meteo['Air_Temperature_celsius'])) /
      (237.3 + (daily_meteo['Air_Temperature_celsius']) ))) / 100)
    daily_meteo['ea'] = ((101.325 * np.exp(13.3185 * (1 - (373.15 / (daily_meteo['Air_Temperature_celsius'] + 273.15))) -
      1.976 * (1 - (373.15 / (daily_meteo['Air_Temperature_celsius'] + 273.15)))**2 -
      0.6445 * (1 - (373.15 / (daily_meteo['Air_Temperature_celsius'] + 273.15)))**3 -
      0.1229 * (1 - (373.15 / (daily_meteo['Air_Temperature_celsius'] + 273.15)))**4)) * daily_meteo['Relative_Humidity_percent']/100)
    daily_meteo['ea'] = (daily_meteo['Relative_Humidity_percent']/100) * 10**(9.28603523 - 2322.37885/(daily_meteo['Air_Temperature_celsius'] + 273.15))
    startDate = pd.to_datetime(daily_meteo.loc[0, 'date'])
    
    ## calibration parameters
    daily_meteo['Shortwave_Radiation_Downwelling_wattPerMeterSquared'] = daily_meteo['Shortwave_Radiation_Downwelling_wattPerMeterSquared'] 
    daily_meteo['Ten_Meter_Elevation_Wind_Speed_meterPerSecond'] = daily_meteo['Ten_Meter_Elevation_Wind_Speed_meterPerSecond'] * windfactor # wind speed multiplier
    
    ## light
    # Package ID: knb-lter-ntl.31.30 Cataloging System:https://pasta.edirepository.org.
    # Data set title: North Temperate Lakes LTER: Secchi Disk Depth; Other Auxiliary Base Crew Sample Data 1981 - current.
    secview0 = pd.read_csv(secchifile)
    secview0['sampledate'] = pd.to_datetime(secview0['sampledate'])
    secview = secview0.loc[secview0['sampledate'] >= startDate]
    if secview['sampledate'].min() >= startDate:
      firstVal = secview.loc[secview['sampledate'] == secview['sampledate'].min(), 'secnview'].values[0]
      firstRow = pd.DataFrame(data={'sampledate': [startDate], 'secnview':[firstVal]})
      secview = pd.concat([firstRow, secview], ignore_index=True)
    secview['dt'] = (secview['sampledate'] - secview['sampledate'][0]).astype('timedelta64[s]') + 1
    secview['kd'] = 1.7 / secview['secnview']
    secview['kd'] = secview.set_index('sampledate')['kd'].interpolate(method="linear").values
    
    return([daily_meteo, secview])
  
def initial_profile(initfile, nx, dx, depth, processed_meteo):
  meteo = processed_meteo
  startDate = meteo['date'].min()
  obs = pd.read_csv(initfile)
  obs['datetime'] = pd.to_datetime(obs['datetime'])
  obs['ditt'] = abs(obs['datetime'] - startDate)
  init_df = obs.loc[obs['ditt'] == obs['ditt'].min()]
  if max(depth) > init_df.Depth_meter.max():
    lastRow = init_df.loc[init_df.Depth_meter == init_df.Depth_meter.max()]
    init_df = pd.concat([init_df, lastRow], ignore_index=True)
    init_df.loc[init_df.index[-1], 'Depth_meter'] = max(depth)
    
  profile_fun = interp1d(init_df.Depth_meter.values, init_df.Water_Temperature_celsius.values)
  out_depths = np.linspace(0, nx*dx, nx) # these aren't actually at the 0, 1, 2, ... values, actually increment by 1.0412; make sure okay
  u = profile_fun(out_depths)
  
  # TODO implement warning about profile vs. met start date
  
  return(u)

def get_hypsography(hypsofile, dx, nx):
  hyps = pd.read_csv(hypsofile)
  out_depths = np.linspace(1, nx*dx, nx)
  area_fun = interp1d(hyps.Depth_meter.values, hyps.Area_meterSquared.values)
  area = area_fun(out_depths)
  area[area == 0] = 1e-2 # TODO: confirm this is correct
  depth = np.linspace(1, nx*dx, nx)
  
  volume = 0.5 * (area[:-1] + area[1:]) * np.diff(depth)
  volume = np.append(volume, 1000)
  
  return([area, depth, volume])
