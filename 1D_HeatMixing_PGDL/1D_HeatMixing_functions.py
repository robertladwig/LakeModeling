import numpy as np
import pandas as pd
from math import pi, exp
from scipy.interpolate import interp1d
import datetime

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

def longwave(cc, sigma, Tair, ea, emissivity, Jlw):  # longwave radiation into
  Tair = Tair + 273.15
  p = (1.33 * ea/Tair)
  Ea = 1.24 * (1 + 0.17 * cc**2) * p**(1/7)
  lw = emissivity * Ea *sigma * Tair**4
  return(lw)

def backscattering(emissivity, sigma, Twater, eps): # backscattering longwave 
  # radiation from the lake
  Twater = Twater + 273.15
  back = -1 * (eps * sigma * (Twater)**4) 
  return(back)

def sensible(p2, B, Tair, Twater, Uw): # convection / sensible heat
  Twater = Twater + 273.15
  Tair = Tair + 273.15
  fu = 4.4 + 1.82 * Uw + 0.26 *(Twater - Tair)
  sensible = -1 * ( p2 * B * fu * (Twater - Tair)) 
  return(sensible)

def latent(Tair, Twater, Uw, p2, pa, ea, RH): # evaporation / latent heat
  Twater = Twater + 273.15
  Tair = Tair + 273.15
  Pressure = pa / 100
  fu = 4.4 + 1.82 * Uw + 0.26 *(Twater - Tair)
  fw = 0.61 * (1 + 10**(-6) * Pressure * (4.5 + 6 * 10**(-5) * Twater**2))
  ew = fw * 10 * ((0.7859+0.03477* Twater)/(1+0.00412* Twater))
  latent = -1* fu * p2 * (ew - ea)# * 1.33) // * 1/6
  return(latent)

def run_thermalmodel(u, startTime, endTime,
  area,
  volume,
  zmax,
  nx,
  dt,
  dx,
  daily_meteo,
  secview,
  ice=False,
  Hi=0,
  iceT=6,
  supercooled=0,
  scheme='explicit',
  kd_light=None,
  denThresh=1e-3,
  reflect=0.3,
  infra=0.7,
  eps=0.97,
  emissivity=0.97,
  sigma=5.67e-8,
  p2=1,
  B=0.61,
  g=9.81,
  Cd=0.0013, # momentum coeff (wind)
  meltP=5,
  dt_iceon_avg=0.8,
  Hgeo=0.1, # geothermal heat
  KEice=1/1000,
  Ice_min=0.1,
  pgdl_mode='off'):
  return(1)
  
  ## linearization of driver data, so model can have dynamic step
  Jsw_fillvals = tuple(daily_meteo.Shortwave_Radiation_Downwelling_wattPerMeterSquared.values[[0, -1]])
  Jsw = interp1d(daily_meteo.dt.values, daily_meteo.Shortwave_Radiation_Downwelling_wattPerMeterSquared.values, kind = "linear", fill_value=Jsw_fillvals, bounds_error=False)
  Jlw_fillvals = tuple(daily_meteo.Longwave_Radiation_Downwelling_wattPerMeterSquared.values[[0,-1]])
  Jlw = interp1d(daily_meteo.dt.values, daily_meteo.Longwave_Radiation_Downwelling_wattPerMeterSquared.values, kind = "linear", fill_value=Jlw_fillvals, bounds_error=False)
  Tair_fillvals = tuple(daily_meteo.Air_Temperature_celsius.values[[0,1]])
  Tair = interp1d(daily_meteo.dt.values, daily_meteo.Air_Temperature_celsius.values, kind = "linear", fill_value=Tair_fillvals, bounds_error=False)
  ea_fillvals = tuple(daily_meteo.ea.values[[0,-1]])
  ea = interp1d(daily_meteo.dt.values, daily_meteo.ea.values, kind = "linear", fill_value=ea_fillvals, bounds_error=False)
  Uw_fillvals = tuple(daily_meteo.Ten_Meter_Elevation_Wind_Speed_meterPerSecond.values[[0, -1]])
  Uw = interp1d(daily_meteo.dt.values, daily_meteo.Ten_Meter_Elevation_Wind_Speed_meterPerSecond.values, kind = "linear", fill_value=Uw_fillvals, bounds_error=False)
  CC_fillvals = tuple(daily_meteo.Cloud_Cover.values[[0,-1]])
  CC = interp1d(daily_meteo.dt.values, daily_meteo.Cloud_Cover.values, kind = "linear", fill_value=CC_fillvals, bounds_error=False)
  Pa_fillvals = tuple(daily_meteo.Surface_Level_Barometric_Pressure_pascal.values[[0,-1]])
  Pa = interp1d(daily_meteo.dt.values, daily_meteo.Surface_Level_Barometric_Pressure_pascal.values, kind = "linear", fill_value=Pa_fillvals, bounds_error=False)
  kd_fillvals = tuple(secview.kd.values[[0,-1]])
  kd = interp1d(secview.dt.values, secview.kd.values, kind = "linear", fill_value=kd_fillvals, bounds_error=False)
  RH_fillvals = tuple(daily_meteo.Relative_Humidity_percent.values[[0,-1]])
  RH = interp1d(daily_meteo.dt.values, daily_meteo.Relative_Humidity_percent.values, kind = "linear", fill_value=RH_fillvals, bounds_error=False)
  
  step_times = np.arange(startTime, endTime, dt)
  nCol = len(step_times)
  um = np.full([nx, nCol], np.nan)
  kzm = np.full([nx, nCol], np.nan)
  n2m = np.full([nx, nCol], np.nan)
  therm_z = np.full([1,nCol], np.nan)
  mix_z = np.full([1,nCol], np.nan)
  Him = np.full([1,nCol], np.nan)
  
  if pgdl_mode == 'on':
    um_diff = np.full([nx, nCol], np.nan)
    um_mix = np.full([nx, nCol], np.nan)
    um_conv = np.full([nx, nCol], np.nan)
    um_ice = np.full([nx, nCol], np.nan)
    n2_pgdl = np.full([nx, nCol], np.nan)
    meteo_pgdl = np.full([9, nCol], np.nan)
  
  if not kd_light is None:
    def kd(n): # using this shortcut for now / testing if it works
      return kd_light
    
  start_time = datetime.datetime.now()
  
  times = np.arange(startTime, endTime, dt)
  for idn, n in enumerate(times):
    un = u
    dens_u_n2 = calc_dens(u)
    time_ind = np.where(times == n)
    print(idn)
    print(n)
    if pgdl_mode == 'on':
      n2 = 9.81/np.mean(dens_u_n2) * (dens_u_n2[1:] - dens_u_n2[:-1])/dx
      n2_pgdl[:,idn] = np.concatenate([n2, np.array([np.nan])])
    
    kz = eddy_diffusivity(dens_u_n2, depth, 9.81, 998.2, ice, area) / 86400
    
    if ice and Tair(n) <= 0:
      kzn = kz
      absorp = 1 - 0.7
      infra = 1 - absorp
    elif (ice & Tair(n) >= 0):
      kzn = kz
      absorp = 1 - 0.3
      infra = 1 - absorp
    elif not ice:
      kzn = kz
      absorp = 1 - reflect
      infra = 1 - absorp
    
    kzm[:,idn] = kzn
      
      
