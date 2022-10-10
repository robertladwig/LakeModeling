import numpy as np
import pandas as pd
import os
from math import pi, exp, sqrt, log, atan
from scipy.interpolate import interp1d
from copy import deepcopy
import datetime
from ancillary_functions import calc_cc, buoyancy_freq, center_buoyancy

## function to calculate density from temperature
def calc_dens(wtemp):
    dens = (999.842594 + (6.793952 * 1e-2 * wtemp) - (9.095290 * 1e-3 *wtemp**2) +
      (1.001685 * 1e-4 * wtemp**3) - (1.120083 * 1e-6* wtemp**4) + 
      (6.536336 * 1e-9 * wtemp**5))
    return dens

## this is our attempt for turbulence closure, estimating eddy diffusivity
def eddy_diffusivity(rho, depth, g, rho_0, ice, area):
    buoy = np.ones(len(depth)) * 7e-5
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
    daily_meteo['Cloud_Cover'] = calc_cc(date = daily_meteo['date'],
                                                airt = daily_meteo['Air_Temperature_celsius'],
                                                relh = daily_meteo['Relative_Humidity_percent'],
                                                swr = daily_meteo['Shortwave_Radiation_Downwelling_wattPerMeterSquared'],
                                                lat = 43, lon = -89.41,
                                                elev = 258)
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
  area[-1] = area[-2] - 1 # TODO: confirm this is correct
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

def PSIM(zeta):
  # Function to compute stability functions for momentum
  if zeta < 0.0:
    X = (1 - 16*zeta)**0.25
    psim = 2*log((1 + X)/2) + log((1 + X*X)/2)-2*atan(X) + pi/2 
  elif zeta > 0.0:
    if zeta > 0.5:
      if zeta > 10.0:
        psim = log(zeta) - 0.76*zeta - 12.093
      else:
        psim = 0.5/(zeta*zeta) - 4.25/zeta - 7.0*log(zeta) - 0.852
    else:
      psim = -5*zeta
  # Stable case
  else:
    psim = 0.0
  return(psim)



def PSITE(zeta):
  # Function to compute stability functions for sensible and latent heat
  if zeta < 0.0:
    X = (1 - 16*zeta)**0.25
    psite = 2*log((1 + X*X)/2)
  elif zeta > 0.0:# Stable case
    if zeta > 0.5:
      if zeta > 10.0:
        psite = log(zeta) - 0.76*zeta - 12.093
      else:
        psite = 0.5/(zeta*zeta) - 4.25/zeta - 7.0*log(zeta) - 0.852
    else: 
      psite = -5*zeta
  else:
    psite = 0.0
  return(psite)

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

def run_thermalmodel(
  u, 
  startTime, 
  endTime,
  area,
  volume,
  depth,
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
  KEice=0,
  Ice_min=0.1,
  pgdl_mode='off'):
    

  
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
  n2m = np.full([(nx-1), nCol], np.nan)
  mix = np.full([1,nCol], np.nan)
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
    un = deepcopy(u)
    dens_u_n2 = calc_dens(u)
    time_ind = np.where(times == n)
    
    if pgdl_mode == 'on':
      n2 = 9.81/np.mean(dens_u_n2) * (dens_u_n2[1:] - dens_u_n2[:-1])/dx
      n2_pgdl[:,idn] = np.concatenate([n2, np.array([np.nan])])
    kz = eddy_diffusivity(dens_u_n2, depth, 9.81, 998.2, ice, area) / 86400

    if ice and Tair(n) <= 0:
      kzn = kz
      absorp = 1 - 0.7
      infra = 1 - absorp
    elif (ice and Tair(n) >= 0):
      kzn = kz
      absorp = 1 - 0.3
      infra = 1 - absorp
    elif not ice:
      kzn = kz
      absorp = 1 - reflect
      infra = 1 - absorp
    kzm[:,idn] = kzn
    
    ## (1) Heat addition
    # surface heat flux
    Q = (absorp * Jsw(n) + longwave(cc = CC(n), sigma = sigma, Tair = Tair(n), ea = ea(n), emissivity = emissivity, Jlw = Jlw(n)) + #longwave(emissivity = emissivity, Jlw = Jlw(n)) +
            backscattering(emissivity = emissivity, sigma = sigma, Twater = un[0], eps = eps) +
            latent(Tair = Tair(n), Twater = un[0], Uw = Uw(n), p2 = p2, pa = Pa(n), ea=ea(n), RH = RH(n)) + 
            sensible(p2 = p2, B = B, Tair = Tair(n), Twater = un[0], Uw = Uw(n)))  
    
    # heat addition over depth
    H =  (1- infra) * (Jsw(n))  * np.exp(-(kd(n) ) * depth)
    Hg = (area[:-1]-area[1:])/dx * Hgeo/(4181 * calc_dens(un[0]))
    Hg = np.append(Hg, Hg.min())
    
    ## (2) DIFFUSION
    if scheme == 'implicit':
        u[0] = (un[0] + 
        (Q * area[0]/(dx)*1/(4184 * calc_dens(un[0]) ) + abs(H[0+1]-H[0]) * area[0]/(dx) * 1/(4184 * calc_dens(un[0]) ) + 
        Hg[0]) * dt/area[0])
      # all layers in between
        for i in range(1,(nx-1)):
            u[i] = (un[i] + (
                abs(H[i+1]-H[i]) * area[i]/(dx) * 1/(4184 * calc_dens(un[i]) ) + Hg[i])* dt/area[i])
      # bottom layer
        u[(nx-1)] = (un[(nx-1)] + abs(H[(nx-1)]-H[(nx-1)-1]) * area[(nx-1)]/(area[(nx-1)]*dx) * 1/(4181 * calc_dens(un[(nx-1)])) +Hg[(nx-1)]/area[(nx-1)]) * dt
        # IMPLEMENTATION OF CRANK-NICHOLSON SCHEME
        j = len(un)
        y = np.zeros((len(un), len(un)))

        alpha = (dt/dx**2) * kzn

        az = - alpha # subdiagonal
        bz = 2 * (1 + alpha) # diagonal
        cz = - alpha # superdiagonal
    
        bz[:, 0] = 1
        az[:, nx-2] = 0
        bz[:, nx-1] = 1
        cz[:, 0] = 0
    
        az = az[:,1:]
        cz = cz[:,:-1]

        y = np.diag_embed(bz, offset=0)+np.diag_embed(az,offset=-1)+np.diag_embed(cz,offset=1) #slightly efficient way of computing the diagonal matrices
        y[:, nx-1, nx-1] = 1
    
        mn = np.zeros_like(un)  
        mn[:, 0] = un[:, 0]
        mn[:,nx-1] = un[:, nx-1]
        
        mn[:, 1:nx-1] = alpha[:, 1:nx-1]*un[:, :nx-2] + 2 * (1 - alpha[:,1:nx-1])*un[:,1:nx-1] + alpha[:,1:nx-1]*un[:,1:nx-1] #is be same as the loop
    
    # DERIVED TEMPERATURE OUTPUT FOR NEXT MODULE
        u = np.linalg.solve(y, mn)
    # TODO: implement / figure out this
    # TODO: what???
    if scheme == 'explicit':
      u[0] = (un[0] + 
        (Q * area[0]/(dx)*1/(4184 * calc_dens(un[0]) ) + abs(H[0+1]-H[0]) * area[0]/(dx) * 1/(4184 * calc_dens(un[0]) ) + 
        Hg[0]) * dt/area[0])
      # all layers in between
      for i in range(1,(nx-1)):
        u[i] = (un[i] + (area[i] * kzn[i] * 1 / dx**2 * (un[i+1] - 2 * un[i] + un[i-1]) +
          abs(H[i+1]-H[i]) * area[i]/(dx) * 1/(4184 * calc_dens(un[i]) ) + Hg[i])* dt/area[i])
      # bottom layer
      u[(nx-1)] = (un[(nx-1)] +
      (abs(H[(nx-1)]-H[(nx-1)-1]) * area[(nx-1)]/(area[(nx-1)]*dx) * 1/(4181 * calc_dens(un[(nx-1)])) +
      Hg[(nx-1)]/area[(nx-1)]) * dt)
                                                           
    if pgdl_mode == 'on':
      um_diff[:, idn] = u
    ## (3) TURBULENT MIXING OF MIXED LAYER
    # the mixed layer depth is determined for each time step by comparing kinetic 
    # energy available from wind and the potential energy required to completely 
    # mix the water column to a given depth
    Zcv = np.sum(depth * area) / sum(area)  # center of volume
    tau = 1.225 * Cd * Uw(n) ** 2 # wind shear is air density times wind velocity 
    if (Uw(n) <= 15):
      c10 = 0.0005 * sqrt(Uw(n))
    else:
      c10 = 0.0026
      
    shear = sqrt((c10 * calc_dens(un[0]))/1.225) *  Uw(n) # shear velocity
    # coefficient times wind velocity squared
    KE = shear *  tau * dt # kinetic energy as function of wind
    
    if ice:
      KE = KE * KEice
    
    maxdep = 0
    for dep in range(0, nx-1):
      if dep == 0:
        PE = (abs(g *   depth[dep] *( depth[dep+1] - Zcv)  *
             # abs(calc_dens(u[dep+1])- calc_dens(u[dep])))
             abs(calc_dens(u[dep+1])- np.mean(calc_dens(u[0])))))
      else:
        PEprior = deepcopy(PE)
        PE = (abs(g *   depth[dep] *( depth[dep+1] - Zcv)  *
            # abs(calc_dens(u[dep+1])- calc_dens(u[dep]))) + PEprior
            abs(calc_dens(u[dep+1])- np.mean(calc_dens(u[0:(dep+1)])))) + PEprior)
            
      if PE > KE:
        maxdep = dep - 1
        break
      elif dep > 0 and PE < KE:
        u[(dep - 1):(dep+1)] = np.sum(u[(dep-1):(dep+1)] * volume[(dep-1):(dep+1)])/np.sum(volume[(dep-1):(dep+1)])
      
      maxdep = dep
      
    mix[0,idn] = KE/PE #append(mix, KE/PE)
    therm_z[0,idn] = depth[maxdep] #append(therm.z, maxdep)
    if pgdl_mode == 'on':
      um_mix[:, idn] = u

    ## (4) DENSITY INSTABILITIES
    # convective overturn: Convective mixing is induced by an unstable density 
    # profile. All groups of water layers where the vertical density profile is 
    # unstable are mixed with the first stable layer below the unstable layer(s) 
    # (i.e., a layer volume weighed means of temperature and other variables are 
    # calculated for the mixed water column). This procedure is continued until 
    # the vertical density profile in the whole water column becomes neutral or stable.
    dens_u = calc_dens(u) 
    diff_dens_u = np.diff(dens_u) 
    diff_dens_u[abs(diff_dens_u) <= denThresh] = 0
    while np.any(diff_dens_u < 0):
      dens_u = calc_dens(u)
      for dep in range(0, nx-1):
        if dens_u[dep+1] < dens_u[dep] and abs(dens_u[dep+1] - dens_u[dep]) >= denThresh:
          u[(dep):(dep+2)] = np.sum(u[(dep):(dep+2)] * volume[(dep):(dep+2)])/np.sum(volume[(dep):(dep+2)])
          break
        
      dens_u = calc_dens(u)
      diff_dens_u = np.diff(dens_u)
      diff_dens_u[abs(diff_dens_u) <= denThresh] = 0
      
    dens_u_n2 = calc_dens(u)
    n2 = 9.81/np.mean(dens_u_n2) * (dens_u_n2[1:] - dens_u_n2[:-1])/dx
    if np.max(n2) > 1e-4:
      max_n2 = depth[np.argmax(n2)]
    else:
      max_n2 = np.max(depth)
    mix_z[0, idn] = max_n2
    if pgdl_mode == 'on':
      um_conv[:, idn] = u
      
      
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
    iceT = iceT * (1 - x) + u[0] * x
    if (iceT <= 0) and Hi < Ice_min:
      if Hi < 0:
        Hi = 1e-5
        
      # if (any(u <= 0) == TRUE){
      supercooled = u < 0
      initEnergy = np.sum((0-u[supercooled])*area[supercooled] * dx * 4.18E6)
      
      if ice == False:
        Hi = Ice_min+(initEnergy/(910*333500))/np.max(area)
      else:
        if (Tair(n) > 0):
          Tice = 0
          Hi = (Hi - max([0, meltP * dt*((absorp*Jsw(n))+(longwave(cc = CC(n), sigma = sigma, Tair = Tair(n), ea = ea(n), emissivity = emissivity, Jlw = Jlw(n)) +
                                                           backscattering(emissivity = emissivity, sigma = sigma, Twater = un[0], eps = eps) +
                                                           latent(Tair = Tair(n), Twater = un[0], Uw = Uw(n ), p2 = p2, pa = Pa(n), ea=ea(n),  RH = RH(n)) + 
                                                           sensible(p2 = p2, B = B, Tair = Tair(n), Twater = un[0], Uw = Uw(n))) )/(1000*333500)]))
        else:
          Tice =  ((1/(10 * Hi)) * 0 +  Tair(n)) / (1 + (1/(10 * Hi))) 
          Hi = max(Ice_min, sqrt(Hi**2 + 2 * 2.1/(910 * 333500)* (0 - Tice) * dt))
      ice = True
      if (Hi > 0):
        u[supercooled] = 0
        u[0] = 0
      Him[0, idn] = Hi
    elif (ice == True and Hi >= Ice_min):
      if (Tair(n) > 0):
        Tice = 0
        Hi = (Hi - max([0, meltP * dt*((absorp*Jsw(n))+(backscattering(emissivity = emissivity, sigma = sigma, Twater = un[0], eps = eps) +
                                                         latent(Tair = Tair(n), Twater = un[0], Uw = Uw(n ), p2 = p2, pa = Pa(n), ea=ea(n*dt),  RH = RH(n)) + 
                                                         sensible(p2 = p2, B = B, Tair = Tair(n), Twater = un[0], Uw = Uw(n))) )/(1000*333500)]))
      else:
        Tice =  ((1/(10 * Hi)) * 0 +  Tair(n*dt)) / (1 + (1/(10 * Hi))) 
        Hi = max(Ice_min, sqrt(Hi**2 + 2 * 2.1/(910 * 333500)* (0 - Tice) * dt))
      
      u[supercooled] = 0
      u[0] = 0
      Him[0, idn] = Hi
    elif ice == True and Hi < Ice_min:
      ice = False
    
    n2m[:,idn] = n2
    um[:,idn] = u
    
    if pgdl_mode == 'on':
      um_ice[:, idn] = u
      meteo_pgdl[0, idn] = Tair(n)
      meteo_pgdl[1, idn] = (longwave(cc = CC(n), sigma = sigma, Tair = Tair(n), ea = ea(n), emissivity = emissivity, Jlw = Jlw(n)) -
        backscattering(emissivity = emissivity, sigma = sigma, Twater = un[0], eps = eps))
      meteo_pgdl[2, idn] = latent(Tair = Tair(n), Twater = un[0], Uw = Uw(n), p2 = p2, pa = Pa(n), ea=ea(n), RH = RH(n))
      meteo_pgdl[3, idn] = sensible(p2 = p2, B = B, Tair = Tair(n), Twater = un[0], Uw = Uw(n))
      meteo_pgdl[4, idn] = Jsw(n)
      meteo_pgdl[5, idn] = kd(n)
      meteo_pgdl[6, idn] = shear
      meteo_pgdl[7, idn] = tau
      meteo_pgdl[8, idn] = np.nanmax(area)

  end_time = datetime.datetime.now()
  print((end_time - start_time))
  
  bf_sim = np.apply_along_axis(center_buoyancy, axis=1, arr = um.T, depths=depth)
  

  df_z_df_sim = pd.DataFrame({'time': times, 'thermoclineDep': bf_sim})

  df_z_df_sim['epi'] = np.nan
  df_z_df_sim['hypo'] = np.nan
  df_z_df_sim['tot'] = np.nan
  df_z_df_sim['stratFlag'] = np.nan
  for j in range(df_z_df_sim.shape[0]):
    if np.isnan(df_z_df_sim.loc[j, 'thermoclineDep']):
      cur_z = 1
      cur_ind = 0
    else:
      cur_z = df_z_df_sim.loc[j, 'thermoclineDep']
      cur_ind = np.max(np.where(depth < cur_z))
      
    df_z_df_sim.loc[j, 'epi'] = np.sum(um[0:(cur_ind + 1), j] * area[0:(cur_ind+1)]) / np.sum(area[0:(cur_ind+1)])
    df_z_df_sim.loc[j, 'hypo'] = np.sum(um[ cur_ind:, j] * area[cur_ind:]) / np.sum(area[cur_ind:])
    df_z_df_sim.loc[j, 'tot'] = np.sum(um[:,j] * area) / np.sum(area)
    if calc_dens(um[-1,j]) - calc_dens(um[0,j]) >= 0.1 and np.mean(um[:,j]) >= 4:
      df_z_df_sim.loc[j, 'stratFlag'] = 1
    else:
      df_z_df_sim.loc[j, 'stratFlag'] = 0
  
  dat = {'temp' : um,
          'diff' : kzm,
          'mixing' : mix,
          'buoyancy' : n2m,
          'icethickness' : Hi,
          'iceflag' : ice,
          'icemovAvg' : iceT,
          'supercooled' : supercooled,
          'mixingdepth' : mix_z,
          'thermoclinedepth' : therm_z,
          'endtime' : endTime, 
          'average' : df_z_df_sim}
  if pgdl_mode == 'on':
    dat = {'temp' : um,
               'diff' : kzm,
               'mixing' : mix,
               'buoyancy' : n2m,
               'icethickness' : Hi,
               'iceflag' : ice,
               'icemovAvg' : iceT,
               'supercooled' : supercooled,
               'mixingdepth' : mix_z,
               'thermoclinedepth' : therm_z,
               'endtime' : endTime, 
               'average' : df_z_df_sim,
               'temp_diff' : um_diff,
               'temp_mix' : um_mix,
               'temp_conv' : um_conv,
               'temp_ice' : um_ice,
               'meteo_input' : meteo_pgdl,
               'buoyancy_pgdl' : n2_pgdl}
  
  return(dat)

