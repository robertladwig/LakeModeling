import datetime
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from scipy.interpolate import interp1d
from math import pi, cos, sin, isinf

## function to calculate density from temperature
def calc_dens(wtemp):
    dens = (999.842594 + (6.793952 * 1e-2 * wtemp) - (9.095290 * 1e-3 *wtemp**2) +
      (1.001685 * 1e-4 * wtemp**3) - (1.120083 * 1e-6* wtemp**4) + 
      (6.536336 * 1e-9 * wtemp**5))
    return dens

## python implementation of gotmtools::calc_cc
def calc_cc(date, airt,  swr, lat, lon, elev,  relh = None, dewt = None,daily = False): 
    if daily == True:
        date = pd.date_range(start=date[0], end=(date.iloc[-1] + timedelta(hours=23)), freq='1H') 
    yday = date.dt.dayofyear.values 
    hour = date.dt.hour.values
    hour[hour == 0] = 24
    std_mer = np.linspace(-90, 90, 13)
    Lsm = std_mer[np.argmin(abs(lon - std_mer))]
    Hsc = 1390
    cd = 0.06
    Rg = 0.045
    theta = lat * pi/180
    r = 1 + 0.017 * np.cos((2 * pi/365) * (186 - yday))
    d = 23.45 * pi/180 * np.cos((2 * pi/365) * (172 - yday))
    dts = (1/15) * (Lsm - lon)
    value = (sin(theta) * np.sin(d))
    value = value/(cos(theta) * np.cos(d))
    tss = (12/pi) * np.arccos(-1*value) + dts + 12
    tsu = -tss + (2 * dts) + 24
    gamma = np.zeros(len(tss))
    dum = np.where(np.logical_and(hour > tsu, hour < tss))
    gamma[dum] = 1
    dum1 = np.where(hour <= 12)
    dum2 = np.where(hour > 12)
    hb1 = pi/12 * (hour - 1 - dts)
    hb1[dum1] = hb1[dum1] + pi
    hb1[dum2] = hb1[dum2] - pi
    hb = hb1
    dum3 = np.where(hb1 > 2 * pi)
    hb[dum3] = hb[dum3] - 2 * pi
    dum4 = np.where(hb1 < 0)
    hb[dum4] = hb[dum4] + 2 * pi
    he1 = pi/12 * (hour - dts)
    he1[dum1] = he1[dum1] + pi
    he1[dum2] = he1[dum2] - pi
    he = he1
    dum3 = np.where(he1 > 2 * pi)
    he[dum3] = he[dum3] - 2 * pi
    dum4 = np.where(he1 < 0)
    he[dum4] = he[dum4] + 2 * pi
    Ho = (Hsc/(r**2) * (sin(theta) * np.sin(d) + 12/pi * cos(theta) * 
        np.cos(d) * (np.sin(he) - np.sin(hb))) * gamma)
    w = (he + hb)/2
    alpha1 = abs(sin(theta) * np.sin(d) + cos(theta) * np.cos(d) * np.cos(w))
    alpha = np.arctan(alpha1/np.sqrt(1 - alpha1**2))
    theta_am1 = ((288 - 0.0065 * elev)/288)**5.256
    theta_am2 = np.sin(alpha) + 0.15 * ((alpha * 180/pi) + 3.855)**(-1.253)
    theta_am = theta_am1/theta_am2
    if dewt is None:
        dewt = (243.04 * (np.log(relh/100) + ((17.625 * airt)/(243.04 + 
            airt)))/(17.625 - np.log(relh/100) - ((17.625 * airt)/(243.04 + 
            airt))))
    if daily == True:
        dewt = np.repeat(dewt, 24)
    Pwc = 0.85 * np.exp(0.11 + 0.0614 * dewt)
    a2 = np.exp(-(0.465 + 0.134 * Pwc) * (0.179 + 0.421 * np.exp(-0.721 * theta_am)) * theta_am)
    a1 = np.exp(-(0.465 + 0.134 * Pwc) * (0.129 + 0.171 * np.exp(-0.88 * theta_am)) * theta_am)
    at = (a2 + 0.5 * (1 - a1 - cd))/(1 - 0.5 * Rg * (1 - a1 - cd))
    Ho = at * Ho
    dum5 = np.where(Ho < 0)
    Ho.iloc[dum5] = 1
    df = pd.DataFrame({'DateTime' : date.values, 'Ho' : Ho.values})
    if daily == True:
        df['DateTime'] = df.DateTime.dt.date
        dfd = df.groupby(['DateTime'])['Ho'].mean()
        df = dfd.reset_index()
    df['swr'] = swr
    df['ccsim'] = np.nan
    
    df.loc[df.Ho >= df.swr, "ccsim"] = ((1 - (df['swr']/df['Ho'])) / 0.65).apply(np.sqrt)

    df.loc[df.ccsim > 1, "ccsim"] = 1
    ccsim_fillvals = tuple(df.ccsim.dropna().values[[0,-1]])
    df['dt'] = (df['DateTime'] - df['DateTime'][0]).astype('timedelta64[s]') + 1
    df_notNA = df.dropna(subset=['ccsim'])
    ccsim_fun = interp1d(df_notNA.dt.values, df_notNA.ccsim.values, kind = "linear", fill_value=ccsim_fillvals, bounds_error=False)
    ccsim = ccsim_fun(df.dt.values)
    
    return(ccsim)

def buoyancy_freq(wtr, depths):
    rhoVar = calc_dens(wtr)
    numDepths = len(depths)
    n2 = 9.81 / rhoVar[:-1] * (rhoVar[1:] - rhoVar[:-1]) / (depths[1:] - depths[:-1])
    n2depths = (depths[1:] + depths[:-1]) / 2 # not returning this; should be okay if grid is regularly spaced
    return(n2)

def center_buoyancy(wtr, depths):
    N2 = buoyancy_freq(wtr, depths)
    dz = depths[1:] - depths[:-1]
    areas = dz * N2
    cent_depths = (depths[1:] + depths[:-1]) / 2 
    
    areas[areas < 0] = 0
    cent_buoyancy = np.sum(cent_depths * areas)/np.sum(areas)
    if isinf(cent_buoyancy):
        return(np.nan)
    else:
        return(cent_buoyancy)
