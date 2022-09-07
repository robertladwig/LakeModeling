## function to calculate density from temperature
calc_dens <-function(wtemp){
  dens = 999.842594 + (6.793952 * 1e-2 * wtemp) - (9.095290 * 1e-3 *wtemp**2) + 
    (1.001685 * 1e-4 * wtemp**3) - (1.120083 * 1e-6* wtemp**4) + 
    (6.536336 * 1e-9 * wtemp**5)
  return(dens)
}

## this is our attempt for turbulence closure, estimating eddy diffusivity
eddy_diffusivity <-function(rho, depth, g, rho_0, ice, area, diffmethod){
  buoy = rep(1, (nx)) * 7e-5
  for (i in seq(1, nx-1)){#range(0, nx - 1):
    buoy[i] = ( abs(rho[i+1] - rho[i]) / (depth[i+1] - depth[i]) * g/rho_0 )
  }
  buoy[nx] = ( abs(rho[nx-1] - rho[nx]) / abs(depth[nx-1] - depth[nx]) * 
                     g/rho_0 )
  
  low_values_flags = buoy < 7e-5  # Where values are low
  buoy[low_values_flags] = 7e-5
  
  if ( diffmethod == 1){
    if (ice){
      ak <- 0.000898
    } else{
      ak <- 0.00706 *( max(area)/1E6)**(0.56)
    }
    
    kz = ak * (buoy)**(-0.43)
  } else if (diffmethod == 2){
    kzmax = 0.410 * sqrt(max(area)/1E6)
    kz = kzmax * ((7.5 * 10^(-5))/buoy)**(0.43)
  }

  return(kz)
}

provide_meteorology <- function(meteofile, secchifile,
                                pressurefile, windfactor = 1,
                                shortwavefactor = 1){
  meteo <-read_csv(meteofile)
  startDate <- meteo$Date[1]
  meteo$Cloud_Cover <- gotmtools::calc_cc(date = as.POSIXct(meteo$Date),
                                          airt = meteo$AirTemp,
                                          relh = meteo$RelHum,
                                          swr = meteo$ShortWave,
                                          lat = 43, lon = -89.41,
                                          elev =  280)
  meteo$lw <- gotmtools::calc_in_lwr(cc = meteo$Cloud_Cover,
                                     airt = meteo$AirTemp,
                                     relh = meteo$RelHum)
  meteo$ea <- (meteo$RelHum/100) * 10^(9.28603523 - 2322.37885/(meteo$AirTemp + 273.15))
  
  meteo$doy <- lubridate::yday(meteo$Date)
  meteo$dt <- as.POSIXct(meteo$Date) - (as.POSIXct(meteo$Date)[1]) + 1
  
  
  light <-read_csv(secchifile) %>% 
    filter(pond == 'B') %>%
    filter(doy >= meteo$doy[1])
  light$dt <- ((light$doy) - ((light$doy)[1]) + 1)* 86400
  light$kd <- 1.7 / light$secchi
  light$kd  <- zoo::na.approx(light$kd)
  
  pressure <-read_csv(pressurefile) %>% 
    mutate(doy = lubridate::yday(valid)) %>%
    filter(doy >= meteo$doy[1])
  pressure$dt <- ((pressure$doy) - ((pressure$doy)[1]) + 1)* 86400
  
  
  meteo$WindSpeed <- meteo$WindSpeed * windfactor
  meteo$ShortWave <- meteo$ShortWave * shortwavefactor
  
  return(list(meteo, light, pressure))
}

initial_profile <- function(initfile, nx, dx, processed_meteo){
  meteo <- processed_meteo
  startDate <- meteo$Date[1]
  obs <- read.csv(initfile)
  obs <- obs %>%
    filter(pond == 'B', site_id == 20) %>%
    mutate('doy' = lubridate::yday(datetime)) 
  init.time <- which(obs$datetime == paste0(startDate,' 12:00:00'))
  obs[init.time,]
  u = approx(obs$temp_depth_m[init.time], obs$temp_c[init.time],
             seq(0, nx * dx, length.out= nx))$y
  return(u)
}

get_hypsography <- function(hypsofile, dx, nx){
  hyps <- read_csv(hypsofile)
  area = approx(hyps$depths,hyps$areas,seq(1,nx*dx, 
                           length.out= nx))$y
  area[which.min(area)] <- 1e-2
  depth = depth= seq(1,nx*dx, length.out = nx)
  volume <- c(rev(diff(pracma::cumtrapz(area, depth))*(-1)),0)
  volume[which(volume == 0)] = min(volume[-which(volume == 0)])
  volume <- rep(0, (length(depth)-1))
  for (p in 1:length(volume)){
    volume[p] <- pracma::trapz(depth[p:(p+1)],area[p:(p+1)])
  }
  volume <- c(volume, 1000)
  return(list(area, depth, volume))
}

get_macrophyte <- function(canopyfile, biomassfile, processed_meteo){
  meteo <- processed_meteo
  
  macro <- read_csv(biomassfile) %>%
    filter(doy >= meteo$doy[1])
  macro$dt <- ((macro$doy) - ((macro$doy)[1]) + 1)* 86400
  macro$dt[1] = 1
  
  canpy <- read_csv(canopyfile) %>%
    filter(doy.f >= meteo$doy[1])
  if (canpy$doy.f[1] > meteo$doy[1]){
    canpy <- rbind(data.frame('doy.f' = meteo$doy[1],
                              'pond' = FALSE,
                              'mean_canopy' = canpy$mean_canopy[1]),
                   canpy)
  }
  canpy$dt <- ((canpy$doy.f) - ((canpy$doy.f)[1]) + 1)* 86400
  canpy$dt[1] = 1

  biomass <- read_csv('bc/biomass.csv')%>%
    filter(doy >= meteo$doy[1])
  if (biomass$doy[1] > meteo$doy[1]){
    biomass <- rbind(data.frame('doy' = meteo$doy[1],
                                'mean_biomass' = biomass$mean_biomass[1],
                                'pond' = FALSE),
                     biomass)
  }
  biomass$dt <- ((biomass$doy) - ((biomass$doy)[1]) + 1)* 86400
  biomass$dt[1] = 1
  
  return(list(canpy, biomass))
}

longwave <- function(cc, sigma, Tair, ea, emissivity, Jlw){  # longwave radiation into
  Tair = Tair + 273.15
  p <- (1.33 * ea/Tair)
  Ea <- 1.24 * (1 + 0.17 * cc**2) * p**(1/7)
  lw <- emissivity * Ea *sigma * Tair**4
  return(lw)
}
backscattering <- function(emissivity, sigma, Twater, eps){ # backscattering longwave 
  # radiation from the lake
  Twater = Twater + 273.15
  back = (eps * sigma * (Twater )^4) 
  return((-1) * back)
}
sensible <- function(p2, B, Tair, Twater, Uw){ # convection / sensible heat
  Twater = Twater + 273.15
  Tair = Tair + 273.15
  fu = 4.4 + 1.82 * Uw + 0.26 *(Twater - Tair)
  sensible <- ( p2 * B * fu * (Twater - Tair)) 
  return((-1) * sensible)
}
latent <- function(Tair, Twater, Uw, p2, pa, ea, RH){ # evaporation / latent heat
  Twater = Twater + 273.15
  Tair = Tair + 273.15
  Pressure = pa / 100
  fu = 4.4 + 1.82 * Uw + 0.26 *(Twater - Tair)
  fw = 0.61 * (1 + 10^(-6) * Pressure * (4.5 + 6 * 10^(-5) * Twater**2))
  ew = fw * 10 * ((0.7859+0.03477* Twater)/(1+0.00412* Twater))
  latent = fu * p2 * (ew - ea)# * 1.33) #* 1/6
  return((-1) * latent)
}

run_thermalmodel <- function(u, startTime, endTime, 
                             ice = FALSE, 
                             Hi = 0, 
                             iceT = 6, 
                             supercooled = 0,
                             scheme = 'explicit', 
                             kd_light = NULL,
                             densThresh = 1e-3, 
                             reflect = 0.3,
                             infra = 0.7, 
                             eps = 0.97, 
                             emissivity = 0.97,
                             sigma = 5.67 * 10^(-8), 
                             p2 = 1, 
                             B = 0.61,
                             g = 9.81,
                             Cd = 0.0013, # momentum coefficient (wind)
                             meltP = 1, 
                             dt_iceon_avg = 0.8,
                             Hgeo = 0.1, # geothermal heat
                             KEice = 1/1000, 
                             Ice_min = 0.1,
                             area, # area
                             depth, # depth
                             volume, # volume
                             zmax, # maximum lake depth
                             nx, # number of layers we will have
                             dt, # 24 hours times 60 min/hour times 60 seconds/min
                             dx,
                             meteo,
                             light,
                             pressure,
                             canpy,
                             biomass,
                             Cdplant = 1,
                             ahat = 0.4,
                             km = 0.04,
                             windfactor = 1,
                             shortwavefactor = 1,
                             diffusionfactor = 1,
                             diffmethod = 1,
                             rho_plant = 998.2){ # spatial step){
  
  ## linearization of driver data, so model can have dynamic step
  Jsw <- approxfun(x = meteo$dt, y = meteo$ShortWave * shortwavefactor, method = "linear", rule = 2)
  Jlw <- approxfun(x = meteo$dt, y = meteo$lw, method = "linear", rule = 2)
  Tair <- approxfun(x = meteo$dt, y = meteo$AirTemp, method = "linear", rule = 2)
  ea <- approxfun(x = meteo$dt, y = meteo$ea, method = "linear", rule = 2)
  Uw <- approxfun(x = meteo$dt, y = meteo$WindSpeed * windfactor, method = "linear", rule = 2)
  CC <- approxfun(x = meteo$dt, y = meteo$Cloud_Cover, method = "linear", rule = 2)
  Pa <- approxfun(x = pressure$dt, y = pressure$bpres_avg*100, method = "linear", rule = 2)
  kd <- approxfun(x = light$dt, y = light$kd, method = "constant", rule = 2)
  macroheight <- approxfun(x = canpy$dt, y = canpy$mean_canopy, method = "linear", rule = 2)
  macrobiomss <- approxfun(x = biomass$dt, y = biomass$mean_biomass, method = "linear", rule = 2)
  
  um <- matrix(NA, ncol =length( seq(startTime, endTime, dt)/dt), nrow = nx)
  kzm <- matrix(NA, ncol = length( seq(startTime, endTime, dt)/dt), nrow = nx)
  n2m <- matrix(NA, ncol = length( seq(startTime, endTime, dt)/dt), nrow = nx)
  mix <- rep(NA, length = length( seq(startTime, endTime, dt)/dt))#(floor(endTime/dt - startTime/dt)))
  therm.z <- rep(NA, length =length( seq(startTime, endTime, dt)/dt))
  mix.z <- rep(NA, length = length( seq(startTime, endTime, dt)/dt))
  Him <- rep(NA, length = length( seq(startTime, endTime, dt)/dt))
  macroz <- rep(NA, length = length( seq(startTime, endTime, dt)/dt))
  
  if (!is.null(kd_light)){
    kd <- approxfun(x = seq(1, endTime, 1), y = rep(kd_light, (endTime)), method = "constant", rule = 2)
  } 
  
  start.time <- Sys.time()
  ## modeling code for vertical 1D mixing and heat transport
  for (n in seq(startTime, endTime, dt)){#1:(floor(endTime/dt - startTime/dt))){  #iterate through time 1:floor(nt/dt)
    
    un = u # prior temperature values
    kz = eddy_diffusivity(calc_dens(un), depth, 9.81, 998.2, ice, area, diffmethod) / 86400 * diffusionfactor
    
    if (ice & Tair(n) <= 0){
      kzn = kz
      absorp = 1 - 0.7
      infra = 1 - absorp
    } else if (ice & Tair(n) >= 0){
      kzn = kz
      absorp = 1 - 0.3
      infra = 1 - absorp
    } else if (!ice) {
      kzn = kz   
      absorp = 1 - reflect# 0.3
      infra = 1 - absorp
    }
    kzm[, match(n, seq(startTime, endTime, dt))] <- kzn
    
    ## (1) Heat addition
    # surface heat flux
    Q <- (absorp * Jsw(n) + longwave(cc = CC(n), sigma = sigma, Tair = Tair(n), ea = ea(n), emissivity = emissivity, Jlw = Jlw(n)) + #longwave(emissivity = emissivity, Jlw = Jlw(n)) +
            backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
            latent(Tair = Tair(n), Twater = un[1], Uw = Uw(n), p2 = p2, pa = Pa(n), ea=ea(n), RH = RH(n)) + 
            sensible(p2 = p2, B = B, Tair = Tair(n), Twater = un[1], Uw = Uw(n)))  
    
    # heat addition over depth
    Hmacrophytes <- seq(dx, zmax, length.out = nx) 
    Hmacrophytes <- ifelse(Hmacrophytes < ((max(Hmacrophytes) - macroheight(n) )), 0, 1) #c(rep(0,2),rep(1,8)) # height of macrophytes (abundance)
    macroz[match(n, seq(startTime, endTime, dt))] = macroheight(n)
    
    P = macrobiomss(n)
    H =  (1- infra) * (Jsw(n))  * #
      exp(-(kd(n) + km * P * Hmacrophytes) *seq(dx,nx*dx,length.out=nx)) 
    
    Hg <- (area-lead(area))/dx * Hgeo/(4181 * calc_dens(un[1])) 
    Hg[which(is.na(Hg))] <- min(Hg, na.rm = TRUE)
    
    # add heat to all layers
    if (scheme == 'implicit'){
      ## (2) DIFFUSION
      # surface layer
      un[1] = un[1] +    Q * area[1]/(dx)*1/(4184 * calc_dens(un[1]) ) * dt/area[1]
      # # bottom layer
      un <- un  + (H * area/(dx) * 1/(4184 * calc_dens(un) ))* dt/area
      
      j <- length(kzn)
      y <- array(0, c(j,j))
      
      # all other layers in between
      # Linearized heat conservation equation matrix (diffusion only)
      az <- (dt/dx**2) * kzn                                         #coefficient for i-1
      cz <- (dt/dx**2) * kzn                #coefficient for i+1
      bz <- 1 + 2 * (dt/dx**2) * kzn                                                         #coefficient for i+1
      #Boundary conditions, surface
      az[1] <- 0
      #cz(1) remains unchanged
      bz[1]<- 1 #+ az[1] + (dt/dx**2) * kzn    
      #Boundary conditions, bottom
      #az(end) remains unchanged
      cz[length(cz)] <- 0
      bz[length(bz)] <- 1 #+ (dt/dx**2) * kzn    + cz[length(cz)]
      y[0 + 1:(j - 1) * (j + 1)] <- -cz[-length(bz)]	# superdiagonal
      y[1 + 0:(j - 1) * (j + 1)] <- bz	# diagonal
      y[2 + 0:(j - 2) * (j + 1)] <- -az[-1] 	# subdiagonal
      
      y[1,2] <- 0#- 2 * (dt/dx**2) * kzn[1]           
      y[nrow(y), (ncol(y)-1)] = 0#-2 * (dt/dx**2) * kzn[ncol(y)]           
      
      u <- solve(y, un)
    }
    
    ## (2) DIFFUSION
    # surface layer
    if (scheme == 'explicit'){
      u[1] = un[1] +
        (Q * area[1]/(dx)*1/(4184 * calc_dens(un[1]) ) +
           abs(H[1+1]-H[1]) * area[1]/(dx) * 1/(4184 * calc_dens(un[1]) ) +
           Hg[1]) * dt/area[1]
      
      # all other layers in between
      for (i in 2:(nx-1)){
        u[i] = un[i] +
          (area[i] * kzn[i] * 1 / dx**2 * (un[i+1] - 2 * un[i] + un[i-1]) +
             abs(H[i+1]-H[i]) * area[i]/(dx) * 1/(4184 * calc_dens(un[i]) ) +
             Hg[i])* dt/area[i]
      }
      # bottom layer
      u[nx] = un[nx] +
        abs(H[nx]-H[nx-1]) * area[nx]/(area[nx]*dx) * 1/(4181 * calc_dens(un[nx]) +
                                                           Hg[nx]/area[nx]) * dt
    }
    
    ## (3) TURBULENT MIXING OF MIXED LAYER
    # the mixed layer depth is determined for each time step by comparing kinetic 
    # energy available from wind and the potential energy required to completely 
    # mix the water column to a given depth
    Zcv <- seq(1, nx) %*% area / sum(area) # center of volume
    tau = 1.225 * Cd * Uw(n)^2 # wind shear is air density times wind velocity 
    if (Uw(n) <= 15) {
      c10 = 0.0005 * sqrt(Uw(n))
    } else {
      c10 = 0.0026
    }
    shear = sqrt((c10 * calc_dens(un[1]))/1.225) *  Uw(n) # shear velocity
    # coefficient times wind velocity squared
    KE = shear *  tau * dt # kinetic energy as function of wind
    
    if (ice){
      KE = KE * KEice
    }
    
    maxdep = 1
    for (dep in 1:(nx-1)){
      if (dep == 1){
        PE = abs(g *  ( seq(1,nx)[dep] - Zcv)  * calc_dens(u[dep]) * dx)
        # PE = abs(g *   seq(1,nx)[dep] *( seq(1,nx)[dep+1] - Zcv)  *
        #            abs(calc_dens(u[dep+1])- mean(calc_dens(u[1:dep]))))
        # DKE = Hmacrophytes[dep]*(calc_dens(u[dep])* ahat * Cdplant) *Uw(n)^3 * dt  *dx
        DKE = Hmacrophytes[dep]*(rho_plant* ahat * Cdplant) *Uw(n)^3 * dt  *dx
      } else {
        PEprior = PE
        PE = abs(g *  ( seq(1,nx)[dep] - Zcv)  * calc_dens(u[dep]) * dx +
                   PEprior)
        # PE = abs(g *   seq(1,nx)[dep] *( seq(1,nx)[dep+1] - Zcv)  *
                   # abs(calc_dens(u[dep+1])- mean(calc_dens(u[1:dep])))) + PEprior
        DKEprior = DKE
        # DKE = Hmacrophytes[dep]*(calc_dens(u[dep]) * ahat * Cdplant) *Uw(n)^3 * dt  *dx + DKEprior # rho_mp
        DKE = Hmacrophytes[dep]*(rho_plant * ahat * Cdplant) *Uw(n)^3 * dt  *dx + DKEprior # rho_mp
        KE = KE - DKE
      }
      if (PE > KE){
        maxdep = dep-1
        break
      } else if (dep>1 & PE < KE ){
        u[(dep-1):dep] = (u[(dep-1):dep] %*% volume[(dep-1):dep])/sum(volume[(dep-1):dep])
      }
      maxdep = dep
    }
    # u[1:maxdep] = (u[1:(maxdep)] %*% volume[1:(maxdep)])/sum(volume[1:(maxdep)]) #mean(u[1:maxdep])
    mix[match(n, seq(startTime, endTime, dt))] <- KE/PE #append(mix, KE/PE)
    therm.z[match(n, seq(startTime, endTime, dt))] <- maxdep #append(therm.z, maxdep)
    
    ## (4) DENSITY INSTABILITIES
    # convective overturn: Convective mixing is induced by an unstable density 
    # profile. All groups of water layers where the vertical density profile is 
    # unstable are mixed with the first stable layer below the unstable layer(s) 
    # (i.e., a layer volume weighed means of temperature and other variables are 
    # calculated for the mixed water column). This procedure is continued until 
    # the vertical density profile in the whole water column becomes neutral or stable.
    dens_u = calc_dens(u) 
    diff_dens_u <- (diff(dens_u)) 
    diff_dens_u[abs(diff(dens_u)) <= densThresh] = 0
    while (any(diff_dens_u < 0)){
      dens_u = calc_dens(u) 
      for (dep in 1:(nx-1)){
        if (dens_u[dep+1] < dens_u[dep] & abs(dens_u[dep+1] - dens_u[dep]) >= densThresh){
          u[dep:(dep+1)] = (u[dep:(dep+1)] %*% volume[dep:(dep+1)])/sum(volume[dep:(dep+1)]) #mean(u[dep:(dep+1)])
          break
        }
      }
      dens_u = calc_dens(u) 
      diff_dens_u <- (diff(dens_u)) 
      diff_dens_u[abs(diff(dens_u)) <= densThresh] = 0
    }
    
    dens_u_n2 = calc_dens(u) 
    n2 <- 9.81/mean(calc_dens(u)) * (lead(dens_u_n2) - lag(dens_u_n2))/dx
    max.n2 <- ifelse(max(n2, na.rm = T) > 1E-4, which.max(n2) * dx, dx * nx)
    mix.z[match(n, seq(startTime, endTime, dt))] <- max.n2
    
    
    
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
    iceT = iceT * (1 - x) + u[1] * x
    if ((iceT <= 0) == TRUE & Hi < Ice_min){
      # if (any(u <= 0) == TRUE){
      supercooled <- which(u < 0)
      initEnergy <- sum((0-u[supercooled])*area[supercooled] * dx * 4.18E6)
      
      if (ice != TRUE) {
        Hi <- Ice_min+(initEnergy/(910*333500))/max(area)
      } else {
        if (Tair(n) > 0){
          Tice <- 0
          Hi = Hi -max(c(0, meltP * dt*((absorp*Jsw(n))+(longwave(cc = CC(n), sigma = sigma, Tair = Tair(n), ea = ea(n), emissivity = emissivity, Jlw = Jlw(n)) +
                                                           backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
                                                           latent(Tair = Tair(n), Twater = un[1], Uw = Uw(n ), p2 = p2, pa = Pa(n), ea=ea(n),  RH = RH(n)) + 
                                                           sensible(p2 = p2, B = B, Tair = Tair(n), Twater = un[1], Uw = Uw(n))) )/(1000*333500)))
        } else {
          Tice <-  ((1/(10 * Hi)) * 0 +  Tair(n)) / (1 + (1/(10 * Hi))) 
          Hi <- min(Ice_min, sqrt(Hi**2 + 2 * 2.1/(910 * 333500)* (0 - Tice) * dt))
        }
      }
      ice = TRUE
      if (Hi > 0){
        u[supercooled] = 0
        u[1] = 0
      }
      Him[ match(n, seq(startTime, endTime, dt))] <- Hi
    } else if (ice == TRUE & Hi >= Ice_min) {
      if (Tair(n) > 0){
        Tice <- 0
        Hi = Hi -max(c(0, meltP * dt*((absorp*Jsw(n))+(backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
                                                         latent(Tair = Tair(n), Twater = un[1], Uw = Uw(n ), p2 = p2, pa = Pa(n), ea=ea(n*dt),  RH = RH(n)) + 
                                                         sensible(p2 = p2, B = B, Tair = Tair(n), Twater = un[1], Uw = Uw(n))) )/(1000*333500))) 
      } else {
        Tice <-  ((1/(10 * Hi)) * 0 +  Tair(n*dt)) / (1 + (1/(10 * Hi))) 
        Hi <- min(Ice_min, sqrt(Hi**2 + 2 * 2.1/(910 * 333500)* (0 - Tice) * dt))
      }
      u[supercooled] = 0
      u[1] = 0
      Him[ match(n, seq(startTime, endTime, dt))] <- Hi
    } else if (ice == TRUE & Hi < Ice_min){
      ice = FALSE 
    }
    
    n2m[, match(n, seq(startTime, endTime, dt))] <- n2
    um[, match(n, seq(startTime, endTime, dt))] <- u
    
    
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  # time =  startDate + seq(1, ncol(um))*dt#/24/3600
  
  df.sim <- data.frame(cbind(seq(startTime,endTime,dt), t(um)) )
  colnames(df.sim) <- c("datetime", as.character(paste0('wtemp.',seq(1,nrow(um))*dx)))
  df.sim$datetime <-   seq(startTime, endTime, dt)#/24/3600
  
  ## averaged responses
  bf.sim <- apply(df.sim[,-1], 1, function(x) rLakeAnalyzer::buoyancy.freq(wtr = x, depths = seq(1,nrow(um))*dx))
  
  bf.sim <- apply(df.sim[,-1], 1, function(x) rLakeAnalyzer::center.buoyancy(wtr = x, depths = seq(1,nrow(um))*dx))
  
  
  # z.bf.sim <- apply(bf.sim,2, function(x) which.max(x))
  
  df.z.df.sim <- data.frame('time' = df.sim$datetime, 'z' = bf.sim)
  
  avg.epi.sim <- NULL
  avg.hyp.sim <- NULL
  avg.tot.sim <- NULL
  for (j in 1:nrow(df.z.df.sim)){
    d = df.sim[,-1]
    if (is.na(df.z.df.sim$z[j])){
      df.z.df.sim$z[j] = 1
    }
    avg.epi.sim <- append(avg.epi.sim,((as.numeric(d[j,1:df.z.df.sim$z[j]], na.rm = T) %*% 
                                          area[1:df.z.df.sim$z[j]] )/ 
                                         sum(area[1:df.z.df.sim$z[j]])))
    avg.hyp.sim <- append(avg.hyp.sim,((as.numeric(d[j,df.z.df.sim$z[j]:ncol(d)], na.rm = T)%*% 
                                          area[df.z.df.sim$z[j]:ncol(d)] )/ 
                                         sum(area[df.z.df.sim$z[j]:ncol(d)])))
    avg.tot.sim <- append(avg.tot.sim,((as.numeric(d[j,1:ncol(d)], na.rm = T)%*% 
                                          area[1:ncol(d)] )/ 
                                         sum(area[1:ncol(d)])))
  }
  
  stratFlag = rep(NA, length = ncol(um))
  for (v in 1:length(stratFlag)){
    stratFlag[v] = ifelse((calc_dens(um[nrow(um),v]) - calc_dens(um[1,v])) >= 0.1 &
                            mean(um[,v]) >= 4, 1, 0)
  }
  
  
  df.avg.sim <- data.frame('time' = df.sim$datetime,
                           'epi' = avg.epi.sim,
                           'hyp' = avg.hyp.sim,
                           'tot' = avg.tot.sim,
                           'stratFlag' = stratFlag,
                           'thermoclineDep' = bf.sim)
  
  return(list('temp'  = um,
              'diff' = kzm,
              'mixing' = mix,
              'buoyancy' = n2m,
              'icethickness' = Hi,
              'iceflag' = ice,
              'icemovAvg' = iceT,
              'supercooled' = supercooled,
              'mixingdepth' = mix.z,
              'thermoclinedepth' = therm.z,
              'macroheight' = macroz,
              'endtime' = endTime,
              'average' = df.avg.sim))
}

wrapper_scales <- function(x, lb, ub){
  y <-  lb+(ub-lb)/(10)*(x)
  return(y)
}
