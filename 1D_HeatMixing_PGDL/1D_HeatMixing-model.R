#' Created on Thu Aug 19 11:29:34 2021
#' 
#' @author: Robert Ladwig
#' @email: rladwig2@wisc.edu

## remove everything from workspace
rm(list = ls())

# set wd to current dir of script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## colors for plotting
library(tidyverse)
library(RColorBrewer)
library(patchwork)

source('1D_HeatMixing_functions.R')

## lake configurations
zmax = 25 # maximum lake depth
nx = 25 # number of layers we will have
dt = 3600  # 24 hours times 60 min/hour times 60 seconds/min
dx = zmax/nx # spatial step

## area and depth values of our lake 
hyps_all <- get_hypsography(hypsofile = 'bc/LakeEnsemblR_bathymetry_standard.csv',
                            dx = dx, nx = nx)

## atmospheric boundary conditions
meteo_all <- provide_meteorology(meteofile = 'bc/LakeEnsemblR_meteo_standard.csv',
                    secchifile = 'bc/light.csv', 
                    windfactor = 0.8)

## here we define our initial profile
u_ini <- initial_profile(initfile = 'bc/obs.txt', nx = nx, dx = dx,
                     depth = hyps_all[[2]],
                     processed_meteo = meteo_all[[1]])

### EXAMPLE RUNS
hydrodynamic_timestep = 24 * 3600 #24/4 * dt
total_runtime <- 365
startingDate <- meteo_all[[1]]$datetime[1]

temp <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
              nrow = nx)
avgtemp <- matrix(NA, ncol = 6,
                nrow = total_runtime * hydrodynamic_timestep/ dt)
temp_diff <- matrix(NA, ncol = (total_runtime * hydrodynamic_timestep/ dt) ,
               nrow = nx)
temp_mix <- matrix(NA, ncol = (total_runtime * hydrodynamic_timestep/ dt) ,
                    nrow = nx)
temp_conv <- matrix(NA, ncol = (total_runtime * hydrodynamic_timestep/ dt) ,
                    nrow = nx)
temp_ice <- matrix(NA, ncol = (total_runtime * hydrodynamic_timestep/ dt) ,
                    nrow = nx)
diff <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
               nrow = nx)
meteo_output <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
               nrow = 9)
buoyancy <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
               nrow = nx)

# do the meteo interpolation all at the start
meteo = get_interp_drivers(meteo_all=meteo_all, total_runtime=total_runtime, 
                           hydrodynamic_timestep=hydrodynamic_timestep, dt=dt, method="integrate")

if (exists('res')) {remove('res')}

for (i in 1:total_runtime){
  if (exists('res')){
    u = res$temp[, ncol(res$temp)]
    startTime = res$endtime
    endTime =  res$endtime + hydrodynamic_timestep -1
    ice = res$icefla
    Hi = res$icethickness 
    iceT = res$icemovAvg
    supercooled = res$supercooled
    kd_light = NULL
    matrix_range = max(1, (startTime/dt)):((endTime/dt)) # this returns floats, not ints, after first round through? seems to cause issue down below in setting up avgtime
    matrix_range_start = max(1, round(startTime/dt) + 1)
    matrix_range_end = round(endTime/dt)
  } else {
    u = u_ini
    startTime = 1
    endTime = hydrodynamic_timestep - 1
    ice = FALSE
    Hi = 0
    iceT = 6
    supercooled = 0
    kd_light = NULL
    matrix_range = max(1, (startTime/dt)):((endTime/dt)+1)
    matrix_range_start = max(1, round(startTime/dt))
    matrix_range_end = round(endTime/dt)
  }
  res <-  run_thermalmodel(u = u, 
                            startTime = startTime, 
                            endTime =  endTime, 
                            ice = ice, 
                            Hi = Hi, 
                            iceT = iceT,
                            supercooled = supercooled,
                            kd_light = kd_light,
                            zmax = zmax,
                            nx = nx,
                            dt = dt,
                            dx = dx,
                            area = hyps_all[[1]], # area
                            depth = hyps_all[[2]], # depth
                            volume = hyps_all[[3]], # volume
                            daily_meteo = meteo[,matrix_range_start:matrix_range_end],
                            # secview = meteo_all[[2]],
                            Cd = 0.0008,
                            pgdl_mode = 'on',
                           scheme = 'implicit')

  temp[, matrix_range_start:matrix_range_end] =  res$temp
  diff[, matrix_range_start:matrix_range_end] =  res$diff
  avgtemp[matrix_range_start:matrix_range_end,] <- as.matrix(res$average)
  temp_diff[, matrix_range_start:matrix_range_end] =  res$temp_diff
  temp_mix[, matrix_range_start:matrix_range_end] =  res$temp_mix
  temp_conv[, matrix_range_start:matrix_range_end] =  res$temp_conv
  temp_ice[, matrix_range_start:matrix_range_end] =  res$temp_ice
  buoyancy[, matrix_range_start:matrix_range_end] =  res$temp
  meteo_output[, matrix_range_start:matrix_range_end] =  res$meteo_input
  buoyancy[, matrix_range_start:matrix_range_end] = res$buoyancy_pgdl
  
  average <- res$average %>%
    mutate(datetime = as.POSIXct(startingDate + time),
           Date = as.Date(datetime, format = "%m/%d/%Y")) %>%
    summarise_all(mean)
  
}


# plotting for checking model output and performance
plot(seq(1, ncol(temp))*dt/24/3600, temp[1,], col = 'red', type = 'l', 
     xlab = 'Time (d)', ylab='Temperature (degC)', ylim=c(-1,35), lwd = 2)
for (i in 2:nx){
  lines(seq(1, ncol(temp))*dt/24/3600, temp[i,], 
        lty = 'dashed',lwd =2)
}

time =  startingDate + seq(1, ncol(temp), 1) * dt
avgtemp = as.data.frame(avgtemp)
colnames(avgtemp) = c('time', 'epi', 'hyp', 'tot', 'stratFlag', 'thermoclineDep')

ggplot(avgtemp) +
  geom_line(aes(time, epi, col = 'epilimnion')) +
  geom_line(aes(time, hyp, col = 'hypolimnion')) +
  geom_line(aes(time, tot, col = 'total')) +
  theme_minimal() 
ggplot(avgtemp) +
  geom_line(aes(time, stratFlag)) +
  theme_minimal() 
ggplot(avgtemp) +
  geom_line(aes(time, thermoclineDep)) +
  geom_point(aes(time, stratFlag, col = as.factor(stratFlag))) +
  scale_y_reverse() +
  theme_minimal() 

df <- data.frame(cbind(time, t(temp)) )
colnames(df) <- c("time", as.character(paste0(seq(1,nrow(temp)))))
m.df <- reshape2::melt(df, "time")
m.df$time <- time

ggplot(m.df, aes((time), as.numeric(variable))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(-2,35),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Temp [degC]')+
  scale_y_reverse() 

# save data for PGDL
df <- data.frame(cbind(time, t(temp)) )
colnames(df) <- c("time", as.character(paste0('tempDegC_total04_',seq(1,nrow(temp)))))
df$time <- time
df <- df %>%
  filter(time >=  '2009-06-04 09:00:00' & time <= '2009-09-01 00:00:00')
write.csv(df, file = 'output/temp_total04.csv', row.names = F)

df <- data.frame(cbind(time, t(temp_diff)) )
colnames(df) <- c("time", as.character(paste0('tempDegC_diff01_',seq(1,nrow(temp)))))
df$time <- time
df <- df %>%
  filter(time >=  '2009-06-04 09:00:00' & time <= '2009-09-01 00:00:00')
write.csv(df, file = 'output/temp_diff01.csv', row.names = F)

df <- data.frame(cbind(time, t(temp_mix)) )
colnames(df) <- c("time", as.character(paste0('tempDegC_mix02_',seq(1,nrow(temp)))))
df$time <- time
df <- df %>%
  filter(time >=  '2009-06-04 09:00:00' & time <= '2009-09-01 00:00:00')
write.csv(df, file = 'output/temp_mix02.csv', row.names = F)

df <- data.frame(cbind(time, t(temp_conv)) )
colnames(df) <- c("time", as.character(paste0('tempDegC_conv03_',seq(1,nrow(temp)))))
df$time <- time
df <- df %>%
  filter(time >=  '2009-06-04 09:00:00' & time <= '2009-09-01 00:00:00')
write.csv(df, file = 'output/temp_conv03.csv', row.names = F)

df <- data.frame(cbind(time, t(diff)) )
colnames(df) <- c("time", as.character(paste0('diffM2s-1_',seq(1,nrow(temp)))))
df$time <- time
df <- df %>%
  filter(time >=  '2009-06-04 09:00:00' & time <= '2009-09-01 00:00:00')
write.csv(df, file = 'output/diff.csv', row.names = F)

df <- data.frame(cbind(time, t(buoyancy)) )
colnames(df) <- c("time", as.character(paste0('n2S-2_',seq(1,nrow(temp)))))
df$time <- time
df <- df %>%
  filter(time >=  '2009-06-04 09:00:00' & time <= '2009-09-01 00:00:00')
write.csv(df, file = 'output/buoyancy.csv', row.names = F)

df <- data.frame(cbind(time, t(meteo_output)) )
colnames(df) <- c("time", "AirTemp_degC", "Longwave_Wm-2",
                  "Latent_Wm-2", "Sensible_Wm-2", "Shortwave_Wm-2",
                  "lightExtinct_m-1","ShearVelocity_mS-1", "ShearStress_Nm-2",
                  "Area_m2")
df$time <- time
df <- df %>%
  filter(time >=  '2009-06-04 09:00:00' & time <= '2009-09-01 00:00:00')
write.csv(df, file = 'output/meteorology_input.csv', row.names = F)



# observed data
# Package ID: knb-lter-ntl.130.29 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: High Frequency Water Temperature Data - Lake  Mendota Buoy 2006 - current.
inUrl2  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/130/29/63d0587cf326e83f57b054bf2ad0f7fe" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")

dt2 <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "sampledate",     
                 "year4",     
                 "month",     
                 "daynum",     
                 "hour",     
                 "depth",     
                 "wtemp",     
                 "flag_wtemp"    ), check.names=TRUE)

unlink(infile2)

# attempting to convert dt2$sampledate dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp2sampledate<-as.Date(dt2$sampledate,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp2sampledate) == length(tmp2sampledate[!is.na(tmp2sampledate)])){dt2$sampledate <- tmp2sampledate } else {print("Date conversion failed for dt2$sampledate. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp2sampledate) 
if (class(dt2$year4)=="factor") dt2$year4 <-as.numeric(levels(dt2$year4))[as.integer(dt2$year4) ]               
if (class(dt2$year4)=="character") dt2$year4 <-as.numeric(dt2$year4)
if (class(dt2$month)=="factor") dt2$month <-as.numeric(levels(dt2$month))[as.integer(dt2$month) ]               
if (class(dt2$month)=="character") dt2$month <-as.numeric(dt2$month)
if (class(dt2$daynum)=="factor") dt2$daynum <-as.numeric(levels(dt2$daynum))[as.integer(dt2$daynum) ]               
if (class(dt2$daynum)=="character") dt2$daynum <-as.numeric(dt2$daynum)
if (class(dt2$depth)=="factor") dt2$depth <-as.numeric(levels(dt2$depth))[as.integer(dt2$depth) ]               
if (class(dt2$depth)=="character") dt2$depth <-as.numeric(dt2$depth)
if (class(dt2$wtemp)=="factor") dt2$wtemp <-as.numeric(levels(dt2$wtemp))[as.integer(dt2$wtemp) ]               
if (class(dt2$wtemp)=="character") dt2$wtemp <-as.numeric(dt2$wtemp)
if (class(dt2$flag_wtemp)!="factor") dt2$flag_wtemp<- as.factor(dt2$flag_wtemp)


dt2$bhour <- ifelse(dt2$hour %/% 100 >= 1, dt2$hour/100, dt2$hour)
dt2$datetime <- as.POSIXct(paste0(dt2$sampledate,' ',dt2$bhour,':00:00'), format = "%Y-%m-%d %H:%M:%S")

obs <- dt2 %>%
  select(datetime, depth, wtemp)

# idx <- na.omit(match(as.POSIXct(df.sim$datetime), as.POSIXct(obs$datetime) ))
idx <- (match(as.POSIXct(obs$datetime), as.POSIXct(df$time) ))

# df.sim <- df.sim[idx, ]
obs <- obs[which(!is.na(idx)), ]

idz <- which(obs$depth %in% seq(0,24,1))
obs = obs[idz,]

obs <- data.frame(obs)
obs$depth <- factor(obs$depth)

wide.obs <- reshape(obs, idvar = "datetime", timevar = "depth", direction = "wide")

m.wide.obs <- reshape2::melt(wide.obs, "datetime")
m.wide.obs$time <- wide.obs$time

ggplot(m.wide.obs, aes((datetime), as.numeric(variable))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(-2,35),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Temp [degC]')+
  scale_y_reverse() 


df <- data.frame(cbind(time, t(temp)) )
colnames(df) <- c("time", as.character(paste0('tempDegC_total04_',seq(1,nrow(temp)))))
df$time <- time
df <- df %>%
  filter(time >=  '2009-06-04 09:00:00' & time <= '2009-09-01 00:00:00')

dat0 = wide.obs[,-1]
dat0[dat0 < 5] = NA
dat = zoo::na.approx(dat0) # apply(as.matrix(wide.obs[,-1]), 1, function(x) zoo::na.approx(x))
# dat2 = do.call(rbind, dat)
# dat2 = matrix(unlist(dat), ncol = 21, byrow = T)
dat3 = apply(as.matrix(dat), 1, function(x) approx(seq(0,20,1),x,seq(1,25,1), method = 'linear', rule=2)$y)
# dat3 = apply(as.matrix(wide.obs[,-1]), 1, function(x) zoo::na.approx(x, x = seq(0,20,1),x_out =seq(1,25,1)))
dat4 = t(dat3)
dat5 = apply(as.matrix(dat4), 2, function(x) approx(wide.obs$datetime - wide.obs$datetime[1], x,
                                                    df$time-df$time[1], method = 'linear', rule=2)$y)

dat.df <- data.frame(cbind(df$time, dat5))
colnames(dat.df) <- c("time", as.character(paste0('tempDegC_total04_',seq(1,nrow(temp)))))
dat.df$time <- df$time
dat.df <- dat.df %>%
  filter(time >=  '2009-06-04 09:00:00' & time <= '2009-09-01 00:00:00')

m.dat.df <- reshape2::melt(dat.df, "time")
m.dat.df$time <- df$time

ggplot(m.dat.df, aes((time), as.numeric(variable))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(-2,35),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Temp [degC]')+
  scale_y_reverse() 

write.csv(dat.df, file = 'output/observed_temp.csv', row.names = F)


library(ncdf4)
library(gotmtools)
# GOTM
# GOTMr::run_gotm('alt_lakeModels/GOTM/')
nc <- nc_open(file.path('alt_lakeModels/', 
                               "GOTM", "output", "output.nc"))
nc_close(nc)
diff <- gotmtools::get_vari(ncdf = file.path('alt_lakeModels/', 
                                  "GOTM", "output", "output.nc"), var = "avh",
                 print = FALSE)
z <- gotmtools::get_vari(ncdf = file.path('alt_lakeModels/', 
                               "GOTM", "output", "output.nc"), var = "z",
              print = FALSE)

# Add in obs depths which are not in depths and less than mean depth
depths <- seq(0, min(z[1, -1]), by = -1 )
obs_dep_neg <- seq(1,nrow(temp)) * (-1)
add_deps <- obs_dep_neg[!(obs_dep_neg %in% depths)]
depths <- c(add_deps, depths)
depths <- depths[order(-depths)]

message("Interpolating GOTM temp to include obs depths... ",
        paste0("[", Sys.time(), "]"))
got <- gotmtools::setmodDepths(diff, z, depths = depths, print = T)
message("Finished interpolating! ",
        paste0("[", Sys.time(), "]"))

got <- maditr::dcast(got, date ~ depths)
got <- got[, c(1, (27:2))]
str_depths <- abs(as.numeric(colnames(got)[2:ncol(got)]))
colnames(got) <- c("datetime", paste("diffM2s-2_", str_depths, sep = ""))

df <- as.data.frame(got) %>%
  filter(datetime >=  '2009-06-04 09:00:00' & datetime <= '2009-09-01 00:00:00')
write.csv(df, file = 'output/diff_gotm.csv', row.names = F)

# Simstrat
SimstratR::run_simstrat('alt_lakeModels/Simstrat/')

diff <- read.table(file.path('alt_lakeModels', "Simstrat", "output", "nuh_out.dat"), header = TRUE,
                   sep = ",", check.names = FALSE)
diff[, 1] <- as.POSIXct(diff[, 1] * 3600 * 24 + 6* 3600, origin = paste0("1995-01-01"))
# In case sub-hourly time steps are used, rounding might be necessary
# diff[, 1] <- round_date(diff[, 1], unit = seconds_to_period(timestep))

# First column datetime, then depth from shallow to deep
diff <- diff[, c(1, ncol(diff):2)]

# Remove columns without any value
diff <- diff[, colSums(is.na(diff)) < nrow(diff)]

# Add in obs depths which are not in depths and less than mean depth
mod_depths <- as.numeric(colnames(diff)[-1])
depths <- seq(0, min(z[1, -1]), by = -1 )
obs_dep_neg <- seq(1,nrow(temp)) * (-1)
add_deps <- obs_dep_neg[!(obs_dep_neg %in% mod_depths)]
depths <- c(add_deps, mod_depths)
depths <- depths[order(-depths)]

if(length(depths) != (ncol(diff) - 1)){
  message("Interpolating Simstrat temp to include obs depths... ",
          paste0("[", Sys.time(), "]"))
  
  
  # Create empty matrix and interpolate to new depths
  wat_mat <- matrix(NA, nrow = nrow(diff), ncol = length(depths))
  for(i in seq_len(nrow(diff))) {
    y <- as.vector(unlist(diff[i, -1]))
    wat_mat[i, ] <- approx(mod_depths, y, depths, rule = 2)$y
  }
  message("Finished interpolating! ",
          paste0("[", Sys.time(), "]"))
  df <- data.frame(wat_mat)
  df$datetime <- diff[, 1]
  df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
  colnames(df) <- c("datetime", paste0("wtr_", abs(depths)))
  diff <- df
}else{
  # Set column headers
  str_depths <- abs(as.numeric(colnames(diff)[2:ncol(diff)]))
  colnames(diff) <- c("datetime", paste0("diffM2s-2_", str_depths))
}

df <- diff %>%
  filter(datetime >=  '2009-06-04 09:00:00' & datetime <= '2009-09-01 00:00:00')
write.csv(df, file = 'output/diff_simstrat.csv', row.names = F)
