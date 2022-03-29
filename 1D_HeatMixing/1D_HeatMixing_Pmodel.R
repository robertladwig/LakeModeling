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
dt = 3600 # 24 hours times 60 min/hour times 60 seconds/min
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
hydrodynamic_timestep = 24 * dt
total_runtime <- 365
startingDate <- meteo_all[[1]]$datetime[1]

temp <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
              nrow = nx)
avgtemp <- matrix(NA, ncol = 6,
                nrow = total_runtime * hydrodynamic_timestep/ dt)
thermoclinedep <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
               nrow = 1)
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
    matrix_range = max(1, (startTime/dt)):((endTime/dt))
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
                            daily_meteo = meteo_all[[1]],
                            secview = meteo_all[[2]],
                            Cd = 0.0008)

  temp[, matrix_range] =  res$temp
  avgtemp[matrix_range,] <- as.matrix(res$average)
  thermoclinedep[, matrix_range] =  res$thermoclinedepth
  
  average <- res$average %>%
    mutate(datetime = as.POSIXct(startingDate + time),
           Date = as.Date(datetime, format = "%m/%d/%Y")) %>%
    # group_by(datetime) %>%
    summarise_all(mean)
  
  ## run C-P-O2 model with input from ''average''
  ## derive kd value and put this in as input for ''run_thermalmodel'', e.g. 
  ## '' kd <- 0.5 '' and then change L 59 to '' kd_light = kd '' 
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
  geom_line(aes(time, stratFlag, col = as.factor(stratFlag))) +
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
  xlab('') +
  scale_y_reverse() 






















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

obs <- data.frame(obs)
obs$depth <- factor(obs$depth)

wide.obs <- reshape(obs, idvar = "datetime", timevar = "depth", direction = "wide")
m.obs <- reshape2::melt(wide.obs, id = 'datetime')
m.obs$datetime <- as.POSIXct(m.obs$datetime)
m.obs$group <- 'obs'

m.obs$variable <-  factor(m.obs$variable, levels=paste0('wtemp.',seq(0,24,1)))

## averaged responses
bf.obs <- apply(wide.obs[,-1], 1, function(x) rLakeAnalyzer::buoyancy.freq(wtr = x, depths = as.numeric(unique(obs$depth))))

z.bf.obs <- apply(bf.obs,2, function(x) which.max(x))

df.z.df.obs <- data.frame('time' = wide.obs$datetime, 'z' = as.numeric(z.bf.obs))
df.z.df.sim <- data.frame('time' = time, 'z' =avgtemp$thermoclineDep)

df.z.df.obs.f = df.z.df.obs %>%
  filter(time >= "2009-01-01 00:00:00" & time <= '2009-12-31 23:00:00')

g.therm <- ggplot() +
  geom_line(data = df.z.df.obs.f,
            aes(time, z, col = 'observed')) +
  geom_line(data = df.z.df.sim,
            aes(time, z, col = 'simulated')) +
  scale_y_reverse() + xlab('') + ylab('Thermocline depth (m)') +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(0.10, 0.15),
        legend.background = element_rect(fill = "white", color = "black")); g.therm

areas = hyps_all[[1]]
areas <- c(max(areas), areas)
area = approx(seq(0,nx,1), areas,
            c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10,
              11, 12, 13, 14, 15, 16, 17, 18, 19,20))$y

avg.epi.obs <- NULL
avg.hyp.obs <- NULL
avg.tot.obs <- NULL
for (j in 1:nrow(df.z.df.obs)){
  d = wide.obs[,-1]
  if (is.na(df.z.df.obs$z[j])){
    df.z.df.obs$z[j] = 1
  }
  if (any(is.na(as.numeric(d[j,1:df.z.df.obs$z[j]], na.rm = T))) == T){
    a = which(is.na(as.numeric(d[j,1:df.z.df.obs$z[j]], na.rm = T)))
    avg.epi.obs <- append(avg.epi.obs,((na.omit(as.numeric(d[j,1:df.z.df.obs$z[j]], na.rm = T)) %*% 
                                          area[1:df.z.df.obs$z[j]][-a] )/ 
                                         sum(area[1:df.z.df.obs$z[j]][-a])))
  } else{
    avg.epi.obs <- append(avg.epi.obs,((na.omit(as.numeric(d[j,1:df.z.df.obs$z[j]], na.rm = T)) %*% 
                                          area[1:df.z.df.obs$z[j]] )/ 
                                         sum(area[1:df.z.df.obs$z[j]])))
  }
 
 
  b =which(is.na( as.numeric(d[j,df.z.df.obs$z[j]:ncol(d)], na.rm = T)))
  avg.hyp.obs <- append(avg.hyp.obs,((na.omit(as.numeric(d[j,df.z.df.obs$z[j]:ncol(d)], na.rm = T))%*% 
                                        area[df.z.df.obs$z[j]:ncol(d)][-b] )/ 
                                       sum(area[df.z.df.obs$z[j]:ncol(d)][-b])))
  avg.tot.obs <- append(avg.tot.obs,((as.numeric(d[j,1:ncol(d)], na.rm = T)%*% 
                                        area[1:ncol(d)] )/ 
                                       sum(area[1:ncol(d)])))
}


df.avg.obs <- data.frame('time' = wide.obs$datetime,
                         'epi' = avg.epi.obs,
                         'hyp' = avg.hyp.obs,
                         'tot' = avg.tot.obs,
                         'type' = 'obs')

df.avg.obs.f = df.avg.obs %>%
  filter(time >= "2009-01-01 00:00:00" & time <= '2009-12-31 23:00:00')

avgtemp$time <- time
g.avg <- ggplot(avgtemp) +
  geom_line(aes(time, epi, col = 'simulated epi'), alpha = 0.7) +
  geom_line(aes(time, hyp, col = 'simulated hyp'), alpha = 0.7) +
  # geom_line(aes(time, tot, col = 'simulated hyp'), alpha = 0.7)  +
  geom_point(data = df.avg.obs.f,
             aes(time, epi, col = 'observed epi'), alpha = 0.7) +
  geom_point(data = df.avg.obs.f,
             aes(time, hyp, col = 'observed hyp'), alpha = 0.7) +

  xlab('') + ylab('Area weighed water temp. (deg C)') +
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = c(0.87, 0.85),
        legend.background = element_rect(fill = "white", color = "black")); g.avg



