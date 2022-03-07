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

source('1D_HeatMixing_functions_IBM.R')

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
diff <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
               nrow = nx)
nind = 1e4
individuals <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
               nrow = nind)
tracers <- matrix(NA, ncol = total_runtime * hydrodynamic_timestep/ dt,
                      nrow = nx)
avgtemp <- matrix(NA, ncol = 6,
                nrow = total_runtime * hydrodynamic_timestep/ dt)
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
    agents = res$agents[, ncol(res$agents)]
    tracer = res$tracer[, ncol(res$tracer)]
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
    agents = c(rep(10,nind))
    tracer = u_ini * 1e-6
    tracer[10] = 100
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
                            Cd = 0.0008,
                            agents = agents,
                           tracer = tracer)

  temp[, matrix_range] =  res$temp
  diff[, matrix_range] =  res$diff
  individuals[, matrix_range] =  res$agents
  tracers[, matrix_range] =  res$tracer
  avgtemp[matrix_range,] <- as.matrix(res$average)
  
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

time =  seq(1, ncol(temp), 1)
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

ggplot(m.df, aes((time), as.numeric(variable))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(-2,35),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Temp [degC]')+
  scale_y_reverse() 

df.tracer <- data.frame(cbind(time, t(tracers)) )
colnames(df.tracer) <- c("time", as.character(paste0(seq(1,nrow(tracers)))))
m.df.tracer <- reshape2::melt(df.tracer, "time")

plot(colSums(tracers))

ggplot(m.df.tracer, aes((time), as.numeric(variable))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(0,100),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Tracer [-]')+
  scale_y_reverse() 

df.ind <- data.frame(cbind(time, t(individuals)) )
colnames(df.ind) <- c("time", as.character(paste0(seq(from = 1e-7,  by =1e-8, length.out = nind))))
m.df.ind <- reshape2::melt(df.ind, "time")

df.diff <- data.frame(cbind(time, t(diff)) )
colnames(df.diff) <- c("time", as.character(paste0(seq(1,nrow(diff)))))
m.df.diff <- reshape2::melt(df.diff, "time")
library(gganimate)
g1<-ggplot(m.df) +
  geom_path(aes(value, as.numeric(as.character(variable)))) +
  theme_minimal()  +xlab('Temp [degC]') +
  ylab('Depth') + 
  scale_y_reverse() 

## vertical temperature profiles
for (i in seq(1,ncol(temp), length.out = 200)){
  n = i
  i = floor(i)
  
  sim = m.df.diff %>% 
    filter(time == time[i]) %>%
    mutate(depth =  as.numeric(as.character(variable)))
  inds = m.df.ind %>% 
    filter(time == time[i]) %>%
    mutate(depth =  as.numeric(as.character(variable)))

  g1<- ggplot() +
    
    # facet_wrap(~ factor(variable, level = c(paste0('wtemp.',seq(0,24,1)))), scales = 'free') +
    geom_point(data = inds, aes(depth, value, col = variable), size = 3) +
    geom_path(data = sim, aes(value, 
                              depth), size = 1.2) +
    xlab('Kz (m2 s-1)') + ylab('depth (m)')+
    scale_y_reverse() +
    ggtitle( time[i]) + 
    labs(col='') +
    xlim(1e-7, 1e-4) +
    theme_light() +
    theme(legend.position = "none") 
  
  g2<-ggplot(inds, aes(y = value)) + 
    geom_density()+ 
    scale_y_reverse(limits=c(nx,0)) +
    theme_light() +
    ylab('')+
    theme(legend.position = "none") 
  
  g=g1 + g2 + plot_layout(widths = c(4, 1)); g
  
  ggsave(paste0('../../animation_ibm/pic_',match(n, seq(1,ncol(temp),length.out=200)),'.png'),
         width = 4, height = 5, units = 'in')
  
}
  

