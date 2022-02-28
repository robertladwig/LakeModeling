## python implementation of gotmtools::calc_cc
def calc_cc(date, airt, relh = None, dewt = None, swr, lat, long, elev, daily = False: 
    if daily == True:
        date = pd.date_range(start=date[0], end=(date[-1] + 23*60*60), freq='1H') # TODO: ended here; need to check
        seq.POSIXt(from = date[1], to = (date[length(date)] + 23 * 60 * 60), by = "1 hour") # delete this after checking above works
    yday <- yday(date)
    hour <- hour(date)
    hour[hour == 0] <- 24
    std.mer = seq(-90, 90, 15)
    Lsm = std.mer[which.min(abs(long - std.mer))]
    Hsc = 1390
    cd = 0.06
    Rg = 0.045
    theta = lat * pi/180
    r = 1 + 0.017 * cos((2 * pi/365) * (186 - yday))
    d = 23.45 * pi/180 * cos((2 * pi/365) * (172 - yday))
    dts = (1/15) * (Lsm - long)
    value = (sin(theta) * sin(d))
    value = value/(cos(theta) * cos(d))
    tss = (12/pi) * acos(-value) + dts + 12
    tsu = -tss + (2 * dts) + 24
    gamma = rep(0, length(tss))
    dum = which(hour > tsu & hour < tss)
    gamma[dum] = 1
    dum1 = which(hour <= 12)
    dum2 = which(hour > 12)
    hb1 = pi/12 * (hour - 1 - dts)
    hb1[dum1] = hb1[dum1] + pi
    hb1[dum2] = hb1[dum2] - pi
    hb = hb1
    dum3 = which(hb1 > 2 * pi)
    hb[dum3] = hb[dum3] - 2 * pi
    dum4 = which(hb1 < 0)
    hb[dum4] = hb[dum4] + 2 * pi
    he1 = pi/12 * (hour - dts)
    he1[dum1] = he1[dum1] + pi
    he1[dum2] = he1[dum2] - pi
    he = he1
    dum3 = which(he1 > 2 * pi)
    he[dum3] = he[dum3] - 2 * pi
    dum4 = which(he1 < 0)
    he[dum4] = he[dum4] + 2 * pi
    Ho = Hsc/(r^2) * (sin(theta) * sin(d) + 12/pi * cos(theta) * 
        cos(d) * (sin(he) - sin(hb))) * gamma
    w = (he + hb)/2
    alpha1 = abs(sin(theta) * sin(d) + cos(theta) * cos(d) * 
        cos(w))
    alpha = atan(alpha1/sqrt(1 - alpha1^2))
    theta_am1 = ((288 - 0.0065 * elev)/288)^5.256
    theta_am2 = sin(alpha) + 0.15 * ((alpha * 180/pi) + 3.855)^(-1.253)
    theta_am = theta_am1/theta_am2
    if dewt is None:
        dewt <- 243.04 * (log(relh/100) + ((17.625 * airt)/(243.04 + 
            airt)))/(17.625 - log(relh/100) - ((17.625 * airt)/(243.04 + 
            airt)))
    if daily == True:
        dewt = rep(dewt, each = 24)
    Pwc = 0.85 * exp(0.11 + 0.0614 * dewt)
    a2 = exp(-(0.465 + 0.134 * Pwc) * (0.179 + 0.421 * exp(-0.721 * 
        theta_am)) * theta_am)
    a1 = exp(-(0.465 + 0.134 * Pwc) * (0.129 + 0.171 * exp(-0.88 * 
        theta_am)) * theta_am)
    at = (a2 + 0.5 * (1 - a1 - cd))/(1 - 0.5 * Rg * (1 - a1 - 
        cd))
    Ho = at * Ho
    dum5 = which(Ho < 0)
    Ho[dum5] = 1
    df = data.frame(DateTime = date, Ho = Ho)
    if (daily == TRUE) {
        df = aggregate(list(Ho = df$Ho), by = list(DateTime = cut(df[, 
            1], "1 day")), mean, na.rm = T)
    }
    df$swr = swr
    df$ccsim <- NA
    for (i in 1:nrow(df)) {
        if (df$Ho[i] < df$swr[i]) {
            df$ccsim[i] <- NaN
        }
        else {
            df$ccsim[i] <- sqrt((1 - (df$swr[i]/df$Ho[i]))/0.65)
        }
    }
    ccsim = df$ccsim
    ccsim[ccsim > 1] <- 1
    sta = min(which(!is.nan(ccsim)))
    stp = max(which(!is.nan(ccsim)))
    ccsim[sta:stp] <- zoo::na.approx(ccsim[sta:stp])
    if (sta != 1) {
        ccsim[1:sta] <- ccsim[sta]
    }
    if (stp != length(ccsim)) {
        ccsim[stp:length(ccsim)] <- ccsim[stp]
    }
    return(ccsim)
}
