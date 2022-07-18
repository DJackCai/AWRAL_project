### Function to create model-equivalent composite to match MODIS data

library(lubridate)
library(signal)

DateSeq <- function(year1, year2){
  dateseq<- seq(as.Date(paste0(year1,"-01-01")), as.Date(paste0(year2,"-12-31")),1)
  return(dateseq)
}

Monthly_Dates <- function(year1, year2) {
  # Extract first day of month 
  firstdates <- unique(ymd(format(dates_seq(year1,year2), "%Y-%m-01")))
  return(firstdates)
}


### create the sequence for dates of MODIS composite ######
MODIS_DATES <- function(year1, year2) {
  dates <- NULL
  yearseq <- seq(year1, year2, 1)
  for (i in 1:length(yearseq)) {
    dates_add <- seq(as.Date(paste0(yearseq[i],"-01-01")),as.Date(paste0(yearseq[i],"-12-31")),8)
    dates <- c(dates,as.character(dates_add))
  }
  return(dates)
}


Sim_LAI_composite <- function(Simdates,tiffDates,LAIsim_ts) {
  # Simdates: time steps in daily simulations
  # tiffDates: dates of product
  # LAIsim_ts: full time series of LAI simulations

  LAI_sim_comp <- NULL
  LAI_tiff_index <- c(which(Simdates%in%tiffDates))

  for (i in 1:length(LAI_tiff_index)) {

    if (i!=length(LAI_tiff_index)) {
      LAI_sim = LAIsim_ts[LAI_tiff_index[i]:(LAI_tiff_index[i+1]-1)]
      LAI_comp = mean(LAI_sim,na.rm=T)
      LAI_sim_comp = c(LAI_sim_comp,LAI_comp)
    } else {
      LAI_sim = LAIsim_ts[ (LAI_tiff_index[i]) : (length(LAIsim_ts)) ]
      LAI_comp = mean(LAI_sim,na.rm=T)
      LAI_sim_comp = c(LAI_sim_comp,LAI_comp)
    }
  }

  return(LAI_sim_comp)
}


      #### 10 April - Smoothed & interpolated leaf area index ##### 

Interp_SG <- function(x,Composite_dates, full_dates, n=11) {
  # x: series to be filtered
  # Composite_dates: dates of composite (assumed the first day)
  # full_dates: complete interpolated series

  if (is.vector(x)) {
    x_filter <- sgolayfilt(x,p=3,n=n)
    x_interp <- spline(x = Composite_dates, y = x_filter, method="natural",xout = full_dates)$y
  } else if (is.matrix(x)) {
  x_filter <- apply(x,2,function(x)
  {sgolayfilt(x,p = 3,n = n)} )

  x_interp <- apply(x_filter,2,function(x) {
    spline(x = Composite_dates, y = x, method="natural",xout = full_dates)$y
  } ) } else if (is.data.frame(x)) {
    x_filter <- apply(x,2,function(x)
    {sgolayfilt(x,p = 3,n = n)} )

    x_interp <- apply(x_filter,2,function(x) {
      spline(x = Composite_dates, y = x, method="natural",xout = full_dates)$y
    } )
  }

  return (x_interp)
}


