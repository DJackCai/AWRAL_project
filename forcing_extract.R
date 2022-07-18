 # forcingpath: has forcing data and static parameters

forcingpath <-"/scratch/os22/dc0105/AWRA/forcings/"
# forcingpath <- "./forcings/" 

maxmin_res <-  function(x)
{
  return((x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)))
}

static_extract <- function(catname,gridnum) {
  STATIC_PAR <- as.matrix(read.csv(paste0(forcingpath,catname,"static.csv")))
  static_grid <- STATIC_PAR[,gridnum]
  static_par_grid <- list(ftree = static_grid[1], S0FC = static_grid[2],
                            SsFC = static_grid[3],  SdFC = static_grid[4],
                            meanPET = static_grid[5])
}

# Set default running period to be 2014-2019 (3 April)

dates_1321 <- seq(as.Date("2013-01-01"),as.Date("2021-12-31"),1)
start_calib = which(dates_1321 == '2014-01-01')   # spin up starts in 2014

end_Vali = which(dates_1321 == '2018-12-31')   
end_SMAP = which(dates_1321 == '2019-12-31')   

        #### Forcing for the specific grid ###### 

grid_forcing_extract <- function(catname, gridnum,start = start_calib,end = end_Vali ) {
  
   # gridnum: string like "grid11" 
   ######### Important: for grids, remove the YYYY-MM-DD columns (1:3) #######
  
	PRECIP <- as.matrix(read.csv(paste0(forcingpath, catname,"rain_day_resampled_201321.csv"),header = T))[start:end,-c(1:3)]
  RAD    = as.matrix(read.csv(paste0(forcingpath,catname, "solar_exposure_day_resampled_201321.csv"),header = T))[start:end,-c(1:3)]
  TMIN   = as.matrix(read.csv(paste0(forcingpath, catname,"temp_min_day_resampled_201321.csv")))[start:end,-c(1:3)]
  TMAX   = as.matrix(read.csv(paste0(forcingpath, catname, "temp_max_day_resampled_201321.csv")))[start:end,-c(1:3)]
  U2 <- as.matrix(read.csv(paste0(forcingpath, catname,"wind_resampled_201321.csv")))[start:end,-c(1:3)]
  
  ## extract for each grid 
  PG_grid1 = PRECIP[, gridnum]
  
  RG_grid1 = RAD[, gridnum]
  RG_grid1[RG_grid1 < 0.1] = 0.1
  RG_grid1 = 11.57*RG_grid1
  
  TA_grid1 = TMIN[,gridnum] + 0.75*(TMAX[,gridnum]-TMIN[,gridnum])
  PE_grid1 = 610.8 * exp(17.27*TMIN[,gridnum]/(237.3+TMIN[,gridnum]))
  
  PAIR_grid1 = 97500.00
  U2_grid1<- U2[,gridnum]
   
  U2_grid1  = U2_grid1 *  (1-(1-0.5)*0.25)/(0.5)
    
  # entire time period of different variables of forcing for the grid 
  all_forcing_grid1 <- list(PG=PG_grid1,RG=RG_grid1,TA=TA_grid1,
                            PE=PE_grid1,U2=U2_grid1)
  
  return(all_forcing_grid1)
  
}


grid_forcing_catchment <- function(catname, N_GRID,start = start_calib,end = end_Vali) {
  # N_GRID: number of grid surrounding and including the catchment 
  # catname: string like "GoulburnRiver_", notice the underscore
  
  PRECIP <- as.matrix(read.csv(paste0(forcingpath, catname,"rain_day_resampled_201321.csv"),header = T))[start:end,]
  PG = PRECIP[,4:(N_GRID+3)]
  RAD    = as.matrix(read.csv(paste0(forcingpath,catname, "solar_exposure_day_resampled_201321.csv"),header = T))[start:end,]
  RG = RAD[,4:(N_GRID+3)]
  RG[RG < 0.1] = 0.1
  RG = 11.57*RG
  TMIN   = as.matrix(read.csv(paste0(forcingpath, catname,"temp_min_day_resampled_201321.csv")))[start:end,]
  TMAX   = as.matrix(read.csv(paste0(forcingpath, catname, "temp_max_day_resampled_201321.csv")))[start:end,]
  TA   = TMIN[,4:(N_GRID+3)] + 0.75*(TMAX[,4:(N_GRID+3)] - TMIN[,4:(N_GRID+3)])
 
  PE   = 610.8 * exp(17.27*TMIN[,4:(N_GRID+3)]/(237.3+TMIN[,4:(N_GRID+3)]))
  PAIR = rep(97500.00,N_GRID)
  
  U2 <- as.matrix(read.csv(paste0(forcingpath, catname,"wind_resampled_201321.csv")))[start:end,]
  U2 <- U2[,4:(N_GRID+3)]
  # convert to effective daytime windspeed
   U2 = (1-(1-0.5)*0.25)/(0.5) * U2
  
  # entire time period of different variables of forcing for the grid 
  
  all_forcing <- list(PG = PG, RG = RG, TA = TA, PE = PE, U2 = U2, PAIR = PAIR)
  return(all_forcing)
}


      #### grid extract from 1994 to 2014  ######

grid_forcing_catchment_longterm <- function(catname, N_GRID,start ,end  ) {
  
  # N_GRID: number of grid surrounding and including the catchment 
  # catname: string like "GoulburnRiver_", notice the underscore
  
  PRECIP <- as.matrix(read.csv(paste0(forcingpath, catname,"rain_day_resampled_19922014.csv"),header = T))[start:end,]
  PG = PRECIP[,4:(N_GRID+3)]
  RAD    = as.matrix(read.csv(paste0(forcingpath,catname, "solar_exposure_day_resampled_19922014.csv"),header = T))[start:end,]
  RG = RAD[,4:(N_GRID+3)]
  RG[RG < 0.1] = 0.1
  RG = 11.57*RG
  TMIN   = as.matrix(read.csv(paste0(forcingpath, catname,"temp_min_day_resampled_19922014.csv")))[start:end,]
  TMAX   = as.matrix(read.csv(paste0(forcingpath, catname, "temp_max_day_resampled_19922014.csv")))[start:end,]
  TA   = TMIN[,4:(N_GRID+3)] + 0.75*(TMAX[,4:(N_GRID+3)] - TMIN[,4:(N_GRID+3)])
  
  PE   = 610.8 * exp(17.27*TMIN[,4:(N_GRID+3)]/(237.3+TMIN[,4:(N_GRID+3)]))
  PAIR = rep(97500.00,N_GRID)
  
  U2 <- as.matrix(read.csv(paste0(forcingpath, catname,"wind_resampled_19922014.csv")))[start:end,]
  U2 <- U2[,4:(N_GRID+3)]
  
  # convert to effective daytime windspeed
  U2 = ((1-(1-0.5)*0.25)/(0.5)) * U2
  
  ### add: gap filling the NA values to make sure the full series can run 
  
  for (i in 1:N_GRID) {
    na_PG_index = which(is.na(PG[,i])==T)
    PG[na_PG_index,i]  = median(PG[,i], na.rm=T)
    na_RG_index = which(is.na(RG[,i])==T)
    RG[na_RG_index,i]  = median(RG[,i], na.rm=T)
    na_TMIN_index = which(is.na(TMIN[,i])==T)
    TMIN[na_TMIN_index,i] = median(TMIN[,i], na.rm=T)
    na_TMAX_index = which(is.na(TMAX[,i])==T)
    TMAX[na_TMAX_index,i] = median(TMAX[,i], na.rm=T)
    na_U2_index = which(is.na(U2[,i])==T)
    U2[na_U2_index,i] = median(U2[,i], na.rm=T)
  }
  
  # entire time period of different variables of forcing for the grid 
  
  all_forcing = list(PG = PG, RG = RG, TA = TA, PE = PE, U2 = U2, PAIR = PAIR)
  return(all_forcing)
  
}



