rm(list=ls())
source("/scratch/os22/dc0105/AWRA/codes/awra_run_SA_funcs.R")
source("/scratch/os22/dc0105/AWRA/codes/LAI_composite.R")
source("/scratch/os22/dc0105/AWRA/codes/forcing_extract.R") 
source("/scratch/os22/dc0105/AWRA/codes/AWRAL_TimeStep_funcs.R")

library(sensitivity)

PROJ_LATLON = '+proj=longlat +datum=WGS84'

forcingpath <-"/scratch/os22/dc0105/AWRA/forcings/"
storepath <-"/scratch/os22/dc0105/AWRA/Outputs/"
streampath <- "/scratch/os22/dc0105/AWRA/Streamdata/"
Sobolpath <- "/scratch/os22/dc0105/AWRA/"

## Accept argument about the catchment name and the grid number 

catstring <- commandArgs(trailingOnly = TRUE)
catnum <- strsplit(as.character(catstring),"-")
catname = (as.character(catnum[[1]][1]))
gridnum = as.numeric(catnum[[1]][2])

      # forcing
dates_1321 <- seq(as.Date("2013-01-01"),as.Date("2021-12-31"),1)
start_SA = which(dates_1321 == '2014-01-01')   # spin up starts in 2014
end_SA = which(dates_1321 == '2018-12-31')

all_forcing_grid1 <- grid_forcing_extract(catname,paste0("grid",gridnum),
                                          start =  start_SA, end = end_SA )

STATIC_PAR <- read.csv(paste0(forcingpath,catname,"static.csv"))

static_grid1 <- STATIC_PAR[,paste0("grid",gridnum)]
static_par_grid1 <- list(ftree = static_grid1[1], S0FC = static_grid1[2],
                         SsFC = static_grid1[3],  SdFC = static_grid1[4],
                         meanPET = static_grid1[5])
                   
  ## Extract LAI data 

dates_1518 = DateSeq(2015,2018)

  ## create dates of composite (monthly here)
LAI_monthly_dates <- as.Date(Monthly_Dates(2015,2018)) # monthly

# LAI_5km_monthly <- as.matrix(read.csv(paste0(forcingpath,catname,"monthlyLAI_1519.csv")))

    ## MOD15A2H 2012-21 monthly data 
LAI_5km_monthly <- as.matrix(read.csv(paste0(forcingpath,catname,"MOD15A2H_monthly1221.csv")))

 # monthly LAI data from 2015 to 2018
# LAI_5km_1519_monthly <- LAI_5km_monthly[1:which(LAI_monthly_dates=="2019-12-01"),]

LAI_5km_1518_monthly <- LAI_5km_monthly[37:84,]
LAI_grid1 = LAI_5km_1518_monthly[,paste0("grid",gridnum)]

     ######  Sobol sensitivity analysis ########

X1_grid1<- hydromad::parameterSets(AWRA_par_ranges,samples = 20000, method="latin.hypercube")
X2_grid1 <- hydromad::parameterSets(AWRA_par_ranges,samples = 20000, method="latin.hypercube")
X_all_grid1 <- paramat_sobol(X1_grid1,X2_grid1)

batchsize <- nrow(X_all_grid1)/40
sobol_response <- c()

for (i in 1:40) {
  sobol_iterate <- list(model = Run_AWRAL_SA_LAI, X1 = X1_grid1,
                        X2 = X2_grid1, nboot = 100, conf = 0.95,
                        X = X_all_grid1[((i-1)*batchsize+1):(i*batchsize), ],
                        call = match.call())
  class(sobol_iterate) <- "sobol2002"
  
  response(sobol_iterate,all_forcing=all_forcing_grid1,
           LAI_data = LAI_grid1,
           static_par = static_par_grid1,
           Simdates = dates_1518,
           tiffDates = LAI_monthly_dates )

  sobol_response <- c(sobol_response,sobol_iterate$y)
  
  write.table(sobol_iterate$y,paste0(storepath,catname,"grid1_R_MOD152AH.txt"),append = T,row.names = F,col.names=F,quote=F)
  if (i%%10==0) {cat(paste0("batch",i),sobol_iterate$y[1],"\n")}
}

sobol_grid1_LAI <- list(model = Run_AWRAL_SA_LAI, X1 = X1_grid1,
                      X2 = X2_grid1, nboot = 200, conf = 0.90,
                      X = X_all_grid1,
                      call = match.call())

class(sobol_grid1_LAI) <- "sobol2002"
sobol_grid1_LAI$y <- sobol_response

# compute TSI & FSI by the tell function #### 

tell(sobol_grid1_LAI)       

names_AWRA_par = names(AWRA_par_ranges)
Sobol_grid1_LAI_result<- data.frame(Par = names_AWRA_par,
                                  TSI = sobol_grid1_LAI$T$original,
                                  FSI = sobol_grid1_LAI$S$original)

write.csv(Sobol_grid1_LAI_result,paste0(Sobolpath,catname,"Sobol_grid1_R_MOD152AH.csv"),quote=F,row.names = F)



