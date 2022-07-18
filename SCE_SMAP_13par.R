source("/scratch/os22/dc0105/AWRA/codes/AWRAL_TimeStep_funcs.R")
library(hydromad)
source("/scratch/os22/dc0105/AWRA/codes/Qsim_Eval.R")
source("/scratch/os22/dc0105/AWRA/codes/Run_AWRAL_Calib.R")
source("/scratch/os22/dc0105/AWRA/codes/forcing_extract.R") 

source("/scratch/os22/dc0105/AWRA/codes/Run_AWRAL_calib_SM_13par_grid.R")

# Accept the argument for grid number 
gridnum <-commandArgs(trailingOnly = TRUE)
gridnum<-as.numeric(c(gridnum)[1])

forcingpath <-"/scratch/os22/dc0105/AWRA/forcings/"
storepath <- "/scratch/os22/dc0105/AWRA/"
 
# catname = "415237"
# catname = "612230"

# catname = "137201"
# catname = "143110"

# catname = "401217"
#  catname = "227237"

# catname = "606195"
#  catname = "419053"
#  catname = "405229"
#  catname = "136301"

# catname = "419035"
# catname = "238235"

# catname = "607144"
# catname = "606001"

# catname = "142001"

# catname = "401230"

# catname = "403217"
# catname = "401217"

# catname = "403222"

# catname = "401212"
# catname = "403232"
#  catname = "401210"
# catname = "405264"
# catname = "403244"
# catname = "408202"


 # catname = "121002"

# catname = "136203"
# catname = "130319"
# catname  = "422319"

# catname = "416305"

# catname = "422338"
# catname  = "206014"
# catname = "110003"
# catname = "229650"
# catname = "402206"

# catname =  "926002"

# catname = "405227"
# catname = "405218"
# catname = "237200"
# catname = "497"
# catname = "226222"
# catname = "206025"

# catname = "116001"
# catname = "116004"
catname = "116006"
# catname = "918003"
# catname = "925001"
# catname = "927001" 
# catname = "224213"
# catname = "181"


## SCE control par
control.SCE = list(trace = 1, ncomplex = 7, maxit =2)
control.SCE = modifyList(list(fnscale = -1), control.SCE) 
if (isTRUE(hydromad.getOption("trace"))){
  control$trace <- 1 } 
bestModel <- NULL
bestFunVal <- Inf*control.SCE$fnscale


### parameter range 
Vc1 = c(0.05,0.6);     Gfrac_max2 = c(0.2,0.6) ; 
Us02 = c(2,10)
Ud01 = c(0.1,7) ;      FsoilEmax1   = c(0.1,0.9) ; FsoilEmax2  = c(0.1,0.9)
S_sls1 = c(0.05,0.45) ; FdrainFC1 = c(0.005,0.1) ; FdrainFC2 = c(0.005,0.1)
#S_sls2 = c(0.03,0.9); ER_frac_ref1 = c(0.02,0.28) ; 
beta1 = c(0.7,8.4) ; beta2 = c(0.7,8.4)

LAImax2 = c(5.5,9)
Tgrow2 = c(40,200)

parlist<-list(Vc1=Vc1, Gfrac_max2=Gfrac_max2,
              Us02=Us02,
              Ud01 = Ud01,
              FsoilEmax1   = FsoilEmax1,
              FsoilEmax2   = FsoilEmax2,
              S_sls1 = S_sls1,  
             # S_sls2 = S_sls2, 
             # ER_frac_ref1 = ER_frac_ref1,
              FdrainFC1 = FdrainFC1,    FdrainFC2 = FdrainFC2,
              beta1 = beta1,     beta2 = beta2,
              LAImax2=LAImax2,    Tgrow2 = Tgrow2)
stopifnot(length(parlist)==13)

lower <-sapply(parlist, min) ; upper <- sapply(parlist, max)
initpars <-sapply(parlist,mean)

 
## do_sce objfun calculation #####

do_sce_SMAP_13par <-function(pars,SMAPts,all_forcing,STATIC_PAR){
  
  # input: par vector (to be updated at each step)
  # all_forcing: list of entire forcing, for the specific grid 
  # SMAPts: SMAPts for the grid 
  # STATIC_PAR: static parameter
  
  update_Vc1 <-pars[1];  update_Gfrac_max2 <- pars[2] 
  update_Us02 <- pars[3];   update_Ud01 <- pars[4];
  update_FsoilEmax1 <- pars[5];  update_FsoilEmax2 <- pars[6]; 
  update_S_sls1 <- pars[7];
  #update_S_sls2 <- pars[7]  ; update_ER_frac_ref1 <- pars[8]
  update_FdrainFC1<- pars[8];   update_FdrainFC2<- pars[9];
  update_beta1 <- pars[10];  update_beta2 <- pars[11]; 
  update_LAImax2 <- pars[12];  update_Tgrow2 <- pars[13]
  
  FORCING = empty_forcing() 
  
  ## Run model - call the run function that takes par vector as inputs
  thisMod <- Run_AWRAL_calib_SM_13par_grid(pars=pars,all_forcing = all_forcing,
                                      STATIC_PAR=STATIC_PAR, N_GRID = 1 )
  
  S0mean <- thisMod$S0mean[-c(1:365),]
  
  SM_sim <- c(S0mean[-1],NA)
  SMAP_obs <- SMAPts
  stopifnot(length(SM_sim)==length(SMAP_obs))
  
  ## compute the objective function ###
  cor_SMAP <- cor(SM_sim,SMAP_obs,use="p")
  thisVal<- cor_SMAP
  
  if (isTRUE(thisVal*control.SCE$fnscale < bestFunVal*control.SCE$fnscale)) {
    bestModel <<- thisMod
    bestFunVal <<- thisVal}
  return(thisVal)  # function to be optimised
}


    ### Run SCE #####

dates_1321 <- seq(as.Date("2013-01-01"),as.Date("2021-12-31"),1)
start_calib = which(dates_1321 == '2014-01-01')   # spin up starts in 2014
end_calib = which(dates_1321 == '2018-12-31')   

STATIC_PAR <- as.matrix(read.csv(paste0(forcingpath,catname,'_resampled_static.csv')))
static_grid <- STATIC_PAR[,paste0("grid",gridnum)]
static_par_grid <- list(ftree = static_grid[1], S0FC = static_grid[2],
                        SsFC = static_grid[3],  SdFC = static_grid[4],
                        meanPET = static_grid[5])

all_forcing_grid <- grid_forcing_extract(catname,paste0("grid",gridnum),
                                         start =  start_calib, end = end_calib )

          ### SMAP data from 2015 to 2018

SMAP1519 <- read.csv(paste0(forcingpath, catname,"_resampled_SMAP1519.csv"))
SMAP_grid_all <- SMAP1519[,paste0("grid",gridnum)]

dates_SMAP <- seq(as.Date("2015-01-01"),as.Date("2019-12-31"),1)
SMAP_end <- which(dates_SMAP == '2018-12-31') 
SMAP_grid <- SMAP_grid_all[1:SMAP_end]


set.seed(258)  
SCE_run <-SCEoptim(do_sce_SMAP_13par, par = initpars,
                   # parameters for running the AWRA and compute objfun
                   all_forcing =  all_forcing_grid, 
                   STATIC_PAR =  static_par_grid,
                   SMAPts =  SMAP_grid,
                   lower = lower, upper = upper, 
                   control = control.SCE)

par_opt <- SCE_run$par
r_opt = SCE_run$value

# store the optimised parameters 

write.table(cbind(paste0("grid",gridnum),t(par_opt)),
            file=paste0(storepath, catname,"grid",gridnum,'_resamp_13parSMAP_1518.txt'),
            row.names=FALSE, quote=F,col.names = c("grid",names(parlist)))

write.table(cbind(paste0("grid", gridnum),t(r_opt)),
            file=paste0(storepath, catname,"grid",gridnum,'_RSMAP_1518.txt'),
            row.names=FALSE, quote=F,col.names = c("grid","Objfun"))

print(paste0("Optimisation for grid",gridnum, " is finished"))
