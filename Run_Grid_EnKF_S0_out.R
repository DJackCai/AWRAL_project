forcingpath <-"/scratch/os22/dc0105/AWRA/forcings/"
storepath <- "/scratch/os22/dc0105/AWRA/"

# forcingpath = "./forcings/"

 ##### Set-up functions #######
source("/scratch/os22/dc0105/AWRA/codes/Qsim_Eval.R")
source("/scratch/os22/dc0105/AWRA/codes/Run_AWRAL_Calib.R")
source("/scratch/os22/dc0105/AWRA/codes/forcing_extract.R")
source("/scratch/os22/dc0105/AWRA/codes/Run_AWRAL_calib_Q_grid.R")
source("/scratch/os22/dc0105/AWRA/codes/AWRAL_TimeStep.R")
source("/scratch/os22/dc0105/AWRA/codes/PARAMS_set.R")
source("/scratch/os22/dc0105/AWRA/codes/Compute_wetness.R")
source("/scratch/os22/dc0105/AWRA/codes/AWRAL_Anci_Funcs.R")


    ### 7.18: use S0 to avoid the issue of wetness (directly MVM) ####


Run_Grid_EnKF_S0_out = function(IF, catname, Ne, gridnum, par_mat, 
                             start_date = "2014-01-01", end_date= "2018-12-31", store_length = 1461, 
                             SMAP_Err = NULL, rescaling = "MVM",
                             Rain_Err = 0.6 , Temp_Err = 2.5 , 
                             Rad_Err = 50) {
  
  ### Compute for a specific IF #####
  
  ## PARAM set & N_GRID ####
  N_GRID = 1
  static_grid <- as.matrix(read.csv(paste0(forcingpath,catname,"_resampled_static.csv")))[,paste0("grid",gridnum)]
  static_par_grid <- list(ftree = static_grid[1], S0FC = static_grid[2],
                          SsFC = static_grid[3],  SdFC = static_grid[4],
                          meanPET = static_grid[5])
  
  dates_1321 = dates_seq(2013,2021)
  start_run = which(dates_1321 == as.Date(start_date))   # spin up starts in 2014
  end_run = which(dates_1321 == as.Date(end_date)) 
  
  ## grid-specific forcing for all time  #####
  all_forcing <- grid_forcing_extract(catname = catname, gridnum = gridnum ,
                                      start = start_run, end= end_run)
  PG <- all_forcing$PG;    RG <- all_forcing$RG
  TA <- all_forcing$TA;    PE <- all_forcing$PE
  U2 <- all_forcing$U2;    PAIR = rep(97500.00,1)
  
  ## Calibrated parameters  ######
  PARAMS = PARAMS_OPT_CATCH(N_GRID = N_GRID, par_mat =  par_mat, 
                            STATIC_PAR = static_par_grid)
  
  meanPET  = static_par_grid$meanPET
  meanP   =  mean(all_forcing$PG)
  humidity = meanP/meanPET
  Sgref = pmax(0.0,  8.15 * meanP^2.34)
  # FdrainFC = pmax(pmin(0.0685*rbind(humidity,humidity)^3.179,0.3),0.005)
  K_gw   = 0.047*humidity^-0.0508
  K_rout = 0.141*meanPET + 0.284
  
  PARAMS$Sgref    = Sgref
  #  PARAMS$FdrainFC = FdrainFC
  PARAMS$K_gw     = K_gw
  PARAMS$K_rout   = K_rout
  
  print("-------- Optimised parameters have been finished -----------")
  
  #### baseline model run - calibrated model #####
  
  STATES = initial_states(PARAMS = PARAMS,N_GRID = N_GRID, meanP = meanP)
  mod_run = Run_AWRAL_calib_Q_grid(pars = par_mat, all_forcing = all_forcing,
                                   STATIC_PAR = static_par_grid, N_GRID = 1)
  S0mean_sim = mod_run$S0mean[-c(1:365),]
  
  
   ### SMAP to S0, via repeated bias correction #####
  SMAP_data  = as.matrix(read.csv(paste0(forcingpath,catname, "_resampled_SMAP1519.csv")))[1:1461,paste0("grid",gridnum)]
  SMAP_shift = dplyr::lag(SMAP_data, n=1)
  
  if (rescaling == "CDF" ) {
    SMAP_scaled  = mapply(function(x,y) {
      SMAP_cdfmatch <- CDFt(y,x,x)
      return(SMAP_cdfmatch$DS) }, as.data.frame(SMAP_shift),
      as.data.frame(S0mean_sim)) } else if (rescaling == "MVM" ) {
        SMAP_scaled  = mapply(function(x,y) { return(Mean_var_match(x,y)) },
                                      as.data.frame(SMAP_shift), as.data.frame(S0mean_sim))
      }
   
  Fc = static_grid[2]
  SMAP_scaled[SMAP_scaled >= Fc] = Fc
  SMAP_scaled[SMAP_scaled < 0.0] = 0.0
  
  ## repeated MVM to avoid remaining truncation bias ###
  
  for (s in 1:10) {
    SMAP_scaled = mapply(function(x,y) { return(Mean_var_match(x,y)) },
                                 as.data.frame(SMAP_scaled), as.data.frame(S0mean_sim))
    SMAP_scaled[SMAP_scaled < 0.0] = 0.0
    SMAP_scaled[SMAP_scaled >= Fc] = Fc
  }
  stopifnot(mean(SMAP_scaled, na.rm=T)-mean(S0mean_sim, na.rm=T) <=0.01 )
  
  print("-------- Pre-processing of SMAP data have been finished -----------")
  
  ## Transformed SMAP error variance after rescaling ####
  
  if (is.null(SMAP_Err)) {
    sd_AWRA = sd(S0mean_sim,na.rm=T)
    sd_SMAP = sd(SMAP_shift,na.rm=T)
    SMAP_Err = 0.04 * (sd_AWRA/sd_SMAP)
  } else {
    # if there is pre-defined error inputs
    stopifnot(ncol(SMAP_shift) == length(SMAP_Err))
    SMAP_Err = SMAP_Err }
  
  ##### EnKF module #####
  
  STATES = initial_states(PARAMS = PARAMS, N_GRID = N_GRID, meanP = meanP)
  NL <- length(PG)  # grid-version 
  FORCING = empty_forcing()
  
  # 1) state variable matrix, evolving at each time ######
  
  #  - forecast state variables
  Xf_S01    = matrix(NA,Ne,N_GRID)  # hru 1 = DR; hru 2 = SR
  Xf_S02    = matrix(NA,Ne,N_GRID)
  Xf_Ss1    = matrix(NA,Ne,N_GRID)
  Xf_Ss2    = matrix(NA,Ne,N_GRID)
  Xf_Sd1    = matrix(NA,Ne,N_GRID)
  Xf_Sd2    = matrix(NA,Ne,N_GRID)
  Xf_Sg     = matrix(NA,Ne,N_GRID)
  Xf_Sr     = matrix(NA,Ne,N_GRID)
  Xf_Mleaf1 = matrix(NA,Ne,N_GRID)
  Xf_Mleaf2 = matrix(NA,Ne,N_GRID)
  w0_f = w0_a = matrix(NA,Ne,N_GRID)
  
  # modelled observation from the ensemble of forecast states
  hXf_w0   = matrix(NA,Ne,N_GRID)
  Xf_S0mean = matrix(NA,Ne,N_GRID)
  QTOT_f = matrix(NA,Ne,N_GRID)
  
  #  analysis state variables
  Xa_S01    = matrix(NA,Ne,N_GRID)
  Xa_S02    = matrix(NA,Ne,N_GRID)
  Xa_Ss1    = matrix(NA,Ne,N_GRID)
  Xa_Ss2    = matrix(NA,Ne,N_GRID)
  Xa_Sd1    = matrix(NA,Ne,N_GRID)
  Xa_Sd2    = matrix(NA,Ne,N_GRID)
  Xa_Sg     = matrix(NA,Ne,N_GRID)
  Xa_Sr     = matrix(NA,Ne,N_GRID)
  Xa_Mleaf1 = matrix(NA,Ne,N_GRID)
  Xa_Mleaf2 = matrix(NA,Ne,N_GRID)
  
  #### 2) storage of all time  ######
  
  ### Things to store for the tuning #### 
  Innov_mean_f = Innov_mean_a = S0_Ensmean = HPfHT_full = matrix(NA,store_length,N_GRID)
  
  # store the ensemble #
  S0_a_Ens = QTOT_mean = NULL 
  
  ## 3) Spin-up period ####
  for (k in 1:NL) {
    
    if (k <=365) {  
      if (k==1) {
        print("-------- Start spin-up period ----------")
        Xa_S01    = matrix(STATES$S0[1,],Ne,N_GRID)
        Xa_S02    = matrix(STATES$S0[2,],Ne,N_GRID)
        Xa_Ss1    = matrix(STATES$Ss[1,],Ne,N_GRID)
        Xa_Ss2    = matrix(STATES$Ss[2,],Ne,N_GRID)
        Xa_Sd1    = matrix(STATES$Sd[1,],Ne,N_GRID)
        Xa_Sd2    = matrix(STATES$Sd[2,],Ne,N_GRID)
        Xa_Sg     = matrix(STATES$Sg,    Ne,N_GRID)
        Xa_Sr     = matrix(STATES$Sr,    Ne,N_GRID)
        Xa_Mleaf1 = matrix(STATES$Mleaf[1,],Ne,N_GRID)
        Xa_Mleaf2 = matrix(STATES$Mleaf[2,],Ne,N_GRID)
      }
      
      for (i in 1:N_GRID) {
        Xa_S01[,i] = pmax(pmin(Xa_S01[,i],PARAMS$S0FC[1,i]),0.01)
        Xa_S02[,i] = pmax(pmin(Xa_S02[,i],PARAMS$S0FC[2,i]),0.01)
        Xa_Ss1[,i] = pmax(pmin(Xa_Ss1[,i],PARAMS$SsFC[1,i]),0.01)
        Xa_Ss2[,i] = pmax(pmin(Xa_Ss2[,i],PARAMS$SsFC[2,i]),0.01)
      }
      
      ### Perturbation of forcing for each ensemble member 
      FORCING$pe   = rbind(PE[k],PE[k])
      FORCING$pair = PAIR
      FORCING$u2   = rbind(U2[k],U2[k])
      
      # perturbation
      Pg = PG[k];   Rg = RG[k];  Ta = TA[k]
      RAIN_L = (1-Rain_Err)*Pg;     RAIN_H = (1+Rain_Err)*Pg
      
      # repeat run for each ensemble member 
      
      for (j in 1:Ne) {
        
        Pens = mapply(runif, n=1, min = RAIN_L, max = RAIN_H)
        FORCING$Pg = rbind(Pens,Pens)      
        
        Ta_ens = mapply(runif, n = 1, min = Ta - Temp_Err, max = Ta + Temp_Err)
        FORCING$Ta = rbind(Ta_ens,Ta_ens)
        
        Rg_ens = mapply(runif, n = 1, min = Rg - Rad_Err, max = Rg + Rad_Err)
        FORCING$Rg = rbind(Rg_ens, Rg_ens)
        
        STATES$S0[1,] = Xa_S01[j,]
        STATES$S0[2,] = Xa_S02[j,]
        STATES$Ss[1,] = Xa_Ss1[j,]
        STATES$Ss[2,] = Xa_Ss2[j,]
        STATES$Sd[1,] = Xa_Sd1[j,]
        STATES$Sd[2,] = Xa_Sd2[j,]
        STATES$Sg     = Xa_Sg[j,]
        STATES$Sr     = Xa_Sr[j,]
        STATES$Mleaf[1,]   = Xa_Mleaf1[j,]
        STATES$Mleaf[2,]   = Xa_Mleaf2[j,]
        
        OUT = AWRAL_TimeStep(FORCING,STATES,PARAMS)
        
        ## forecast state for the ensemble 
        Xf_S01[j,] = OUT$S0[1,]
        Xf_S02[j,] = OUT$S0[2,]
        Xf_Ss1[j,] = OUT$Ss[1,]
        Xf_Ss2[j,] = OUT$Ss[2,]
        Xf_Sd1[j,] = OUT$Sd[1,]
        Xf_Sd2[j,] = OUT$Sd[2,]
        Xf_Sg[j,]  = OUT$Sg
        Xf_Sr[j,]  = OUT$Sr
        Xf_Mleaf1[j,] = OUT$Mleaf[1,]
        Xf_Mleaf2[j,] = OUT$Mleaf[2,]
        QTOT_f[j,]    = OUT$Qtot       # flux, not getting updated 
        Xf_S0mean[j,] = OUT$S0mean
        
        stopifnot(length(which(is.na(Xf_S01[j,])==T)) == 0)
        stopifnot(length(which(is.na(Xf_S02[j,])==T)) == 0)
        stopifnot(length(which(is.na(Xf_Ss1[j,])==T)) == 0)
        stopifnot(length(which(is.na(Xf_Ss2[j,])==T)) == 0)
      }
      
      for (i in 1:N_GRID) {
        Xf_S01[,i] = pmax(pmin(Xf_S01[,i],PARAMS$S0FC[1,i]),0.01)
        Xf_S02[,i] = pmax(pmin(Xf_S02[,i],PARAMS$S0FC[2,i]),0.01)
        Xf_Ss1[,i] = pmax(pmin(Xf_Ss1[,i],PARAMS$SsFC[1,i]),0.01)
        Xf_Ss2[,i] = pmax(pmin(Xf_Ss2[,i],PARAMS$SsFC[2,i]),0.01)
        Xf_Sd1[,i] = pmax(pmin(Xf_Sd1[,i],PARAMS$SdFC[1,i]),0.01)
        Xf_Sd2[,i] = pmax(pmin(Xf_Sd2[,i],PARAMS$SdFC[2,i]),0.01)
      }
      Xa_S01 = Xf_S01;   Xa_S02 = Xf_S02
      Xa_Ss1 = Xf_Ss1;   Xa_Ss2 = Xf_Ss2
      Xa_Sd1 = Xf_Sd1;   Xa_Sd2 = Xf_Sd2
      Xa_Sg  = Xf_Sg;    Xa_Sr  = Xf_Sr
      Xa_Mleaf1 = Xf_Mleaf1;      Xa_Mleaf2 = Xf_Mleaf2
      
    } else {
      
      #### 4) ASSIMILATION PERIOD #######  
      for (i in 1:N_GRID) {
        Xa_S01[,i] = pmax(pmin(Xa_S01[,i],PARAMS$S0FC[1,i]),0.01)
        Xa_S02[,i] = pmax(pmin(Xa_S02[,i],PARAMS$S0FC[2,i]),0.01)
        Xa_Ss1[,i] = pmax(pmin(Xa_Ss1[,i],PARAMS$SsFC[1,i]),0.01)
        Xa_Ss2[,i] = pmax(pmin(Xa_Ss2[,i],PARAMS$SsFC[2,i]),0.01)
      }
      ### perturbation  ####
      
      FORCING$pe   = rbind(PE[k],PE[k])
      FORCING$pair = PAIR
      FORCING$u2   = rbind(U2[k],U2[k])
    
      Pg = PG[k];     Rg = RG[k];    Ta = TA[k]
      RAIN_L = (1- Rain_Err)*Pg;     RAIN_H = (1 + Rain_Err) *Pg
      
      ## Ensemble forecast ####
      
      for (j in 1:Ne) {
        
        Pens = mapply(runif, n=1, min = RAIN_L, max = RAIN_H)
        FORCING$Pg = rbind(Pens,Pens)     
        Ta_ens = mapply(runif, n = 1, min = Ta - Temp_Err , max = Ta + Temp_Err)
        FORCING$Ta = rbind(Ta_ens,Ta_ens)
        Rg_ens = mapply(runif, n=1, min = Rg - Rad_Err, max = Rg + Rad_Err)
        FORCING$Rg = rbind(Rg_ens, Rg_ens)
        
        ### Update the analysis state after assimilation ###

        STATES$S0[1,] = Xa_S01[j,]
        STATES$S0[2,] = Xa_S02[j,]
        STATES$Ss[1,] = Xa_Ss1[j,]
        STATES$Ss[2,] = Xa_Ss2[j,]
        STATES$Sd[1,] = Xa_Sd1[j,]
        STATES$Sd[2,] = Xa_Sd2[j,]
        STATES$Sg     = Xa_Sg[j,]
        STATES$Sr     = Xa_Sr[j,]
        STATES$Mleaf[1,]   = Xa_Mleaf1[j,]
        STATES$Mleaf[2,]   = Xa_Mleaf2[j,]
        
        OUT = AWRAL_TimeStep(FORCING,STATES,PARAMS)
        
        ## forecast state for the ensemble 
        Xf_S01[j,] = OUT$S0[1,]
        Xf_S02[j,] = OUT$S0[2,]
        Xf_Ss1[j,] = OUT$Ss[1,]
        Xf_Ss2[j,] = OUT$Ss[2,]
        Xf_Sd1[j,] = OUT$Sd[1,]
        Xf_Sd2[j,] = OUT$Sd[2,]
        Xf_Sg[j,]  = OUT$Sg
        Xf_Sr[j,]  = OUT$Sr
        Xf_Mleaf1[j,] = OUT$Mleaf[1,]
        Xf_Mleaf2[j,] = OUT$Mleaf[2,]
        QTOT_f[j,]    = OUT$Qtot       # flux, not getting updated 
        Xf_S0mean[j,] = OUT$S0mean  #### Ensemble S0 (H(x))

        stopifnot(length(which(is.na(Xf_Ss1[j,])==T)) == 0)
        stopifnot(length(which(is.na(Xf_Ss2[j,])==T)) == 0)
      }
      
      for (i in 1:N_GRID) {
        Xf_S01[,i] = pmax(pmin(Xf_S01[,i],PARAMS$S0FC[1,i]),0.01)
        Xf_S02[,i] = pmax(pmin(Xf_S02[,i],PARAMS$S0FC[2,i]),0.01)
        Xf_Ss1[,i] = pmax(pmin(Xf_Ss1[,i],PARAMS$SsFC[1,i]),0.01)
        Xf_Ss2[,i] = pmax(pmin(Xf_Ss2[,i],PARAMS$SsFC[2,i]),0.01)
        Xf_Sd1[,i] = pmax(pmin(Xf_Sd1[,i],PARAMS$SdFC[1,i]),0.01)
        Xf_Sd2[,i] = pmax(pmin(Xf_Sd2[,i],PARAMS$SdFC[2,i]),0.01)
        hXf_w0[,i] = pmax(pmin(hXf_w0[,i],1),0)
      }
      
      S0mean_f = sweep(Xf_S01,MARGIN =  2, PARAMS$Fhru[1,], "*") + 
        sweep(Xf_S02,MARGIN =  2, PARAMS$Fhru[2,], "*")
      Ssmean_f = sweep(Xf_Ss1,MARGIN =  2, PARAMS$Fhru[1,], "*") + 
        sweep(Xf_Ss2,MARGIN =  2, PARAMS$Fhru[2,], "*")
      
      ####  5) Grid specific Kalman gain calculation, using the evolved ensemble ######
      
      ### Grid by grid inflation of ensemble variance 
      
      Dev_S01 = (Xf_S01 - matrix(colMeans(Xf_S01, na.rm = T),Ne,N_GRID,byrow=T))
      Dev_S02 = (Xf_S02 - matrix(colMeans(Xf_S02, na.rm = T),Ne,N_GRID,byrow=T))
      
      Xf_S01 = sweep(Dev_S01,2, IF, "*") + matrix(colMeans(Xf_S01, na.rm = T), Ne, N_GRID, byrow=T)
      Xf_S02 = sweep(Dev_S02,2, IF, "*") + matrix(colMeans(Xf_S02, na.rm = T), Ne, N_GRID, byrow=T)
      
      for (i in 1:N_GRID) {
        Xf_S01[,i] = pmax(pmin(Xf_S01[,i], PARAMS$S0FC[1,i]), 0.01)
        Xf_S02[,i] = pmax(pmin(Xf_S02[,i], PARAMS$S0FC[2,i]), 0.01)
      }
      
      ##### HPfHT: variance of H(x)  (i.e. S0mean) ######
      
      HPfHT = apply(Xf_S0mean,2,function(x) {var(x, na.rm=T)})     # Model variance 
      SMAP_OBS =  SMAP_scaled[k-365]
      
      for (i in 1:N_GRID) {
        
        ### gain calculation: grid specific ####
        
        PfHT = matrix(NA,4,1)
        PfHT[1] = cov(Xf_S01[,i], Xf_S0mean[,i], use = 'pairwise.complete.obs')
        PfHT[2] = cov(Xf_S02[,i], Xf_S0mean[,i], use = 'pairwise.complete.obs')
        PfHT[3] = cov(Xf_Ss1[,i], Xf_S0mean[,i], use = 'pairwise.complete.obs')
        PfHT[4] = cov(Xf_Ss2[,i], Xf_S0mean[,i], use = 'pairwise.complete.obs')
        
        if (is.na(SMAP_OBS[i])) {
          K   = rep(0,4)
          OBS = rep(-9999,Ne)  ## avoid NA in states  
          
        } else {
          OBS = rnorm(Ne, mean = SMAP_OBS[i], sd = SMAP_Err[i])
          K = PfHT / (HPfHT[i] + (SMAP_Err[i]**2))   # change to the space of error variance (4th May) 
          OBS[OBS<0] = 0   }
        
        # 6) state updating for all ensemble members #######
        
        Xa_S01[,i] = Xf_S01[,i] + K[1] * (OBS - Xf_S0mean[,i])
        Xa_S02[,i] = Xf_S02[,i] + K[2] * (OBS - Xf_S0mean[,i])
        Xa_Ss1[,i] = Xf_Ss1[,i] + K[3] * (OBS - Xf_S0mean[,i])
        Xa_Ss2[,i] = Xf_Ss2[,i] + K[4] * (OBS - Xf_S0mean[,i])
        
        Innov_grid = OBS - Xf_S0mean[,i]
        Innov_mean_f[k-365,i] = mean(Innov_grid, na.rm=T)
      }
      
      for (i in 1:N_GRID) {
        Xa_S01[,i] = pmax(pmin(Xa_S01[,i],PARAMS$S0FC[1,i]),0.01)
        Xa_S02[,i] = pmax(pmin(Xa_S02[,i],PARAMS$S0FC[2,i]),0.01)
        Xa_Ss1[,i] = pmax(pmin(Xa_Ss1[,i],PARAMS$SsFC[1,i]),0.01)
        Xa_Ss2[,i] = pmax(pmin(Xa_Ss2[,i],PARAMS$SsFC[2,i]),0.01)
      }
      
      Xa_Sd1 = Xf_Sd1;    Xa_Sd2 = Xf_Sd2
      Xa_Sg  = Xf_Sg;     Xa_Sr  = Xf_Sr
      Xa_Mleaf1 = Xf_Mleaf1;      Xa_Mleaf2 = Xf_Mleaf2
      
      ### Grid analysed S0 and Ss #####
      
      S0mean_a = sweep(Xa_S01, MARGIN =  2, PARAMS$Fhru[1,], "*") + 
        sweep(Xa_S02, MARGIN =  2, PARAMS$Fhru[2,], "*")
      Ssmean_a = sweep(Xa_Ss1, MARGIN =  2, PARAMS$Fhru[1,], "*") + 
        sweep(Xa_Ss2, MARGIN =  2, PARAMS$Fhru[2,], "*")
      
      ###### 7) save ensemble statistics from  evolving matrix  #######
      
      # Ensemble mean of streamflow ##
      QTOT_mean = c(QTOT_mean, mean(QTOT_f, na.rm=T))
      for (i in 1:N_GRID) {
      
      Innov_mean_a[k-365, i] = mean(OBS -S0mean_a[,i] , na.rm=T)
      HPfHT_full[k-365, i] = HPfHT
      S0_a_Ens = rbind(S0_a_Ens, S0mean_a[,i])
      S0_Ensmean[k-365, i] = mean(S0mean_a[,i], na.rm=T)
      
        }
      
      if (k%%800 ==0 ) {print(paste0("Ensemble run for time ", k, " is finished.")) } 
    }
  }

  ## Output 
  return(list(
    SMAP_scaled  = SMAP_scaled, SMAP_Err = SMAP_Err,  QTOT_mean = QTOT_mean, 
    HPfHT = HPfHT_full, S0_Ensmean = S0_Ensmean, 
    Innov_f = Innov_mean_f, Innov_a = Innov_mean_a, 
    S0_a_Ens = S0_a_Ens
  ))
  
}

