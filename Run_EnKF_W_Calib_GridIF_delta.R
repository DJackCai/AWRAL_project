source("awra_run_func.R")
source("Forcing_extracts.R")
source("PARAMS_set.R")
source("AWRAL_TimeStep.R")
source("Run_AWRAL_Calib_Q.R")
source("Wetness_scaling.R")

##### Assimilation of SMAP soil moisture data into the AWRA-L model 

## Reference: Renzullo, L. J., van Dijk, A. I. J. M., Perraud, J.-M., Collins, D., Henderson, B., Jin, H., Smith, A. B. and McJannet, D. L., 2014. 
## Continental satellite soil moisture data assimilation improves root-zone moisture analysis for water resources assessment, Journal of Hydrology, 519: 2747-2762. 


Run_AWRAL_EnKF_W_Calib_GridIF <- function(catname, IF_vec, Ne, par_mat, 
                                          start_date = "2014-01-01", end_date= "2018-12-31", store_length = 1461, 
                                          SMAP_Err = NULL, rescaling = "MVM", 
                                          Rain_Err = 0.6 , Temp_Err = 2.5 , 
                                          Rad_Err = 50) {
  
  # catname : catchment name (catchment ID)
  # Ne: number of ensemble members 
  # start_date: start date of model run, here, the start date of the spin-up period
  # end_date: end date of model run
  # SMAP_Err: error variance of SMAP data, can supply our own, otherwise use difference estimator 
  # Rain_Err : multiplicative perturbation to rainfall 
  # Temp_Err : temperature additive perturbation
  # Rad_Err : solar radiation additive perturbation
  # par_mat: optimised parameter matrix 
  # Updated - IF_vec: grid-specific inflation factor for optimal assimilation
  
  
  ## PARAM set & N_GRID ####
  STATIC_PAR <- as.matrix(read.csv(paste0("./forcings/",catname,"_resampled_static.csv")))
  N_GRID = ncol(STATIC_PAR)
  
  dates_1321 = dates_seq(2013,2021)
  start_run = which(dates_1321 == as.Date(start_date))   # spin up starts in 2014
  end_run = which(dates_1321 == as.Date(end_date)) 
  
  ## forcing for all time  #####
  all_forcing <- grid_forcing_catchment(catname = catname, N_GRID = N_GRID ,
                                        start = start_run, end= end_run)
  PG <- all_forcing$PG
  RG <- all_forcing$RG
  TA <- all_forcing$TA
  PE <- all_forcing$PE
  U2 <- all_forcing$U2
  PAIR = rep(97500.00,N_GRID)
  
  ## Calibrated parameters  ######
  PARAMS = PARAMS_OPT_CATCH(N_GRID = N_GRID, par_mat =  par_mat, 
                            STATIC_PAR = STATIC_PAR)
  meanPET  = as.vector(STATIC_PAR[5,])
  meanP   = colMeans(all_forcing$PG)
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
  
  mod_run = Run_AWRAL_calib_Q(pars = par_mat, all_forcing = all_forcing, STATIC_PAR = STATIC_PAR, N_GRID =N_GRID  )
  
  ##### Convert SMAP observation to wetness from 2015-2018 #####
  
  SMAP_data  = as.matrix(read.csv(paste0("./forcings/",catname, "_resampled_SMAP1519.csv")))[1:1461,paste0("grid",1:N_GRID)]
  S0mean_sim = mod_run$S0mean[-c(1:365),]
  
  # Converted to relative wetness # 
  SMAP_wetness = apply(SMAP_data, 2, function(x) {
    SMAP_wt = quantile(x, prob=0.02, na.rm=T)
    SMAP_fc = quantile(x, prob=0.98, na.rm=T)
    return( (x-SMAP_wt)/(SMAP_fc-SMAP_wt))
  })
  ## shift to AWRA time 
  SMAP_wetness_shift = apply(SMAP_wetness,2, function(x) {return(c(NA,x[-length(x)]))} )
  SMAP_wetness_shift[SMAP_wetness_shift < 0.0] = 0.0
  SMAP_wetness_shift[SMAP_wetness_shift > 1.0] = 1.0
  
  ## Rescaled SMAP wetness to AWRA wetness to make sure the difference is 0 on average 
  
  AWRA_S0max = apply(S0mean_sim,2,quantile,prob=0.98,na.rm=T)
  AWRA_S0min = apply(S0mean_sim,2,quantile,prob=0.02,na.rm=T)
  AWRA_wetness = apply(S0mean_sim,2,function(x) {
    AWRA_wt = quantile(x, prob=0.02, na.rm=T)
    AWRA_fc = quantile(x, prob=0.98, na.rm=T)
    return( (x-AWRA_wt)/(AWRA_fc-AWRA_wt))
  })
  AWRA_wetness[AWRA_wetness < 0.0] = 0.0
  AWRA_wetness[AWRA_wetness > 1.0] = 1.0
  
  if (rescaling == "CDF" ) {
    SMAP_wetness_scaled  = mapply(function(x,y) { 
      SMAP_cdfmatch <- CDFt(y,x,x)
      return(SMAP_cdfmatch$DS) }, as.data.frame(SMAP_wetness_shift), 
      as.data.frame(AWRA_wetness)) } else if (rescaling == "MVM" ) {
        SMAP_wetness_scaled  = mapply(function(x,y) { return(Mean_var_match(x,y)) },
                                      as.data.frame(SMAP_wetness_shift), as.data.frame(AWRA_wetness))
      }
  
  SMAP_wetness_scaled[SMAP_wetness_scaled < 0.0] = 0.0
  SMAP_wetness_scaled[SMAP_wetness_scaled > 1.0] = 1.0
  
  ## SMAP wetness error variance after rescaling
  
  if (is.null(SMAP_Err)) {
    # SMAP_Err = apply(SMAP_wetness_scaled,2,Diff_error_estim) 
    sd_w_AWRA = apply(AWRA_wetness,2,sd,na.rm=T)
    sd_w_SMAP = apply(SMAP_wetness_scaled,2,sd,na.rm=T)
    SMAP_maxmin =  apply(SMAP_data,2,quantile,0.98,na.rm=T) - apply(SMAP_data,2,quantile,0.02,na.rm=T)
    SMAP_Err = 0.04 * (sd_w_AWRA/sd_w_SMAP) / SMAP_maxmin
  } else {
    # if there is pre-defined error inputs
    stopifnot(ncol(SMAP_wetness_scaled) == length(SMAP_Err))
    SMAP_Err = SMAP_Err }
  
  print("-------- Pre-processing of SMAP data have been finished -----------")
  
           ########## EnKF module #########
  
  STATES = STATES_unpert = initial_states(PARAMS = PARAMS,N_GRID = N_GRID, meanP = meanP)
  NL =  dim(PG)[1] 
  FORCING = FORCING_Unpert = empty_forcing()
  
  
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
  QTOT_f = ET_f = FSAT_f = QR_f = Qg_f = SE_f =  matrix(NA,Ne,N_GRID)
  
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
  S01_median_f = S01_median_a = S02_median_f = S02_median_a = matrix(NA,store_length,N_GRID)
  Ss1_median_f = Ss1_median_a = Ss2_median_f = Ss2_median_a = matrix(NA,store_length,N_GRID)
  
  S0_min_f =   S0_max_f = S0_75_f = S0_25_f = S0_median_f = S0_mean_f = matrix(NA,store_length,N_GRID)
  Ss_min_f =   Ss_max_f = Ss_75_f = Ss_25_f = Ss_median_f = Ss_mean_f = matrix(NA,store_length,N_GRID)
  
  w0_min_f =   w0_max_f = w0_75_f = w0_25_f = w0_median_f = w0_mean_f = matrix(NA,store_length,N_GRID)
  w0_min_a =   w0_max_a = w0_75_a = w0_25_a = w0_median_a = w0_mean_a = matrix(NA,store_length,N_GRID)
  
  S0_min_a =   S0_max_a = S0_75_a = S0_25_a = S0_median_a = S0_mean_a = matrix(NA,store_length,N_GRID)
  Ss_min_a =   Ss_max_a = Ss_75_a = Ss_25_a = Ss_median_a = Ss_mean_a = matrix(NA,store_length,N_GRID)
  QTOT_median =  QTOT_mean   = QTOT_75 = QTOT_25 =  QTOT_max = QTOT_min = matrix(NA,store_length,N_GRID)
  
  ET_median =  ET_mean   = ET_75 = ET_25 =  ET_max = ET_min = matrix(NA,store_length,N_GRID)
  SG_median =  SG_mean   = SG_75 = SG_25 =  SG_max = SG_min = matrix(NA,store_length,N_GRID)
  FSAT_median =  FSAT_mean   = FSAT_75 = FSAT_25 = matrix(NA,store_length,N_GRID)
  
  S01_ens10 = S01_ens20 = S01_ens30  = matrix(NA,store_length,N_GRID)
  S02_ens10 = S02_ens20 = S02_ens30  = matrix(NA,store_length,N_GRID)
  Kgain_10 = Kgain_20 = Kgain_30 =  matrix(NA,store_length,4)
  
  QR_mean = QR_median = QR_75 = QR_25 = matrix(NA,store_length,N_GRID)
  QG_mean = QG_median = QG_75 = QG_25 = matrix(NA,store_length,N_GRID)
  
  SE_mean = SE_median = SE_75 = SE_25 = matrix(NA,store_length,N_GRID)
  
  Innov_mean = Innov_median = Innov_75 = Innov_25 = matrix(NA,store_length,N_GRID)
  
  BIAS_S01 = BIAS_S02 = BIAS_Ss1 = BIAS_Ss2 = matrix(NA,store_length,N_GRID)
  
  
  listnames = 1:N_GRID
  w0_ensemble  = sapply(listnames, function(x) NULL)
  
  Q_grid_ensemble = sapply(listnames, function(x) NULL)
  Qtot_ensemble = NULL
  
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
      ## unperturbed run 
      FORCING_Unpert$Pg = rbind(Pg,Pg)
      FORCING_Unpert$Rg = rbind(Rg,Rg)
      FORCING_Unpert$Ta = rbind(Ta,Ta)
      
      STATES_unpert$S0    = rbind(apply(Xa_S01,2,mean),apply(Xa_S02,2,mean))
      STATES_unpert$Ss    = rbind(apply(Xa_Ss1,2,mean),apply(Xa_Ss2,2,mean))
      STATES_unpert$Sd    = rbind(apply(Xa_Sd1,2,mean),apply(Xa_Sd2,2,mean))
      STATES_unpert$Sg    = apply(Xa_Sg,2,mean)
      STATES_unpert$Sr    = apply(Xa_Sr,2,mean)
      STATES_unpert$Mleaf = rbind(apply(Xa_Mleaf1,2,mean),apply(Xa_Mleaf2,2,mean))
      
      OUT_unpert = AWRAL_TimeStep(FORCING_Unpert,STATES_unpert,PARAMS)
      
      ### Perturbation of forcing for each ensemble member 
      FORCING$pe   = rbind(PE[k,],PE[k,])
      FORCING$pair = PAIR
      FORCING$u2   = rbind(U2[k,],U2[k,])
      # perturbation
      Pg = PG[k,]; Rg = RG[k,];  Ta = TA[k,]
      RAIN_L = (1-Rain_Err)*Pg;     RAIN_H = (1+Rain_Err)*Pg
      
      # repeat run for each ensemble member 
      
      for (j in 1:Ne) {
        
        Pens = mapply(runif, n=1, min = RAIN_L, max = RAIN_H)
        FORCING$Pg = rbind(Pens,Pens)      
        
        Ta_ens = mapply(runif, n = 1, min = Ta - Temp_Err, max = Ta + Temp_Err)
        FORCING$Ta = rbind(Ta_ens,Ta_ens)
        
        Rg_ens = mapply(runif, n=1, min = Rg - Rad_Err, max = Rg + Rad_Err)
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
        ET_f[j, ]     = OUT$ETmean
        FSAT_f[j, ] = OUT$fsat
        
        stopifnot(length(OUT$S0mean)==length(AWRA_S0min))
        hXf_w0[j, ] = (OUT$S0mean-AWRA_S0min)/(AWRA_S0max-AWRA_S0min)
        # hXf_w0[j, ] = sweep(sweep(OUT$S0mean,2,AWRA_S0min,"-"), 2, AWRA_S0max - AWRA_S0min, "/")
        hXf_w0[hXf_w0 < 0.0] = 0.0
        hXf_w0[hXf_w0 > 1.0] = 1.0
        
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
      
      ## Bias correction for forecast model states ####
      
      Bias_S01 = apply(Xf_S01,2,mean) - OUT_unpert$S0[1,]
      Bias_S02 = apply(Xf_S02,2,mean) - OUT_unpert$S0[2,]
      Bias_Ss1 = apply(Xf_Ss1,2,mean) - OUT_unpert$Ss[1,]
      Bias_Ss2 = apply(Xf_Ss2,2,mean) - OUT_unpert$Ss[2,]
      
      Bias_Sd1 = apply(Xf_Sd1,2,mean) - OUT_unpert$Sd[1,]
      Bias_Sd2 = apply(Xf_Sd2,2,mean) - OUT_unpert$Sd[2,]
      
      Xf_S01_Correct = sweep(Xf_S01, 2, Bias_S01, "-")
      Xf_S02_Correct = sweep(Xf_S02, 2, Bias_S02, "-")
      Xf_Ss1_Correct = sweep(Xf_Ss1, 2, Bias_Ss1, "-")
      Xf_Ss2_Correct = sweep(Xf_Ss2, 2, Bias_Ss2, "-")
      Xf_Sd1_Correct = sweep(Xf_Sd1, 2, Bias_Sd1, "-")
      Xf_Sd2_Correct = sweep(Xf_Sd2, 2, Bias_Sd2, "-")
      for (i in 1:N_GRID) {
        Xf_S01_Correct[,i] = pmax(pmin(Xf_S01_Correct[,i],PARAMS$S0FC[1,i]),0.01)
        Xf_S02_Correct[,i] = pmax(pmin(Xf_S02_Correct[,i],PARAMS$S0FC[2,i]),0.01)
        Xf_Ss1_Correct[,i] = pmax(pmin(Xf_Ss1_Correct[,i],PARAMS$SsFC[1,i]),0.01)
        Xf_Ss2_Correct[,i] = pmax(pmin(Xf_Ss2_Correct[,i],PARAMS$SsFC[2,i]),0.01)
        Xf_Sd1_Correct[,i] = pmax(pmin(Xf_Sd1_Correct[,i],PARAMS$SdFC[1,i]),0.1)
        Xf_Sd2_Correct[,i] = pmax(pmin(Xf_Sd2_Correct[,i],PARAMS$SdFC[2,i]),0.1)
      }
      Xa_S01 = Xf_S01_Correct;   Xa_S02 = Xf_S02_Correct
      Xa_Ss1 = Xf_Ss1_Correct;   Xa_Ss2 = Xf_Ss2_Correct
      Xa_Sd1 = Xf_Sd1_Correct;   Xa_Sd2 = Xf_Sd2_Correct
      Xa_Sg  = Xf_Sg  ;  Xa_Sr  = Xf_Sr
      Xa_Mleaf1 = Xf_Mleaf1;      Xa_Mleaf2 = Xf_Mleaf2
      stopifnot(max(apply(Xf_S01_Correct,2,mean)  - OUT_unpert$S0[1,]) < 0.001)
      
    } else {
      
      #### 4) ASSIMILATION PERIOD #######  
      for (i in 1:N_GRID) {
        Xa_S01[,i] = pmax(pmin(Xa_S01[,i],PARAMS$S0FC[1,i]),0.01)
        Xa_S02[,i] = pmax(pmin(Xa_S02[,i],PARAMS$S0FC[2,i]),0.01)
        Xa_Ss1[,i] = pmax(pmin(Xa_Ss1[,i],PARAMS$SsFC[1,i]),0.01)
        Xa_Ss2[,i] = pmax(pmin(Xa_Ss2[,i],PARAMS$SsFC[2,i]),0.01)
      }
      
      ## unperturbed run ###
      FORCING$pe   =  FORCING_Unpert$pe = rbind(PE[k,],PE[k,])
      FORCING$pair =  FORCING_Unpert$pair = PAIR
      FORCING$u2   =  FORCING_Unpert$u2  = rbind(U2[k,],U2[k,])
      FORCING_Unpert$Pg = rbind(Pg,Pg)
      FORCING_Unpert$Rg = rbind(Rg,Rg)
      FORCING_Unpert$Ta = rbind(Ta,Ta)
      
      STATES_unpert$S0    = rbind(apply(Xa_S01,2,mean),apply(Xa_S02,2,mean))
      STATES_unpert$Ss    = rbind(apply(Xa_Ss1,2,mean),apply(Xa_Ss2,2,mean))
      STATES_unpert$Sd    = rbind(apply(Xa_Sd1,2,mean),apply(Xa_Sd2,2,mean))
      STATES_unpert$Sg    = apply(Xa_Sg,2,mean)
      STATES_unpert$Sr    = apply(Xa_Sr,2,mean)
      STATES_unpert$Mleaf = rbind(apply(Xa_Mleaf1,2,mean),apply(Xa_Mleaf2,2,mean))
      
      OUT_unpert = AWRAL_TimeStep(FORCING_Unpert,STATES_unpert,PARAMS)
    
      ### perturbed run  ####
      Pg = PG[k,];     Rg = RG[k,];    Ta = TA[k,]
      RAIN_L = (1- Rain_Err)*Pg;     RAIN_H = (1 + Rain_Err) *Pg
      
      ## ensemble forecast, parallel run for each grid ####
      
      for (j in 1:Ne) {
        
        Pens = mapply(runif, n=1, min = RAIN_L, max = RAIN_H)
        FORCING$Pg = rbind(Pens,Pens)     
        Ta_ens = mapply(runif, n = 1, min = Ta - Temp_Err , max = Ta + Temp_Err)
        FORCING$Ta = rbind(Ta_ens,Ta_ens)
        Rg_ens = mapply(runif, n=1, min = Rg - Rad_Err, max = Rg + Rad_Err)
        FORCING$Rg = rbind(Rg_ens, Rg_ens)
        
        ### Update the analysis state after assimilation
        
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
        ET_f[j, ]     = c(1,1) %*% (OUT$ETtot * PARAMS$Fhru) 
        FSAT_f[j, ] = OUT$fsat
        QR_f[j, ]  =  c(1,1) %*%  (OUT$QR * PARAMS$Fhru)
        Qg_f[j,]   =  OUT$Qg
        SE_f[j,]  =  c(1,1) %*%  (OUT$SE * PARAMS$Fhru)
        
        stopifnot(length(OUT$S0mean)==length(AWRA_S0min))
        hXf_w0[j, ] = (OUT$S0mean-AWRA_S0min)/(AWRA_S0max-AWRA_S0min)
        
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
        hXf_w0[,i] = pmax(pmin(hXf_w0[,i],1),0)
      }
      
      ### Bias correction for forecast model states  ######
      
      Bias_S01 = apply(Xf_S01,2,mean) - OUT_unpert$S0[1,]
      Bias_S02 = apply(Xf_S02,2,mean) - OUT_unpert$S0[2,]
      Bias_Ss1 = apply(Xf_Ss1,2,mean) - OUT_unpert$Ss[1,]
      Bias_Ss2 = apply(Xf_Ss2,2,mean) - OUT_unpert$Ss[2,]
      
      Bias_Sd1 = apply(Xf_Sd1,2,mean) - OUT_unpert$Sd[1,]
      Bias_Sd2 = apply(Xf_Sd2,2,mean) - OUT_unpert$Sd[2,]
      
      Xf_S01_Correct = sweep(Xf_S01, 2, Bias_S01, "-")
      Xf_S02_Correct = sweep(Xf_S02, 2, Bias_S02, "-")
      Xf_Ss1_Correct = sweep(Xf_Ss1, 2, Bias_Ss1, "-")
      Xf_Ss2_Correct = sweep(Xf_Ss2, 2, Bias_Ss2, "-")
      Xf_Sd1_Correct = sweep(Xf_Sd1, 2, Bias_Sd1, "-")
      Xf_Sd2_Correct = sweep(Xf_Sd2, 2, Bias_Sd2, "-")
      
      BIAS_S01[k-365,] = Bias_S01
      BIAS_S02[k-365,] = Bias_S02
      BIAS_Ss1[k-365,] = Bias_Ss1
      BIAS_Ss2[k-365,] = Bias_Ss2
      BIAS_QTOT[k-365, ] = Bias_QTOT
      stopifnot(max(apply(Xf_S01_Correct,2,mean)  - OUT_unpert$S0[1,]) < 0.001)
      
      ### correct to SM range 
      for (i in 1:N_GRID) {
        Xf_S01_Correct[,i] = pmax(pmin(Xf_S01_Correct[,i],PARAMS$S0FC[1,i]),0.01)
        Xf_S02_Correct[,i] = pmax(pmin(Xf_S02_Correct[,i],PARAMS$S0FC[2,i]),0.01)
        Xf_Ss1_Correct[,i] = pmax(pmin(Xf_Ss1_Correct[,i],PARAMS$SsFC[1,i]),0.01)
        Xf_Ss2_Correct[,i] = pmax(pmin(Xf_Ss2_Correct[,i],PARAMS$SsFC[2,i]),0.01)
        Xf_Sd1_Correct[,i] = pmax(pmin(Xf_Sd1_Correct[,i],PARAMS$SdFC[1,i]),0.01)
        Xf_Sd2_Correct[,i] = pmax(pmin(Xf_Sd2_Correct[,i],PARAMS$SdFC[2,i]),0.01)
      }
      
      ## grid S0 and Ss after bias correction
      S0mean_f = sweep(Xf_S01_Correct, MARGIN =  2, PARAMS$Fhru[1,], "*") + 
        sweep(Xf_S02_Correct,MARGIN =  2, PARAMS$Fhru[2,], "*")
      Ssmean_f = sweep(Xf_Ss1_Correct, MARGIN =  2, PARAMS$Fhru[1,], "*") + 
        sweep(Xf_Ss2_Correct,MARGIN =  2, PARAMS$Fhru[2,], "*")
      
      hXf_w0_Correct = sweep(sweep(S0mean_f,2,AWRA_S0min,"-"), 2, AWRA_S0max - AWRA_S0min, "/")
      hXf_w0_Correct[hXf_w0_Correct>1] = 1
      hXf_w0_Correct[hXf_w0_Correct<0] = 0
      
      if (k%%200 ==0 ) {print(paste0("Bias correction for time ", k, " is finished.")) }
    
      ####  5) Grid specific Kalman gain calculation, using the evolved ensemble ######
      
      ### Grid by grid inflation of ensemble variance
      Dev_S01 = (Xf_S01_Correct - matrix(colMeans(Xf_S01_Correct, na.rm = T),Ne,N_GRID,byrow=T))
      Dev_S02 = (Xf_S02_Correct - matrix(colMeans(Xf_S02_Correct, na.rm = T),Ne,N_GRID,byrow=T))
      
      Xf_S01_Correct = sweep(Dev_S01,2, IF_vec, "*") + matrix(colMeans(Xf_S01_Correct), Ne, N_GRID, byrow=T)
      Xf_S02_Correct = sweep(Dev_S02,2, IF_vec, "*") + matrix(colMeans(Xf_S02_Correct), Ne, N_GRID, byrow=T)
      
      for (i in 1:N_GRID) {
        Xf_S01_Correct[,i] = pmax(pmin(Xf_S01_Correct[,i], PARAMS$S0FC[1,i]), 0.01)
        Xf_S02_Correct[,i] = pmax(pmin(Xf_S02_Correct[,i], PARAMS$S0FC[2,i]), 0.01)
      }
      
      ##### HPfHT: now is the variance of H(x) ######
      HPfHT = apply(hXf_w0_Correct,2,function(x) {var(x, na.rm=T)})     # Model variance 
      SMAP_OBS =  SMAP_wetness_scaled[k-365,]
      
      for (i in 1:N_GRID) {
        
        ### gain calculation: grid specific ####
        PfHT = matrix(NA,4,1)
        PfHT[1] = cov(Xf_S01_Correct[,i], hXf_w0_Correct[,i], use = 'pairwise.complete.obs')
        PfHT[2] = cov(Xf_S02_Correct[,i], hXf_w0_Correct[,i], use = 'pairwise.complete.obs')
        PfHT[3] = cov(Xf_Ss1_Correct[,i], hXf_w0_Correct[,i], use = 'pairwise.complete.obs')
        PfHT[4] = cov(Xf_Ss2_Correct[,i], hXf_w0_Correct[,i], use = 'pairwise.complete.obs')
        
        if (is.na(SMAP_OBS[i])) {
          K   = rep(0,4)
          OBS = rep(-9999.,Ne)
          
        } else {
          OBS = rnorm(Ne, mean = SMAP_OBS[i], sd = SMAP_Err[i])
          K = PfHT / (HPfHT[i] + (SMAP_Err[i]**2))   # change to the space of error variance (4th May) 
          OBS[OBS<0] = 0   }
        
        # store the gain 
        if (i == 10) { Kgain_10[k-365,] = K } else if (i == 20) {
          Kgain_20[k-365,] = K    } else if (i == 30) {
            Kgain_30[k-365,] = K    } 
        
        # 6) state updating for all ensemble members #######
        Xa_S01[,i] = Xf_S01_Correct[,i] + K[1] * (OBS - hXf_w0_Correct[,i])
        Xa_S02[,i] = Xf_S02_Correct[,i] + K[2] * (OBS - hXf_w0_Correct[,i])
        Xa_Ss1[,i] = Xf_Ss1_Correct[,i] + K[3] * (OBS - hXf_w0_Correct[,i])
        Xa_Ss2[,i] = Xf_Ss2_Correct[,i] + K[4] * (OBS - hXf_w0_Correct[,i])
        
        Innov_grid = OBS - hXf_w0_Correct[,i]
        Innov_mean[k-365,i] = mean(OBS - hXf_w0_Correct[,i], na.rm=T)
        Innov_median[k-365,i] = median(OBS - hXf_w0_Correct[,i], na.rm=T)
        Innov_75[k-365,i] = quantile(OBS - hXf_w0_Correct[,i], 0.75, na.rm=T)
        Innov_25[k-365,i] = quantile(OBS - hXf_w0_Correct[,i], 0.25, na.rm=T)
        
      }
      
      for (i in 1:N_GRID) {
        Xa_S01[,i] = pmax(pmin(Xa_S01[,i],PARAMS$S0FC[1,i]),0.01)
        Xa_S02[,i] = pmax(pmin(Xa_S02[,i],PARAMS$S0FC[2,i]),0.01)
        Xa_Ss1[,i] = pmax(pmin(Xa_Ss1[,i],PARAMS$SsFC[1,i]),0.01)
        Xa_Ss2[,i] = pmax(pmin(Xa_Ss2[,i],PARAMS$SsFC[2,i]),0.01)
      }
      
      S0mean_a = sweep(Xa_S01, MARGIN =  2, PARAMS$Fhru[1,], "*") + 
        sweep(Xa_S02, MARGIN =  2, PARAMS$Fhru[2,], "*")
      Ssmean_a = sweep(Xa_Ss1, MARGIN =  2, PARAMS$Fhru[1,], "*") + 
        sweep(Xa_Ss2, MARGIN =  2, PARAMS$Fhru[2,], "*")
      
      Xa_Sd1 = Xf_Sd1;   Xa_Sd2 = Xf_Sd2
      Xa_Sg  = Xf_Sg;    Xa_Sr  = Xf_Sr
      Xa_Mleaf1 = Xf_Mleaf1;      Xa_Mleaf2 = Xf_Mleaf2
      
      ###### 7) save ensemble statistics from the evolving matrix  #######
      
      S01_median_f[k-365,] = apply(Xf_S01_Correct,2,function(x) {median(x,na.rm=T)} )
      S01_median_a[k-365,] = apply(Xa_S01,2,function(x) {median(x,na.rm=T)})
      S02_median_f[k-365, ] = apply(Xf_S02_Correct,2,function(x) {median(x,na.rm=T)})
      S02_median_a[k-365, ] = apply(Xa_S02,2,function(x) {median(x,na.rm=T)}) 
      Ss1_median_f[k-365, ] = apply(Xf_Ss1_Correct,2,function(x) {median(x,na.rm=T)})
      Ss1_median_a[k-365,]  = apply(Xa_Ss1,2,function(x) {median(x,na.rm=T)})
      Ss2_median_f[k-365, ] = apply(Xf_Ss2_Correct,2,function(x) {median(x,na.rm=T)})
      Ss2_median_a[k-365,]  = apply(Xa_Ss2,2,function(x) {median(x,na.rm=T)})
      
      S0_min_f[k-365,] = apply(S0mean_f,2,function(x) {min(x,na.rm=T)} )
      S0_max_f[k-365,] = apply(S0mean_f,2,function(x) {max(x,na.rm=T)})
      S0_75_f[k-365, ] = apply(S0mean_f,2,function(x) {quantile(x,prob = 0.75,na.rm=T)})
      S0_25_f[k-365, ] = apply(S0mean_f,2,function(x) {quantile(x,prob = 0.25,na.rm=T)}) 
      S0_median_f[k-365, ] = apply(S0mean_f,2,function(x) {median(x,na.rm=T)})
      S0_mean_f[k-365,]  = apply(S0mean_f,2,function(x) {mean(x,na.rm=T)})
      
      Ss_min_f[k-365,] = apply(Ssmean_f,2,function(x) {min(x,na.rm=T)} )
      Ss_max_f[k-365,] = apply(Ssmean_f,2,function(x) {max(x,na.rm=T)})
      Ss_75_f[k-365, ] = apply(Ssmean_f,2,function(x) {quantile(x,prob = 0.75,na.rm=T)})
      Ss_25_f[k-365, ] = apply(Ssmean_f,2,function(x) {quantile(x,prob = 0.25,na.rm=T)}) 
      Ss_median_f[k-365, ] = apply(Ssmean_f,2,function(x) {median(x,na.rm=T)})
      Ss_mean_f[k-365,]  = apply(Ssmean_f,2,function(x) {mean(x,na.rm=T)})
      
      S0_min_a[k-365,] = apply(S0mean_a,2,function(x) {min(x,na.rm=T)} )
      S0_max_a[k-365,] = apply(S0mean_a,2,function(x) {max(x,na.rm=T)})
      S0_75_a[k-365, ] = apply(S0mean_a,2,function(x) {quantile(x,prob = 0.75,na.rm=T)})
      S0_25_a[k-365, ] = apply(S0mean_a,2,function(x) {quantile(x,prob = 0.25,na.rm=T)}) 
      S0_median_a[k-365, ] = apply(S0mean_a,2,function(x) {median(x,na.rm=T)})
      S0_mean_a[k-365,]  = apply(S0mean_a,2,function(x) {mean(x,na.rm=T)})
      
      Ss_min_a[k-365,] = apply(Ssmean_a,2,function(x) {min(x,na.rm=T)} )
      Ss_max_a[k-365,] = apply(Ssmean_a,2,function(x) {max(x,na.rm=T)})
      Ss_75_a[k-365, ] = apply(Ssmean_a,2,function(x) {quantile(x,prob = 0.75,na.rm=T)})
      Ss_25_a[k-365, ] = apply(Ssmean_a,2,function(x) {quantile(x,prob = 0.25,na.rm=T)}) 
      Ss_median_a[k-365, ] = apply(Ssmean_a,2,function(x) {median(x,na.rm=T)})
      Ss_mean_a[k-365,]  = apply(Ssmean_a,2,function(x) {mean(x,na.rm=T)})
      
      # S01_ens10[k-365,] = Xa_S01[10,]  # ensemble member 10 for all grids
      # S01_ens20[k-365,] = Xa_S01[20,]  # ensemble member 20 for all grids
      # S01_ens30[k-365,] = Xa_S01[30,]
      # 
      # S02_ens10[k-365,] = Xa_S02[10,]  # ensemble member 10 for all grids
      # S02_ens20[k-365,] = Xa_S02[20,]  # ensemble member 20 for all grids
      # S02_ens30[k-365,] = Xa_S02[30,]
      
      QTOT_median[k-365, ] = apply(QTOT_f, 2, function(x) {median(x,na.rm=T)})
      QTOT_mean[k-365, ]  = apply(QTOT_f, 2, function(x) {mean(x,na.rm=T)})
      QTOT_75[k-365,]     = apply(QTOT_f, 2, function(x) {quantile(x,prob = 0.75,na.rm=T)})
      QTOT_25[k-365,]     = apply(QTOT_f, 2, function(x) {quantile(x,prob = 0.25,na.rm=T)})
      QTOT_min[k-365,] = apply(QTOT_f,2,function(x) {min(x,na.rm=T)} )
      QTOT_max[k-365,] = apply(QTOT_f,2,function(x) {max(x,na.rm=T)})
      
      QG_mean[k-365,] = apply(Qg_f, 2, function(x) {mean(x,na.rm=T)})
      QG_median[k-365,] = apply(Qg_f, 2, function(x) {median(x,na.rm=T)})
      QG_75[k-365,] = apply(Qg_f, 2, quantile,prob=0.75,na.rm=T)
      QG_25[k-365,] = apply(Qg_f, 2, quantile,prob=0.25,na.rm=T)
      
      QR_mean[k-365,] = apply(QR_f, 2, function(x) {mean(x,na.rm=T)})
      QR_median[k-365,] = apply(QR_f, 2, function(x) {median(x,na.rm=T)})
      QR_75[k-365,] = apply(QR_f, 2, quantile,prob=0.75,na.rm=T)
      QR_25[k-365,] = apply(QR_f, 2, quantile,prob=0.25,na.rm=T)
      
      ET_median[k-365, ] = apply(ET_f, 2, function(x) {median(x,na.rm=T)})
      ET_mean[k-365, ]  = apply(ET_f, 2, function(x) {mean(x,na.rm=T)})
      ET_75[k-365,]     = apply(ET_f, 2, function(x) {quantile(x,prob = 0.75,na.rm=T)})
      ET_25[k-365,]     = apply(ET_f, 2, function(x) {quantile(x,prob = 0.25,na.rm=T)})
      # ET_min[k-365,] = apply(ET_f,2,function(x) {min(x,na.rm=T)} )
      # ET_max[k-365,] = apply(ET_f,2,function(x) {max(x,na.rm=T)})
      
      SG_median[k-365, ] = apply(Xa_Sg, 2, function(x) {median(x,na.rm=T)})
      SG_mean[k-365, ]  = apply(Xa_Sg, 2, function(x) {mean(x,na.rm=T)})
      SG_75[k-365,]     = apply(Xa_Sg, 2, function(x) {quantile(x,prob = 0.75,na.rm=T)})
      SG_25[k-365,]     = apply(Xa_Sg, 2, function(x) {quantile(x,prob = 0.25,na.rm=T)})
      
      FSAT_median[k-365, ] = apply(FSAT_f, 2, function(x) {median(x,na.rm=T)})
      FSAT_mean[k-365, ]  = apply(FSAT_f, 2, function(x) {mean(x,na.rm=T)})
      FSAT_75[k-365,]     = apply(FSAT_f, 2, function(x) {quantile(x,prob = 0.75,na.rm=T)})
      FSAT_25[k-365,]     = apply(FSAT_f, 2, function(x) {quantile(x,prob = 0.25,na.rm=T)})
      
      
      
      SE_median[k-365, ] = apply(SE_f, 2, function(x) {median(x,na.rm=T)})
      SE_mean[k-365, ]  = apply(SE_f, 2, function(x) {mean(x,na.rm=T)})
      SE_75[k-365,]     = apply(SE_f, 2, function(x) {quantile(x,prob = 0.75,na.rm=T)})
      SE_25[k-365,]     = apply(SE_f, 2, function(x) {quantile(x,prob = 0.25,na.rm=T)})
      
      # wetness 
      # w0_f =  matrix(PARAMS$Fhru[1,] / PARAMS$S0FC[1,],Ne,N_GRID,byrow=T) * Xf_S01 +
      #  matrix(PARAMS$Fhru[2,] / PARAMS$S0FC[2,],Ne,N_GRID,byrow=T) * Xf_S02
      
      w0_f = hXf_w0_Correct
      w0_min_f[k-365,] = apply(w0_f,2,function(x) {min(x,na.rm=T)} )
      w0_max_f[k-365,] = apply(w0_f,2,function(x) {max(x,na.rm=T)})
      w0_75_f[k-365, ] = apply(w0_f,2,function(x) {quantile(x,prob = 0.75,na.rm=T)})
      w0_25_f[k-365, ] = apply(w0_f,2,function(x) {quantile(x,prob = 0.25,na.rm=T)}) 
      w0_median_f[k-365, ] = apply(w0_f,2,function(x) {median(x,na.rm=T)})
      w0_mean_f[k-365,]  = apply(w0_f,2,function(x) {mean(x,na.rm=T)})
      
      ## analysed wetness for this time step
      # max - min rescaling, make sure it's column wise operation 
      w0_a = sweep(sweep(S0mean_a,2,AWRA_S0min,"-"), 2, AWRA_S0max - AWRA_S0min, "/")
      w0_a[w0_a < 0] = 0
      w0_a[w0_a > 1] = 1
   
      w0_min_a[k-365,] = apply(w0_a, 2,function(x) {min(x,na.rm=T)} )
      w0_max_a[k-365,] = apply(w0_a, 2,function(x) {max(x,na.rm=T)})
      w0_75_a[k-365, ] = apply(w0_a, 2,function(x) {quantile(x,prob = 0.75,na.rm=T)})
      w0_25_a[k-365, ] = apply(w0_a, 2,function(x) {quantile(x,prob = 0.25,na.rm=T)}) 
      w0_median_a[k-365, ] = apply(w0_a, 2,function(x) {median(x,na.rm=T)})
      w0_mean_a[k-365,]  = apply(w0_a, 2,function(x) {mean(x,na.rm=T)})
      
      ## store ensemble member of wetness for rank histogram computation ##
      
      for (i in 1:N_GRID) {
        w0_ensemble[[i]] = rbind(w0_ensemble[[i]], w0_a[,i])
        Q_grid_ensemble[[i]] = rbind(Q_grid_ensemble[[i]],QTOT_f[,i])
      }
      
      ## catchment-wide QTOT
      QTOT_catch = apply(QTOT_f,1,mean)
      Qtot_ensemble = rbind(Qtot_ensemble, QTOT_catch)
      
      if (k%%400 ==0 ) {print(paste0("Ensemble run for time ", k, " is finished.")) } 
    }
  }
  
  return(list(
    SMAP_Err = SMAP_Err, 
    S01_median_f = S01_median_f, S01_median_a = S01_median_a,
    S02_median_f  = S02_median_f, S02_median_a  = S02_median_a,
    Ss1_median_f = Ss1_median_f, Ss1_median_a = Ss1_median_a,
    Ss2_median_f  = Ss2_median_f, Ss2_median_a  = Ss2_median_a,
    
    S0_min_f = S0_min_f, S0_max_f = S0_max_f, S0_75_f = S0_75_f, S0_25_f = S0_25_f,
    S0_median_f = S0_median_f, S0_mean_f = S0_mean_f, 
    
    Ss_min_f = Ss_min_f,  Ss_max_f = Ss_max_f, Ss_75_f = Ss_75_f,
    Ss_median_f = Ss_median_f, Ss_mean_f = Ss_mean_f ,
    
    S0_min_a = S0_min_a, S0_max_a =S0_max_a,
    S0_75_a = S0_75_a,  S0_25_a = S0_25_a,
    S0_median_a = S0_median_a, S0_mean_a = S0_mean_a,
    
    Ss_min_a = Ss_min_a, Ss_max_a = Ss_max_a, 
    Ss_75_a = Ss_75_a,  Ss_25_a = Ss_25_a, 
    Ss_median_a = Ss_median_a, Ss_mean_a = Ss_mean_a,
    
    QTOT_median = QTOT_median, QTOT_mean =QTOT_mean ,
    QTOT_75 = QTOT_75, QTOT_25 =QTOT_25 ,
    QTOT_max = QTOT_max, QTOT_min = QTOT_min, 
    
    ET_median = ET_median, ET_mean =ET_mean ,
    ET_75 = ET_75, ET_25 =ET_25 ,
    # ET_max = ET_max, ET_min = ET_min, 
    
    SG_median = SG_median, SG_mean = SG_mean ,
    SG_75 = SG_75, SG_25 = SG_25 ,
    
    FSAT_median = FSAT_median, FSAT_mean  = FSAT_mean, 
    FSAT_75 = FSAT_75,   FSAT_25 = FSAT_25 ,
    
    QG_mean = QG_mean, QG_median = QG_median, 
    QG_75 = QG_75,  QG_25 = QG_25,
    QR_mean = QR_mean, QR_median = QR_median,
    QR_75 = QR_75 ,  QR_25 = QR_25,
    SE_mean = SE_mean, SE_median = SE_median,
    SE_75 = SE_75 ,  SE_25 = SE_25,
    
    w0_min_f = w0_min_f, w0_max_f = w0_max_f, w0_75_f = w0_75_f, w0_25_f = w0_25_f,
    w0_median_f = w0_median_f, w0_mean_f = w0_mean_f, 
    w0_min_a = w0_min_a, w0_max_a = w0_max_a, w0_75_a = w0_75_a, w0_25_a = w0_25_a,
    w0_median_a = w0_median_a, w0_mean_a = w0_mean_a, 
    
    # S01_ens10 = S01_ens10 , S01_ens20 = S01_ens20, S01_ens30 = S01_ens30, 
    # S02_ens10 = S02_ens10 , S02_ens20 = S02_ens20, S01_ens30 = S02_ens30, 
    Kgain_10 = Kgain_10, Kgain_20 = Kgain_20, Kgain_30 = Kgain_30,
    w0_ensemble = w0_ensemble, Qtot_ensemble = Qtot_ensemble,
    Q_grid_ensemble = Q_grid_ensemble,
    mod_run = mod_run,
    SMAP_w_scaled =  SMAP_wetness_scaled,
    AWRA_OL_w = AWRA_wetness,
    # innovation 
    Innov_mean = Innov_mean ,  Innov_median = Innov_median,
    Innov_75 = Innov_75, Innov_25 = Innov_25,
    BIAS_S01 = BIAS_S01, BIAS_S02 = BIAS_S02, 
    BIAS_Ss1 = BIAS_Ss1, BIAS_SS2 = BIAS_SS2,
  ))
  
}
