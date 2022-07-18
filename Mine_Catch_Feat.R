### Functions that helps extract the catchment characteristics #### 
dates_1321 = dates_seq(2013,2021)
start_calib = which(dates_1321 == '2014-01-01')   # spin up starts in 2014
end_calib = which(dates_1321 == '2017-12-31')   
start_Vali = which(dates_1321 == '2018-01-01')   # spin up starts in 2014
end_Vali = which(dates_1321 == '2018-12-31')   

dates_9214 = dates_seq(1992,2014)
start_calib_9414 = which(dates_9214 == '1994-01-01')   # spin up starts in 2014
end_calib_9414 = which(dates_9214 == '2014-12-31')   

library(EcoHydRology); library(imputeTS)


annual_sum = function(full_ts,year1,year2) {
  years = format(DateSeq(year1,year2),"%Y")
  stopifnot(length(full_ts) == length(years))
  full_df = data.frame(Year = years, value = full_ts)
  annual_data = full_df%>%group_by(Year)%>%
    summarise(annual = sum(value,na.rm=T))
  return(annual_data$annual) 
  
}

mean_annual = function(full_ts,year1,year2) {
  years = format(DateSeq(year1,year2),"%Y")
  stopifnot(length(full_ts) == length(years))
  full_df = data.frame(Year = years, value = full_ts)
  annual_data = full_df%>%group_by(Year)%>%
    summarise(annual = sum(value,na.rm=T))
  mean_annual = mean(annual_data$annual)
  
  return( mean_annual)
  
}



catch_mining_framework = function(catid, cormethod = "spearman", 
                                  Q_thres = 0.005, Rain_thres = 0.01) {
  
  ### Collate a set of catchment characteristics that help 
  ### anticipate the data assimilation efficacy
  
  ## Read in shp data, ancillary data,streamflow data and SMAP data  
  STATIC_PAR <- as.matrix(read.csv(paste0("./forcings/",catid,"_resampled_static.csv")))
  N_GRID <- ncol(STATIC_PAR)
  shp_info = site_shp(id =  catid)
  area = shp_info$area_km2
  all_forcing_9414 <- grid_forcing_catchment_longterm(catname = catid,N_GRID = N_GRID ,
                                                      start = start_calib_9414, end= end_calib_9414)
  all_forcing_1418 <- grid_forcing_catchment(catname = catid,N_GRID = N_GRID ,
                                             start = 366, end= end_Vali)
  all_forcing_1518 <- grid_forcing_catchment(catname = catid,N_GRID = N_GRID ,
                                             start = 731, end= end_Vali)
  
  Qval_9414 = read.zoo(paste0("./Streamdata/Qobs_",catid,"_19942014.csv"),header=T,
                       sep=",",index.column = 1)$Q%>%as.vector()
  Qval_1518 = read.zoo(paste0("./Streamdata/Qobs_",catid,"_1518.csv"),header=T,
                       sep=",",index.column = 1)$Q%>%as.vector()
  
  SMAP_dat_site = as.matrix(read.csv(paste0("./forcings/",catid,"_resampled_SMAP1519.csv")))[,paste0("grid",1:N_GRID)]
  SMAP_dat_shift = apply(SMAP_dat_site,2, dplyr::lag, n= 1)
  SMAP_dat_shift_1518 = SMAP_dat_shift[1:1461,]
  
  ##### STAGE 1: features related to catchment runoff mechanisms ##### 
  
  ## 1) catchment humidity ###
  PG_9414 = all_forcing_9414$PG;  PG_1518 = all_forcing_1518$PG
  PG_9418 = rbind(PG_9414, PG_1518)
  meanPET  = as.vector(STATIC_PAR[5,])
  meanP   = colMeans(PG_9418)
  grid_humidity = meanP/meanPET  # long-term humidity at each grid 
  catch_humidity = mean(grid_humidity)
  
  ## 2) ftree #### 
  grid_ftree = STATIC_PAR[1,]
  catch_ftree = mean(grid_ftree)
  
  ## 3) mean annual total rain and streamflow #### 
  
  ## mean annual rainfall for the catchment (avg over grids)
  # grid_annual_rain = apply(PG_9418, 2, annual_sum, year1 = 1994, year2 = 2018)
  # catch_annual_rain = apply(grid_annual_rain,1,mean,na.rm=T)  ## annual rain ts
  
  grid_mean_annual_rain = apply(PG_9418, 2, mean_annual, year1 = 1994, year2 = 2018)
  catch_annual_meanP = mean(grid_mean_annual_rain)
  
  Qval_9418 = c(Qval_9414,Qval_1518)
  catch_annual_Q = annual_sum(Qval_9418,year1 = 1994, year2 = 2018)
  catch_meanQ = mean(catch_annual_Q)
  
  ### 4) Runoff coefficient (average daily ROC), subset fot specific days ### 
  
  minP_9418 = apply(PG_9418,1,min)
  P_Q_index = which(Qval_9418 >= Q_thres & minP_9418 > Rain_thres)
  
  PG_9418_subset = PG_9418[P_Q_index,]
  Qval_subset = matrix(Qval_9418[P_Q_index],length(P_Q_index),1)
  
  ROC_daily_9418_grids = 1/(sweep(PG_9418_subset, MARGIN = 1, Qval_subset, FUN = "/"))
  ROC_daily_9418_catch = apply(ROC_daily_9418_grids,1,mean,na.rm=T)
  ROC_daily_mean = mean(ROC_daily_9418_catch,na.rm=T)
  
  ### ROC for 2015-2018 
  
  minP_1518 = apply(PG_1518,1,min)
  P_Q_index_1518 = which(Qval_1518 >= Q_thres & minP_1518 > Rain_thres)
  PG_1518_sub = PG_1518[P_Q_index_1518, ]
  Qval_1518_sub = matrix(Qval_1518[P_Q_index_1518],length(P_Q_index_1518),1)
  
  ROC_daily_1518_grids = 1/(sweep(PG_1518_sub, MARGIN = 1, Qval_1518_sub, FUN = "/"))
  ROC_daily_1518_catch = apply(ROC_daily_1518_grids,1,mean,na.rm=T)
  
  #### 5) antecedent SM and Q coupling (Crow et al., 2017, 2018) #### 
  
  
  ## 22 June: don't subset the data 
  SMAP_lag1518 = apply(SMAP_dat_shift_1518,2,dplyr::lag,n=1)
  
  cor_SMAP_Q = apply(SMAP_dat_shift_1518, 2, function(x) { cor(x,Qval_1518,use='p',method = cormethod)   })
  lag1_SMAP_Q = apply(SMAP_lag1518,2,function(x) {cor(x,Qval_1518, use="p",method = cormethod)})
  
  catch_cor_SMAPQ = mean(cor_SMAP_Q,na.rm=T); catch_lag1_SMAPQ = mean(lag1_SMAP_Q,na.rm=T)
  
  cor_SMAP_ROC = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(SMAP_dat_shift_1518),
                        as.data.frame(ROC_daily_1518_grids))
  lag1_SMAP_ROC = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(SMAP_lag1518),
                         as.data.frame(ROC_daily_1518_grids))
  
  catch_cor_SMAPROC = mean(cor_SMAP_ROC, na.rm=T)
  catch_lag1_SMAPROC = mean(lag1_SMAP_ROC, na.rm=T)
  
  ### 6) Autocorrelation of Q - don't leave out #### 
  
  Q_fill = Qval_9418
  Q_fill[-c(which(Qval_9418>0.01))] = NA 
  
  Q_fill_1518 = Qval_1518
  Q_fill_1518[-c(which(Qval_1518>0.01))] = NA 
  
  acf1_Q_9418 = acf(Q_fill,lag.max = 1,plot = F,na.action = na.pass)[["acf"]][2]
  # acf1_Q_1518 = acf(Q_fill_1518,lag.max = 1,plot = F,na.action = na.pass)[["acf"]][2]
  
  #### 7) control of rainfall on Q ##### 
  # PG_9418_lag1 = apply(PG_9418,2,dplyr::lag,n=1)
  # PG_9418_lag2 = apply(PG_9418,2,dplyr::lag,n=2)
  # PG_9418_lag3 = apply(PG_9418,2,dplyr::lag,n=3)
  # 
  # Rain_sub = PG_9418[P_Q_index,]
  # Rain_sub_lag1 = PG_9418_lag1[P_Q_index,]
  # Rain_sub_lag2 = PG_9418_lag2[P_Q_index,]
  # Rain_sub_lag3 = PG_9418_lag3[P_Q_index,]
  # 
  # cor_Rain_Q  = apply(Rain_sub,2,function(x) {cor(x, Qval_subset, use="p",method = cormethod)})
  # lag1_Rain_Q = apply(Rain_sub_lag1,2,function(x) {cor(x,Qval_subset, use="p",method = cormethod)})
  # lag2_Rain_Q = apply(Rain_sub_lag2,2,function(x) {cor(x,Qval_subset, use="p",method = cormethod)})
  # lag3_Rain_Q = apply(Rain_sub_lag3,2,function(x) {cor(x,Qval_subset, use="p",method = cormethod)})
  # 
  # catch_cor_RQ = mean(cor_Rain_Q, na.rm=T)
  # catch_lag1_RQ = mean(lag1_Rain_Q,na.rm=T)
  # catch_lag2_RQ = mean(lag2_Rain_Q,na.rm=T)
  # catch_lag3_RQ = mean(lag3_Rain_Q,na.rm=T)
  # 
  
  ### 7) Estimated baseflow index after separation #### 
  
  Qval_9418_full = Qval_9418%>%na_interpolation(.,option = "spline")
  bf_sep = BaseflowSeparation(streamflow = as.vector(Qval_9418_full), 
                              filter_parameter = 0.925, passes = 3)
  
  BFI_9418 = sum(bf_sep$bt)/sum(Qval_9418_full)
  
  Qval_1518_full = Qval_1518%>%na_interpolation(.,option = "spline")
  bf_sep_1518 = BaseflowSeparation(streamflow = as.vector(Qval_1518_full), 
                                   filter_parameter = 0.925, passes = 3)
  BFI_1518 = sum(bf_sep_1518$bt)/sum(Qval_1518_full)
  
  #### STAGE 2: Features related to model representation and coupling #### 
  
  ### 1) connection between SM states and streamflow ##### 
  par_opt_Q = as.matrix(read.table(paste0(catid,"_Qopt_Q_Resamp_19942014.txt"),header=T))
  
  fit_mod_9414 <- Run_AWRAL_calib_Q(par_opt_Q, all_forcing = all_forcing_9414, 
                                    STATIC_PAR =STATIC_PAR, N_GRID)
  fit_mod_1518 <- Run_AWRAL_calib_Q(par_opt_Q, all_forcing = all_forcing_1418,  
                                    STATIC_PAR =STATIC_PAR, N_GRID)
  
  ## lag SM states #### 
  S0mean_9414 = fit_mod_9414$S0mean;     Ssmean_9414 = fit_mod_9414$Ssmean
  S0mean_1518 = fit_mod_1518$S0mean[-c(1:365), ]
  Ssmean_1518 = fit_mod_1518$Ssmean[-c(1:365), ]
  S0mean_9418 = rbind(S0mean_9414, S0mean_1518 )
  Ssmean_9418 = rbind(Ssmean_9414, Ssmean_1518)
  SG_9414 = fit_mod_9414$SG
  SG_1518 = fit_mod_1518$SG[-c(1:365),]
  SG_9418 = rbind(SG_9414, SG_1518)
  
  ## Compute lag series 
  S0mean_9418_lag = apply(S0mean_9418,2,dplyr::lag, n= 1)
  Ssmean_9418_lag = apply(Ssmean_9418,2,dplyr::lag, n= 1)
  S0mean_1518_lag = apply(S0mean_1518,2,dplyr::lag, n= 1)
  Ssmean_1518_lag = apply(Ssmean_1518,2,dplyr::lag, n= 1)
  
  Qtot_9414 = fit_mod_9414$QTOT;   
  Qtot_1518 = fit_mod_1518$QTOT[-c(1:365),]
  Qtot_9418 = rbind(Qtot_9414, Qtot_1518)
  Qsim_9414 = apply(fit_mod_9414$QTOT, 1, mean);  Qsim_1518 = apply(fit_mod_1518$QTOT[-c(1:365),], 1, mean)
  Qsim_9418 = c(Qsim_9414, Qsim_1518)
  
  P_Q_index_sim = which(Qsim_9418 >= Q_thres & minP_9418 > Rain_thres)
  
  Qsim_9418_sub = Qsim_9418[P_Q_index_sim]
  Qtot_9418_sim_sub = Qtot_9418[P_Q_index_sim, ]
  
  S0mean_9418_sub = S0mean_9418[P_Q_index_sim,]
  S0mean_9418_sub_lag = S0mean_9418_lag[P_Q_index_sim,]
  
  Ssmean_9418_sub = Ssmean_9418[P_Q_index_sim,]
  Ssmean_9418_sub_lag = Ssmean_9418_lag[P_Q_index_sim,]
  
  cor_S0_Q_9418 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(S0mean_9418_sub),
                         as.data.frame(Qtot_9418_sim_sub))
  cor_Ss_Q_9418 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(Ssmean_9418_sub),
                         as.data.frame(Qtot_9418_sim_sub))
  lag_S0_Q_9418 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(S0mean_9418_sub_lag),
                         as.data.frame(Qtot_9418_sim_sub))
  lag_Ss_Q_9418 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(Ssmean_9418_sub_lag),
                         as.data.frame(Qtot_9418_sim_sub))
  
  ### coupling from 2015-2018 ### 
  P_Q_index_sim_1518 = which(Qsim_1518 >= Q_thres & minP_1518 > Rain_thres)
  
  Qtot_1518_sim_sub = Qtot_1518[P_Q_index_sim_1518, ]
  # Qsim_1518_sub = Qsim_1518[P_Q_index_sim_1518]
  
  S0mean_1518_sub = S0mean_1518[P_Q_index_sim_1518,]
  S0mean_1518_sub_lag = S0mean_1518_lag[P_Q_index_sim_1518,]
  Ssmean_1518_sub = Ssmean_1518[P_Q_index_sim_1518,]
  Ssmean_1518_sub_lag = Ssmean_1518_lag[P_Q_index_sim_1518, ]
  
  cor_S0_Q_1518 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(S0mean_1518_sub),
                         as.data.frame(Qtot_1518_sim_sub))
  cor_Ss_Q_1518 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(Ssmean_1518_sub),
                         as.data.frame(Qtot_1518_sim_sub))
  lag_S0_Q_1518 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(S0mean_1518_sub_lag),
                         as.data.frame(Qtot_1518_sim_sub))
  lag_Ss_Q_1518 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(Ssmean_1518_sub_lag),
                         as.data.frame(Qtot_1518_sim_sub))
  
  ## catchment average for coupling 
  
  catch_cor_S0Q_9418 = mean(cor_S0_Q_9418)
  catch_cor_SsQ_9418 = mean(cor_Ss_Q_9418)
  catch_lag_S0Q_9418 = mean(lag_S0_Q_9418)
  catch_lag_SsQ_9418 = mean(lag_Ss_Q_9418)
  
  catch_cor_S0Q_1518 = mean(cor_S0_Q_1518)
  catch_cor_SsQ_1518 = mean(cor_Ss_Q_1518)
  catch_lag_S0Q_1518 = mean(lag_S0_Q_1518)
  catch_lag_SsQ_1518 = mean(lag_Ss_Q_1518)
  
  
  ### Coupling between SM but with ROC #####
  
  # PG_9418_sim_sub  = PG_9418[P_Q_index_sim, ]
  # 
  # ROC_9418_sim = mapply(function(x,y) {return(x/y) }, as.data.frame(Qtot_9418_sim_sub),
  #                       as.data.frame(PG_9418_sim_sub))
  # 
  # PG_1518_sim_sub  = PG_1518[P_Q_index_sim_1518, ]
  # ROC_1518_sim = mapply(function(x,y) {return(x/y) }, as.data.frame(Qtot_1518_sim_sub),
  #                       as.data.frame(PG_1518_sim_sub))
  # 
  # cor_S0_ROC_9418 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(S0mean_9418_sub),
  #                          as.data.frame(ROC_9418_sim))
  # cor_Ss_ROC_9418 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(Ssmean_9418_sub),
  #                          as.data.frame(ROC_9418_sim))
  # lag_S0_ROC_9418 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(S0mean_9418_sub_lag),
  #                          as.data.frame(ROC_9418_sim))
  # lag_Ss_ROC_9418 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(Ssmean_9418_sub_lag),
  #                          as.data.frame(ROC_9418_sim))
  # 
  # cor_S0_ROC_1518 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(S0mean_1518_sub),
  #                          as.data.frame(ROC_1518_sim))
  # cor_Ss_ROC_1518 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(Ssmean_1518_sub),
  #                          as.data.frame(ROC_1518_sim))
  # lag_S0_ROC_1518 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(S0mean_1518_sub_lag),
  #                          as.data.frame(ROC_1518_sim))
  # lag_Ss_ROC_1518 = mapply(function(x,y) {cor(x,y,use='p',method = cormethod)}, as.data.frame(Ssmean_1518_sub_lag),
  #                          as.data.frame(ROC_1518_sim))
  # 
  # catch_cor_S0ROC_9418 = mean(cor_S0_ROC_9418)
  # catch_cor_SsROC_9418 = mean(cor_Ss_ROC_9418)
  # catch_lag_S0ROC_9418 = mean(lag_S0_ROC_9418)
  # catch_lag_SsROC_9418 = mean(lag_Ss_ROC_9418)
  # 
  # catch_cor_S0ROC_1518 = mean(cor_S0_ROC_1518)
  # catch_cor_SsROC_1518 = mean(cor_Ss_ROC_1518)
  # catch_lag_S0ROC_1518 = mean(lag_S0_ROC_1518)
  # catch_lag_SsROC_1518 = mean(lag_Ss_ROC_1518)
  # 
  
  #### 2) importance of baseflow in runoff generation process ##### 
  
  QR_1518  = fit_mod_1518$QR[-c(1:365), ]
  QG_1518  = fit_mod_1518$QG[-c(1:365), ]
  QR_9418 = rbind(fit_mod_9414$QR, QR_1518)
  QG_9418 = rbind(fit_mod_9414$QG, QG_1518)
  
  Rof_9418 = QR_9418 + QG_9418
  Rof_1518 = QR_1518 + QG_1518
  
  ### Subset the runoff days 
  ROf_days_9418 = which(apply(Rof_9418,1,min) > 0)
  ROf_days_1518 = which(apply(Rof_1518,1,min) > 0)
  
  QR_9418_Qdays = QR_9418[ROf_days_9418,]
  QG_9418_Qdays = QG_9418[ROf_days_9418,]
  
  QR_1518_Qdays = QR_1518[ROf_days_1518,]
  QG_1518_Qdays = QG_1518[ROf_days_1518,]
  
  Rof_9418_Qdays = Rof_9418[ROf_days_9418,]
  Rof_1518_Qdays = Rof_1518[ROf_days_1518,]
  
  QG_prop_9418 = mapply(function(x,y) {sum(x)/sum(y)}, as.data.frame(QG_9418_Qdays),
                        as.data.frame(Rof_9418_Qdays))
  QG_prop_1518 = mapply(function(x,y) {sum(x)/sum(y)}, as.data.frame(QG_1518_Qdays),
                        as.data.frame(Rof_1518_Qdays))
  
  catch_QG_prop_9418  = mean(QG_prop_9418 , na.rm=T)
  catch_QG_prop_1518  = mean(QG_prop_1518 , na.rm=T)
  
  ### 3) importance of saturation excess vs infiltration excess (Crow et al., 2017) #### 
  
  QR_days_9418  =  which(apply(QR_9418,1,min) > 0) 
  QR_days_1518  =  which(apply(QR_1518,1,min) > 0) 
  
  SE_9414= fit_mod_9414$SE
  SE_1518  = fit_mod_1518$SE[-c(1:365), ]
  SE_9418 = rbind(SE_9414, SE_1518)
  
  SE_9418_QR = SE_9418[QR_days_9418,]
  SE_1518_QR= SE_1518[QR_days_1518,]
  
  QR_9418_QR= QR_9418[QR_days_9418,]
  QR_1518_QR= QR_1518[QR_days_1518,]
  
  #  SE_prop_9418 = mapply(function(x,y) {sum(x)/sum(y)}, as.data.frame(SE_9418_Qdays),
  #                       as.data.frame(QR_9418_Qdays))
  #  SE_prop_1518 = mapply(function(x,y) {sum(x)/sum(y)}, as.data.frame(SE_1518_Qdays),
  #                       as.data.frame(QR_1518_Qdays))
  #  catch_SE_prop_9418 = mean(SE_prop_9418 ,na.rm=T)
  #  catch_SE_prop_1518 = mean(SE_prop_1518,na.rm=T)
  
  SE_prop_9418 = mapply(function(x,y) {x/y}, as.data.frame(SE_9418_QR),
                        as.data.frame(QR_9418_QR))
  SE_prop_1518 = mapply(function(x,y) {x/y}, as.data.frame(SE_1518_QR),
                        as.data.frame(QR_1518_QR))
  SE_prop_grids_9418 = apply(SE_prop_9418,2,median ,na.rm=T)
  catch_SE_prop_9418 = mean(SE_prop_grids_9418,na.rm=T)
  SE_prop_grids_1518 = apply(SE_prop_1518,2,median,na.rm=T)
  catch_SE_prop_1518 = mean(SE_prop_grids_1518,na.rm=T)
  
  ### 4) vertical coupling between layer 1 & 2 with Sg #### 
  
  cor_SsSG_9418 = mapply(function(x,y) { return(cor(x,y,use = "p",method = cormethod)) }, 
                         as.data.frame(Ssmean_9418), as.data.frame(SG_9418))
  cor_SsSG_1518 = mapply(function(x,y) { return(cor(x,y,use = "p",method = cormethod)) }, 
                         as.data.frame(Ssmean_1518), as.data.frame(SG_1518))
  
  coup_SsSG_9418 =  mean(cor_SsSG_9418, na.rm=T)
  coup_SsQG_1518 =  mean(cor_SsSG_1518, na.rm=T)
  
  catch_summary = data.frame(
    ## basic features
    humidity = catch_humidity, ftree =catch_ftree, area = area, 
    meanP = catch_annual_meanP, meanQ = catch_meanQ,
    ROC = ROC_daily_mean, 
    ### catchment wetness vs streamflow coupling
    R_SMAPQ = catch_cor_SMAPQ, lag1_SMAP_Q = catch_lag1_SMAPQ,
    R_SMAPROC = catch_cor_SMAPROC, lag1_SMAPROC = catch_lag1_SMAPROC, 
    
    ### Routing impact + forcing control 
    # acf1_Q = acf1_Q_9418, acf1_Q_1518 = acf1_Q_1518, 
    # 
    # R_Qrain = catch_cor_RQ, lag1_Qrain = catch_lag1_RQ, 
    # lag2_Qrain = catch_lag2_RQ,lag3_Qrain = catch_lag3_RQ, 
    # 
    
    #### baseflow contribution 
    BFI_9418 = BFI_9418, BFI_1518 = BFI_1518, 
    
    ## Modelled connection between SM state and Q
    R_S0Qsim = catch_cor_S0Q_9418,     R_S0Qsim_1518 = catch_cor_S0Q_1518,
    R_SsQsim = catch_cor_SsQ_9418,     R_SsQsim_1518 = catch_cor_SsQ_1518,
    lag1_S0Qsim = catch_lag_S0Q_9418,  lag1_S0Qsim_1518 = catch_lag_S0Q_1518,
    lag1_SsQsim = catch_lag_SsQ_9418,  lag1_SsQsim_1518 = catch_lag_SsQ_1518,
    
    # R_S0ROC = catch_cor_S0ROC_9418,     R_S0ROC_1518 = catch_cor_S0ROC_1518,
    # R_SsROC = catch_cor_SsROC_9418,     R_SsROC_1518 = catch_cor_SsROC_1518,
    # lag1_S0ROC = catch_lag_S0ROC_9418,  lag1_S0ROC_1518 = catch_lag_S0ROC_1518,
    # lag1_SsROC = catch_lag_SsROC_9418,  lag1_SsROC_1518 = catch_lag_SsROC_1518,
    # 
    ## importance of baseflow and saturation excess (correctable)
    QG_prop = catch_QG_prop_9418, QG_prop_1518 = catch_QG_prop_1518,
    SE_prop = catch_SE_prop_9418,  SE_prop_1518 = catch_SE_prop_1518,
    
    ## vertical coupling
    R_SsSG = coup_SsSG_9418, R_SsSG_1518 = coup_SsQG_1518 )
  
}






