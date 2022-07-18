## Parameter to set default parameter sets 
### Parameter set with static parameter replaced option (16 March)

PARAMS_set <- function(II,STATIC_PAR =NULL){
  ## II: the 2 x NGRID "template" for each grid and two HRUs
  ## STATIC_PAR: list object of static parameters 
  # ftree, 3 parameters for SzFC, and meanPET
  
  if (is.null(STATIC_PAR)) {
    PARAMS = list(
      Nhru        = 2,
      Fhru        = II * c(0.5,0.5),
      # LAI 
      SLA         = II * c(3,10),
      # FRACTIONAL VEGETATION COVER 
      LAIref      = II * c(2.5,1.4),
      # saturated fraction
      Sgref       = II * c(20,20),
      # storage at field capacity - wetness 
      S0FC        = II * c(30,30),
      SsFC        = II * c(200,200),
      SdFC        = II * c(1000,1000),
      # Fraction of day - Calculate E0 and other effective forcing
      fday        = II * c(0.5,0.5),    # 
      # Maximum transpiration
      Vc          = II * c(0.35,0.65),
      # shortwave radiation balance
      alb_dry     = II * c(0.26,0.26),
      alb_wet     = II * c(0.16,0.16),
      w0ref_alb   = II * c(0.3,0.3),
      
      # Groudn heat flux parameter
      Gfrac_max   = II * c(0.3,0.3),     
      fvegref_G   = II * c(0.22,0.22),
      
      # Aerodynamic conductance 
      hveg	    = II * c(10,0.5),
      # Root uptake 
      Us0         = II * c(6,6),   ## physiological maximum root water update 
      Ud0         = II * c(4,0),
      wslimU      = II * c(0.3,0.3),    # uptake limiting relative water content
      wdlimU      = II * c(0.3,0.3),
      # Canopy conductance 
      cGsmax      = II * c(0.03,0.03),
      # soil evaporation
      FsoilEmax   = II * c(0.2,0.5),
      w0limE      = II * c(0.85,0.85),
      # Open water evaporation
      FwaterE     = II * c(0.7,0.7), 
      
      # Rainfall interception evaporation
      S_sls       = II * c(0.1,0.1),
      ER_frac_ref = II * c(0.2,0.05),
      
      # Surface runoff
      InitLoss    = II * c(5,5),
      PrefR       = II * c(150,150),
      
      # Drainage (S0,Ss,Sd)
      FdrainFC    = II * c(0.029,0.029),
      beta        = II * c(4.5,4.5),
      # Capillary rise 
      Fgw_conn    = II * c(1,1),
      # Routing for groundwater and streamflow discharge 
      K_gw        = 0.06,      # gets replaced through spatialising step below
      K_rout      = 0.77,      # gets replaced through spatialising step below
      
      # VEGETATION ADJUSTMENT FOR EQUILIBRIUM BIOMASS
      LAImax      = II * c(8,8),
      Tgrow       = II * c(1000,150),
      Tsenc       = II * c(60,10)) }     else {
        
        ### If static parameters for all grids are supplied ####
        
        if (is.list(STATIC_PAR)) {
          ftree_val <- STATIC_PAR$ftree
          S0FC_val <- STATIC_PAR$S0FC
          SsFC_val <- STATIC_PAR$SsFC
          SdFC_val <- STATIC_PAR$SdFC   } else if (is.matrix(STATIC_PAR)) {
            
            ftree_val <- STATIC_PAR[1,]
            S0FC_val <- STATIC_PAR[2,]
            SsFC_val <- STATIC_PAR[3,]
            SdFC_val <- STATIC_PAR[4,]
          }
        
        PARAMS = list(
          Nhru        = 2,
          Fhru        = rbind(ftree_val,1-ftree_val), 
          # LAI 
          SLA         = II * c(3,10),
          # FRACTIONAL VEGETATION COVER 
          LAIref      = II * c(2.5,1.4),
          # saturated fraction
          Sgref       = II * c(20,20),
          # storage at field capacity - wetness 
          S0FC        = rbind(S0FC_val,S0FC_val),
          SsFC        = rbind(SsFC_val,SsFC_val),
          SdFC        = rbind(SdFC_val,SdFC_val),
          # Fraction of day - Calculate E0 and other effective forcing
          fday        = II * c(0.5,0.5),    # 
          # Maximum transpiration
          Vc          = II * c(0.35,0.65),
          # shortwave radiation balance
          alb_dry     = II * c(0.26,0.26),
          alb_wet     = II * c(0.16,0.16),
          w0ref_alb   = II * c(0.3,0.3),
          
          # Groudn heat flux parameter
          Gfrac_max   = II * c(0.3,0.3),     
          fvegref_G   = II * c(0.22,0.22),
          
          # Aerodynamic conductance 
          hveg	    = II * c(10,0.5),
          # Root uptake 
          Us0         = II * c(6,6),   ## physiological maximum root water update 
          Ud0         = II * c(4,0),
          wslimU      = II * c(0.3,0.3),    # uptake limiting relative water content
          wdlimU      = II * c(0.3,0.3),
          # Canopy conductance 
          cGsmax      = II * c(0.03,0.03),
          # soil evaporation
          FsoilEmax   = II * c(0.2,0.5),
          w0limE      = II * c(0.85,0.85),
          # Open water evaporation
          FwaterE     = II * c(0.7,0.7), 
          
          # Rainfall interception evaporation
          S_sls       = II * c(0.1,0.1),
          ER_frac_ref = II * c(0.2,0.05),
          
          # Surface runoff
          InitLoss    = II * c(5,5),
          PrefR       = II * c(150,150),
          
          # Drainage (S0,Ss,Sd)
          FdrainFC    = II * c(0.029,0.029),
          beta        = II * c(4.5,4.5),
          # Capillary rise 
          Fgw_conn    = II * c(1,1),
          # Routing for groundwater and streamflow discharge 
          K_gw        = 0.06,      # gets replaced through spatialising step below
          K_rout      = 0.77,      # gets replaced through spatialising step below
          
          # VEGETATION ADJUSTMENT FOR EQUILIBRIUM BIOMASS
          LAImax      = II * c(8,8),
          Tgrow       = II * c(1000,150),
          Tsenc       = II * c(60,10) )
      }
  return(PARAMS)
  
}

 empty_forcing <- function() {
  FORCING  = list( Pg   = c(),
                   Rg   = c(),
                   Ta   = c(),
                   pe   = c(),
                   pair = c(),
                   u2   = c() )
  return(FORCING)
}

 
 

PARAMS_OPT_GRID <- function(N_GRID, par_mat, STATIC_PAR) {
  # par_mat: N * p matrix of the grid-specific optimised parameter
  # grid-specific: STATIC_PAR needs to be a list ###
  
  II = rbind(rep(1,N_GRID), rep(1,N_GRID))
  Init_PARAMS <- PARAMS_set(II,STATIC_PAR = STATIC_PAR)
  
  #####  Update the initial parameters by the pars to be optimised  #####
  
  Init_PARAMS$Vc <- rbind(par_mat[,1],II[2,] * 0.65)
  
  Init_PARAMS$Gfrac_max <-rbind(II[1,]*0.3,par_mat[,2]) 
  
  Init_PARAMS$Us0 =   rbind(II[1,]*6,par_mat[,3])     
  Init_PARAMS$Ud0 = rbind(par_mat[,4],   II[2,] * 0)   
  Init_PARAMS$FsoilEmax = rbind(par_mat[,5], par_mat[,6])
  Init_PARAMS$S_sls =   rbind(par_mat[,7], II[2,] *0.1)
  Init_PARAMS$FdrainFC = rbind(par_mat[,8], par_mat[,9])
  Init_PARAMS$beta <- rbind(par_mat[,10], par_mat[,11])
  Init_PARAMS$LAImax <- rbind(II[1,] * 8,   par_mat[,12])   
  Init_PARAMS$Tgrow <- rbind(II[1,] * 1000,   par_mat[,13])    
  
  PARAMS = Init_PARAMS
  return(PARAMS)
  
}

### Common optimised parameter for the entire catchment 

PARAMS_OPT_CATCH = function(N_GRID, par_mat, STATIC_PAR) {
  # par_mat: here, is a common parameter vector applied to all grids in the catchment 
  
  II = as.matrix(rbind(rep(1,N_GRID),rep(1,N_GRID)))
  Init_PARAMS <- PARAMS_set(II,STATIC_PAR = STATIC_PAR)
  
  #####  Update the initial parameters by the pars to be optimised  #####
  
  ###### in the Fhru format
  
  Init_PARAMS$Vc <- II*c(par_mat[1],0.65)
  Init_PARAMS$Gfrac_max <- II*c(0.3,par_mat[2])
  
  Init_PARAMS$Us0 = II*c(6,par_mat[3])
  Init_PARAMS$Ud0 = II*c(par_mat[4],0)
  Init_PARAMS$FsoilEmax = II * c(par_mat[5],par_mat[6])
  Init_PARAMS$S_sls = II * c(par_mat[7],0.1)
  # Init_PARAMS$ER_frac_ref = II*c(pars[9],0.05)
  Init_PARAMS$FdrainFC = II * c(par_mat[8],par_mat[9])
  Init_PARAMS$beta = II * c(par_mat[10],par_mat[11])
  Init_PARAMS$LAImax = II * c(8,par_mat[12])
  Init_PARAMS$Tgrow = II * c(1000,par_mat[13])
  
  PARAMS = Init_PARAMS
  return(PARAMS)
  
}

