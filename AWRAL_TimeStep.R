AWRAL_TimeStep = function(forcing, state, par) {
  #% Australian Water Resources Assessment Landscape (AWRA-L) model
  #% This script contains the AWRA-L time step model
  
  ### "state" contains the state variables (updated at each time step)
  ### "forcing" contains the forcing data
  ### "out" contains the output variables
  
  # 18 June: where it gets updated - store infiltration excess and saturation excess 
  # to improve understanding of hydrologic behaviour
  
  #Assign forcing variables
  
  Pg   = forcing$Pg       # Daily precip (mm)
  Rg   = forcing$Rg       # Daytime average SW Radiation (W m**-2)
  Ta   = forcing$Ta       # Air temp (deg C)
  pe   = forcing$pe       # Actual vapour pressure (Pa)
  pair = forcing$pair     # Air pressure (Pa)
  u2   = forcing$u2       # 2-m Wind speed (m s**-1)
  
  #Assign state variables
  
  # should give Ss and S0 first  
  S0    = state$S0
  Ss    = state$Ss
  Sd    = state$Sd
  
  S0t   = state$S0[1,] # topsoil storage under trees
  S0g   = state$S0[2,] # topsoil storage under grass
  
  Sst   = state$Ss[1,] # shallow soil storage under trees
  Ssg   = state$Ss[2,] # shallow soil storage under grass
  
  Sdt   = state$Sd[1,] # deep soil storage under trees
  Sdg   = state$Sd[2,] # deep soil storage under grass
  Sg    = state$Sg
  Sr    = state$Sr
  Mleaf = state$Mleaf
  
  # Assign Parameters
  Nhru   = par$Nhru	#Number of HRU units
  Fhru   = par$Fhru	
  SLA    = par$SLA
  LAIref = par$LAIref
  Sgref  = par$Sgref
  S0FC   = par$S0FC
  SsFC   = par$SsFC
  SdFC   = par$SdFC
  fday   = par$fday	# currently static ... to be updated 
  Vc     = par$Vc
  alb_dry   = par$alb_dry
  alb_wet   = par$alb_wet 	# (3.2)
  w0ref_alb = par$w0ref_alb    	# (3.2)
  Gfrac_max = par$Gfrac_max    	# (3.5)
  fvegref_G = par$fvegref_G    	# (3.5)
  hveg   = par$hveg        	# (3.7)
  Us0    = par$Us0         	# (2.3)
  Ud0    = par$Ud0         	# (2.3)
  wslimU = par$wslimU
  wdlimU = par$wdlimU
  cGsmax = par$cGsmax
  FsoilEmax   = par$FsoilEmax    	# (4.5)
  w0limE      = par$w0limE     	# (4.5)
  FwaterE     = par$FwaterE   	# (4.7)
  S_sls       = par$S_sls  	# (4.2)
  ER_frac_ref = par$ER_frac_ref  	# (4.2)
  InitLoss    = par$InitLoss     	# (2.2)
  PrefR       = par$PrefR      	# (2.2)
  FdrainFC    = par$FdrainFC     	# (2.4)
  beta        = par$beta    	# (2.4)
  Fgw_conn    = par$Fgw_conn     	# (2.6)
  K_gw        = par$K_gw  	# (2.5)
  K_rout      = par$K_rout	# (2.7)
  LAImax      = par$LAImax	# (5.5)
  Tgrow       = par$Tgrow		# (5.4)
  Tsenc       = par$Tsenc		# (5.4)
  
  # Diagnostic equations
  LAI   = SLA * Mleaf    		# (5.3)
  fveg  = 1. - exp(-LAI / LAIref)	# (5.3)
  fsoil = 1. - fveg
  
  ### S0 needs to be a vector -- loaded first
  
  w0 = S0 / S0FC     		# (2.1)
  ws = Ss / SsFC     		# (2.1)
  wd = Sd / SdFC     		# (2.1)
  
  # Spatialise catchment fractions
  
  fwater = rbind(pmin.int(0.005,0.007 * Sr^0.75),pmin.int(0.005,0.007 * Sr^0.75))
  
  # Sgref suits for unified value ##  
  fsat   = rbind(pmin.int(1, pmax.int(pmin.int(0.005, 0.007 * Sr^0.75), Sg / Sgref)),
                 pmin.int(1, pmax.int(pmin.int(0.005, 0.007 * Sr^0.75), Sg / Sgref)))
  Sghru  = rbind(Sg,Sg)
  
  # CALCULATION OF PET
  # Conversions and coefficients (3.1)
  pes 		= 610.8 * exp(17.27 * Ta /(237.3 + Ta))     # saturated vapour pressure
  fRH 		= pe / pes;     # relative air humidity
  cRE 		= 0.03449 + 4.27e-5 * Ta
  Caero 		= fday * 0.176 * (1 + Ta / 209.1) * (pair - 0.417 * pe)*(1 - fRH)
  keps 		= 1.4e-3 * ((Ta / 187) * (Ta / 187) + Ta / 107 + 1) * (6.36 * pair + pe) / pes
  Rgeff 		= Rg / fday
  
  # shortwave radiation balance (3.2)
  alb_veg 	= 0.452 * Vc
  alb_soil	= alb_wet + (alb_dry - alb_wet) * exp(-w0 / w0ref_alb)
  alb 		= fveg * alb_veg + fsoil * alb_soil
  RSn 		= (1 - alb) * Rgeff
  
  # long wave radiation balance (3.3 to 3.5)
  StefBolz 	= 5.67e-8
  Tkelv 	 	= Ta + 273.16
  RLin	   	= (0.65 * (pe / Tkelv)^0.14) * StefBolz * Tkelv * Tkelv * Tkelv * Tkelv #(3.3)
  RLout 	 	= 1 * StefBolz * Tkelv * Tkelv * Tkelv * Tkelv                          #(3.4)
  RLn      	= RLin - RLout
  fGR 		= Gfrac_max * (1. - exp(-fsoil / fvegref_G))  		                #(3.5)
  
  # ENERGY AVAILABLE FOR SENSIBLE & GROUND HEAT FLUX
  Rneff 	 	= (RSn + RLn) * (1. - fGR)
  
  # Aerodynamic conductance (3.7)
  fh       = log(813. / hveg - 5.45)
  ku2      = 0.305 / (fh * (fh + 2.3))
  ga       = ku2*u2
  
  #### Finally: calculate Potential evaporation
  kalpha      = 1. + Caero * ga / Rneff
  E0          = cRE * (1. / (1. + keps)) * kalpha * Rneff * fday
  E0[E0<0]  = 0.0 ## prevents E0 from becoming negative
  
  # CALCULATION OF ET FLUXES AND ROOT WATER UPTAKE
  # Root water uptake constraint (4.4)
  Usmax       = Us0 * pmin.int(ws / wslimU, 1.) # issue 24 Feb 
  Udmax       = Ud0 * pmin.int(wd / wdlimU, 1.)
  U0          = pmax.int(Usmax, Udmax)  # invalid input type here?? - becomes a vector 
  
  # Maximum transpiration (4.3)
  Gsmax       = cGsmax * Vc
  gs          = fveg * Gsmax
  ft          = 1 / (1 + (keps / (1 + keps)) * ga / gs)
  Etmax       = ft * E0
  
  # Actual transpiration (4.1)
  Et          = pmin.int(U0, Etmax)
  # Root water uptake distribution (2.3)
  Us          = pmin.int( (Usmax / (Usmax + Udmax + 1.e-2)) * Et, Ss - 1e-2 )
  Ud          = pmin.int( (Udmax / (Usmax + Udmax + 1.e-2)) * Et, Sd - 1e-2)
  
  Et          = Us + Ud      # to ensure mass balance
  
  Et[Et<0]    = 0.0
  ### here it requires a matrix format from pmin - 
  # because pmin.int returns a vector 
  #Et          = matrix(Us + Ud,nrow=2,byrow=T)   
  
  # Soil evaporation (4.5)
  fsoilE      = FsoilEmax * pmin.int(w0 / w0limE, 1.)
  Es          = (1 - fsat) * fsoilE * (E0 - Et)
  
  # Groundwater evaporation (4.6)
  Eg          = (fsat - fwater) * FsoilEmax * (E0 - Et)
  
  # Open water evaporation (4.7)
  Er          = fwater * FwaterE * (E0 - Et)
  
  # Rainfall interception evaporation (4.2)
  Sveg        = S_sls * LAI
  fER         = ER_frac_ref * fveg
  Pwet        = -log(1 - fER / fveg) * Sveg / fER
  Ei          = (Pg < Pwet) * fveg * Pg + 
    (Pg >= Pwet) * (fveg * Pwet + fER * (Pg - Pwet))
  
  # CALCULATION OF WATER BALANCES
  # surface water fluxes (2.2)
  Pn     	    = Pg - Ei - InitLoss
  Pn[Pn<0] = 0.0
  Rhof        = (1 - fsat) * ( Pn / (Pn + PrefR) ) * Pn   # infiltration excess part 
  Rsof        = fsat * Pn                                 # saturation excess part 
  QR          = Rhof + Rsof
  I           = Pg - Ei - QR
  
  # SOIL WATER BALANCES (2.1 & 2.4)
  # Topsoil water balance (S0)
  S0          = S0  + I - Es
  S0[S0 < 0] = 0.0
  SzFC        = S0FC
  Sz          = S0
  wz          = pmax.int(Sz,1.e-2) / SzFC
  fD          = (wz > 1) * pmax.int(FdrainFC, 1 - 1 / wz) + 
    (wz <= 1) * FdrainFC * exp(beta * (wz - 1))
  Dz          = pmin.int(fD * Sz, Sz - 1e-2)
  Dz[Dz < 0] = 0.0
  D0          = Dz
  S0          = S0 - D0
  
  # Shallow root zone water balance (Ss)
  Ss          = Ss + D0 - Us
  Ss[Ss < 0] = 0.0
  SzFC        = SsFC
  Sz          = Ss
  wz          = pmax.int(Sz, 1.e-2) / SzFC
  fD          = (wz > 1) * pmax.int(FdrainFC, 1 - 1 / wz) + (wz <= 1) * FdrainFC * exp(beta * (wz - 1))
  Dz          = pmin.int(fD * Sz, Sz - 1e-2)
  Dz[Dz < 0] = 0.0
  Ds          = Dz
  Ss          = Ss - Ds
  
  # Deep root zone water balance (Sd) (2.6)
  Sd          = Sd + Ds - Ud
  Sd[Sd < 0] = 0.0
  SzFC        = SdFC
  Sz          = Sd
  wz          = pmax.int(Sz,1.e-2) / SzFC
  fD          = (wz > 1) * pmax.int(FdrainFC, 1 - 1 / wz) + (wz <= 1) * FdrainFC * exp(beta * (wz - 1))
  Dz          = pmin.int(fD * Sz, Sz - 1e-2)
  Dz[Dz < 0] = 0.0
  Dd          = Dz
  Sd          = Sd - Dd
  Y           = pmin.int(Fgw_conn * pmax.int(wdlimU * SdFC - Sd, 0.0), Sghru)
  Sd          = Sd + Y
  
  
  # CATCHMENT WATER BALANCE
  # Groundwater store water balance (Sg) (2.5)
  
  ## Note: no HRU separation here - use colSums
  NetGf       =   colSums(Fhru * (Dd - Eg - Y))
  Sg          =   Sg + NetGf
  Sg[Sg<0] = 0.0
  Qg          =   pmin.int(Sg, (1 - exp(-K_gw)) * Sg)
  Sg          =   Sg - Qg
  
  # Surface water store water balance (Sr) (2.7)
  Sr          =   Sr +  colSums(Fhru * (QR - Er)) + Qg
  Sr[Sr<0] = 0.0
  Qtot        =   pmin.int(Sr, (1 - exp(-K_rout)) * Sr)
  Sr          =   Sr - Qtot
  
  # VEGETATION ADJUSTMENT (5) - update the vegetation biomass at next step
  fveq        =   (1. / pmax.int((E0 / U0) - 1, 1e-3)) * (keps / (1 + keps)) * (ga / Gsmax)
  fvmax       =   1. - exp(-LAImax / LAIref)
  fveq        =   pmin.int(fveq, fvmax)
  dMleaf      =   -log(1. - fveq) * LAIref / SLA-Mleaf
  Mleafnet    =  (dMleaf>0) * (dMleaf / Tgrow) + (dMleaf < 0) * dMleaf / Tsenc 
  Mleaf       =   Mleaf + Mleafnet
  
  ### Updating diagnostic variables!!
  LAI         =   SLA * Mleaf    		# (5.3)
  fveg        =   1. - exp(-LAI / LAIref) # (5.3)
  fsoil       =   1. - fveg
  w0          =   S0 / S0FC     # (2.1)
  ws          =   Ss / SsFC     # (2.1)
  wd          =   Sd / SdFC     # (2.1)
  
  #### ASSIGN OUTPUT VARIABLES ####
  
  
  ETtot      = Es + Eg + Er+ Ei+ Et
  out = list( 
    ### fluxes
    #E0      =   colSums(Fhru * E0),
    #Ee      =   colSums(Fhru * (Es + Eg + Er + Ei)),
    Et      =   colSums(Fhru * Et),
    Ei      =   colSums(Fhru * Ei),
    #Etot    =   colSums(Fhru * (Es + Eg + Er + Ei)) + colSums(Fhru * Et),
    Qtot    =   Qtot,
    Qg      =   Qg,
    QR  =  QR, 
    gwflux  =   NetGf,
    SE = Rsof,      # saturation excess  Rsof = fsat * Pn                                
    IE = Rhof ,     # infiltration excess: Rhof  = (1 - fsat) * ( Pn / (Pn + PrefR) ) * Pn
    
    #	D       =   colSums(Fhru * Dd),
    
    # HRU specific drainage
    Eg      =   colSums(Fhru * Eg),
    D1    = Dd[1],
    D2    = Dd[2],
    Et1   = Et[1],
    Et2   = Et[2],
    ## Evapotranspiration flux 
    ETtot      = Es + Eg + Er+ Ei+ Et,
    #	ET1   = (Es + Eg + Er+ Ei+ Et)[1],#ETtot[1],
    #	ET2   = (Es + Eg + Er+ Ei+ Et)[2],#ETtot[2],
    
    # states
    # 1) water balance state 
    S0      =    S0,
    Ss      =    Ss,
    Sd      =    Sd,
    Sg      =    Sg,
    Sr     =    Sr,
    #	Stot    =   S0 + Ss + Sd + Sg + Sr + colSums(Fhru * Mleaf * 4),  
    
    # 2) vegetation states 
    
    Mleaf   =   Mleaf,
    Mleaf_grid =  colSums(Fhru *Mleaf ),
    LAI    =   colSums(Fhru * LAI),  # LAI for the grid
    
    fveg    =   colSums(Fhru * fveg),
    
    ## 3) Satellite equivalent states
    #	albedo  =   colSums(Fhru * alb),
    EVI     =   colSums(Fhru * (Vc * fveg + 0.07)),      # assume 0.07 is EVI for bare soil
    fsat    =   colSums(Fhru * fsat),
    wunsat  =   colSums(Fhru * w0)  , 
    S0mean =     colSums(Fhru * S0),   
    ETmean = colSums(Fhru * ETtot)
  )
  return(out)
}