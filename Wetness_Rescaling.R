### Scaling satellite/modelled SM estimates to wetness
### And bias correction to allow optimal assimilation

Compute_wetness = function(Original_ts, prob1 = 0.02, prob2= 0.98) {
  
  # Original_ts: can be either SMAP, AWRA or other models 
  
  if (is.vector(Original_ts)) {
    SM_wt = quantile(Original_ts, prob=prob1, na.rm=T)
    SM_fc = quantile(Original_ts, prob=prob2, na.rm=T)
    wetness = (Original_ts-SM_wt)/(SM_fc-SM_wt) } else if (is.data.frame(Original_ts) |is.matrix(Original_ts)){
      
      ## if matrix, compute for each column 
      wetness = apply(Original_ts, 2, function(x) {
        SM_wt = quantile(x, prob=prob1, na.rm=T)
        SM_fc = quantile(x, prob=prob2, na.rm=T)
        return( (x-SM_wt)/(SM_fc-SM_wt))
      }) 
    }
  return(wetness)
} 

### Bias Correction ####  

### Two common ways to correct bias in satellite data before data assimilation

## 1) Mean and variance matching 
## 2) CDF matching 

Mean_var_match <- function (original, target) {
  sd_target <- sd(target,na.rm=T)
  sd_orig   <- sd(original, na.rm=T)
  mean_orig <- mean(original,na.rm=T)
  mean_target <- mean(target, na.rm=T)
  return(  (sd_target/sd_orig)*(original-mean_orig) + mean_target  )
  
}

## vectorised version for multiple grids ###### 

## mapply: column by column operation for a matrix 

library(CDFt)  ## package that can help with CDF matching 

Obs_rescaling  <- function(OBS_df,MOD_df, method = "CDF") {
  ## df1: observation (NL * NGRID)
  ## df2: model space (NL * NGRID)
  if (method == "MVM") {
    rescaled_data  = mapply(function(x,y) { return(Mean_var_match(original = x,target = y)) },
                            as.data.frame(OBS_df), as.data.frame(MOD_df))
  } else if (method == "CDF") {
    
    rescaled_data  = mapply(function(x,y) { 
      cdfmatch_data <- CDFt(y,x,x)
      return(cdfmatch_data$DS) }, as.data.frame(OBS_df), 
      as.data.frame(MOD_df))
    
  }
  
  return(rescaled_data)
  
}
