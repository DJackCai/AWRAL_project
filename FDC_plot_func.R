FDC_plot_func <-function(num_quantiles,obs, OL, DA, mainstring )
{ 
  
  quantiles <- ppoints(num_quantiles)  ### pick percentiles to approximate the cdf
  quantile_OL <- quantile(OL,probs=quantiles, na.rm=T)
  quantile_DA <- quantile(DA,probs=quantiles, na.rm=T)
  quantile_observed <- quantile(obs,probs=quantiles, na.rm=T)
  suppressWarnings(plot(1-quantiles,quantile_observed,type="l",log = "y",
                        ylab="Log Streamflow (mm/day)",
                        xlab="Probability of exceedance",lwd=1.2,
                        main = mainstring, xlim = c(0,1)))
  suppressWarnings(lines(1-quantiles,quantile_OL,log="y",col="blue",lwd=1.2))
  suppressWarnings(lines(1-quantiles,quantile_DA,log="y",col="red",lwd=1.2))
  legend("topright", 
         c("Observed","OL","EnKF"),
         lwd=c(1.2,1.2),lty=c(1,1),col=c("black","blue","red"),bty="n")   }


## Example (9/6/2022): Flow duration curve analysis to see why 
## DA reduces overall bias but bad at other metrics 


FDC_plot_func(500, obs = Qobs_403222$Q[-c(1:97)], 
              OL = Qtot_403222_Qcalib_v9414[-c(1:97)], 
              DA = Qtot_403222_QEnKF[-c(1:97)], 
              mainstring = "Flow duration curve at Morse Ck (403222)")
