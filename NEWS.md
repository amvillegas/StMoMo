StMoMo v0.4.1.9000
----------------------------------------------------------------
* Adjust random walk code to work with missing values correctly


StMoMo v0.4.1
----------------------------------------------------------------
* Add Citation file for new JSS publication
* Correct minor bugs


StMoMo v0.4.0 
----------------------------------------------------------------
* Upgrade plotting function of forecast models (plot.forStMoMo)
  to include nice fan charts similar to those in the forecast
  package 
* Change inteface of fit function (fit.StMoMo) to allow 
  simiplified calling using newly defined StMoMoData class
* Add the possibility of using general arima models (beyond the
  random walk with drift) for the forecasting and simulation of
  the period indexes
* Add new estimation function for the Renshaw and Haberman model
  implementing the approximate constraint Hunt and Villegas (2015)

StMoMo v0.3.1
----------------------------------------------------------------
* Change residual plotting function (plot.resStMoMo) to use a
  more effective divergent colour palette for colourmaps.
* Change the default plotting limits in residual plotting 
  function to have a symmetric scale
* Fix formatting issues in logLik, AIC and BIC functions for 
  models fitted with gnm  


StMoMo v0.3.0
----------------------------------------------------------------
* Add optional parameters kt.lookback and gc.lookback to 
  functions forecast.fitStMoMo, simulate.fitStMoMo and 
  simulate.bootStMoMo. These new parameters allow the forecast
  and simulation of mortality models using only a limited history
  of the period and cohort indexes
  
* Improve documentation of some functions

* Update package vignette
  
* Add citation file


StMoMo v0.2.0
----------------------------------------------------------------

* First stable version of StMoMo released to CRAN