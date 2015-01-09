StMoMo
========================================================

StMoMo (Sthocastic Mortality Modelling) is an R package providing functions to specify 
and fit stochastic mortality models including the Lee-Carter models, the CBD model 
the APC model. The package also includes tools for analysing the goodness of fit 
of the models and performing mortality projections and simulations

The package is still under development, but the current version provides has
many of the final functionalities. For instance you can you StMoMo to fit, 
forecast, and simulate a Lee-Carte model with the code:

 ```R
  #Fit a Lee-Carter model 
  library(StMoMo)
  LC <- lc()
  LCfit<-fit(LC, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext,
            ages = EWMaleData$ages, years = EWMaleData$years)
  plot(LCfit)
  plot(forecast(LCfit))
  
  #Simulate trajectories
  LCsim <- simulate(LCfit)
  par(mfrow=c(1, 2))
  plot(LCfit$years, LCfit$kt[1, ], xlim = range(LCfit$years, LCsim$kt.s$years),
       ylim = range(LCfit$kt, LCsim$kt.s$sim), type = "l",
       xlab = "year", ylab = "kt",
       main = "Lee-Carter: Simulated paths of the period index kt")
  matlines(LCsim$kt.s$years, LCsim$kt.s$sim[1, , ], type = "l", lty = 1)
  
  plot(LCfit$years, (LCfit$Dxt / LCfit$Ext)["65", ],
       xlim = range(LCfit$years, LCsim$years),
       ylim = range((LCfit$Dxt / LCfit$Ext)["65", ], LCsim$rates["65", , ]),
       type = "l", xlab = "year", ylab = "rate",
       main = "Lee-Carter: Simulated mortality rates at age 65")
  matlines(LCsim$years, LCsim$rates["65", , ], type = "l", lty = 1)

````

# Installing:

To instal the latest development version: 

```R
    install.packages("devtools")
    devtools::install_github("amvillegas/StMoMo")
````

# Main function:

### Definition and fitting of mortality models:
    * StMoMo: Creata a General Stochastic Mortality Model
    * lca: Create a Lee-Carter model
    * cbd: Create a CBD mortality model
    * apc: Create an Age-Period-Cohort model
    * fit.StMoMo: fit a stochastic mortality model

### Analysing and visualising parameters and goodness-of-fit:
    * plot.fitStMoMo: Plot the fitted parameter of a model
    * plot.resStMoMo: Plots of the residuals of mortality model

### Forecasting and simulation of mortality:
    * forecast.fitStMoMo: Forecast a mortality model
    * simulate.fitStMoMo: Simulate a mortality model
    * bootstrap.fitStMoMo: Bootstrap a mortality model
    

If you're interested in the package feel free to email andresmauriciovillegas@gmail.com
or track development at http://github.com/amvillegas/StMoMo
