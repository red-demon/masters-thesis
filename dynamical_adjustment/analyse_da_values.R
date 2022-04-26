figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/dynamical_adjustment/Atlantic")
setwd(figpath)
values=read.csv("dynamical_adjustment_values.csv")
typeof(values)

values=values[order(values$df.cor.pred.scenario.anom),]
diff.rmse=values$df.rmse.residuals.ensmean-values$df.rmse.scenario.ensmean
diff.stn=values$df.signal.to.noise.dyn.adj-values$df.signal.to.noise.scenario
plot(values$df.cor.pred.scenario.anom,diff.rmse,type="l")
plot(values$df.cor.pred.scenario.anom,diff.stn,type="l")

