ErrorValCalcFun <- function(para) {
  aSosval <- para[1]
  bSosval <- para[2]
  
  # Call Lissajous_Plots_GenerationFun with parameters
  Lissajous_Plots_GenerationFun(c(aSosval, bSosval))
  
  # Call CreateMultiPanelPlotViscousFun to calculate error
  error2 <- CreateMultiPanelPlotViscousFun()
  
  return(error2)
}