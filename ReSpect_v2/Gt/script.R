library(pracma)

source("~/GitHub/RRespect/ReSpect_v2/Gt/contSpec.r")
source("~/GitHub/RRespect/ReSpect_v2/Gt/SetParameters.r")
source("~/GitHub/RRespect/ReSpect_v2/Gt/GetExpData.r")
source("~/GitHub/RRespect/ReSpect_v2/Gt/InitializeH.r")
source("~/GitHub/RRespect/ReSpect_v2/Gt/LevenMarq.r")
source("~/GitHub/RRespect/ReSpect_v2/Gt/LCurve.r")


setwd("~/GitHub/RRespect/ReSpect_v2/Gt")

par=SetParameters()
H=contSpec()

source("~/GitHub/RRespect/ReSpect_v2/Gt/DiscSpec.r")
source("~/GitHub/RRespect/ReSpect_v2/Gt/GetWeights.r")
source("~/GitHub/RRespect/ReSpect_v2/Gt/GridDensity.r")
source("~/GitHub/RRespect/ReSpect_v2/Gt/MaxwellModes.r")
source("~/GitHub/RRespect/ReSpect_v2/Gt/PlotMaxwellModes.r")
source("~/GitHub/RRespect/ReSpect_v2/Gt/SetParameters.r")

