library(pracma)

source("I:/giorgioluciano.github.io/RRespect/ReSpect_v2/Gt/contSpec.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/SetParameters.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/GetExpData.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/InitializeH.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/LevenMarq.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/LCurve.r")


setwd("I:/giorgioluciano.github.io//GitHub/RRespect/ReSpect_v2/Gt")

par=SetParameters()
H=contSpec()

source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/DiscSpec.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/GetWeights.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/GridDensity.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/MaxwellModes.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/PlotMaxwellModes.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/SetParameters.r")

