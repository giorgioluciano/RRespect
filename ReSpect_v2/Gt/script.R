library(pracma)

source("I:/giorgioluciano.github.io/RRespect/ReSpect_v2/Gstar/contSpec.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gstar/SetParameters.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gstar/GetExpData.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gstar/InitializeH.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gstar/LevenMarq.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gstar/LCurve.r")


#setwd("I:/giorgioluciano.github.io//GitHub/RRespect/ReSpect_v2/Gt")

par=SetParameters()
H=contSpec()

source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/DiscSpec.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/GetWeights.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/GridDensity.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/MaxwellModes.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/PlotMaxwellModes.r")
source("I:/giorgioluciano.github.io//RRespect/ReSpect_v2/Gt/SetParameters.r")

