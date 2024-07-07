setwd("I:/giorgioluciano.github.io/RRespect/ReSpect_v2/Gstar")
#r_files <- list.files(pattern = "\\.r$")
#lapply(r_files, source)

source("I:/giorgioluciano.github.io/RRespect/ReSpect_v2/Gstar/setParameters.r")
source("I:/giorgioluciano.github.io/RRespect/ReSpect_v2/Gstar/contSpec.r")
source("I:/giorgioluciano.github.io/RRespect/ReSpect_v2/Gstar/GetExpData.r")
source("I:/giorgioluciano.github.io/RRespect/ReSpect_v2/Gstar/InitializeH.r")
source("I:/giorgioluciano.github.io/RRespect/ReSpect_v2/Gstar/LevenMarq.r")

par =SetParameters()
H <- contSpec()

#source("InitializeH.r")


# kernelD <- function(H, w, s) {
#   ns <- length(s)
#   hs <- numeric(ns)
#   hs[1] <- 0.5 * log(s[2] / s[1])
#   hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
#   hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
#   
#   n <- length(w)
#   
#   ws <- outer(w, s, "*")
#   ws2 <- ws^2
#   Hsuper <- matrix(exp(H), nrow = 2 * n, ncol = ns, byrow = TRUE) * rep(hs, each = 2 * n)
#   
#   DK <- rbind(ws2 / (1 + ws2), ws / (1 + ws2)) * Hsuper
#   
#   return(DK)
# }
# 
# 
# 
# 
# kernel <- function(H, w, s) {
#   ns <- length(s)
#   hs <- numeric(ns)
#   
#   # Uses trapezoidal rule for integration
#   hs[1] <- 0.5 * log(s[2] / s[1])
#   hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
#   
#   hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
#   
#   # Create meshgrid equivalent
#   ws <- outer(w, s, "*")
#   ws2 <- ws^2
#   
#   # Calculate K
#   K <- c((ws2 / (1 + ws2)) %*% (hs * exp(H)), (ws / (1 + ws2)) %*% (hs * exp(H)))
#   
#   return(K)
# }

