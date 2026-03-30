
# GetWeights_time.R
# Time-domain version: kernel is exp(-T/S)  [n x ns]

library(pracma)
source("common.R")

GetWeights <- function(H, t, s, wb) {
  ns <- length(s);  n <- length(t)

  hs <- numeric(ns)
  hs[1]        <- 0.5 * log(s[2] / s[1])
  hs[ns]       <- 0.5 * log(s[ns] / s[ns-1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns-2)]))

  res <- meshgrid(s, t)
  S   <- res$X;  T <- res$Y
  kern <- exp(-T / S)                              # n x ns

  wij <- kern %*% diag(hs * exp(H))               # n x ns
  K   <- kern %*% (hs * exp(H))                   # n x 1

  for (i in 1:n) wij[i, ] <- wij[i, ] / K[i]

  wt <- colSums(wij)

  # normalise + blend
  wt <- wt / pracma::trapz(log(s), wt)
  wt <- (1 - wb) * wt + wb * mean(wt) * rep(1, ns)
  return(wt)
}
