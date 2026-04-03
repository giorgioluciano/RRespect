
# GetWeights.R


source("common.R")

GetWeights <- function(H, w, s, wb) {
  ns <- length(s)
  n  <- length(w)

  hs <- numeric(ns)
  hs[1]        <- 0.5 * log(s[2] / s[1])
  hs[ns]       <- 0.5 * log(s[ns] / s[ns-1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns-2)]))

  res <- meshgrid(s, w)
  S   <- res$X;  W <- res$Y
  ws  <- S * W;  ws2 <- ws^2

  kern <- rbind(ws2 / (1 + ws2), ws / (1 + ws2))  # 2n x ns
  wij  <- kern %*% diag(hs * exp(H))               # 2n x ns
  K    <- kern %*% (hs * exp(H))                    # 2n x 1

  for (i in 1:(2*n)) wij[i, ] <- wij[i, ] / K[i]

  wt <- colSums(wij)

  # normalise (trapz over log(s))
  wt <- wt / pracma::trapz(log(s), wt)

  # blend with uniform distribution
  wt <- (1 - wb) * wt + wb * mean(wt) * rep(1, ns)

  return(wt)
}
