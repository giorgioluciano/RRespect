# gt_getweights.R — spectral weight computation, time domain (internal)
# Ported from Gt/GetWeights.R; no library() / source() calls.

.gtGetWeights <- function(H, t, s, wb) {
  ns <- length(s);  n <- length(t)

  hs <- numeric(ns)
  hs[1]        <- 0.5 * log(s[2] / s[1])
  hs[ns]       <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))

  kern  <- exp(-outer(t, 1/s))                   # n x ns: K[i,j] = exp(-t[i]/s[j])
  scale <- hs * exp(H)                            # ns-vector
  K     <- as.vector(kern %*% scale)              # n-vector
  wij   <- sweep(kern, 2, scale, "*")             # column-scale (no ns×ns diag alloc)
  wij   <- sweep(wij,  1, K,    "/")             # row-normalise (no loop)

  wt <- colSums(wij)

  wt <- wt / pracma::trapz(log(s), wt)
  wt <- (1 - wb) * wt + wb * mean(wt) * rep(1, ns)
  wt
}
