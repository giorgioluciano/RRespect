
# common_time.R
# Utilities for time-domain (G(t)) relaxation spectrum extraction
# Translated from common.py (pyReSpect-time)

library(pracma)

# Kernel matrix: n x ns   K_(i,j) = hs_j * exp(-t_i / s_j)
getKernMat <- function(s, t) {
  ns <- length(s)
  hs <- numeric(ns)
  hs[1]        <- 0.5 * log(s[2] / s[1])
  hs[ns]       <- 0.5 * log(s[ns] / s[ns-1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns-2)]))

  res <- meshgrid(s, t)   # X = s (cols), Y = t (rows)
  S   <- res$X;  T <- res$Y
  return(exp(-T / S) * matrix(hs, nrow = length(t), ncol = ns, byrow = TRUE))
}

# K(H)(t) = kernMat %*% exp(H)  [+ G0 if plateau]
kernel_prestore <- function(H, kernMat, G0 = NULL) {
  Kh <- as.vector(kernMat %*% exp(H))
  if (!is.null(G0)) Kh <- Kh + G0
  return(Kh)
}

# Read G(t) data: 2-column [t  Gt]  or 3-column [t  Gt  wt]
GetExpData <- function(fname) {
  data <- read.table(fname, header = FALSE)
  cols <- ncol(data)
  to  <- data[,1];  Gto <- data[,2]
  if (cols >= 3) wGo <- data[,3]

  # remove duplicates
  idx <- !duplicated(to)
  to <- to[idx];  Gto <- Gto[idx]
  if (cols >= 3) wGo <- wGo[idx]

  if (cols == 2) {
    # interpolate onto 100 equi-log points
    t  <- exp(seq(log(min(to)), log(max(to)), length.out = 100))
    Gt <- approx(to, Gto, xout = t, rule = 2)$y
    return(list(t = t, Gt = Gt, wexp = rep(1, length(t))))
  } else {
    return(list(t = to, Gt = Gto, wexp = wGo))
  }
}
