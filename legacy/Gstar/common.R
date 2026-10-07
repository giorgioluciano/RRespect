# common.R
# Traduzione fedele di common.py (pyReSpect-freq)

library(pracma)  # per meshgrid

# ── GetExpData ────────────────────────────────────────────────────────────────
GetExpData <- function(fname) {
  data <- read.table(fname, header = FALSE)
  cols <- ncol(data)

  wo   <- data[, 1]
  Gpo  <- data[, 2]
  Gppo <- data[, 3]

  if (cols > 3) {
    wG1 <- data[, 4]
    wG2 <- data[, 5]
  }

  # rimuovi duplicati (come np.unique)
  idx  <- !duplicated(wo)
  wo   <- wo[idx]; Gpo <- Gpo[idx]; Gppo <- Gppo[idx]
  if (cols > 3) { wG1 <- wG1[idx]; wG2 <- wG2[idx] }

  if (cols == 3) {
    # interpola su 100 punti log-spaced
    w   <- exp(seq(log(min(wo)), log(max(wo)), length.out = 100))
    fp  <- approxfun(wo, Gpo,  rule = 2)
    fpp <- approxfun(wo, Gppo, rule = 2)
    Gp  <- fp(w)
    Gpp <- fpp(w)
    Gst <- c(Gp, Gpp)
    return(list(w = w, Gst = Gst, wexp = rep(1, length(Gst))))
  } else {
    return(list(w = wo, Gst = c(Gpo, Gppo), wexp = c(wG1, wG2)))
  }
}

# ── getKernMat ────────────────────────────────────────────────────────────────
getKernMat <- function(s, w) {
  ns  <- length(s)
  hsv <- numeric(ns)
  hsv[1]        <- 0.5 * log(s[2] / s[1])
  hsv[ns]       <- 0.5 * log(s[ns] / s[ns-1])
  hsv[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns-2)]))

  res <- meshgrid(s, w)
  S   <- res$X
  W   <- res$Y
  ws  <- S * W
  ws2 <- ws^2

  kernMat <- rbind(ws2 / (1 + ws2), ws / (1 + ws2)) *
             matrix(hsv, nrow = 2*length(w), ncol = ns, byrow = TRUE)
  return(kernMat)
}

# ── kernel_prestore ───────────────────────────────────────────────────────────
kernel_prestore <- function(H, kernMat, G0 = NULL) {
  Kh <- as.vector(kernMat %*% exp(H))
  if (!is.null(G0)) {
    n   <- nrow(kernMat) / 2
    G0v <- c(rep(G0, n), rep(0, n))
    return(Kh + G0v)
  }
  return(Kh)
}
