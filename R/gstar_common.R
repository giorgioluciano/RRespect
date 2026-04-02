# gstar_common.R — frequency-domain kernel utilities (internal package functions)
# Ported from Gstar/common.R; no library() / source() calls.

# Read G*(w) data: 3-column [w  G'  G''] or 5-column [w  G'  G''  wt1  wt2]
.gstarGetExpData <- function(fname) {
  data <- utils::read.table(fname, header = FALSE)
  cols <- ncol(data)

  wo   <- data[, 1]
  Gpo  <- data[, 2]
  Gppo <- data[, 3]

  if (cols > 3) {
    wG1 <- data[, 4]
    wG2 <- data[, 5]
  }

  idx  <- !duplicated(wo)
  wo   <- wo[idx];  Gpo <- Gpo[idx];  Gppo <- Gppo[idx]
  if (cols > 3) { wG1 <- wG1[idx]; wG2 <- wG2[idx] }

  if (cols == 3) {
    w   <- exp(seq(log(min(wo)), log(max(wo)), length.out = 100))
    fp  <- stats::approxfun(wo, Gpo,  rule = 2)
    fpp <- stats::approxfun(wo, Gppo, rule = 2)
    Gst <- c(fp(w), fpp(w))
    return(list(w = w, Gst = Gst, wexp = rep(1, length(Gst))))
  }
  list(w = wo, Gst = c(Gpo, Gppo), wexp = c(wG1, wG2))
}

# Kernel matrix: 2n x ns   rows 1:n → G', rows (n+1):2n → G''
.gstarGetKernMat <- function(s, w) {
  ns  <- length(s)
  hsv <- numeric(ns)
  hsv[1]        <- 0.5 * log(s[2] / s[1])
  hsv[ns]       <- 0.5 * log(s[ns] / s[ns - 1])
  hsv[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))

  ws  <- outer(w, s)             # n×ns: ws[i,j] = w[i]*s[j]
  ws2 <- ws^2

  rbind(ws2 / (1 + ws2), ws / (1 + ws2)) *
    matrix(hsv, nrow = 2 * length(w), ncol = ns, byrow = TRUE)
}

# K(H)(w) = kernMat %*% exp(H)  [+ G0 on G' rows if plateau]
.gstarKernelPrestore <- function(H, kernMat, G0 = NULL) {
  Kh <- as.vector(kernMat %*% exp(H))
  if (!is.null(G0)) {
    n   <- nrow(kernMat) / 2
    G0v <- c(rep(G0, n), rep(0, n))
    return(Kh + G0v)
  }
  Kh
}
