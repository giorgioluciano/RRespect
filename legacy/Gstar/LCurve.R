# LCurve.R
# L-curve bayesiana per la scelta ottimale di lambda

# Source con path relativo alla posizione di questo file
this_dir <- dirname(sys.frame(1)$ofile)
source(file.path(this_dir, "LevenMarq.R"), encoding = "UTF-8")

lcurve <- function(par, s, Gstar, kernMat, H0, G0 = NULL) {

  lams <- exp(seq(log(par$lam_min), log(par$lam_max), length.out = par$lam_pts))
  n    <- length(Gstar) / 2

  res_norm <- numeric(par$lam_pts)
  reg_norm <- numeric(par$lam_pts)
  H_list   <- vector("list", par$lam_pts)
  G0_list  <- vector("list", par$lam_pts)

  for (i in seq_along(lams)) {
    fit <- levenberg_marquardt(par, s, Gstar, kernMat, H0, lams[i], G0)
    H_list[[i]] <- fit$H
    if (!is.null(G0)) G0_list[[i]] <- fit$G0

    Gpred <- kernMat %*% fit$H
    if (!is.null(G0)) Gpred[1:n] <- Gpred[1:n] + fit$G0

    res_norm[i] <- log(sum((Gpred - Gstar)^2))
    reg_norm[i] <- log(sum(fit$H^2))
  }

  # Curvatura della L-curve (differenze finite secondo ordine)
  x  <- res_norm
  y  <- reg_norm
  dx <- diff(x);  dy <- diff(y)
  ddx <- diff(dx); ddy <- diff(dy)
  idx <- seq_along(ddx)
  kappa <- abs(dx[idx] * ddy - dy[idx] * ddx) /
           (dx[idx]^2 + dy[idx]^2)^1.5

  best_i <- which.max(kappa) + 1

  if (par$verbose)
    cat(sprintf("L-curve: lambda ottimale = %.4e (indice %d)\n",
                lams[best_i], best_i))

  out <- list(lam = lams[best_i], H = H_list[[best_i]])
  if (!is.null(G0)) out$G0 <- G0_list[[best_i]]
  return(out)
}
