# 7/2023: allowing an optional weight column in the input data file
#         improving encapsulation of functions

# Help to find continuous spectrum
# March 2019 major update:
# (i)   added plateau modulus G0 (also in pyReSpect-time) calculation
# (ii)  following Hansen Bayesian interpretation of Tikhonov to extract p(lambda)
# (iii) simplifying lcurve (starting from high lambda to low)
# (iv)  changing definition of rho2 and eta2 (no longer dividing by 1/n and 1/nl)

# 7/2023: allowing an optional weight column in the input data file + encapsulation of private functions
# 3/2019: adding G0 support stored as last gi in transport
#
#
# 12/15/2018: 
# (*) Introducing NNLS optimization of previous optimal solution
#     - can make deltaBaseWeightDist : 0.2 or 0.25 [coarser from 0.05]
#     - wrote new routines: FineTuneSolution, res_tG (for vector of residuals)
# 

initializeDiscSpec <- function(par) {
  # read input; initialize parameters
  if (par$verbose) {
    cat('\n(*) Start\n(*) Loading Data Files: ...', par$GexpFile, '...\n')
  }
  
  # Read experimental data
  data <- GetExpData(par$GexpFile)
  t <- data$t
  Gexp <- data$Gt
  wexp <- data$wG
  
  # Read the continuous spectrum
  fNameH <- 'output/H.dat'
  spectrum <- read.table(fNameH, header = TRUE)
  s <- spectrum[, 1]
  H <- spectrum[, 2]
  
  n <- length(t)
  
  # range of N scanned
  Nmin <- max(floor(0.5 * log10(max(t)/min(t))), 2)
  Nmax <- min(floor(3.0 * log10(max(t)/min(t))), n/4)
  
  if (par$MaxNumModes > 0) {
    Nmax <- min(Nmax, par$MaxNumModes)
  }
  
  Nv <- seq(Nmin, Nmax, by = 1)
  
  # Estimate Error Weight from Continuous Curve Fit
  kernMat <- getKernMat(s, t)
  kernMat <-t(kernMat)
  
  if (par$plateau) {
    f <- file(fNameH, open = "r")
    first_line <- scan(f, what = "", nlines = 1)
    close(f)
    if (length(first_line) > 2) {
      G0 <- as.numeric(first_line[length(first_line)])
    } else {
      stop("Problem reading G0 from H.dat; Plateau = True")
    }
    Gc <- kernel_prestore(H, kernMat, G0)
  } else {
    Gc <- kernel_prestore(H, kernMat)
  }
  
  Cerror <- 1 / sd(wexp * (Gc / Gexp - 1))
  
  return(list(t = t, Gexp = Gexp, wexp = wexp, s = s, H = H, Nv = Nv, Gc = Gc, Cerror = Cerror))
}


MaxwellModes <- function(z, t, Gexp, wexp, isPlateau) {
  tau <- exp(z)
  
  g_res <- nnLLS(t, tau, Gexp, wexp, isPlateau)
  g <- g_res$g
  error <- g_res$error
  condKp <- g_res$condKp
  
  if (isPlateau) {
    izero <- which(g[-length(g)] / max(g[-length(g)]) < 1e-7)
  } else {
    izero <- which(g / max(g) < 1e-7)
  }
  
  if (length(izero) > 0) {
    tau <- tau[-izero]
    g <- g[-izero]
  }
  
  return(list(g = g, tau = tau, error = error, condKp = condKp))
}

library(nnls)
nnLLS <- function(t, tau, Gexp, wexp, isPlateau) {
  K <- exp(-outer(t, tau, "/"))

  
  if (isPlateau) {
    K <- cbind(K, rep(1, length(Gexp)))
  }
  
  Kp <- K * (wexp / Gexp)
  condKp <- kappa(Kp, exact=TRUE)
  g <- nnls(Kp, wexp)$x
  
  GtM <- K %*% g
  error <- sum((wexp * (GtM / Gexp - 1))^2)
  
  return(list(g = g, error = error, condKp = condKp))
}

GetWeights <- function(H, t, s, wb) {
  ns <- length(s)
  n <- length(t)
  
  hs <- numeric(ns)
  wt <- (hs)
  
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns-1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns] / s[1:(ns-2)]))
  
  kern <- exp(-outer(t, s, "/"))
  wij <- kern %*% diag(hs * exp(H))
  K <- kern %*% (hs * exp(H))
  
  for (i in 1:n) {
    wij[i, ] <- wij[i, ] / K[i]
  }
  
  for (j in 0:(ns-1)) {
    wt[j+1] <- sum(wij[, j+1])
  }
  
  wt <- wt / trapz(log(s), wt)
  print(paste("wb:", wb))
  print(paste("mean(wt):", mean(wt)))
  print((1 - wb) * wt)
  print((wb * mean(wt)) * rep(1, length(wt)))
  wt <- (1 - wb) * wt + (wb * mean(wt)) * rep(1, length(wt))
  
  return(wt)
}
#z <- GridDensity(log(s), wt, Nopt)
GridDensity <- function(x, px, N) {
  npts <- 100
  xi <- seq(min(x), max(x), length.out = npts)
  fint <- splinefun(x, px, method = "natural")
  pint <- fint(xi)
  ci <- cumtrapz(xi, pint)
  pint <- pint / max(ci)
  ci <- ci / max(ci)
  
  alfa <- 1 / (N - 1)
  zij <- numeric(N + 1)
  z <- numeric(N)
  
  z <- numeric(N)
  z[1] <- min(x)
  z[N] <- max(x)
  
  beta <- seq(0.5, N - 0.5, by = 1) * alfa
  zij[1] <- z[1]
  zij[N + 1] <- z[N]
  fint_ci <- splinefun(ci, xi, method = "fmm")
  
  zij[2:N] <- fint_ci(beta)
  h <- diff(zij)
  
  beta <- seq(1, N - 1, by = 1) * alfa
  if (N > 2) {
    beta <- seq(1, N - 1, by = 1) * alfa
    z[2:(N-1)] <- fint_ci(beta)
  } # per questioni di indice in confronto a python dobbiamo aggiungere un if se N = 2
  
  return(list(z = z, h = h))
}



library(minpack.lm)

mergeModes_magic <- function(g, tau, imode) {
  # Definizione della funzione di costo
  costFcn_magic <- function(par, g, tau, imode) {
    gn   <- par[1]
    taun <- par[2]
    
    g1   <- g[imode]
    g2   <- g[imode + 1]
    tau1 <- tau[imode]
    tau2 <- tau[imode + 1]
    
    tmin <- min(tau1, tau2) / 10
    tmax <- max(tau1, tau2) * 10
    
    # Definizione della funzione per l'integrale
    normKern_magic <- function(t) {
      Gn <- gn * exp(-t / taun)
      Go <- g1 * exp(-t / tau1) + g2 * exp(-t / tau2)
      return((Gn / Go - 1)^2)
    }
    
    return(integrate(normKern_magic, tmin, tmax)$value)
  }
  
  # Stima iniziale
  iniGuess <- c(g[imode] + g[imode + 1], 0.5 * (tau[imode] + tau[imode + 1]))
  
  # Minimizzazione della funzione di costo
  res <- optim(iniGuess, costFcn_magic, g = g, tau = tau, imode = imode)
  
  # Aggiorna tau con il nuovo valore ottimizzato
  newtau <- tau[-(imode + 1)]
  newtau[imode] <- res$par[2]
  
  return(newtau)
}

FineTuneSolution <- function(tau, t, Gexp, wexp, isPlateau, estimateError = FALSE) {
  success <- FALSE
  
  # Definizione della funzione per i residui
  res_tG <- function(tau, texp, Gexp, wexp, isPlateau) {
    g <- nnLLS(texp, tau, Gexp, wexp, isPlateau)$g
    Gmodel <- rep(0, length(texp))
    
    for (j in 1:length(tau)) {
      Gmodel <- Gmodel + g[j] * exp(-texp / tau[j])
    }
    
    if (isPlateau) {
      Gmodel <- Gmodel + g[length(g)]
    }
    
    residual <- wexp * (Gmodel / Gexp - 1)
    
    return(residual)
  }
  
  tryCatch({
    # Minimizzazione non lineare
    control = nls.lm.control(ftol = 1e-10, ptol = 1e-10, gtol = 1e-10, maxiter = 500, nprint = 1)
    res <- nls.lm(par = tau, fn = res_tG, texp = t, Gexp = Gexp, wexp = wexp, isPlateau = isPlateau,control=control)
    tau <- res$par
    tau0 <- tau
    
    # Stima dell'errore
    if (estimateError) {
      J <- res$jacobian
      cov <- solve(t(J) %*% J) * mean(res$fvec^2)
      dtau <- sqrt(diag(cov))
    }
    
    success <- TRUE
  }, error = function(e) {
    
  })
  
  # Risoluzione finale con i modi di Maxwell
  result <- MaxwellModes(log(tau), t, Gexp, wexp, isPlateau)
  g <- result$g
  tau <- result$tau
  
  if (estimateError && success) {
    if (length(tau) < length(tau0)) {
      nkill <- 0
      for (i in 1:length(tau0)) {
        if (min(abs(tau0[i] - tau)) > 1e-12 * tau0[i]) {
          dtau <- dtau[-(i - nkill)]
          nkill <- nkill + 1
        }
      }
    }
    return(list(g = g, tau = tau, dtau = dtau))
  } else if (estimateError) {
    return(list(g = g, tau = tau, dtau = rep(-1, length(tau))))
  } else {
    return(list(g = g, tau = tau))
  }
}


getDiscSpecMagic <- function(par) {
  # Inizializza i dati
  data <- initializeDiscSpec(par)
  t <- data$t
  Gexp <- data$Gexp
  wexp <- data$wexp
  s <- data$s
  H <- data$H
  Nv <- data$Nv
  Gc <- data$Gc
  Cerror <- data$Cerror
  
  n <- length(t)
  ns <- length(s)
  npts <- length(Nv)
  
  # Range di wtBaseDist
  wtBase <- par$deltaBaseWeightDist * seq(1, 1/par$deltaBaseWeightDist-1)
  AICbst <- rep(0, length(wtBase))
  Nbst <- rep(0, length(wtBase))
  nzNbst <- rep(0, length(wtBase))  # Numero di modi non nulli
  
 
  for (ib in seq_along(wtBase)) {
    wb <- wtBase[ib]
    print(paste("ib:", ib, "wb:", wb))
    # Calcola i pesi
    wt <- GetWeights(H, t, s, wb)
    
    ev <- rep(0, npts)
    nzNv <- rep(0, npts)
    
    # Ciclo su Nv
    for (i in seq_along(Nv)) {
      N <- Nv[i]
      print(i)
      z <- GridDensity(log(s), wt, N)
      hz =z$h
      z <- z$z
      result <- MaxwellModes(z, t, Gexp, wexp, par$plateau)
      g <- result$g
      tau <- result$tau
      ev[i] <- result$error
      nzNv[i] <- length(g)
    }
    
    # Calcola AIC e salva i migliori risultati per questo wb
    AIC <- 2 * Nv + 2 * Cerror * ev
    AICbst[ib] <- min(AIC)
    Nbst[ib] <- Nv[which.min(AIC)]
    nzNbst[ib] <- nzNv[which.min(AIC)]
  }
  
  # Miglior impostazione globale di wb e Nopt
  Nopt <- as.integer(Nbst[which.min(AICbst)])
  wbopt <- wtBase[which.min(AICbst)]
  
  # Ricalcola i migliori dati e ottimizza
  wt <- GetWeights(H, t, s, wbopt)
  z <- GridDensity(log(s), wt, Nopt)
  h <-  z$h
  z <- z$z
  
  
  result <- MaxwellModes(z, t, Gexp, wexp, par$plateau)
  g <- result$g
  tau <- result$tau
  error <- result$error
  
  fineTuned <- FineTuneSolution(tau, t, Gexp, wexp, par$plateau, estimateError = TRUE)
  g <- fineTuned$g
  tau <- fineTuned$tau
  dtau <- fineTuned$dtau
  
  # Controlla se i modi sono abbastanza vicini per essere uniti
  if (length(tau) > 1) {
    indx <- order(tau)
    tau <- tau[indx]
    tauSpacing <- tau[-1] / tau[-length(tau)]
    itry <- 0
    
    if (par$plateau) {
      g[-1] <- g[indx]
    } else {
      g <- g[indx]
    }
    
    while (min(tauSpacing) < par$minTauSpacing && itry < 3) {
      imode <- which.min(tauSpacing)
      tau <- mergeModes_magic(g, tau, imode)
      fineTuned <- FineTuneSolution(tau, t, Gexp, wexp, par$plateau, estimateError = TRUE)
      g <- fineTuned$g
      tau <- fineTuned$tau
      dtau <- fineTuned$dtau
      tauSpacing <- tau[-1] / tau[-length(tau)]
      itry <- itry + 1
    }
    
    if (par$plateau) {
      G0 <- g[length(g)]
      g <- g[-length(g)]
    }
  }
  
  if (par$verbose) {
    cat("(*) Numero ottimale di nodi =", length(g), "\n")
  }
  
  # Plot
  if (par$plotting) {
    plot(wtBase, AICbst, type = "l", main = "AIC", xlab = "baseDistWt", log = "y")
    abline(v = wbopt, col = "gray")
    plot(tau, g, type = "b", log = "xy", main = "Discrete vs Continuous Spectrum", xlab = expression(tau), ylab = expression(g))
    lines(s, exp(H), col = "red")
    legend("bottomright", legend = c("Discrete", "Continuous"), col = c("black", "red"), lty = 1)
  }
  
  # Stampa i risultati
  if (par$verbose) {
    cat("(*) Numero ottimale di nodi =", length(g), "\n")
    if (par$plateau) {
      cat("(*) Modulo del Plateau:", G0, "\n")
      write.table(cbind(g, tau, dtau), file = "output/dmodes.dat", col.names = c("g", "tau", "dtau"))
    } else {
      write.table(cbind(g, tau, dtau), file = "output/dmodes.dat", col.names = c("g", "tau", "dtau"))
    }
  }
  
  return(list(Nopt = Nopt, g = g, tau = tau, error = error))
}

# Eseguire il programma principale
#if (interactive()) {
#  par <- readInput("inp.dat")
#  result <- getDiscSpecMagic(par)
#}


