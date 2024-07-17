# Function: discSpec(par)
#
# Uses the continuous relaxation spectrum extracted using contSpec()
# to determine an approximate discrete approximation.
#
# Input: Communicated by the datastructure "par"
#
#        fGstFile = name of file that contains G(t) in 2 columns [t Gt]
#                   default: 'Gt.dat' is assumed.
#        verbose  = 1, then prints onscreen messages, and prints datafiles
#        plotting = 1, then plots to stdio.
#
#        prune    = 1, then tries to kill modes with -ve g(i)
#  
#        Nopt     = optional argument, if you want to get a spectrum for a
#                   specified number of modes. If absent it will use some
#                   heuristic algorithm to figure out an optimum.
#        
# In addition par.BaseDistWt and par.condWt control the blending of the
# flat profile, and the weight of the "condition number" in determining
# the optimal Nopt. Both are set to 0.5 by default.
#
#
# Output: Nopt    = optimum number of discrete modes
#         g tau   = spectrum
#         error   = error norm of the discrete fit
#        
#         dmodes.dat : Prints the [g tau] for the particular Nopt
#         Nopt.dat   : If Nopt not supplied, then optimum [N error(N) cond(N)]
#         Gfitd.dat  : The discrete G(t) for Nopt [t Gt] 

discSpec <- function(par) {
  
  # Add appropriate subdirectories to search path
  # This is not directly translatable to R, assumed it's handled outside the function
  
  # Load global settings if par is not provided
  if (missing(par)) {
    par <- SetParameters()  # Load in global settings
  }
  
  # Print loading message if verbose is enabled
  if (par$verbose) {
    cat('\n(*) Start\n(*) Loading Data Files: ...')
  }
  
  # Read data files
  data <- ReadData(par$GstFile, 'output/H.dat')
  t <- data$t
  Gt <- data$Gt
  s <- data$s
  H <- data$H
  
  n <- length(t)
  ns <- length(s)
  
  # Find the distribution of nodes you need
  wt <- GetWeights(H, t, s)
  wt <- wt / trapz(log(s), wt)
  wt <- (1 - par$BaseDistWt) * wt + (par$BaseDistWt * mean(wt)) * rep(1, length(wt))
  
  # Try different N: number of Maxwell modes
  Nmax <- min(floor(3 * log(max(t) / min(t))), n / 4)  # maximum N
  Nmin <- max(floor(0.5 * log10(max(t) / min(t))), 3)  # minimum N
  
  Nv <- seq(Nmin, Nmax, by = 1)
  npts <- length(Nv)
  
  ev <- numeric(npts)
  condN <- numeric(npts)
  
  for (i in seq_len(npts)) {
    N <- Nv[i]
    data_tau <- GridDensity(log(s), wt, N)
    z <- data_tau$z
    hz <- data_tau$h
    
    data_gtau <- MaxwellModes(z, t, Gt, par$prune)
    g <- data_gtau$g
    tau <- data_gtau$tau
    ev[i] <- data_gtau$error
    condN[i] <- data_gtau$condKp
  }
  
  
  # Determine optimal number of modes (Nopt)
  if (par$Nopt > 0) {
    Nopt <- par$Nopt
    
    if (par$verbose) {
      cat(sprintf('\n(*) Using %d number of discrete modes\n', Nopt))
    }
  } else {
    emin <- min(ev)
    condNmin <- min(condN)
    
    cost <- (1 - par$condWt) * (ev - emin)^2 + par$condWt * (log(condN / condNmin))^2
    print(cost)
    idx_min <- which.min(cost)
    Nopt <- Nv[idx_min]
    
    if (par$verbose) {
      cat(sprintf('\n(*) Number of optimum nodes = %d\n', Nopt))
    }
  }
  
  # Get data for optimal Nopt
  data_tau <- GridDensity(log(s), wt, Nopt)
  z <- data_tau$z
  
  data_gtau <- MaxwellModes(z, t, Gt, par$prune)
  g <- data_gtau$g
  tau <- data_gtau$tau
  error <- data_gtau$error
  
  # Plotting
  if (par$plotting) {
    par(mfrow = c(2, 1))
    
    plot(tau, g, type = 'o', xlab = 'tau', ylab = 'g')
    
    PlotMaxwellModes(g, tau, t, Gt)
    
    K <- kernel(H, t, s)
    lines(t, K, col = 'red', lwd = 2)
  }
  
  # Printing
  if (par$verbose) {
    f1 <- file('output/dmodes.dat', 'w')
    cat(sprintf('(*) Condition number of matrix equation: %e\n', data_gtau$condKp))
    cat('\n\t\tModes\n\t\t-----\n\n')
    cat('i \t    g(i) \t    tau(i)\n')
    cat('---------------------------------------\n')
    for (i in seq_along(g)) {
      cat(sprintf('%d \t %9.5e \t %9.5e\n', i, g[i], tau[i]))
      cat(f1, sprintf('%d \t %9.5e \t %9.5e\n', i, g[i], tau[i]))
    }
    cat('\n')
    close(f1)
    
    if (par$Nopt == 0) {
      f2 <- file('output/Nopt.dat', 'w')
      for (i in seq_along(Nv)) {
        cat(f2, sprintf('%e\t%e\t%e\n', Nv[i], ev[i], condN[i]))
      }
      close(f2)
    }
    
    f3 <- file('output/Gfitd.dat', 'w')
    K <- kernel(H, t, s)
    for (i in seq_along(t)) {
      cat(f3, sprintf('%e\t%e\n', t[i], K[i]))
    }
    close(f3)
  }
  
  # Return outputs
  return(list(Nopt = Nopt, g = g, tau = tau, error = error))
}


ReadData <- function(fNameGst, fNameH) {
  # Read input data
  
  
  # Leggere i dati s e H dal file fNameH
  data <- read.table(fNameH)  # Se i dati sono in formato tabella senza intestazioni
  s <- data[, 1]  # Prima colonna come vettore s
  H <- data[, 2]  # Seconda colonna come vettore H
  
  # Chiamata alla funzione GetExpData con il primo file
  result <- GetExpData(fNameGst)
  
  # Estrazione dei risultati
  t <- result$t
  Gt <- result$Gt
  
  # Ritorna i risultati
  return(list(t = t, Gt = Gt, s = s, H = H))
}


