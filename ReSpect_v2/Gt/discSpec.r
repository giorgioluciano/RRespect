# Function: discSpec
#
# Uses the continuous relaxation spectrum extracted using contSpec()
# to determine an approximate discrete approximation.
#
# Input: Communicated by the datastructure "par"
#
#        fGstFile = name of file that contains G*(w) in 3 columns [w Gp Gpp]
#                   default: 'Gst.dat' is assumed.
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
#         [g tau] = spectrum
#         error   = error norm of the discrete fit
#        
#         dmodes.dat : Prints the [g tau] for the particular Nopt
#         Nopt.dat   : If Nopt not supplied, then optimum [N error(N) cond(N)]
#         Gfitd.dat  : The discrete G* for Nopt [w Gp Gpp] 
#

discSpec <- function(par = NULL) {
  if (is.null(par)) {
    par <- SetParameters()  # Load in global settings
  }
  
  if (par$verbose) {
    cat("\n(*) Start\n(*) Loading Data Files: ...\n")
  }
  
  data <- ReadData(par$GstFile, 'output/H.dat')
  t <- data$t
  Gt <- data$Gt
  s <- data$s
  H <- data$H
  
  n <- length(t)
  ns <- length(s)
  
  # Find the distribution of nodes you need
  wt <- GetWeights(H, t, s)
  wt <- wt / pracma::trapz(log(s), wt)
  wt <- (1 - par$BaseDistWt) * wt + (par$BaseDistWt * mean(wt)) * rep(1, length(wt))
  
  # Try different N: number of Maxwell modes
  Nmax <- min(floor(3 * log(max(t) / min(t))), n / 4)
  Nmin <- max(floor(0.5 * log10(max(t) / min(t))), 3)
  
  npts  = Nmax - Nmin + 1;
  Nv <- Nmin:Nmax
  ev <- numeric(length(Nv))
  condN <- numeric(length(Nv))
  
  for (i in 1:length(npts)) {
    N <- Nv[i]
    grid <- GridDensity(log(s), wt, N) #check 
    z <- grid$z
    hz <- grid$hz
    
    maxwell_modes <- MaxwellModes(z, t, Gt, par$prune)
    
    g <- maxwell_modes$g
    tau <- maxwell_modes$tau
    ev[i] <- maxwell_modes$error
    condN[i] <- maxwell_modes$cond
  }
  
  
  # Use supplied number of modes or
  if (par$Nopt > 0) {
    Nopt <- par$Nopt
    if (par$verbose) {
      cat(sprintf("\n(*) Using %d number of discrete modes\n", Nopt))
    }
  } else {
    emin <- min(ev)
    condNmin <- min(condN)
    cost <- (1 - par$condWt) * (ev - emin)^2 + par$condWt * (log(condN / condNmin))^2
    Nopt <- Nv[which.min(cost)]
    if (par$verbose) {
      cat(sprintf("\n(*) Number of optimum nodes = %d\n", Nopt))
    }
  }
  
  # Send the best data-set stats
  grid <- GridDensity(log(s), wt, Nopt)
  z <- grid$z
  hz <- grid$hz
  maxwell_modes <- MaxwellModes(z, t, Gt, par$prune)
  g <- maxwell_modes$g
  tau <- maxwell_modes$tau
  error <- maxwell_modes$error
  
  # Some Plotting
  if (par$plotting) {
   
    
    plot(tau, g, type = "o", log = "xy", xlab = "tau", ylab = "g", main = "discrete spectrum")
    
    # [AX,H1,H2] = plotyy(Nv,ev,Nv,log10(condN),'plot');
    # set(AX(1));
    # set(AX(2));
    # set(get(AX(1),'Ylabel'),'String','Error');
    # set(get(AX(2),'Ylabel'),'String','log(cond)');
    # xlabel('number of modes, N')
    # set(H1,'LineStyle','-');
    # set(H2,'LineStyle','--');
    # hold on
    # plot(Nopt,ev[which(Nv==Nopt)],'r*')
    # hold off
    
    PlotMaxwellModes(g, tau, t, Gt)
    # hold on 
    # K = kernel(H,t,s);
    # loglog(w,K[1:n],'r-','LineWidth',2); hold on
    # loglog(w,K[n+1:2*n],'r-','LineWidth',2);
    #xlab('t')
    #ylab('G(t)')
    # hold off
  }
  
  # Some Printing
  if (par$verbose) {
    f1 <- file('output/dmodes.dat', 'w')
    cat(sprintf("(*) Condition number of matrix equation: %e\n", maxwell_modes$cond), file = f1)
    
    cat("\n\t\tModes\n\t\t-----\n\n", file = f1)
    cat("i \t    g(i) \t    tau(i)\n", file = f1)
    cat("---------------------------------------\n", file = f1)
    for (i in 1:length(g)) {
      cat(sprintf('%d \t %9.5e \t %9.5e\n', i, g[i], tau[i]), file = f1)
    }
    cat("\n", file = f1)
    close(f1)
    
    if (par$Nopt == 0) {
      f2 <- file('output/Nopt.dat', 'w')
      for (i in 1:length(Nv)) {
        cat(sprintf('%e\t%e\t%e\n', Nv[i], ev[i], condN[i]), file = f2)
      }
      close(f2)
    }
    
    f3 <- file('output/Gfitd.dat', 'w')
    K <- kernel(H, t, s)
    for (i in 1:n) {
      cat(sprintf('%e\t%e\t%e\n',t[i], K[i], K[n + i]), file = f3)
    }
    close(f3)
  }
  
  return(list(Nopt = Nopt, g = g, tau = tau, error = error))
}

# Read Data Files
ReadData <- function(fNameGst, fNameH) {
  # Leggi i dati di input
  exp_data <- GetExpData(fNameGst)
  t <- exp_data$t
  Gt <- exp_data$Gt
  
  # Leggi lo spettro continuo
  data <- read.table(fNameH, header = FALSE)
  s <- data[, 1]
  H <- data[, 2]
  
  return(list(t = t, Gt = Gt, s = s, H = H))
}

