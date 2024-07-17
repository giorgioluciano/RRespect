library(pracma)
contSpec <- function(par) {
  
   # Carica le impostazioni globali se par non è fornito
  if (missing(par)) {
    par <- SetParameters()  # Carica le impostazioni globali
  }
  
  # Stampa il messaggio di caricamento se verbose è abilitato
  if (par$verbose) {
    cat('\n(*) Start\n(*) Loading Data File:', par$GstFile)
  }
  
  # Carica i dati sperimentali
  data_exp <- GetExpData(par$GstFile)
  t <- data_exp$t
  Gexp <- data_exp$Gt
  
  if (par$verbose) {
    cat('done\n(*) Initial Set up...')
  }
  
  # Imposta alcune variabili interne
  n <- length(t)
  ns <- par$ns  # discretizzazione di 'tau'
  
  tmin <- t[1]
  tmax <- t[n]
  
  if (par$FreqEnd == 1) {
    smin <- exp(-pi/2) * tmin
    smax <- exp(pi/2) * tmax
  } else if (par$FreqEnd == 2) {
    smin <- tmin
    smax <- tmax
  } else if (par$FreqEnd == 3) {
    smin <- exp(+pi/2) * tmin
    smax <- exp(-pi/2) * tmax
  } else {
    stop("Valore di FreqEnd non valido")
  }
  
  hs <- (smax / smin)^(1 / (ns - 1))
  s <- smin * hs^(0:(ns - 1))
  
  Hgs <- InitializeH(Gexp, t, s)
  
  # Trova il lambda ottimale con 'lcurve'
  if (par$verbose) {
    cat('done\n(*) Building the L-curve, ...')
  }
  
  if (par$lamC == 0) {
    result_lcurve <- lcurve(Gexp, Hgs, t, s, par$SmFacLam)
    lamC <- result_lcurve$lamC
    lam <- result_lcurve$lam
    rho <- result_lcurve$rho
    eta <- result_lcurve$eta
  } else {
    lamC <- par$lamC
  }
  
  
  if (par$verbose) {
    cat(sprintf('%e\n(*) Extracting the continuous spectrum, ...', lamC))
  }
  
  H <- LevenMarq(lamC, Gexp, Hgs, t, s)
  
  # Stampa alcuni file di dati
  if (par$verbose) {
    cat('done\n(*) Writing and Printing, ...')
    
    print(result_lcurve$rho_eta)
    
    if (par$lamC == 0) {
      write.table(data.frame(lam = lam, rho = rho, eta = eta), file = "output/rho-eta.dat", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    
    write.table(data.frame(s = s, H = H), file = 'output/H.dat', row.names = FALSE, col.names = FALSE)
    
    K <- kernel(H, t, s)
    write.table(data.frame(t = t, Gfit = K), file = 'output/Gfit.dat', row.names = FALSE, col.names = FALSE)
  }
  
  # Grafici
  if (par$plotting) {
    par(mfrow = c(2, 1))
    
    plot(s, H, type = 'o', xlab = 's', ylab = 'H(s)')
    title('H')
    
    K <- kernel(H, t, s)
    plot(t, Gexp, type = 'o', xlab = 't', ylab = 'Gt(exp)', main = 'G(t)')
    lines(t, K, col = 'black')
    legend('topright', legend = c('Gt(fit)'), col = 'black', lty = 1, cex = 0.8)
  }
  
  if (par$verbose) {
    cat('done\n(*) End\n')
  }
  
  # Restituisci lo spettro adattato H e opzionalmente lamC
  return(list(H = H, lamC = lamC))
}
