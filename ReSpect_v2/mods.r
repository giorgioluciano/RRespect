# Funzione comune per calcolare lo stress
calculate_stress <- function(time, PreSum, Integrand, gi, taui, tstep) {
  dim <- length(time)
  tau12 <- numeric(dim)
  
  for (i in 1:dim) {
    time2 <- seq(time[1], time[i], by = tstep)
    PreSumn <- PreSum(time2, time[i])
    integrand2 <- sapply(time2, function(tp) Integrand(tp, time[i], PreSumn[which(time2 == tp)], gi, taui))
    
    if (i == 1) {
      tau12[i] <- 0
    } else if (i == 2) {
      tau12[i] <- (integrand2[1] + integrand2[2])/2 * tstep
    } else if (i %% 2 == 1) {
      tau12[i] <- SimpsonOdd(integrand2, tstep)
    } else {
      tau12[i] <- SimpsonOdd(integrand2[1:(length(integrand2)-3)], tstep) + 
        SimpsonThreeEight(integrand2[(length(integrand2)-3):length(integrand2)], tstep)
    }
  }
  
  return(tau12)
}

# Funzione comune per il calcolo delle armoniche
calculate_harmonics <- function(tau12t2, Fs, Amp) {
  L <- length(tau12t2)
  Y <- fft(tau12t2)
  
  P2 <- abs(Y/L)
  P1 <- P2[1:(L/2+1)]
  P1[2:(length(P1)-1)] <- 2 * P1[2:(length(P1)-1)]
  
  f <- 2*pi*Fs*(0:(L/2))/L
  harmInd <- seq(1, 61, by=2)
  
  harmonics <- approx(f, P1, xout=harmInd)$y / Amp
  
  P2r <- Re(Y/L)
  P1r <- P2r[1:(L/2+1)]
  P1r[2:(length(P1r)-1)] <- 2 * P1r[2:(length(P1r)-1)]
  harmonicsR <- approx(f, P1r, xout=harmInd)$y / Amp
  
  P2i <- Im(Y/L)
  P1i <- P2i[1:(L/2+1)]
  P1i[2:(length(P1i)-1)] <- 2 * P1i[2:(length(P1i)-1)]
  harmonicsI <- -approx(f, P1i, xout=harmInd)$y / Amp
  
  return(list(harmonicsR = harmonicsR, harmonicsI = harmonicsI))
}

# Funzione base KBKZ
KBKZ_base <- function(para, PointsPerPeriod, NoPeriods, CutPeriods, etas, gi, taui, PreSum) {
  omega <- para[1]
  Amp <- para[2]
  
  mint <- 0
  maxt <- 2*pi/omega*NoPeriods
  tstep <- 2*pi/omega/PointsPerPeriod
  time <- seq(mint, maxt, by=tstep)
  
  Integrand <- function(tp, t, PreSumn, gi, taui) {
    PreSumn * sum(gi/taui * exp(-(t-tp)/taui))
  }
  
  tau12 <- calculate_stress(time, PreSum, Integrand, gi, taui, tstep)
  
  tau12s <- etas*Amp*omega*cos(omega*time)
  tau12t <- tau12 + tau12s
  
  time2 <- time[(PointsPerPeriod*CutPeriods+1):length(time)]
  tau12t2 <- tau12t[(PointsPerPeriod*CutPeriods+1):length(tau12t)]
  
  Fs <- length(time2) / (2*pi*(NoPeriods-CutPeriods))
  
  harmonics <- calculate_harmonics(tau12t2, Fs, Amp)
  
  strain <- Amp * sin(omega * time2)
  strainRt <- Amp * omega * cos(omega * time2)
  
  return(list(time2 = time2, strain = strain, strainRt = strainRt, 
              tau12t2 = tau12t2, harmonicsR = harmonics$harmonicsR, 
              harmonicsI = harmonics$harmonicsI))
}

# Modelli specifici
KBKZ_Laun <- function(para, PointsPerPeriod, NoPeriods, CutPeriods, etas) {
  result <- taugcal(0.256, 9963/7537*4775.94, 3.162, 0.001, 15)
  taui <- result$taui
  gi <- result$gi
  
  f <- para[3]
  n1 <- para[4]
  n2 <- para[5]
  
  PreSum <- function(tp, t) {
    Amp <- para[2]
    omega <- para[1]
    Amp * (sin(omega*t) - sin(omega*tp)) * 
      (f*exp(-n1*abs(Amp*(sin(omega*t) - sin(omega*tp)))) + 
       (1-f)*exp(-n2*abs(Amp*(sin(omega*t) - sin(omega*tp)))))
  }
  
  KBKZ_base(para, PointsPerPeriod, NoPeriods, CutPeriods, etas, gi, taui, PreSum)
}

KBKZ_Lissajous <- function(para, PointsPerPeriod, NoPeriods, CutPeriods, etas) {
  result <- taugcal(0.256, 9963/7537*4775.94, 3.162, 0.001, 15)
  taui <- result$taui
  gi <- result$gi
  
  a <- para[3]
  
  PreSum <- function(tp, t) {
    Amp <- para[2]
    omega <- para[1]
    Amp * (sin(omega*t) - sin(omega*tp)) / (1 + a*Amp^2 * (sin(omega*t) - sin(omega*tp))^2)
  }
  
  KBKZ_base(para, PointsPerPeriod, NoPeriods, CutPeriods, etas, gi, taui, PreSum)
}

KBKZ_SoskeyLissajous <- function(para, PointsPerPeriod, NoPeriods, CutPeriods, etas) {
  gi <- c(177464, 54980, 35272.4, 22493.7, 16039.9, 13572.9, 9004.74, 
          6926.17, 527.529, 2.7341)
  taui <- c(0.00020788, 0.00288704, 0.0234895, 0.163525, 1.03726, 
            5.24526, 15.6349, 35.5925, 167.211, 4810.48)
  
  aSos <- para[3]
  bSos <- para[4]
  
  PreSum <- function(tp, t) {
    Amp <- para[2]
    omega <- para[1]
    Amp * (sin(omega*t) - sin(omega*tp)) * 
      (1 / (1 + aSos * (abs(Amp * (sin(omega*t) - sin(omega*tp))))^bSos))
  }
  
  KBKZ_base(para, PointsPerPeriod, NoPeriods, CutPeriods, etas, gi, taui, PreSum)
}

# Funzione principale
run_KBKZ_model <- function(model_type, para) {
  PointsPerPeriod <- 125
  NoPeriods <- 200
  CutPeriods <- 100
  etas <- 0.01

  if (model_type == "Laun") {
    result <- KBKZ_Laun(para, PointsPerPeriod, NoPeriods, CutPeriods, etas)
  } else if (model_type == "Lissajous") {
    result <- KBKZ_Lissajous(para, PointsPerPeriod, NoPeriods, CutPeriods, etas)
  } else if (model_type == "SoskeyLissajous") {
    result <- KBKZ_SoskeyLissajous(para, PointsPerPeriod, NoPeriods, CutPeriods, etas)
  } else {
    stop("Invalid model type")
  }
  
  return(result)
}

# Funzioni di supporto
SimpsonOdd <- function(integrand, h) {
  (integrand[1] + 4*sum(integrand[seq(2, length(integrand)-1, by=2)]) + 
   2*sum(integrand[seq(3, length(integrand)-2, by=2)]) + 
   integrand[length(integrand)]) * h/3
}

SimpsonThreeEight <- function(integrand, h) {
  3*h/8 * (integrand[1] + 3*integrand[2] + 3*integrand[3] + integrand[4])
}

# Esempi di utilizzo
para_laun <- c(1, 0.1, 0.6658, 10.35, 278.9)
result_laun <- run_KBKZ_model("Laun", para_laun)
print(result_laun)

para_lissajous <- c(1, 0.1, 0.5)
result_lissajous <- run_KBKZ_model("Lissajous", para_lissajous)
print(result_lissajous)

para_soskey <- c(1, 0.1, 0.5, 0.7)
result_soskey <- run_KBKZ_model("SoskeyLissajous", para_soskey)
print(result_soskey)