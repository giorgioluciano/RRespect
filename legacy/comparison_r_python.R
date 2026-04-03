# ================================================
# ReSpectR - COMPARISON SCRIPT COMPLETO
# Discrete + Continuous, Time + Frequency
# COPIA → SALVA → source("compare_all.R")
# ================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)

cat("🚀 ReSpectR - Comparison R vs Python\n\n")

# ===== PATH BASE (modifica qui se necessario) =====
base_R_time   <- "C:/test/Gt/output_R_time"
base_PY_time  <- "C:/test/pyReSpect-time-master/output_PY_time"
base_R_freq   <- "C:/test/Gstar/output_R_freq"
base_PY_freq  <- "C:/test/pyReSpect-freq-master/output_PY_freq"

tests_time <- paste0("test", 1:7)
tests_freq <- c("test1","test1n","test1u","test2","test3","test4","test5","test6","test7")

# ================================================
# FUNZIONI LETTURA
# ================================================

read_dmodes <- function(path) {
  f <- file.path(path, "output", "dmodes.dat")
  if (!file.exists(f)) return(NULL)
  d <- read.table(f, header=FALSE)
  colnames(d) <- c("tau", "G_mode")[1:ncol(d)]
  d$Nopt <- nrow(d)
  return(d)
}

read_gfit_time <- function(path) {
  f <- file.path(path, "output", "Gfit.dat")
  if (!file.exists(f)) return(NULL)
  d <- read.table(f, header=FALSE, col.names=c("time", "G"))
  return(d)
}

read_gfit_freq <- function(path) {
  f <- file.path(path, "output", "Gfit.dat")
  if (!file.exists(f)) return(NULL)
  d <- read.table(f, header=FALSE, col.names=c("omega", "Gprime", "Gpp"))
  return(d)
}

read_H <- function(path) {
  f <- file.path(path, "output", "H.dat")
  if (!file.exists(f)) return(NULL)
  d <- read.table(f, header=FALSE, col.names=c("tau", "H"))
  return(d)
}

# ================================================
# 1. DISCRETE TIME: Nopt + g_rmse
# ================================================
cat("── 1. Discrete TIME comparison ──\n")

disc_time <- data.frame()
for (tt in tests_time) {
  path_r  <- file.path(base_R_time,  tt)
  path_py <- file.path(base_PY_time, tt)
  
  dm_r  <- read_dmodes(path_r)
  dm_py <- read_dmodes(path_py)
  
  if (is.null(dm_r) | is.null(dm_py)) {
    cat(sprintf("  SKIP %s: file mancante\n", tt)); next
  }
  
  nopt_r  <- nrow(dm_r)
  nopt_py <- nrow(dm_py)
  
  # g_rmse solo se stessa dimensione
  if (nopt_r == nopt_py) {
    g_rmse <- sqrt(mean((dm_r$G_mode - dm_py$G_mode)^2))
    tau_rel <- mean(abs(dm_r$tau - dm_py$tau) / abs(dm_py$tau + 1e-20))
  } else {
    g_rmse <- NA; tau_rel <- NA
  }
  
  disc_time <- rbind(disc_time, data.frame(
    test=tt, Nopt_R=nopt_r, Nopt_PY=nopt_py,
    delta=nopt_r-nopt_py, g_rmse=g_rmse, tau_reldiff=tau_rel
  ))
}
print(disc_time)
write.csv(disc_time, "compare_disc_time.csv", row.names=FALSE)

# ================================================
# 2. DISCRETE FREQ: Nopt
# ================================================
cat("\n── 2. Discrete FREQ comparison ──\n")

disc_freq <- data.frame()
for (tt in tests_freq) {
  path_r  <- file.path(base_R_freq,  tt)
  path_py <- file.path(base_PY_freq, tt)
  
  dm_r  <- read_dmodes(path_r)
  dm_py <- read_dmodes(path_py)
  
  if (is.null(dm_r) | is.null(dm_py)) {
    cat(sprintf("  SKIP %s: file mancante\n", tt)); next
  }
  
  nopt_r  <- nrow(dm_r)
  nopt_py <- nrow(dm_py)
  
  if (nopt_r == nopt_py) {
    g_rmse <- sqrt(mean((dm_r$G_mode - dm_py$G_mode)^2))
  } else {
    g_rmse <- NA
  }
  
  disc_freq <- rbind(disc_freq, data.frame(
    test=tt, Nopt_R=nopt_r, Nopt_PY=nopt_py,
    delta=nopt_r-nopt_py, g_rmse=g_rmse
  ))
}
print(disc_freq)
write.csv(disc_freq, "compare_disc_freq.csv", row.names=FALSE)

# ================================================
# 3. CONTINUOUS TIME: H(tau) RMSE
# ================================================
cat("\n── 3. Continuous TIME: H(tau) comparison ──\n")

cont_time <- data.frame()
for (tt in tests_time) {
  path_r  <- file.path(base_R_time,  tt)
  path_py <- file.path(base_PY_time, tt)
  
  h_r  <- read_H(path_r)
  h_py <- read_H(path_py)
  
  if (is.null(h_r) | is.null(h_py)) {
    cat(sprintf("  SKIP %s: H.dat mancante\n", tt)); next
  }
  
  # Interpola su griglia comune
  tau_common <- h_py$tau
  h_r_interp <- approx(h_r$tau, h_r$H, xout=tau_common, rule=2)$y
  
  H_rmse    <- sqrt(mean((h_r_interp - h_py$H)^2, na.rm=TRUE))
  H_maxdiff <- max(abs(h_r_interp - h_py$H), na.rm=TRUE)
  H_reldiff <- mean(abs(h_r_interp - h_py$H) / (abs(h_py$H) + 1e-20), na.rm=TRUE)
  
  cont_time <- rbind(cont_time, data.frame(
    test=tt, H_rmse=H_rmse, H_maxdiff=H_maxdiff, H_reldiff=H_reldiff
  ))
}
print(cont_time)
write.csv(cont_time, "compare_cont_time.csv", row.names=FALSE)

# ================================================
# 4. CONTINUOUS FREQ: H(tau) RMSE
# ================================================
cat("\n── 4. Continuous FREQ: H(tau) comparison ──\n")

cont_freq <- data.frame()
for (tt in tests_freq) {
  path_r  <- file.path(base_R_freq,  tt)
  path_py <- file.path(base_PY_freq, tt)
  
  h_r  <- read_H(path_r)
  h_py <- read_H(path_py)
  
  if (is.null(h_r) | is.null(h_py)) {
    cat(sprintf("  SKIP %s: H.dat mancante\n", tt)); next
  }
  
  tau_common <- h_py$tau
  h_r_interp <- approx(h_r$tau, h_r$H, xout=tau_common, rule=2)$y
  
  H_rmse    <- sqrt(mean((h_r_interp - h_py$H)^2, na.rm=TRUE))
  H_maxdiff <- max(abs(h_r_interp - h_py$H), na.rm=TRUE)
  H_reldiff <- mean(abs(h_r_interp - h_py$H) / (abs(h_py$H) + 1e-20), na.rm=TRUE)
  
  cont_freq <- rbind(cont_freq, data.frame(
    test=tt, H_rmse=H_rmse, H_maxdiff=H_maxdiff, H_reldiff=H_reldiff
  ))
}
print(cont_freq)
write.csv(cont_freq, "compare_cont_freq.csv", row.names=FALSE)

# ================================================
# 5. PLOT Gfit TIME: test1 e test7
# ================================================
cat("\n── 5. Plot G(t): test1 e test7 ──\n")

for (tt in c("1","7")) {
  r_path  <- file.path(base_R_time,  paste0("test",tt))
  py_path <- file.path(base_PY_time, paste0("test",tt))
  
  gt_r  <- read_gfit_time(r_path)
  gt_py <- read_gfit_time(py_path)
  
  if (is.null(gt_r)|is.null(gt_py)) { cat(sprintf("  SKIP test%s\n",tt)); next }
  
  gt <- full_join(gt_r, gt_py, by="time", suffix=c("_R","_PY"))
  
  p <- ggplot(gt, aes(x=time)) +
    geom_line(aes(y=G_PY, color="Python ref"), linewidth=1.2) +
    geom_line(aes(y=G_R,  color="R"),          linewidth=1.0, linetype="22") +
    scale_x_log10(labels=scales::scientific) +
    scale_color_manual(values=c("Python ref"="#000000","R"="#01696f")) +
    labs(title=paste0("test",tt,": G(t) reconstruction"),
         x="time t [s]", y="G(t)") +
    theme_minimal(base_size=12) +
    theme(legend.position="bottom", legend.title=element_blank())
  
  ggsave(paste0("gt_time_test",tt,".png"), p, width=9, height=6, dpi=300, bg="white")
  cat(sprintf("  ✅ gt_time_test%s.png\n", tt))
}

# ================================================
# 6. PLOT Gfit FREQ: test1 e test7
# ================================================
cat("\n── 6. Plot G'(ω): test1 e test7 ──\n")

for (tt in c("1","7")) {
  r_path  <- file.path(base_R_freq,  paste0("test",tt))
  py_path <- file.path(base_PY_freq, paste0("test",tt))
  
  gf_r  <- read_gfit_freq(r_path)
  gf_py <- read_gfit_freq(py_path)
  
  if (is.null(gf_r)|is.null(gf_py)) { cat(sprintf("  SKIP test%s\n",tt)); next }
  
  gf <- full_join(gf_r, gf_py, by="omega", suffix=c("_R","_PY"))
  
  p <- ggplot(gf, aes(x=omega)) +
    geom_line(aes(y=Gprime_PY, color="Python ref"), linewidth=1.2) +
    geom_line(aes(y=Gprime_R,  color="R"),          linewidth=1.0, linetype="22") +
    scale_x_log10(labels=scales::scientific) +
    scale_y_log10(labels=scales::scientific) +
    scale_color_manual(values=c("Python ref"="#000000","R"="#01696f")) +
    labs(title=paste0("test",tt,": G'(\u03c9) reconstruction"),
         x=expression(omega~"[rad/s]"), y=expression(G*"'("*omega*")")) +
    theme_minimal(base_size=12) +
    theme(legend.position="bottom", legend.title=element_blank())
  
  ggsave(paste0("gfreq_test",tt,".png"), p, width=9, height=6, dpi=300, bg="white")
  cat(sprintf("  ✅ gfreq_test%s.png\n", tt))
}

# ================================================
# SOMMARIO FINALE
# ================================================
cat("\n🎉 CONFRONTO COMPLETO:\n")
cat(sprintf("Disc TIME: %d/%d perfetti\n", sum(disc_time$delta==0), nrow(disc_time)))
cat(sprintf("Disc FREQ: %d/%d perfetti\n", sum(disc_freq$delta==0), nrow(disc_freq)))
cat(sprintf("Cont TIME: H_rmse medio = %.2e\n", mean(cont_time$H_rmse, na.rm=TRUE)))
cat(sprintf("Cont FREQ: H_rmse medio = %.2e\n", mean(cont_freq$H_rmse, na.rm=TRUE)))

cat("\n✅ FILE SALVATI:\n")
cat("  compare_disc_time.csv\n")
cat("  compare_disc_freq.csv\n")
cat("  compare_cont_time.csv\n")
cat("  compare_cont_freq.csv\n")
cat("  gt_time_test1.png / gt_time_test7.png\n")
cat("  gfreq_test1.png   / gfreq_test7.png\n")