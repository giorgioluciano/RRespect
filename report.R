library(ggplot2)
library(dplyr)

# ===== FUNZIONE GENERICA =====
plot_gfreq <- function(test_num) {
  r_path  <- paste0("C:/test/Gstar/output_R_freq/test", test_num, "/output/Gfit.dat")
  py_path <- paste0("C:/test/pyReSpect-freq-master/output_PY_freq/test", test_num, "/output/Gfit.dat")
  
  gf_r  <- read.table(r_path,  header=FALSE, col.names=c("omega", "Gprime_R",  "Gpp_R"))
  gf_py <- read.table(py_path, header=FALSE, col.names=c("omega", "Gprime_PY", "Gpp_PY"))
  
  gf <- full_join(gf_r, gf_py, by="omega")
  
  p <- ggplot(gf, aes(x=omega)) +
    geom_line(aes(y=Gprime_PY, color="Python ref"), linewidth=1.2, alpha=0.9) +
    geom_line(aes(y=Gprime_R,  color="R"),          linewidth=1.0, linetype="22", alpha=0.9) +
    scale_x_log10(labels=scales::scientific) +
    scale_y_log10(labels=scales::scientific) +
    scale_color_manual(values=c("Python ref"="#000000", "R"="#01696f")) +
    labs(title=paste0("Frequency test", test_num, ": G'(\u03c9) reconstruction"),
         x=expression(omega~"[rad/s]"), y=expression(G*"'("*omega*")")) +
    theme_minimal(base_size=12) +
    theme(legend.position="bottom", legend.title=element_blank())
  
  outfile <- paste0("gfreq_test", test_num, ".png")
  ggsave(outfile, p, width=9, height=6, dpi=300, bg="white")
  cat("✅ Salvato:", outfile, "\n")
}

# ===== FUNZIONE TIME (dash sottile) =====
plot_gt_time <- function(test_num) {
  r_path  <- paste0("C:/test/Gt/output_R_time/test", test_num, "/output/Gfit.dat")
  py_path <- paste0("C:/test/pyReSpect-time-master/output_PY_time/test", test_num, "/output/Gfit.dat")
  
  gt_r  <- read.table(r_path,  header=FALSE, col.names=c("time", "G_R"))
  gt_py <- read.table(py_path, header=FALSE, col.names=c("time", "G_PY"))
  
  gt <- full_join(gt_r, gt_py, by="time")
  
  p <- ggplot(gt, aes(x=time)) +
    geom_line(aes(y=G_PY, color="Python ref"), linewidth=1.2, alpha=0.9) +
    geom_line(aes(y=G_R,  color="R"),          linewidth=1.0, linetype="22", alpha=0.9) +
    scale_x_log10(labels=scales::scientific) +
    scale_color_manual(values=c("Python ref"="#000000", "R"="#01696f")) +
    labs(title=paste0("Time test", test_num, ": G(t) reconstruction"),
         x="time t [s]", y="G(t)") +
    theme_minimal(base_size=12) +
    theme(legend.position="bottom", legend.title=element_blank())
  
  outfile <- paste0("gt_time_test", test_num, ".png")
  ggsave(outfile, p, width=9, height=6, dpi=300, bg="white")
  cat("✅ Salvato:", outfile, "\n")
}

# ===== ESEGUI TUTTO =====
cat("🎨 Generando 4 plot...\n")
plot_gt_time("1")
plot_gt_time("7")
plot_gfreq("1")
plot_gfreq("7")

cat("\n✅ PRONTI PER OVERLEAF:\n")
cat("   gt_time_test1.png\n")
cat("   gt_time_test7.png\n")
cat("   gfreq_test1.png\n")
cat("   gfreq_test7.png\n")
