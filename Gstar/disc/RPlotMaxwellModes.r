# Load necessary library for plotting
library(ggplot2)

# Function to plot experimental dynamic moduli and dynamic moduli obtained from the CRS
# Input: g, t = spectrum
#        w = n*1 vector contains frequencies
#        Gp, Gpp = experimental data

PlotMaxwellModes <- function(g, t, w, Gp, Gpp) {
  N <- length(g)
  
  ws <- outer(t, w, "*")
  ws2 <- ws^2
  
  GpM <- ws2 / (1 + ws2) %*% g
  GppM <- ws / (1 + ws2) %*% g
  
  # Convert data to a data frame for ggplot2
  data <- data.frame(
    w = rep(w, 2),
    value = c(Gp, Gpp),
    type = rep(c("Gp", "Gpp"), each = length(w)),
    fit = rep(c("Experimental", "Fit"), each = length(w))
  )
  
  fit_data <- data.frame(
    w = rep(w, 2),
    value = c(GpM, GppM),
    type = rep(c("Gp", "Gpp"), each = length(w)),
    fit = rep(c("Experimental", "Fit"), each = length(w))
  )
  
  combined_data <- rbind(data, fit_data)
  
  # Plot using ggplot2
  ggplot(combined_data, aes(x = w, y = value, color = fit)) +
    geom_point(data = subset(combined_data, fit == "Experimental")) +
    geom_line(data = subset(combined_data, fit == "Fit"), size = 1) +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~ type, scales = "free_y") +
    theme_minimal() +
    labs(
      title = "Dynamic Moduli Comparison",
      x = "Frequency (w)",
      y = "Moduli (Gp, Gpp)"
    )
}

# Example usage
g <- runif(10)  # Example spectrum
t <- runif(10)  # Example relaxation times
w <- seq(0.1, 10, length.out = 20)  # Example frequencies
Gp <- sin(w)  # Example Gp data
Gpp <- cos(w)  # Example Gpp data

PlotMaxwellModes(g, t, w, Gp, Gpp)

