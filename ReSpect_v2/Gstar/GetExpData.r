# Function: GetExpData
#
# Reads in the experimental data from the input file
#
# Input:  fname = name of file that contains G*(w) in 3 columns [w Gp Gpp]
#         Space it evenly on a log scale
#
# Output: A n*1 vector "w", and a 2n*1 vector Gexp = [Gp; Gpp]

GetExpData <- function(fname) {
  # Read the data from the file
  data <- read.table(fname, header = FALSE)
  
  # Separate the data into three columns: w, Gp, Gpp
  wo <- data[, 1]
  Gpo <- data[, 2]
  Gppo <- data[, 3]
  
  # Remove any repeated frequency values
  unique_data <- unique(data)
  wo <- unique_data[, 1]
  Gpo <- unique_data[, 2]
  Gppo <- unique_data[, 3]
  
  # Space it evenly on a log scale
  w <-  10^seq(log10(min(wo)), log10(max(wo)), length.out = length(wo))
  Gp <- approx(wo, Gpo, xout = w, method = "linear", rule = 2)$y
  Gpp <- approx(wo, Gppo, xout = w, method = "linear", rule = 2)$y
  
  # Use supersmoother to clean it up further. This may be optional
  Gp <- supsmu(w, Gp)$y
  Gpp <- supsmu(w, Gpp)$y
  
  # Combine Gp and Gpp into a single vector
  Gexp <- c(Gp, Gpp)
  
  return(list(w = w, Gexp = Gexp))
}


# Example usage
#fname <- "Gst.dat"  # replace with your actual file path
#result <- GetExpData(fname)
#w <- result$w
#Gexp <- result$Gexp

#print(w)
#print(Gexp)
