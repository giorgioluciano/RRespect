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
  
  # Separate the data into three columns: w, Gt
  to <- data[, 1]
  Gto <- data[, 2]
  
  
  # Remove any repeated frequency values
  to_unique <- unique(to)
  
  i <- match(to_unique, to)
  
  Gto <- Gto[i]
  to <- to_unique
  
  # Space it evenly on a log scale
  
  t <-  10^seq(log10(min(to)), log10(max(to)), length.out = length(to))
  Gt <- approx(to, Gto, xout = t, method = "linear", rule = 2)$y
  
  return(list(t = t, Gt=Gt))
}


# Example usage
#fname <- "Gst.dat"  # replace with your actual file path
#result <- GetExpData(fname)
#w <- result$w
#Gexp <- result$Gexp

#print(w)
#print(Gexp)
