# Function: SetParameters
# Set parameters for calculation
# Returns a list 'par' containing various parameters

SetParameters <- function() {
  # Initialize parameters
  par <- list(
    verbose = 1,         # Printing to screen and files is ON (=1) or OFF (0)
    plotting = 1,        # Plotting functions are ON (=1) or OFF (0)
    ns = 100,            # Number of grid points to represent the continuous spectrum
    lamC = 0,            # Specify lambda_C instead of using the one inferred from the L-curve
    SmFacLam = 0,        # Smoothing Factor for controlling lambda_C
    GstFile = 'Gt.dat',  # Default filename for Gp and Gpp data
    FreqEnd = 1,         # Treatment of frequency window ends
    Nopt = 0,            # Specify number of discrete modes (>0), 0 to determine automatically
    prune = 1,           # Avoid modes with negative weights (=1) or don't care (0)
    BaseDistWt = 0.5,    # Parameter controls spacing of tau_i
    condWt = 0.5         # Parameter controls relative importance of error and conditioning
  )
  
  # Check consistency of parameters
  CheckConsistency(par)
  
  return(par)
}

# Function: CheckConsistency
# Checks consistency of parameters in 'par'
# Throws an error if parameters are out of valid range

CheckConsistency <- function(par) {
  if (par$BaseDistWt < 0 || par$BaseDistWt > 1) {
    stop('par$BaseDistWt needs to be between 0 and 1')
  }
  
  if (par$condWt < 0 || par$condWt > 1) {
    stop('par$condWt needs to be between 0 and 1')
  }
}

# Example usage:
parameters <- SetParameters()
print(parameters)
