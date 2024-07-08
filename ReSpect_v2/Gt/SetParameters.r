# Define the function SetParameters
SetParameters <- function() {
  par <- list()
  
  # ===================================== #
  # O V E R A L L   P A R A M E T E R S #
  # ===================================== #
  # Printing to screen and files is ON (=1) or OFF (0)
  par$verbose <- 1
  
  # Plotting functions are ON (=1) or OFF (0)
  par$plotting <- 1
  
  # ======================================= #
  # C O N T I N U O U S   S P E C T R U M #
  # ======================================= #
  # Number of grid points to represent the continuous spectrum
  par$ns <- 100
  
  # Specify lambda_C instead of using the one inferred from the L-curve
  # If 0, then use the L-curve method to determine lambda.
  par$lamC <- 0
  
  # Smoothing Factor: Indirect way of controlling lambda_C, relative to
  # the one inferred from L-curve.
  # Set between -1 (lowest lambda explored) and 1 (highest lambda explored);
  # When set to 0, using lambda_C determined from the L-curve
  par$SmFacLam <- 0
  
  # Default filename for Gp and Gpp data
  # It should contain G*(w) in 3 columns [w Gp Gpp]
  par$GstFile <- 'Gst.dat'
  
  # Treatment of frequency window ends:
  #  = 3 : t = 1/w - strict condition
  #  = 2 : t = 1/w
  #  = 1 : t = 1/w + lenient condition
  par$FreqEnd <- 1
  
  # ==================================== #
  # D I S C R E T E    S P E C T R U M #
  # ==================================== #
  # Specify number of discrete modes (>0). If set equal to zero, then 
  # determine automatically, using internal algorithm
  par$Nopt <- 0
  
  # Avoid modes with -ve weights (=1), or don't care (0)
  par$prune <- 1
  
  # Parameter controls spacing of tau_i 
  # = 0 => completely determined by H(tau)
  # = 1 => tau_i are equispaced
  # Set between 0 and 1. 
  par$BaseDistWt <- 0.5
  
  # Parameter controls relative importance of error and conditioning in
  # determining Nopt (operative only if par.Nopt = 0, else it is irrelevant)
  # = 0 => completely determined by error
  # = 1 => completely determined by condition number
  # Set between 0 and 1
  par$condWt <- 0.5
  
  # Check for parameter consistency
  CheckConsistency(par)
  
  return(par)
}

# Define the function CheckConsistency
CheckConsistency <- function(par) {
  if (par$BaseDistWt < 0 || par$BaseDistWt > 1) {
    stop('par$BaseDistWt needs to be between 0 and 1')
  }
  
  if (par$condWt < 0 || par$condWt > 1) {
    stop('par$condWt needs to be between 0 and 1')
  }
}

# Example usage
#parameters <- SetParameters()
#print(parameters)
