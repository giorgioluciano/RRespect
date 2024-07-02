# R already has implemented Friedman's SuperSmoother
# Description
# Smooth the (x, y) values by Friedman's ‘super smoother’.

# Usage
# supsmu(x, y, wt =, span = "cv", periodic = FALSE, bass = 0, trace = FALSE)
# Arguments
# x	
# x values for smoothing

# y	
# y values for smoothing

# wt	
# case weights, by default all equal

# span	
# the fraction of the observations in the span of the running lines smoother, or "cv" to choose this by leave-one-out cross-validation.

# periodic	
# if TRUE, the x values are assumed to be in [0, 1] and of period 1.

# bass	
# controls the smoothness of the fitted curve. Values of up to 10 indicate increasing smoothness.

# trace	
# logical, if true, prints one line of info “per spar”, notably useful for "cv".

# Details
# supsmu is a running lines smoother which chooses between three spans for the lines. The running lines smoothers are symmetric, with k/2 data points each side of the predicted point, and values of k as 0.5 * n, 0.2 * n and 0.05 * n, where n is the number of data points. If span is specified, a single smoother with span span * n is used.

# The best of the three smoothers is chosen by cross-validation for each prediction. The best spans are then smoothed by a running lines smoother and the final prediction chosen by linear interpolation.

# The FORTRAN code says: “For small samples (n < 40) or if there are substantial serial correlations between observations close in x-value, then a pre-specified fixed span smoother (span > 0) should be used. Reasonable span values are 0.2 to 0.4.”

# Cases with non-finite values of x, y or wt are dropped, with a warning.

# Value
# A list with components

# x	
# the input values in increasing order with duplicates removed.

# y	
# the corresponding y values on the fitted curve.

# References
# Friedman, J. H. (1984) SMART User's Guide. Laboratory for Computational Statistics, Stanford University Technical Report No. 1.
# Friedman, J. H. (1984) A variable span scatterplot smoother. Laboratory for Computational Statistics, Stanford University Technical Report No. 5.