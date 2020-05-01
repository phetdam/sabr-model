# contains functions for calibrating the CEV forward model to market implied vol
# data by least squares fitting of the market volatility smiles. most of its
# functions correspond to those in sabr.R, but are for the CEV forward instead.
#
# Changelog:
#
# 08-21-2019
#
# changed calibration method from simulated annealing to nelder-mead.
#
# 08-20-2019
#
# added cevf_ivol_obj and cevf_ivol_apobj, as well as fit_cevf. almost done.
#
# 08-19-2019
#
# initial creation. added function names and cevf_ivol.

# function names
fit_cevf__n <- "fit_cevf"
plot_cevf__n <- "plot_cevf"
cevf_ivol_obj__n <- "cevf_ivol_obj"
cevf_ivol_apobj__n <- "cevf_ivol_apobj"

# magnitude of the barrier parameter used in objective functions
cevf_lambda <- 10 ^ (-10)

# cev forward implied volatility approximation from hagan and woodward's paper.
#
# parameters:
#
# f        forward price
# K        strike price [vector] (same units as F)
# Tx       expiration time in years
# sigma    constant scaling the diffusion coefficient of the forward
# beta     nonnegative constant skewing the diffusion of the forward
cevf_ivol <- function(f, K, Tx, sigma, beta) {
  if (sigma <= 0 || beta < 0) { return(NULL) }
  fav <- (f + K) / 2
  cbeta <- 1 - beta
  return(sigma / fav ^ cbeta *
           (1 + cbeta * (2 + beta) / 24 * ((f - K) / fav) ^ 2 + (cbeta ^ 2) / 24
            * sigma ^ 2 / fav ^ (2 * cbeta) * Tx))
}

# objective function for tuning cev forward parameters to black ivols. returns
# mean sum of squared errors as optimization criteria, and uses interior point
# barrier functions to constrain the search space. allows one to weight data
# points non-uniformly. for calibrating boith sigma and beta.
#
# parameters:
#
# X       vector of parameters to tune. order:
#             [1] sigma
#             [2] beta
# f       forward price
# K       strike price [vector] (same units as F)
# ivol    implied volatilites [vector] derived from market prices
# Tx      expiration time in years
# wts     weights for the ivol observations, default NULL. if NULL, then the
#         resulting weight vector will simply be a vector of 1s. all weights
#         should be nonnegative.
cevf_ivol_obj <- function(X, f, K, ivol, Tx, wts = NULL) {
  # the following code is essentially the same as that in sabr.R for the
  # sabr_ivol_obj function.
  # length of K and ivol
  NK <- length(K)
  Niv <- length(ivol)
  # if length of K, ivol don't match, print error and exit
  if (NK != Niv) {
    cat(cevf_ivol_obj__n, ": error: strike and implied vol vectors must have ", 
        "equal length\n", sep = "")
    return(NULL)
  }
  # get vector length
  N <- NK
  # if wts is NULL, set it to 1
  if (is.null(wts) == TRUE) {
    wts <- 1
  }
  # else if wts is not NULL and has a different length, print error and exit
  else if (length(wts) != N) {
    cat(cevf_ivol_obj__n, ": error: wts vector of incorrect length\n", sep = "")
    return(NULL)
  }
  return(sum(wts * (ivol - cevf_ivol(f, K, Tx, X[1], X[2])) ^ 2) / N
         + cevf_lambda * (-log(X[1]) - log(X[2])))
}

# objective function for tuning cev forward parameters to black ivols. returns
# mean sum of squared errors as optimization criteria, and uses interior point
# barrier functions to constrain the search space. allows one to weight data
# points non-uniformly. for calibrating only sigma, with a priori beta.
#
# parameters:
#
# X       vector of parameters to tune. only contains sigma.
#             [1] sigma
# beta    a priori choice of skew. for illustrative purposes usually.
# f       forward price
# K       strike price [vector] (same units as F)
# ivol    implied volatilites [vector] derived from market prices
# Tx      expiration time in years
# wts     weights for the ivol observations, default NULL. if NULL, then the
#         resulting weight vector will simply be a vector of 1s. all weights
#         should be nonnegative.
cevf_ivol_apobj <- function(X, beta, f, K, ivol, Tx, wts) {
  # the following code is essentially the same as in the above function
  # length of K and ivol
  NK <- length(K)
  Niv <- length(ivol)
  # if length of K, ivol don't match, print error and exit
  if (NK != Niv) {
    cat(cevf_ivol_apobj__n, ": error: strike and implied vol vectors must have ", 
        "equal length\n", sep = "")
    return(NULL)
  }
  # get vector length
  N <- NK
  # if wts is NULL, set it to 1
  if (is.null(wts) == TRUE) {
    wts <- 1
  }
  # else if wts is not NULL and has a different length, print error and exit
  else if (length(wts) != N) {
    cat(cevf_ivol_apobj__n, ": error: wts vector of incorrect length\n", sep = "")
    return(NULL)
  }
  return(sum(wts * (ivol - cevf_ivol(f, K, Tx, X[1], beta)) ^ 2) / N - 
           cevf_lambda * log(X[1])
  )
}

# using the optim function, calibrates the cev forward implied vol formula using
# nelder-mead. the initial guess for sigma will always be the implied volatility
# of the strike closest to f. returns the output of optim(); i.e. par (estimates),
# value (error), counts (iterations), convergence (0 for convergence, anything
# else is error), and message. by default, will choose strikes that are 20% 
# above/below the strike (i.e. 0.8 <= K / atmf <= 1.2), unless the parameter
# k_range is set to something else. use npoints to specify sampling. will use
# cevf_ivol_obj as objective function by default, unless beta is specified with
# an a priori value, in which case cevf_ivol_apobj is used instead.
#
# note: if beta is specified a priori, the brent method is used to solve for
# sigma, which has an upper bound of 10 ^ 16.
#
# parameters:
#
# M        smile data list; assumed to contain a column for strikes and a
#          column for the implied volatility associated with each strike.
# f        value of the ATM forward for this smile
# Tx       expiration date for this option smile in years
# beta     optional value of beta for calibration, beta >= 0. default NULL, in
#          which case beta will also be calibrated to the data.
# scale    optional scaling parameter NEEDED to make strikes in units of the ATM
#          forward price; default is 1. set to < 1 to scale down, > 1 to scale up.
#          note: example uses e-mini options, which require scale = 10 ^ (-2).
# k_range  default 0.2. indicates the range of strikes relative to f to choose;
#          i.e. for all K, 1 - k_range <= K / f <= 1 + k_range; 0 < k_range <= 1.
# npoints  the number of sampling points to use, default NULL. will be more or less
#          equidistantly selected from the desired strike range. if NULL is
#          passed, then all elements will be used. we recommend 20.
# wts      weights for the ivol observations, default NULL. if NULL, then the
#          resulting weight vector will simply be a vector of 1s. all weights
#          should be nonnegative.
# klab     optional string indicating which column contains strike prices, if
#          the labeling is different from the default, which is "strike".
# ilab     optional string indicating column of ivols to use; default "ivol"
# guess    optional initial guess for beta. default sigma is value closest to
#          ATM ivol, so no need to guess; guess for beta is 0 default.
# tol      optional absolute convergence tolerance; default 10 ^ (-8)
fit_cevf <- function(M, f, Tx, beta = NULL, scale = 1, k_range = 0.2,
                     npoints = NULL, wts = NULL, klab = "strike", ilab = "ivol",
                     guess = 0.5, tol = 10 ^ (-8)) {
  # if k_range <= 0 or > 1, print error and exit
  if (k_range <= 0 || k_range > 1) {
    cat(fit_cevf__n, ": error: k_range must be in interval (0, 1]\n", sep = "")
    return(NULL)
  }
  # if sample makes no sense, print error and exit
  if (is.null(npoints) == FALSE && npoints < 1) {
    cat(fit_cevf__n, ": error: must sample at least one point\n", sep = "")
    return(NULL)
  }
  # if beta is not NULL, check for sanity
  #if (is.null(beta) == FALSE && beta < 0) {
  #  cat(fit_cevf__n, ": error: beta must be nonnegative\n", sep = "")
  #  return(NULL)
  #}
  # first get the vector of strikes and vector of implied vols, in the range
  # specified by k_range, and scale strikes
  kvec <- get(klab, M) * scale
  ivec <- get(ilab, M)
  ivec <- ivec[kvec >= (1 - k_range) * f & kvec <= (1 + k_range) * f]
  kvec <- kvec[kvec >= (1 - k_range) * f & kvec <= (1 + k_range) * f]
  # sigma is the implied volatility corresponding to the strike nearest (but >=)
  # than the ATM forward value (strikes in kvec already scaled)
  sigma <- ivec[kvec == min(kvec[kvec >= f])]
  # if npoints is NULL or exceeds the length of kvec, use all data
  if (is.null(npoints) == TRUE || npoints > length(kvec)) {
    npoints = length(kvec)
  }
  # if wts is NULL, do nothing
  if (is.null(wts) == TRUE) { }
  # else if wts is not NULL and has a different length than npoints, error
  else if (length(wts) != npoints) {
    cat(fit_cevf__n, ": error: wts vector of incorrect length\n", sep = "")
    return(NULL)
  }
  # get mask for only npoints of strike/ivol pairs (from ivol_util.R). note
  # that if npoints == length(kvec), the vector will be all TRUE
  mask <- get_strike_mask(kvec, npoints)
  # mask kvec and ivec
  kvec <- kvec[mask]
  ivec <- ivec[mask]
  # holds results from optim; if beta is a priori then will contain calibrated
  # sigma with a priori beta added to $par, but if beta is not specified then
  # both sigma and beta will be present.
  ores <- NULL
  # if beta is NULL, then use cevf_ivol_obj to calibrate beta
  if (is.null(beta) == TRUE) {
    beta <- guess
    ores <- optim(c(sigma, beta), fn = cevf_ivol_obj, f = f, K= kvec, ivol = ivec,
                  Tx = Tx, wts = wts, method = "Nelder-Mead", control = list(abstol = tol))
  }
  # else calibrate only sigma, and stick beta into ores$par
  else {
    ores <- optim(c(sigma), fn = cevf_ivol_apobj, beta = beta, f = f, K = kvec,
                  ivol = ivec, Tx = Tx, wts = wts, method = "Brent", lower = 0,
                  upper = 10 ^ 16, control = list(abstol = tol))
    ores$par <- c(ores$par, beta)
  }
  # return the list of results from optim
  return(ores)
}

# plots the cevf implied volatility against a set of strikes and implied vols.
#
# parameters:
#
# Ks       vector of strikes
# ivs      vector of implied vols, must be same length as K
# f        value of ATM forward; must be in same units as vector of strikes
# Tx       option expiration date (years)
# pars     calibrated cevf parameters. vector of
#           [1] sigma
#           [2] beta
# k_range  range of strikes to be plotted, default 0.2. same as in fit_sabr, so
#          one must choose 0 <= k_range <= 1.
# scale    optional parameter that adjusts units of strikes to implied vols.
# npoints  optional number of strike/ivol pairs to sample from within the range
#          of strike/ivol pairs determined by k_range, default NULL (do not do
#          any sampling). we recommend 20.
# xlab     optional label for x axis
# ylab     optional label for y axis
# main     optional title label; default just displays parameters and expiration
# cex      optional scaling for title, axis marks, axis labels, ticks; default 1
# lwd      optional line width, default 2
plot_cevf <- function(Ks, ivs, f, Tx, pars, k_range = 0.2, scale = 1,
                      npoints = NULL, xlab = "strike / ATM forward",
                      ylab = "implied vol", main = NULL, cex = 1, lwd = 2) {
  # quick sanity checks
  if (length(Ks) != length(ivs)) {
    cat(plot_cevf__n, ": error: must be equal number of strikes and implied ", 
        "vols\n", sep = "")
    return(NULL)
  }
  if (length(pars) != 2) {
    cat(plot_cevf__n, ": error: parameter vector must be length 2\n", sep = "")
    return(NULL)
  }
  if (is.null(npoints) == FALSE && npoints < 1) {
    cat(plot_cevf__n, ": error: must sample at least one point\n", sep = "")
    return(NULL)
  }
  # select strike range (order is important!)
  ivs <- ivs[Ks >= (1 - k_range) * f & Ks <= (1 + k_range) * f]
  Ks <- Ks[Ks >= (1 - k_range) * f & Ks <= (1 + k_range) * f]
  # if npoints is NULL or greater than the number of points, use all data
  if (is.null(npoints) == TRUE || npoints > length(Ks)) {
    npoints = length(Ks)
  }
  # get mask and mask both Ks and ivs
  mask <- get_strike_mask(Ks, npoints)
  Ks <- Ks[mask]
  ivs <- ivs[mask]
  # if main is NULL, the title will just be the parameters and the maturity date
  if (is.null(main)) {
    # can use greek letters by using the unicode representation
    main = paste("\u03c3 = ", round(pars[1], digits = 3), ", \u03b2 = ",
                 round(pars[2], digits = 3), sep = "")
  }
  # matrix plot; plot implied vols as points, sabr vols as interpolated line
  matplot(Ks * scale / f, cbind(ivs, cevf_ivol(f, Ks * scale, Tx, pars[1],
                                               pars[2])), pch = 18, 
          col = c("steelblue", "plum"), lwd = lwd, lty = 1, type = c("p", "l"),
          xlab = xlab, ylab = ylab, main = main, font.main = 1, cex.main = cex,
          cex.lab = cex, cex.axis = cex, cex = cex)
}