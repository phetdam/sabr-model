# contains functions for calibrating the SABR model to market data by least
# squares fitting of the market volatility curve/surface.
#
# Changelog:
#
# 08-20-2019
#
# updated description for fit_sabr. removed defunct abstol, prev_err and cur_err
# local variables, since we are just doing global optimization.
#
# 08-18-2019
#
# updated sabr_ivol_obj and fit_sabr to allow one to also pass an optional vector
# of weights for the implied vol observations, and added some error checking to
# both the functions.
#
# 08-15-2019
#
# edited fit_sabr to allow the sampling of strikes during fitting and plotting.
#
# 08-07-2019
#
# changed lw to lwd in plot_sabr. lw is not a graphical parameter.
#
# 07-17-2019
#
# changed fit_sabr default ivol arg from "call_ivol" to just "ivol". changed
# the default k_range in the fitting function, and added interior point boundary
# functions to the objective function to enforce strict bounds. got rid of
# redundant section in calibration function; we are already using global. added
# the plot_sabr function to make it easier to compare and graph fits.
#
# 07-14-2019
#
# initial creation. migrated over sabr_ivol, sabr_ivol_obj, fit_sabr, sabr_call,
# and the incomplete sabr_put functions from ivol_util.R

# function names
fit_sabr__n <- "fit_sabr"
sabr_put__n <- "sabr_put"
plot_sabr__n <- "plot_sabr"
sabr_ivol_obj__n <- "sabr_ivol_obj"

# sabr implied volatility approximation. original from hagan's paper.
#
# parameters:
#
# f        forward price
# K        strike price [vector] (same units as F)
# Tx       expiration time in years
# alpha    current level of volatility (implied), shifts smile
#          for initial calibration, use atm market ivol initially [vector]
# beta     controls skew and shifts curve
# rho      controls correlation between F and vol process; control skew
# nu       vol of vol process, controls smile curvature
sabr_ivol <- function(f, K, Tx, alpha, beta, rho, nu) {
  if (alpha <= 0 || f == 0 || K == 0|| beta < 0 || nu < 0) { return(NULL) }
  logfdivK <- log(f / K)
  fmulK <- f * K
  z <- nu / alpha * fmulK ^ (0.5 - 0.5 * beta) * logfdivK
  chi <- log((sqrt(1 - 2 * rho * z  + z ^ 2) + z - rho) / (1 - rho))
  # second big half of sabr ivol formula, which includes Tx term
  ihalf <- (1 + ((1 - beta) ^ 2 / 24 * alpha ^ 2 / fmulK ^ (1 - beta) + 
                   0.25 * rho * beta * nu * alpha /  fmulK ^ (0.5 - 0.5 * beta) + 
                   (2 - 3 * rho ^ 2) * nu ^ 2 / 24) * Tx)
  # if f ~ K, then we will get NaN, so return ATM formula. we need to use
  # ifelse to handle vector conditions, else return standard sabr formula/
  ivol <- ifelse(abs(f - K) < 10 ^ (-4), 
                 (alpha / f ^ (1 - beta) * ihalf), 
                 alpha * z / (chi * fmulK ^ (0.5 - 0.5 * beta) * 
                                (1 + (1 - beta) ^ 2 / 24 * logfdivK ^ 2 + 
                                   (1 - beta) ^ 4 / 1920 * logfdivK ^ 4)) * ihalf
  )
  return(ivol)
}

# objective function for tuning sabr parameters to black ivol derived from the
# market prices. returns mean sum of squared errors as optimization criteria, and
# uses interior point barrier functions to limit the optimization space explicitly.
# allows one to weight implied vol observations differently as well.
#
# parameters:
#
# X       vector of parameters to tune. order:
#             [1] alpha
#             [2] rho
#             [3] nu
# beta    predetermined; controls skew and shifts curve
# f       forward price
# K       strike price [vector] (same units as F)
# ivol    implied volatilites [vector] derived from market prices
# Tx      expiration time in years
# wts     weights for the ivol observations, default NULL. if NULL, then the
#         resulting weight vector will simply be a vector of 1s. all weights
#         should be nonnegative.
sabr_ivol_obj <- function(X, beta, f, K, ivol, Tx, wts = NULL) {
  # barrier parameter (just one for simplicity)
  lambda <- 10 ^ (-10)
  # length of K and ivol
  NK <- length(K)
  Niv <- length(ivol)
  # if length of K, ivol don't match, print error and exit
  if (NK != Niv) {
    cat(sabr_ivol_obj__n, ": error: strike and implied vol vectors must have ", 
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
    cat(sabr_ivol_obj__n, ": error: wts vector of incorrect length\n", sep = "")
    return(NULL)
  }
  # a vector we used for the SPX 09-2019 data to forcibly fit the curvature
  #wts <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1000, 1000, 10, 100, 1000, 1000, 2000, 1000)
  return(sum(wts * (ivol - sabr_ivol(f, K, Tx, X[1], beta, X[2], X[3])) ^ 2) / N +
           lambda * (-log(X[1]) - log(X[2] + 1) - log(1 - X[2]) - log(X[3]))
         )
}

# using the optim function, calibrates the sabr implied volatility formula using
# nelder-mead. the initial guess for alpha will always be the implied volatility
# of the strike closest to f. returns the output of optim(); i.e. par (estimates),
# value (error), counts (iterations), convergence (0 for convergence, anything
# else is error), and message. by default, will choose strikes that are 20% 
# above/below the strike (i.e. 0.8 <= K / atmf <= 1.2), unless the parameter
# k_range is set to something else. use npoints to specify sampling.
#
# parameters:
#
# M        smile data list; assumed to contain a column for strikes and a
#          column for the implied volatility associated with each strike.
# f        value of the ATM forward for this smile
# Tx       expiration date for this option smile in years
# beta     selected value of beta for calibration; 0 <= beta <= 1
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
# guess    optional two-element vector of initial guesses for rho and nu. default
#          is 0.5 and 0.5 for both rho and nu. initial guess for alpha is always
#          the implied volatility of the strike closest to f. index:
#          [1] rho   [2] nu
# tol      optional absolute convergence tolerance; default 10 ^ (-8)
fit_sabr <- function(M, f, Tx, beta, scale = 1, k_range = 0.2, npoints = NULL,
                     wts = NULL, klab = "strike", ilab = "ivol",
                     guess = c(0.5, 0.9), tol = 10 ^ (-8)) {
  # guess must be two elements long
  if (length(guess) != 2) {
    cat(fit_sabr__n, ": error: initial guesses for rho, nu only\n", sep = "")
    return(NULL)
  }
  # if k_range <= 0 or > 1, print error and exit
  if (k_range <= 0 || k_range > 1) {
    cat(fit_sabr__n, ": error: k_range must be in interval (0, 1]\n", sep = "")
    return(NULL)
  }
  # if sample makes no sense, print error and exit
  if (is.null(npoints) == FALSE && npoints < 1) {
    cat(fit_sabr__n, ": error: must sample at least one point\n", sep = "")
    return(NULL)
  }
  # prev, current error
  p_err <- Inf
  c_err <- Inf
  # initial guesses for alpha, rho, and nu
  # first get the vector of strikes and vector of implied vols, in the range
  # specified by k_range, and scale strikes
  kvec <- get(klab, M) * scale
  ivec <- get(ilab, M)
  ivec <- ivec[kvec >= (1 - k_range) * f & kvec <= (1 + k_range) * f]
  kvec <- kvec[kvec >= (1 - k_range) * f & kvec <= (1 + k_range) * f]
  # alpha is the implied volatility corresponding to the strike nearest (but >=)
  # than the ATM forward value (strikes in kvec already scaled)
  alpha <- ivec[kvec == min(kvec[kvec >= f])]
  # if npoints is NULL or exceeds the length of kvec, use all data
  if (is.null(npoints) == TRUE || npoints > length(kvec)) {
    npoints = length(kvec)
  }
  # if wts is NULL, do nothing
  if (is.null(wts) == TRUE) { }
  # else if wts is not NULL and has a different length than npoints, error
  else if (length(wts) != npoints) {
    cat(fit_sabr__n, ": error: wts vector of incorrect length\n", sep = "")
    return(NULL)
  }
  # get mask for only npoints of strike/ivol pairs (from ivol_util.R). note
  # that if npoints == length(kvec), the vector will be all TRUE
  mask <- get_strike_mask(kvec, npoints)
  # mask kvec and ivec
  kvec <- kvec[mask]
  ivec <- ivec[mask]
  # set guesses for rho and nu
  rho <- guess[1]
  nu <- guess[2]
  # will hold results from optim; global so no need to do anything else
  ores <- optim(c(alpha, rho, nu), fn = sabr_ivol_obj, beta = beta, f = f,
                K = kvec, ivol = ivec, Tx = Tx, wts = wts,
                method = "Nelder-Mead", control = list(abstol = tol))
  # return the list of results from optim
  return(ores)
}

# given a calibrated sabr model, use black_call to return the price of an call.
# note that alpha, rho, nu in the vector pars should have been calibrated, and
# that beta should have been chosen a priori.
#
# parameters:
#
# pars    calibrated sabr parameters alpha, rho, nu
#             [1] alpha
#             [2] rho
#             [3] nu
# beta    a priori choice for skew
#
# see black_call for documentation on remaining parameters
sabr_call <- function(pars, beta, f, K, Tx, dfactor = Dexp, dfargs = .Dexp_r) {
  ivol <- sabr_ivol(f, K, Tx, pars[1], beta, pars[2], pars[3])
  return(black_call(f, K, ivol, Tx, dfactor = dfactor, dfargs = dfargs))
}

# plots the sabr implied volatility against a set of strikes and implied vols.
# good to produce a nice graph to compare the two.
#
# parameters:
#
# Ks       vector of strikes
# ivs      vector of implied vols, must be same length as K
# f        value of ATM forward; must be in same units as vector of strikes
# Tx       option expiration date (years)
# pars     calibrated sabr parameters. vector of
#           [1] alpha
#           [2] rho
#           [3] nu
# beta     choice of beta (skew) during calibration
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
plot_sabr <- function(Ks, ivs, f, Tx, pars, beta, k_range = 0.2, scale = 1,
                      npoints = NULL, xlab = "strike / ATM forward",
                      ylab = "implied vol", main = NULL, cex = 1, lwd = 2) {
  # quick sanity checks
  if (length(Ks) != length(ivs)) {
    cat(plot_sabr__n, ": error: must be equal number of strikes and implied ", 
        "vols\n", sep = "")
    return(NULL)
  }
  if (length(pars) != 3) {
    cat(plot_sabr__n, ": error: parameter vector must be length 3\n", sep = "")
    return(NULL)
  }
  if (beta < 0 || beta > 1) {
    cat(plot_sabr__n, ": error: illegal value of beta; 0 <= beta <= 1\n", sep = "")
    return(NULL)
  }
  if (is.null(npoints) == FALSE && npoints < 1) {
    cat(plot_sabr__n, ": error: must sample at least one point\n", sep = "")
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
    main = paste("\u03b1 = ", round(pars[1], digits = 3), ", \u03b2 = ", beta,
                 ", \u03c1 = ", round(pars[2], digits = 3), ", \u03bd = ", round(
                   pars[3], digits = 3), sep = "")
  }
  # matrix plot; plot implied vols as points, sabr vols as interpolated line
  matplot(Ks * scale / f, cbind(ivs, sabr_ivol(f, Ks * scale, Tx, pars[1], beta,
                                               pars[2], pars[3])), pch = 18, 
          col = c("steelblue", "magenta"), lwd = lwd, lty = 1, type = c("p", "l"),
          xlab = xlab, ylab = ylab, main = main, font.main = 1, cex.main = cex,
          cex.lab = cex, cex.axis = cex, cex = cex)
}


# given a calibrated sabr model, return the corresponding put price, using the
# longstaff-schwartz least-squares monte carlo algorithm. if flavor is european,
# black_put will be called to return the price, and any parameters
# used to control monte carlo will be ignored. if flavor is american, monte
# carlo simulation will be conducted. NOT FUNCTIONAL
# step_scale
# changes number of steps and size of dt, as number of steps == step_scale * Tx,
# floored. use t-interval for wider tail.
sabr_put <- function(pars, beta, f, K, Tx, dfactor = Dexp, dfargs = .Dexp_r,
                     flavor = "european", sig_lvl = 0.05, step_scale = 365,
                     iter_num = 400) {
  if (length(pars) != 3 || Tx <= 0) { return(NULL) }
  if (flavor != "european" && flavor != "american") {
    cat(sabr_put__n, ": error: option flavor must be european or american\n",
        sep = "")
    return(NULL)
  }
  # if the option is european, just return black_call after calculating ivol
  if (flavor == "european") {
    ivol <- sabr_ivol(f, K, Tx, pars[1], beta, pars[2], pars[3])
    return(black_put(f, K, ivol, Tx, dfactor = dfactor, dfargs = dfargs))
  }
  # else the option is american, so we need to do monte carlo pricing. we control
  # the run length automatically by using sequential control.
  # current forward price and stochastic volatility
  Ft <- NULL
  vol_t <- NULL
  # give the parameters in pars names
  alpha <- pars[1]
  rho <- pars[2]
  nu <- pars[3]
  # "counter" rho, square root 1 - rho ^2 factor
  crho <- sqrt(1 - rho ^ 2)
  # number of steps in path
  nsteps <- floor(step_scale * Tx)
  # time delta and its square root
  dt <- 1 / step_scale
  dt_sqrt <- sqrt(dt)
  # for the number of iterations determined, generate the forward price matrix.
  # note that the forward price matrix holds exercise dates t1 ... t_nsteps
  # which corresponds to the nsteps time intervals that Tx has been split into.
  # therefore, we can simply discount by referencing index (since R 1-indexes)
  paths <- matrix(nrow = iter_num, ncol = nsteps)
  # for the number of specified paths
  for (i in 1:iter_num) {
    # set forward price and vol
    F_t <- f
    vol_t <- alpha
    # generate path
    for (j in 1:nsteps) {
      # generate differential wiener process increments
      dW_1 <- rnorm(1, mean = 0, sd = dt_sqrt)
      dW_2 <- rho * dW_1 + crho * rnorm(1, mean = 0, sd = dt_sqrt)
      # update forward price and volatility, and save forward price to path.
      # the absolute value is just in case the volatility becomes negative.
      F_t <- F_t + vol_t * (F_t ^ beta) * dW_1
      vol_t <- abs(vol_t * (1 + nu * dW_2))
      # save in index i, j
      paths[i, j] <- F_t
    }
  }
  # convert last column of forward prices at k to cash flows
  paths[, nsteps] <- pmax(0, K - paths[, nsteps])
  #print(paths[, 1:5])
  # after generating all the paths, we need to regress each kth column on the
  # respective k - 1th column. we use lm (QR decomposition) to get parameters
  for (c in nsteps:2) {
    # linear model for conditional expectation of DISCOUNTED cash flow at k
    # conditional on k - 1th state variable price. we first collect only the
    # indices of the ITM forward prices at k - 1, as these are the only ones
    # that have exercise value anyways.
    # get all the forward prices at k - 1 to make mask (we will need later too)
    F_all <- paths[, c - 1]
    # boolean mask for ITM forward paths only, length iter_num (number of rows)
    # note that this mask is from the k - 1 prices, not the cash flows at k
    itm_mask <- K - F_all > 0
    # get ITM forwards only at k - 1, set paths[itm_mask, c - 1] to be just cash
    # flows, and get discounted cash flows at k (current). we need all the
    # forwards in order to get a boolean mask of the correct length when
    # comparing estimated discounted cash flows against the cash flows on the
    #ITM cash flow paths.
    F_itm <- paths[itm_mask, c - 1]
    # set the OTM forward prices to discounted continuation value
    paths[!itm_mask, c - 1] <- dfactor(0, dt, dfargs) * paths[!itm_mask, c]
    # if there are no ITM cash flows, don't regress
    if (length(F_itm) > 0) {
      # else calculate exercise cash flows and continuation values
      paths[itm_mask, c - 1] <- pmax(0, K - paths[itm_mask, c - 1])
      dcfs <- dfactor(0, dt, dfargs) * paths[itm_mask, c]
      # build model using masked forward prices and DISCOUNTED cash flows. may 
      # decide to later replace this with a call to optim or nls instead.
      dlambda_model <- lm(dcfs ~ I(F_itm ^ 2) + F_itm)
      # retrieve parameters (debug?)
      dl_params <- dlambda_model$coefficients
      # get new vector of estimated dcfs (not sure of parameter order, check later)
      dcfs_hat <- predict.lm(dlambda_model, newdata = data.frame(F_itm = F_itm))
      # mask for ITM paths that will not be exercised
      nx_itm_mask <- paths[itm_mask, c - 1] < dcfs_hat
      paths[nx_itm_mask, c - 1] <- dfactor(0, dt, dfargs) * paths[nx_itm_mask, c]
    }
    # set any cash flow in the paths[itm_mask, c - 1] less than the ITM paths in
    # dcfs_hat to 0
    #print(itm_mask)
    #print(paths[, (c - 1):c])
    #print(dcfs_hat)
    
    #paths[!itm_mask, c - 1] <- dfactor(0, dt, dfargs) * paths[!itm_mask, c]
    #paths[itm_mask, c - 1][paths[itm_mask, c - 1] < paths[itm_mask, c]] <- 0
    #print(paths[, (c - 1):c])
  }
  #print(paths[, 1:5])
  #print(paths[, (nsteps - 5):nsteps])
  # after creating the cash flow matrix, start from the first time period for
  # each path i, find the first exercise date j of p where a positive cash flow
  # exists, and discount it j periods back. average all of these cash flows.
  # variable to hold all of the discounted cash flows and to get average from
  #print(paths[, 1:5])
  put_price <- 0
  #for (i in 1:iter_num) {
  #  for (j in 1:nsteps) {
  #    if (paths[i, j] > 0) {
  #      put_price <- put_price + dfactor(0, j * dt, dfargs) * paths[i, j]
  #      break
  #    }
  #  }
  #}
  #return(put_price / iter_num)
  return(sum(dfactor(0, dt, dfargs) * paths[, 1]) / iter_num)
}