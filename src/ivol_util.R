# contains functions for retrieving black-76 call and put prices, obtaining the
# black (lognormal) implied vol, and plotting the implied volatility surface.
#
# dependencies:
#
# rgl - for graphing the interactive volatility surface
#
# Changelog:
#
# 08-15-2019
#
# updated ATM points and time to expiration for both gas and index futures. also
# updated the comments on the ATM points, times to expiration, and update time.
# added very useful utility function that can take a vector of strikes and
# return a boolean mask that allows one to select some number of strikes from
# the original vector that are as equidistant from each other as possible.
#
# 08-14-2019
#
# realized i could have just downloaded implied volatility data from the cme
# website using their options settlement tool. i am assuming that these are
# black vols since these are calculated off of options on futures.
#
# 07-21-2019
#
# realized that somehow i had used american index options instead of european
# ones. the CME spx futures standard options are american style; only the end-
# of-month monthly options are european style. this change has been made.
#
# 07-16-2019
#
# updated plot_ivol_surface such that it plots color gradients for surfaces,
# for both the persp and rgl plots. made note about coercing strike bounds
# when plotting with plot_ivol_surface, and updated it to plot with regards to
# percentage max strike in ALL the data used, not just the displayed data.
# also changed description of the file since it is not just about sabr. changed
# the values of Tx and atm forward prices to create two sets of vectors for
# both sets of data that will be used. added term structure of discount rates.
#
# 07-14-2019
#
# changed name of file to ivol_util.R. moved all sabr-related functions to the
# new sabr.R file. moved over sabr_ivol, sabr_ivol_obj, fit_sabr, sabr_call,
# and the incomplete sabr_put functions to the sabr.R file. added how to remove
# NA values from a matrix and how to write to a csv without row names to the
# section with helpful bits of code. added to useful code bit section how to
# retrieve implied volatility, put it with price data, and rearrange col names.
#
# 07-10-2019
#
# added comment on correct viewing angle. tried to fix problems with longstaff-
# schwartz implementation. still not sure what i am doing wrong.
#
# 07-01-2019
#
# did some extra work on the monte carlo put pricing function. need to add a
# check for when none of the paths are ITM, else lm will throw "0 non-NA"error.
#
# 06-29-2019
#
# start to implement longstaff-schwartz algorithm for pricing puts.
#
# 06-26-2019
#
# checked that the call prices generated are close enough to that of the market
# price (they are). changed the default confidence interval width for the put
# pricing function to 0.05 to speed up termination, and added option to pass in
# max iterations. we have problems with the put price being too expensive for
# what it should be. corrected error in using the wrong discount factor for the
# european put and sabr_call.
#
# 06-23-2019
#
# finally got around to writing the monte carlo put pricing function. uses
# sequential run length control to determine when to terminate, but is still
# only a vanilla pricer. will add antithetic path generation to see if that
# improves anything at all. 
#
# 06-15-2019
#
# turns out i was slicing the range of strikes incorrectly; i.e. i should have
# sliced using the original vectors of strikes for implied vol and then strikes.
# slicing the strikes first and then re-slicing for implied vol gives nonsense.
# completed the "adaptive" calibration method, with some modifications.
#
# 06-13-2019
#
# added first draft of the "adaptive" calibration method. turns out that
# global optimization methods work decently well actually (can't really tell
# any performance drop for smile/skew calibration). however, there are problems
# with choice of initial parameters; what's happening?
#
# 06-12-2019
#
# note: implied volatility surface is valid only for calls and European puts
# only, as American puts command an early exercise premium (greater value than
# their intrinsic value, unlike European puts), so the volatility surface
# produced from American puts is completely different ITM and OTM. also when
# calibrating, calibrate using ATM options, as options too far into and out of
# the money have stale prices and low volume.
#
# 06-09-2019
#
# added another useful line of code to reuse; way to select dataframe entries
# relative to where their corresponding strikes are to the ATM forward price.
#
# 06-02-2019
#
# successfully fit the sabr model to the volatility smile. fit noticeably worse
# the more deep in the money or out of the money option implied vols are added.
# best fit seems to be around 0.8 - 1.2 times the ATM forward. started work on
# a function to iteratively optimize parameters to fit the model as close as
# possible to the volatility smile.
#
# 05-31-2019
#
# finally finished plot_ivol_surface. it has grown into a monster of a function,
# supporting both persp() and persp3d() plotting, with smart interpolation and
# sampling of points to produce a nice wireframe vol surface. i cannot believe
# this one function (which should have been trivial!) took this long. added
# analytical formula for the black european put. updated the "useful lines of
# code to reuse" section to change plotting just call smiles to both call and
# put smiles at the same time.
#
# 05-29-2019
#
# mostly complete with implied volatility surface plotting function. added
# string constants for all the functions that require names.
#
# 05-28-2019
#
# started work on function to plot an entire implied volatility surface.
# modified some comments on the implied volatility finding function.
#
# 05-27-2019
#
# initial creation. added objective function for ivol calibration, the
# analytical black-76 model for options on futures, and a discount function
# with constant discount rate. added sabr ivol function. added constants
# for different expirations and ATM forward prices (see below). added function
# to return implied volatilies for options given the market price and necessary
# black formula parameters.

# function names
black_ivol__n <- "black_ivol"
plot_ivol_surface__n <- "plot_ivol_surface"
get_strike_mas__n <- "get_strike_mask"

# ATM forward prices by expiry and expiry dates. for the european S&P index
# futures options, they on the last business day of the contract month. for the
# henry hub natural gas futures options, they expire on the day before the third
# last business day in the month BEFORE the expiration month. the underlying
# futures contract expires on that third to last business day. 
#
# ATM forwards correspond to the underlying futures contract that expires in
# the same month as the option written on it. note that when plotting the implied
# volatility surface, it is recommended to multiply the implied vols from the
# natural gas options by 100 or 1000 since the strikes are often decimals.
#
# dates from which expiry is calculated: 08/13/2019
#
# expiry months for all futures contracts:
# 09/2019, 10/2019, 11/2019, 12/2019, 01/2020

# at the money forwards and expirations for the index options
spx_atmf <- c(2932.75, 2932.75, 2932.75, 2934.25, 2934.25)
spx_Tx <- c(48/365, 79/365, 108/365, 140/365, 171/365)

# at the money forwards and expirations for the natural gas options
hh_atmf <- c(2.159, 2.232, 2.417, 2.53, 2.503)
hh_Tx <- c(43/365, 76/365, 104/365, 135/365, 168/365)

# constant instantaneous discount rate, default argument for args (was 0.03)
.Dexp_r <- c(0.02)

# term structure of risk-free treasury rates since 07/12/2019
.Dr_1m <- 0.0216
.Dr_2m <- 0.0218
.Dr_3m <- 0.0214
.Dr_6m <- 0.0207
.Dr_1y <- 0.0196
.Dr_2y <- 0.0184

# "typical" constant rate exponential discount function.
#
# parameters:
#
# Tx     expiration time in years
# tc     current time, also in years
# args   list of arguments (in this case, just discount rate (1 = 100%))
Dexp <- function(tc, Tx, args = .Dexp_r) {
  return(exp(-args[1] * (Tx - tc)))
}

# black model for european call options on futures.
#
# parameters:
#
# f        today's forward price
# K        strike price [vector]
# sigma    "real" (or equivalent implied) volatility [vector]
# Tx       expiration time in years (scalr)
# dfactor  optional discounting function; must take varargs, where first and
#          second arguments are current time and expiration time, and last is
#          a vector of parameters for the discount function (default Dexp)
# dfargs   optional vector of parameters for dfactor, default c(0.03)
#
# ex. using annual risk-free discount rate of 3%
# > black_call(2000, 2000, 0.35, 1, dfargs = c(log(1 + 0.03)))
# [1] 269.7483
black_call <- function(f, K, sigma, Tx, dfactor = Dexp, dfargs = .Dexp_r) {
  if (f == 0 || K == 0) { return(NULL) }
  if (sigma < 0 || Tx < 0) { return(NULL)}
  d1 <- (log(f / K) + 0.5 * sigma ^ 2 * Tx) / (sigma * sqrt(Tx))
  d2 <- d1 - sigma * sqrt(Tx)
  return(dfactor(0, Tx, dfargs) * (f * pnorm(d1) - K * pnorm(d2)))
}

# black model for european put options on futures.
#
# parameters are the same as in the above function for black calls.
black_put <- function(f, K, sigma, Tx, dfactor = Dexp, dfargs = .Dexp_r) {
  if (f == 0 || K == 0) { return(NULL) }
  if (sigma < 0 || Tx < 0) { return(NULL)}
  d1 <- (log(f / K) + 0.5 * sigma ^ 2 * Tx) / (sigma * sqrt(Tx))
  d2 <- d1 - sigma * sqrt(Tx)
  return(dfactor(0, Tx, dfargs) * (K * pnorm(-d2) - f * pnorm(-d1)))
}

# objective function for finding ivol from market price using black's model.
# returns squared errors as optimization criteria.
#
# do NOT pass a price vector; use individually for each price, forward, strike.
# optimization was done with the optimize() method.
#
# parameters:
#
# ivol     implied volatility parameter to be optimized
# price    market price of the call option
# is_put   optional boolean if for European put options ONLY; default FALSE
# 
# other parameters are the same as with black_call, but with no vectors
black_ivol_obj <- function(ivol, price, f, K, Tx, is_put = FALSE, 
                           dfactor = Dexp, dfargs = .Dexp_r) {
  # which black formula to use (call by default)
  black <- black_call
  # if we specified that this is a put option, use black_put
  if (is_put == TRUE) {
    black <- black_put
  }
  # return squared error between black model price and actual price
  return((price - black(f, K, ivol, Tx, dfactor = dfactor, 
                                  dfargs = dfargs)) ^ 2)
}

# returns a 2-column matrix of implied volatilies and the error between the
# market price and the black option price given the implied volatility, given a
# vector of market prices for options, the ATM forward price, a vector of
# strikes, the expiration date, and optional discount factor and arguments for
# the discount factor function.
#
# uses the optimize() function to find implied volatility, and will dynamically
# increase the upper search bound if the implied volatility estimate does not
# minimize black_ivol_obj() within the specified tolerance.
#
# reminder: strikes must be in units of the ATM forward price, not vice versa!
#
# parameters:
#
# prices     vector of call or put option prices
# f          current ATM forward price
# K          vector of call or put option strikes
# Tx         expiration date (year)
#
# dfactor, dfargs are optional discount factor arguments; read above comments
# for documentation.
#
# is_put     optional boolean to specify is we are extracting implied vols from
#            put prices instead of call prices. default FALSE, use ONLY if the
#            put prices in question are European puts!
# tol        optional tolerance for calculating implied volatility
# ubound     optional starting upper bound for implied volatility
# ilim       optional max number of iterations to increase upper search bound
#            before the algorithm deems the ivol estimate "good enough"
black_ivol <- function(prices, f, K, Tx, dfactor = Dexp, dfargs = .Dexp_r, 
                       is_put = FALSE, tol = 10 ^ (-6), ubound = 10, ilim = 5) {
  # if the price and strike vectors are of differing lengths, print error
  if (length(prices) != length(K)) {
    cat(black_ivol__n, ": error: price and strike vectors of differing ",
        "lengths\n", sep = "")
    return(NULL)
  }
  # lower bound, upper bound
  lb <- 0
  ub <- ubound
  # max number of iterations to increase upper search bound
  il <- ilim
  # number of prices/strikes
  N <- length(prices)
  # vector of 0s for ivol and errors
  ivol <- rep(0, N)
  err <- rep(0, N)
  # for each element, find its implied vol using optimize() and black_ivol_obj()
  for (i in 1:N) {
    # results from previous optimize() call
    res <- NULL
    # while there is no previous result (first iteration) or while error > tol
    # AND il > 0 (i.e. we can keep expanding the search upper bound)
    while ((is.null(res) == TRUE) || (res$objective > tol && il > 0)) {
      res <- optimize(black_ivol_obj, c(lb, ub), prices[i], f, K[i], Tx, 
                      dfactor = dfactor, dfargs = dfargs, is_put = is_put, 
                      tol = tol)
      # if the objective is greater than tol
      if (res$objective > tol) {
        # set lower bound to old upper bound, and crudely increase upper bound
        lb <- ub
        ub <- ub * 10   # this is very crude
        # decrement il
        il <- il - 1
      }
    }
    # reset il, lower bound, and upper bound
    il <- ilim
    lb <- 0
    ub <- ubound
    # set ivol[i] and err[i]
    ivol[i] <- res$minimum
    err[i] <- res$objective
  }
  # return matrix of ivol and err
  return(cbind(ivol, err))
}

# plot that generates the entire volatility surface. takes in multiple lists,
# each of which are assumed to have labeled columns, with consistent labeling
# across all the matrices. each matrix must have at the very least a column
# with implied volatilities and a corresponding column of strikes. the function
# will get the largest and smallest possible strikes from all strikes across
# all lists, and will then create a matrix of values where rows correspond to
# the strike  - (lowest_strike) + 1 (due to 1-indexing), and columns correspond
# to different maturities. returns named list of the matrix to send for vol
# plotting, and the vector of maturites (in order of columns) to use for scaling
# the expiration plotting axis. will plot the surface by default.
#
# built-in approx() is used for linear interpolation, persp() is used for
# plotting the standard surface, while persp3d() from rgl is used for the
# interactive surface (if desired).
#
# parameters:
#
# X           list of lists, each of which are assumed to have labeled columns
#             and contain the appropriate column labels for ivols and strikes.
# Tx          vector of expirations; must match the length of X.
# scale       optional integral value to scale expirations. default is 1. strikes
#             are automatically scaled to [0, 1]
# same_k      optional boolean indicating whether or not to plot only the strikes
#             that are shared among all expirations (give or take a few points)
# make_plot   optional boolean indicating whether to generate a plot or not.
# colors      optional vector of colors to be passed to persp() only. recommended
#             color vector is default. use the colors() command to see a list
#             of possible fill colors. note that using a single color will get
#             rid of the gradient effect.
# theta, phi  viewing angle parameters to pass to persp() only. theta controls
#             the rotation of the plot, while phi controls the viewing angle.
#             default recommended; common choices for viewing the surface are
#             (listed as {theta, phi}) {30, 5}, {40, 5}, {-30, 10}, {140, 13}
# opengl      optional boolean indicating whether to use rgl to generate the
#             interactive version of the surface; default FALSE.
# klab        optional string indicating which column contains strike prices,
#             if the labeling is different from the default.
# ilab        optional string indicating which column of ivols to use.
# xlab        optional string for the x label of the plot
# ylab        optional string for the y label of the plot
# zlab        optional string for the z label of the plot
# main        optional string for the plot title
plot_ivol_surface <- function(X, Tx, scale = 1, same_k = FALSE, make_plot = TRUE,
                              colors = c("skyblue1", "cyan", "green", "yellow",
                              "orange", "red"), theta = 30, phi = 5, opengl = FALSE,
                              klab = "strike", ilab = "ivol", 
                              xlab = "% max strike", ylab = "expiration",
                              zlab = "implied vol", main = "") {
  # if X is not a list, print error
  if (typeof(X) != "list") {
    cat(plot_ivol_surface, ": error: X must be of type list\n", sep = "")
    return (NULL)
  }
  # if the length of X and Tx do not match, print error
  if (length(X) != length(Tx)) {
    cat(plot_ivol_surface__n, ": error: length of X (", length(X), ") and Tx (", 
        length(Tx), ") must match\n", sep = "")
    return(NULL)
  }
  # number of expirations
  Nx <- length(Tx)
  # get the range's maximum and minimum strikes and global max
  kmin <- Inf
  kmax <- -Inf
  gkmax <- -Inf
  # if same_k is TRUE, we need the largest min and smallest max strikes
  if (same_k == TRUE) {
    # need to flip sign since we have changed the comparison method
    kmin <- -Inf
    kmax <- Inf
    for (i in 1:Nx) {
      kmin <- max(kmin, min(get(klab, X[[i]])))
      kmax <- min(kmax, max(get(klab, X[[i]])))
    }
  }
  # else just get global min and max
  else {
    # need to use double bracketing for lists
    for (i in 1:Nx) {
      # the get() function allows us to return the named object
      kmin <- min(kmin, get(klab, X[[i]]))
      kmax <- max(kmax, get(klab, X[[i]]))
    }
  }
  # get global max
  for (i in 1:Nx) {
    gkmax <- max(gkmax, get(klab, X[[i]]))
  }
  # if any of these values are illegal, print error and exit
  if ((kmin < 0) || (kmax < 0) || (kmax - kmin < 0)) {
    cat(plot_ivol_surface__n, ": error: illegal kmin (", kmin, ") and kmax (", 
        kmax,")\n", sep = "")
    return(NULL)
  }
  # make matrix of kmax - kmin + 1 rows and Nx columns, and for each expiration
  # (column) populate with each vector of interpolated implied vols. there will
  # be NAs because not each vector is of the same length; each will start at
  # the index of the minimum strike for that expiry - kmin + 1. the comparison
  # for selecting ony strikes that are >= kmin and <= kmax only has an effect
  # when same_k == TRUE
  ivol_surf <- matrix(nrow = kmax - kmin + 1, ncol = Nx)
  for (i in 1:Nx) {
    # get vector of strikes and vector of ivols, where strikes >= kmin, <= kmax
    kvec <- get(klab, X[[i]])
    ivec <- get(ilab, X[[i]])
    # use one ampersand for element-wise logical and; select implied volatilities
    # first before modifying kvec
    ivec <- ivec[kvec >= kmin & kvec <= kmax]
    kvec <- kvec[kvec >= kmin & kvec <= kmax]
    # result from linear interpolation
    int_iv <- approx(kvec, ivec, n = max(kvec) - min(kvec) + 1)
    # insert into ivol_surf column i at int_iv$x - kmin + 1
    ivol_surf[int_iv$x - kmin + 1, i] <- int_iv$y
  }
  
  # plot, if make_plot is TRUE, using [0, 1] x axis, scaled expiration vector
  # as y axis, and plot with either persp() or persp3d() based on value of the
  # opengl boolean.
  if (make_plot == TRUE) {
    
    #persp(y = scale * Tx, z = ivol_surf, theta = 35, phi = 35, xlab = xlab,
    #      ylab = ylab, zlab = zlab, main = main, ticktype = "detailed", 
    #      cex.axis = 0.7, cex.lab = 0.8, font.main = 1, border = NA, col = "magenta")
    #rgl.material(shininess = 20, col = "magenta", lwd = 0)
    
    # number of rows in ivol_surf
    ivol_n <- nrow(ivol_surf)
    # number of facets along x axis
    face_n <- 2 * log(ivol_n, base = 2)
    # if opengl is TRUE, plot with rgl
    if (opengl == TRUE) {
      # open new device so existing device focus is not overridden
      open3d()
      # get facet edge points and calculate the midpoints by using diagonals; no need
      # to calculate facets since rgl uses points only
      ivol_fps <- ivol_surf[seq(1, ivol_n, len = face_n), ]
      # calculate breaks for buckets of implied vols
      breaks <- hist(ivol_fps, plot = FALSE)$breaks
      # get colors for lowest breaks to highest breaks, cutting by breaks ivol_fmids
      # into different sector buckets, labeling with the colors produced from colors
      col_grad <- colorRampPalette(colors)(length(breaks) - 1)
      factors <- cut(ivol_fps, breaks = breaks, labels = col_grad)
      # plot with rgl; select only face_n points to make the surface appear like
      # an interpolated wireframe surface. surface itself is gradient unlit. use xlim
      # to plot only the relevant x points (as defined by the x sequences). divide
      # by the global max to correctly scale strikes
      persp3d(x = seq(kmin, kmax, len = face_n) / gkmax, 
              xlim = c(kmin, kmax) / gkmax, y = scale * Tx, 
              z = ivol_surf[seq(1, ivol_n, len = face_n), ], 
              xlab = xlab, ylab = ylab, zlab = zlab, cex = 0.8, col = col_grad[factors],
              ambient = "black", emission = "black", lit = FALSE,
              back = "fill", front = "fill", 
              smooth = TRUE)
      # add wireframe surface, using the same face_n points. although lit = FALSE
      # there is still lighting effect, why?
      surface3d(x = seq(kmin, kmax, len = face_n) / gkmax, y = scale * Tx, 
                z = ivol_surf[seq(1, ivol_n, len = face_n), ],
                front = "lines", back = "lines", lit = FALSE)
    }
    # else use persp() r default, with color gradient along z axis
    else {
      # calculate midpoint of each facet on implied vol surface
      # get facet edge points and calculate the midpoints by using diagonals
      ivol_fps <- ivol_surf[seq(1, ivol_n, len = face_n), ]
      ivol_fmids <- NULL
      # picture: uses the corners marked by 'x's
      #              ...
      #  i - 1  |___________|       
      #  i      |          x|
      #  i + 1  |x__________|
      for (i in 1:(Nx - 1)) {
        ivol_fmids <- c(ivol_fmids, (ivol_fps[-1, i] + ivol_fps[-face_n, i + 1]) / 2)
      }
      # calculate breaks for buckets of implied vols
      breaks <- hist(ivol_fmids, plot = FALSE)$breaks
      # get colors for lowest breaks to highest breaks, cutting by breaks ivol_fmids
      # into different sector buckets, labeling with the colors produced from colors
      col_grad <- colorRampPalette(colors)(length(breaks) - 1)
      factors <- cut(ivol_fmids, breaks = breaks, labels = col_grad)
      # plot using percentage of max strike, with scaled expirations and factors
      persp(x = seq(kmin, kmax, len = face_n) / gkmax, y = scale * Tx,
           z = ivol_surf[seq(1, ivol_n, len = face_n), ],
           xlab = xlab, ylab = ylab, zlab = zlab, cex.axis = 0.7, cex.lab = 0.8,
           ticktype = "detailed", col = col_grad[factors], theta = theta, phi = phi)
    }
  }
  # return matrix of interpolated implied volatilities
  return(ivol_surf)
  #return(ivol_surf_df)
}

# utility function. for a given vector of floats, sorted in ascending order,
# for a specified natural number npoints, will return a boolean mask for the
# original vector that will select npoints values from the original vector that
# are as equidistant to each other as possible. used for strike sampling, as
# sometimes there are more strikes clustered in some ranges than in others.
#
# parameters:
#
# Ks       vector of strikes, sorted in ascending order.
# npoints  number of TRUEs you want in the mask; i.e. number of equidistant
#          strikes you would like to get from the original Ks
get_strike_mask <- function(Ks, npoints) {
  # sanity check; Ks must not be NULL
  if (is.null(Ks) == TRUE) {
    cat(get_strike_mask__n, ": error: expected numeric vector, received NULL\n",
        sep = "")
    return(NULL)
  }
  if (npoints < 1) {
    cat(get_strike_mask__n, ": error: npoints must be a natural number\n",
        sep = "")
    return(NULL)
  }
  # if npoints > length of Ks, then simply return vector of TRUE of length(Ks)
  if (npoints > length(Ks)) {
    return(rep(TRUE, times = length(Ks)))
  }
  # create mask of all FALSE, length of Ks
  mask <- rep(FALSE, times = length(Ks))
  # we choose strikes that are most equidistant from each other, including the
  # maximum and minimum strikes. we use a loop and multiple min/max calls.
  # hypothetical range of length npoints perfectly equidistant strikes
  eqks <- seq(from = min(Ks), to = max(Ks), length.out = npoints)
  # for each one of the equidistant strikes, find the strike that is closest
  # to it and extract it with a boolean mask, and or that mask into mask
  for (i in 1:npoints) {
    # mask to select smallest value greater than eqks[i]
    up_mask <- Ks == min(Ks[Ks >= eqks[i]])
    # mask to select greatest value smaller than eqks[i]
    down_mask <- Ks == max(Ks[Ks <= eqks[i]])
    # calculate absolute differences between up_mask value and down_mask value
    up_diff <- abs(Ks[up_mask] - eqks[i])
    down_diff <- abs(Ks[down_mask] - eqks[i])
    # if up_diff <= down_diff, or up_mask with mask
    if (up_diff <= down_diff) {
      mask <- mask | up_mask
    }
    # else or it with down_mask
    else {
      mask <- mask | down_mask
    }
  }
  # return the mask
  return(mask)
}

# useful lines of code to reuse:
#
# [1] remove all NA rows from matrix
#
# a <- a[complete.cases(a), ]
#
# [2] write matrix to .csv file without writing row names as well
#
# write.csv(df, "file.csv", row.names = FALSE)
#
# [3] getting implied volatility from a price column, cbinding output to matrix,
#     rearranging columns, and changeing column names
#
# ivol <- black_ivol(price_vector, atmf, strike_vector, Tx (years))
# df <- cbind(df, ivol[, 1])
# df <- cbind(df, ivol[, 2])
# # way of rearranging columns is to index with a out-of-order numerical vector
# # where each integer corresponds to a column number
# df <- df[, c(1, 4, 5, 2, 3)]
# colnames(df) <- c( ... column names in desired order ...)
#
# [3] compare black prices to actual market call prices for a given expiry
#     (suppose 06/2019) to a fractional power to make the disparity more
#     apparent. the black model prices are clearly out of alignment, and
#     make very clear the existence of a non-flat volatility surface.
# 
#matplot(df$strike / (atmf[1] * 100), (cbind(df$call_prior_settle,
#        df1$put_prior_settle, black_call(atmf[1], df$strike / 100, 0.3, Tx[1]),
#        black_put(atmf[1], df$strike / 100, 0.3, Tx[1]))) ^ 0.2, type = "l",
#        lty = 1, col = c(4, 2, "skyblue", "pink"), xlab = "strike / ATM forward",
#        ylab = "(implied vol) ^ 0.2", 
#        main = "Market price vs. Black model price for 06/2019 expiation")
#
# [4] plot the implied volatility smile for a given expiry (suppose 06/2019),
#     using dots to show the exact points of ivol at ATM forward / strike. 
#     font.main = 1 changes the bold title font to more normal-looking font.
#     we opt to standardize by ATM f to make clear where the ATM point is and
#     to make reading the x-axis easier. compares puts and calls.
#
#matplot(df$strike / (atmf[1] * 100), cbind(df$call_ivol, df$put_ivol), 
#        pch = 18, col = c("steelblue", "palevioletred"), 
#        xlab = "strike / ATM forward", ylab = "implied vol", 
#        main = "Implied volatility for 06/2019 expiration", font.main = 1)
#
#     if drawing lines instead of points, remember to use lty = 1 with type = "l",
#     so that all the lines drawn are solid lines.
#     although c("steelblue", "palevioletred") looks nice, when plotting lines,
#     use default "blue" and "red" instead, so c(4, 2)
#
#     we can plot a vertical line to show where the ATM forward price is, and
#     add a legend to indicate which set of lines/dots is which option type. bty
#     gets rid of the [ugly] legend border.
#
#abline(v = 1, col = 2)
#legend("topright", legend = c("call_ivol", "put_ivol"),
#       col = c("steelblue", "palevioletred"), cex = 0.8, pch = 18, bty = "n")
#
# [5] select entries from a dataframe where the strike / 100 (to be in index
#     points like the ATM forward price) is within a certain range of the
#     corresponding ATM forward price. suppose we are using option data for the
#     06/2019 expiration. 0 <= L < 1 < U set upper and lower bounds for strikes
#
#df2[df2$strike / 100 >= L * atmf[2] & df2$strike / 100 <= U * atmf[2], ]
#
# [6] plot specific portion of strikes for ivol surface. choose a smaller frame
#     in the list and set same_k to TRUE to force it to be the boundary.
#
#plot_ivol_surface(list(df1, df2, df3, df4_atm), Tx, same_k = TRUE,
#                  make_plot = TRUE, theta = 140, phi = 13, opengl = FALSE)

