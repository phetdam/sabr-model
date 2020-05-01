# calibrates the model to implied vol smile of European options on E-mini S&P
# 500 index futures expiring in December 2019, with time to expiration
# calculated from 08-13-2019.
#
# important note: recommended to run in R interpreter since Rscript cannot
# convert the unicode characters in graph and does not open display device.
# so in R interpreter, type source("./spx_test.R") if file is in current dir.
#
# Changelog:
#
# 05-01-2020
#
# changed data directory, as data file is not in current directory, and all
# directory references as file has been moved to the demos directory.
#
# 08-20-2019
#
# initial creation. requires both ivol_util.R, sabr.R, and the data file
# e-mini_smile_12-2019.csv to all be in the same directory.

# load required files
source("../src/ivol_util.R")
source("../src/sabr.R")
# read data from file
spx_data <- read.csv("../data/e-mini_smile_12-2019.csv")
# choose a priori beta, we will choose beta = 1 for demo purposes
beta <- 1
# calibrate and get parameter list; spx_atmf and spx_Tx are vectors of ATM forward
# and expiration data for the index futures, indexed by expiration. i usually
# select only ATM range 0.8 - 1.2 of forward, and sample a few points. the
# fit_sabr function can be found in sabr.R
pars <- fit_sabr(spx_data, spx_atmf[4], spx_Tx[4], beta, k_range = 0.2, npoints = 20)
# use the plot_sabr.R function to make a nice plot. cex controls point/label size.
plot_sabr(spx_data$strike, spx_data$ivol, spx_atmf[4], spx_Tx[4], pars$par, beta,
          k_range = 0.2, cex = 1.5, npoints = 20)