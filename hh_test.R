# calibrates the model to implied vol smile of European options on Henry Hub
# natural gas futures expiring in November 2019, with time to expiration
# calculated from 08-13-2019.
#
# important note: recommended to run in R interpreter since Rscript cannot
# convert the unicode characters in graph and does not open display device.
# so in R interpreter, type source("./hh_test.R") if file is in current dir.
#
# Changelog:
#
# 08-20-2019
#
# initial creation. requires both ivol_util.R, sabr.R, and the data file
# hh_smile_11-2019.csv to all be in the same directory.

# load required files
source("./ivol_util.R")
source("./sabr.R")
# read data from file
hh_data <- read.csv("./hh_smile_11-2019.csv")
# choose a priori beta, we will choose beta = 1 / 2 for demo purposes
beta <- 1 / 2
# calibrate and get parameter list; hh_atmf and hh_Tx are vectors of ATM forward
# and expiration data for natural gas futures, indexed by expiration. i usually
# select only ATM range 0.8 - 1.2 of forward, and sample a few points. the
# fit_sabr function can be found in sabr.R
pars <- fit_sabr(hh_data, hh_atmf[3], hh_Tx[3], beta, k_range = 0.2, npoints = 20)
# use the plot_sabr.R function to make a nice plot. cex controls point/label size.
plot_sabr(hh_data$strike, hh_data$ivol, hh_atmf[3], hh_Tx[3], pars$par, beta,
          k_range = 0.2, cex = 1.5, npoints = 20)