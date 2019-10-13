# sabr_model

![./banner.png] (./banner.png)

by Derek Huang

_last updated on: 10-13-2019_  
_file created on: 10-13-2019_

The sabr_model repo contains some selected R code, plots of fits to implied volatility smiles, and surface plots of implied volatility surfaces. Conducted under a FAST grant received from the NYU Dean's Undergraduate Research Fund, the materials here are part of the research I produced during the summer of 2019. Only recently did I decide to release the better parts of the material onto GitHub, to serve as a testament to my first "real" attempt at pricing options. Before, I had implemented the Cox binomial model in Python and compared its shortcomings to market prices, but I consider my implementation to be primitive and honestly not my best work. The code in this repository is not production quality, but is at least prototyping quality. I think most of it should work on versions of R >= 3.6.1.

## Source

The main R code used for the project, located in the top-level directory.

* __cevf.R:__ Contains functions for calibrating the CEV forward model to market implied volatility data by least squares fitting of the market smiles. Contains a fitting function, objective function for the fitting function, implied volatility approximation for the CEV forward model, and a function to plot the model's fit on a smile.

* __ivol_util.R:__ Contains implied volatility related utility functions, as well as data on ATM forward levels times to expiration, and discount curve. Contains a discounting function, unused function for extracting Black volatilities, a function to plot an implied vol surface, functions to return Black call/put prices, and a useful strike sampling function. The file's end also contains several blocks of commented-out code that I commonly used and ended up recording so I wouldn't forget them.

* __sabr.R:__ Contains functions for calibrating the SABR model to market implied volatility data by least squares fitting. Contains the same functions as in cevf.R but modified for SABR.

__hh_test.R__ and __spx_test.R__ are sample testing scripts that calibrate a model (specify within the script) to a particular smile (specify in the script). As noted in their heading comments, please run them from the R interpreter by source()ing them, as graphics problems seem to occur when RScript is used. Contains same functions as in cevf.R, but for SABR instead.

## Directories

### data

Contains the implied volatility smile data for the models to fit. Implied volatility data  was downloaded from CME's QuikStrike Option Settlement Tool, and copied by hand to the .csv files in the directory. Smile data is from European options on E-mini S&P 500 futures and options on Henry Hub natural gas futures. Expirations are monthly from September 2019 to January 2020.

### \*\_fits

Contains .png files depicting sample fits of CEV forward and SABR model fits to the implied volatility data. _hh\_\*_ directories contain fits to implied vol data from European options on Henry Hub natural gas futures, while _spx\_\*_ directories contain fits to implied vol data from European options on E-mini S&P 500 futures. _\*cevf\*_ indicates the CEV forward model is being fitted, while _\*sabr\*_ indicates the SABR model is being fitted.

## Images

The two .png files in the top-level directory are plots of implied vol surfaces, one for the Henry hub data (vol\_surface\_hh\_shared\_crop.png), and one for the E-mini S&P 500 data (vol\_surface\_spx\_shared_crop.png). 