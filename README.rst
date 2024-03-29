.. README.rst

   last updated on: 2022-02-04
   file created on: 2019-10-13

sabr-model
==========


.. image:: ./banner.png
   :alt: ./banner.png

A repo with some selected R code, plots of fits to implied volatility smiles,
and surface plots of implied volatility surfaces, conducted under a FAST grant
received from the NYU Dean's Undergraduate Research Fund.

The materials here are part of the research I produced during the summer
of 2019. In 2020 I decided to release the better parts of the material
onto GitHub, to serve as a testament to my first "real" attempt at pricing
options. Before, I had implemented the Cox binomial model in Python and
compared its shortcomings to market prices, but I consider my implementation to
be primitive and honestly not my best work. The code in this repository is not
production quality, but is at least prototyping quality. Most of the code
should work on versions of R >= 3.6.1.

Includes PDF slides accompanying my presentation video for NYU's 2020
Undergraduate Research Conference. Due to the COVID-19 outbreak in the United
States, the conference was switched to a virtual format.

Directories
-----------

``data``
~~~~~~~~

Contains the implied vol smile data for the models to fit. Implied vol data was
downloaded from CME's QuikStrike Option Settlement Tool and copied by hand to
the ``.csv`` files in the directory. Smile data is from European options on
E-mini S&P 500 futures and options on Henry Hub natural gas futures.

Expirations are monthly from September 2019 to January 2020.

``demos``
~~~~~~~~~

Contains ``hh_test.R`` and ``spx_test.R``, short demos that calibrate a model
(specify within the script) to a particular smile (specify in the script). As
noted in their heading comments, please run them from the R interpreter by
``source()``\ ing them, as graphics problems seem to occur when using RScript.
Currently configured to fit SABR.

``*_fits``
~~~~~~~~~~

Contains ``.png`` files with sample fits of CEV forward and SABR model fits to
implied vol data. ``hh_*`` directories contain fits to implied vol data from
European options on Henry Hub natural gas futures, while ``spx_*`` directories
contain fits to implied vol data from European options on E-mini S&P 500
futures. ``*_cevf_*`` indicates the CEV forward model is being fitted, while
``*_sabr_*`` indicates the SABR model is being fitted. Although not perfect,
SABR fits are pretty good for the most part, while the CEV forward model fits
tend to be much worse.

``presentation``
~~~~~~~~~~~~~~~~

Contains PDF slides for accompanying my presentation video for NYU's 2020
Undergraduate Research Conference.

``src``
~~~~~~~

Contains the main R code used for the project formerly located in the top-level
directory. File descriptions below.

``cevf.R``
   Contains functions for calibrating the CEV forward model to market implied
   vol data by least squares fitting of the market smiles. Contains a
   fitting function, objective function for the fitting function, implied
   vol approximation for the CEV forward model, and a function to plot the
   model's fit on a smile.

``ivol_util.R``
   Contains implied vol related utility functions, as well as data on ATM
   forward levels, times to expiration, and discount curve. Contains a
   discounting function, unused function for extracting Black volatilities, a
   function to plot an implied vol surface, functions to return Black call/put
   prices, and a useful strike sampling function. The end of the file contains
   several commented blocks code that I have kept as reference.

``sabr.R``
   Contains functions for calibrating the SABR model to market implied vol data
   by least squares fitting. Contains the same functions as in ``cevf.R`` but
   modified for SABR.

``vol_surfaces``
~~~~~~~~~~~~~~~~

Contains files ``vol_surface_hh_shared_crop.png`` and
``vol_surface_spx_shared_crop.png``, which are plots of implied vol surfaces
constructed from the Henry Hub and E-mini S&P 500 data respectively.