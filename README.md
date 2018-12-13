# SFG
Python scripts for SFG (Sum frequency generation) spectral data analysis, including custom extensions
 to cover the analysis and plotting of other spectra (IR, Raman, UV) and Langmuir compression isotherms
as well as automatized analysis of sample data gathered during a cruise in the Baltic Sea(Physical Chemistry)

## Modules

# scm

scm (abbreviation for SessinControllManager) is the basic class for interactive and IDE-based communication with the
database backend. It provides plotting functionality, systematic fetching of experimental data based on match criteria
and access to analytic tools.

# sfg

The sfg module contains the necessary classes for the analysis of Sum Frequency Generation spectroscopic raw data. It
handles normalization, holds values for the reference laser intensities, performs integration, smoothing and peak
picking.

# langmuir

This module provides a class to deal with experimental data of Langmuir trough measurements. It can provide values like
the maximum surface pressure, the lift-off point and the surface elasticity.

# gasex

The gasex module contains classes which simplify the interaction with the big and diverse data set obtained during the
GasEx cruise.

# spectools
An experimental module containing classes for the visualization and analysis of UV and infrared spectral data.

# texbackend

This experimental module allows the user to transform the output of other classes in this package directly into LaTeX
source code, including figures for plots.

# sfg_eval

This module provides a variety of plotting functions for lists of SFG spectra objects.

# gasex_eval

This module holds the functions that are used for the visualization of the sample data related to the GasEx cruise.

# lt_gui and new_guy

Auxiliary modules used for the Qt interface of the SessionControlManager (interactive plotting).