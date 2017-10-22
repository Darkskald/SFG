# SFG
Python scripts for SFG (Sum frequency generation) spectral data analysis (Physical Chemistry)

## Parts

### Classes
The Classes file contains the class library for the spectral data handling. Classes are provided for 
loading raw data, processing according to the usual procedure for SFG data, plotting and analysis of 
the data. This Classes are under development and currently at a very early stage.

### IPyInterpreter
This file containts the implementation of the IpyInterpreter Class. This class is meant to provide 
access to the SFG spectra database with using the Ipython command line and simplify data handling
by calling classes und functions from the rest of the module. When using this package, this class 
is the central instance of interaction for the user.

### IR
Useful tools to handle experimental IR and Raman data and compare them with SFG results

### UV
Tool to plot experimental UV/Vis spectra

### detailed_analysis
A section where functions for more detailed analysis of the spectra (peak picking, integration...)
are established and tested before final implementation.

### functions
This file is a container for useful code snippets which simplify chemistry-related and data analysis
related calculations
