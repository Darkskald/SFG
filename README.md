# SFG
Python scripts for SFG (Sum frequency generation) spectral data analysis, including custom extensions
to cover the analysis and plotting of other spectra (IR, Raman, UV) and Langmuir compression isotherms
as well as automatized analysis of sample data gathered during a research cruise in the Baltic Sea. It relies on a 
*sqlite* database containing all the experimental raw data and metainformation. This means that the import functionality
converting the raw data to the SQL table are logically separated from each other.

## Modules

### newport

The new import module collecting all metainformation and raw data and writing it into the SQL table used by the rest of
the program. This module is meant to be kept in the same directory as the raw data on the local drive of the author, but 
for reasons of convenience it is stored in the same GitHub repository.

### spectrum

The spectrum module provides a way to work with spectrum objects derived from a common abstract base 
class. This approach ensures that, independent of the individual kind of the spectrum , x and y values and their 
corresponding units (provided as TeX-strings) can be visualized conveniently. The spectra classes implement methods like
 normalization, smoothing and ASCII export to operate on their data conveniently.
 
 ### novel_gasex
 
 This reimplementation of the old analysis code for the data obtained during the cruise in June and September 2018. It
 is based on the Pandas library for Data Science and is organized in a vertical hierarchy: The data of the single 
 measurements (SFG, Surface tTension, Langmuir Isotherms) is mapped to the corresponding **Sample** objects. All samples
 taken at a single cruise station are then mapped to a **Station** object. For each of this Station objects, average values
 for different parameters are calculated and added to a new dataframe which then can be used in interactive mode (eg. 
 Ipython or Jupyter Notebook) to perform analysis on the cruise data. To simplify the comparison between both the
 cruises, the class **Cruise** collects all Stations belonging to them and calculates the desired values cruise-wise.
 
 ## gasex_analysis
 
 A helper module to provide plotting functions for the data obtained due to the **novel_gasex** functionality. 
 
 ## legacy code
 
 The rest of the Python scripts in this repository is legacy code from earlier versions of the data analysis routines
 and will be converted to fit into the new framework in the near future.