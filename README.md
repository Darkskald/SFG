# SFG
Python scripts for SFG (Sum frequency generation) spectral data analysis, including custom extensions
to cover the analysis and plotting of other spectra (IR, Raman, UV) and Langmuir compression isotherms
as well as automatized analysis of sample data gathered during a research cruise in the Baltic Sea. It relies on a 
*sqlite* database maintained by ORM (sqlalchemy) by containing all the experimental raw data and metainformation. This means that the import functionality
converting the raw data to the SQL table are logically separated from each other.

## Modules

### importer

The importer module collects all raw data and converts them in a form were they can be used easily by the database
handling classes to be persisted via ORM.

### orm

This module provides the necessary classes to interact with the measurement database via ORM.

### spectrum

The spectrum module provides a way to work with spectrum objects derived from a common abstract base 
class. This approach ensures that, independent of the individual kind of the spectrum , x and y values and their 
corresponding units (provided as TeX-strings) can be visualized conveniently. The spectra classes implement methods like
 normalization, smoothing and ASCII export to operate on their data conveniently.
 
 ### gasex
 
 This scripts were written to analyze the data obtained during the cruise in June and September 2018. It
 is based on the Pandas library for Data Science and is organized in a vertical hierarchy: The data of the single 
 measurements (SFG, Surface tTension, Langmuir Isotherms) is mapped to the corresponding **Sample** objects. All samples
 taken at a single cruise station are then mapped to a **Station** object. For each of this Station objects, average values
 for different parameters are calculated and added to a new dataframe which then can be used in interactive mode (eg. 
 Ipython or Jupyter Notebook) to perform analysis on the cruise data. 
 
 