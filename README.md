# SFG

This package provides utilities to handle the results of measurement data acquired in the context 
of investigation of the Sea Surface Microlayer (SML). It deals with Vibrational Sum Frequency
Generation (VSFG) spectra, Langmuir compression isotherms and surface tension measurements. 

## spectrum

The spectrum sub-package contains tools to deal with spectra, implemented as classes with a unified
set of common properties, in a general way and with additional functionality with respect to VSFG
spectra and Langmuir compression isotherms.

## orm

The orm subpackage provides classes and functions to organize the data as a relational database
(in this case sqlite) to ensure a maximum of separation between analysis of the data and the initial
import. The idea behind this division is the comfort of having just a single portable database file
containing all results and metadata in a very accessible way rather than having to take care of 
dozens of subdirectories and hundreds of files. The package employs sqlalchemy as object-relational
mapper to have a maximum of abstraction from the actual SQL interaction.

## natural samples 

This sub-package bundles custom functionality to deal with data obtained during a time series 
surfactant observation and two consecutive cruises, targeting the interference between gas exchange
and surfactant presence, in the Baltic Sea.
