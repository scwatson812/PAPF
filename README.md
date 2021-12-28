# RKAD
Contains code which implements Ripley's K function for areal data as described in 'A Generalization of Ripley's K function for the Detection of Spatial Clustering in Areal Data'

The file Calculate Test Statistic.R generates data on a square grid and computes the test statistic based on Ripley's K function for areal data.  
This file can be modified by replacing the square grid with a shapefile of the users choice.

The file Simulate the Null Distribution of the Test Statistics.R performs Monte Carlo simulations to estimate the null distribution of the test statistic.  
Multiple instances of this code can be run in parallel to improve computation time.
