# RKAD
Contains code which implements Ripley's K function for areal data as described in 'A Generalization of Ripley's K function for the Detection of Spatial Clustering in Areal Data'.

The file Calculate m(r i y) for observed data.R generates data on a square grid.  For each observed unite and each specified radius, the proportion of the circle centered at the centriod of that unit with the given radius is calculated.  
This file can be modified by replacing the square grid with a shapefile of the users choice.

The file Simulate the Null Distribution.R performs Monte Carlo simulations to estimate the null distribution of the test statistic.  
However, the test statistic itself is not output, but rather for each simulation, the proportion of the circle centered at the centriod of that unit with the given radius.  
Multiple instances of this code can be run in parallel to improve computation time.

The file Perform the Hypothesis Test.R takes the output from the other two files, computes the test statistics, and performs the hypothesis test.
