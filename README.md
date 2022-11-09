# PAPF
Contains code which implements the positive area proportion function (PAPF) hypothesis test for the detection of clustering and/or dispersion in binary areal data as described in 'A Hypothesis Test for Detecting Distance-Specific Clustering and Dispersion in Areal Data', by Stella Self, Anna Overby, Anja Zgodic, David White, Alexander McLain and Caitlin Dyckman.

The file 'Calculate m_i(r y) for observed data.R' generates data on a square grid.  For each observed unit and each specified radius, the proportion of the circle centered at the centriod of that unit with the given radius is calculated.  
This file can be modified by replacing the square grid with a shapefile of the users choice.

The file 'Simulate the Null Distribution.R' performs Monte Carlo simulations to estimate the null distribution of the test statistic using a square grid.
For each simulation, the proportions of the circle centered at the centriod of each observed unit with each radius are returned.  
Multiple instances of this code can be run in parallel to improve computation time.

The file 'Perform the Hypothesis Test.R' takes the output from the other two files, computes the test statistic, and performs the hypothesis test.
