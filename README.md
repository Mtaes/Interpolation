# Interpolation

Program lets you compare polynomial interpolation using Lagrange polynomials and cubic spline interpolation.

Data is loaded from csv file. Sample file data1.csv contains two columns; first is "x", second "y".

Program accepts two arguments: path to csv file and resolution of points used in interpolation, e.g. value of 2 means that only every second point will be used.

Results of polynomial interpolation are saved in "resultsLagrange.txt" file and results of cubic spline interpolation are saved in "resultsSplineInterpolation.txt". First column of both files contains values of "x", second values of "y" from original file and third values of "y" obtained through interpolation.