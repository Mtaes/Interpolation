#pragma once
#define BUFFER_SIZE 100

//Contains informations about points
struct Informations {
	double  *x;
	double *fx;
	int count;
};

//Reads data from CSV file and returns it in Informations structure
Informations readFromCSV(char *fileName);

//Returns subset of points from original data
//wsp - resolution of subset (e.g value 2 - take every second point)
Informations preparePoints(int wsp, double *x, double *fx, int numberOfPoints);

//Use Lagrange polynomials
double lagrange(double searchX, Informations points);

//Use Cubic Spline interpolation
double cubicSplineInterpolation(double searchX, Informations points);