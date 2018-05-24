#include "interpolation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


Informations readFromCSV(char *fileName) {
	Informations result;
	int numberOfPoints = 512;
	char delim = ',';
	char buffer[BUFFER_SIZE];
	char *number;
	int i = 0;

	FILE *fs = fopen(fileName, "r");
	if (!fs) {
		printf("The file could not be opened.\n");
		result.count = NULL;
		result.x = nullptr;
		result.fx = nullptr;
		return result;
	}

	result.x = (double*)malloc(sizeof(double) * numberOfPoints);
	result.fx = (double*)malloc(sizeof(double) * numberOfPoints);

	while (fgets(buffer, BUFFER_SIZE, fs) != NULL)
	{
		if (i == numberOfPoints) {
			numberOfPoints *= 2;
			result.x = (double *)realloc(result.x, sizeof(double) * numberOfPoints);
			result.fx = (double *)realloc(result.fx, sizeof(double) * numberOfPoints);
		}

		number = strtok(buffer, &delim);
		if (number != NULL) {
			result.x[i] = atof(number);
		}
		number = strtok(NULL, &delim);
		if (number != NULL) {
			result.fx[i] = atof(number);
		}
		i++;
	}

	result.count = i;

	fclose(fs);
	return result;
}

Informations preparePoints(int wsp, double *x, double *fx, int numberOfPoints) {
	Informations values;
	int limit = (numberOfPoints / wsp);
	if (numberOfPoints%wsp != 0) {
		limit++;
	}
	values.count = limit;
	values.x = (double *)malloc(sizeof(double)*(limit + 1));
	values.fx = (double *)malloc(sizeof(double)*(limit + 1));

	for (int i = 0; i < limit; i++) {
		values.x[i] = x[i * wsp];
		values.fx[i] = fx[i * wsp];
	}

	if (values.x[limit - 1] != x[numberOfPoints - 1] || values.fx[limit - 1] != fx[numberOfPoints - 1]) {
		values.x[limit] = x[numberOfPoints - 1];
		values.fx[limit] = fx[numberOfPoints - 1];
		values.count++;
		limit++;
	}

	return values;
}

double lagrange(double searchX, Informations points) {
	double sum = 0;
	double denominator;
	double numerator;

	for (int i = 0; i < points.count; i++) {
		numerator = 1;
		denominator = 1;
		for (int j = 0; j < points.count; j++) {
			if (i != j) {
				denominator *= (points.x[i] - points.x[j]);
				numerator *= (searchX - points.x[j]);
			}
		}
		sum += ((points.fx[i] * numerator) / denominator);
	}

	return sum;
};

double *calculateB(Informations points) {
	double *b = (double *)malloc(sizeof(double) * points.count);
	double temp;
	for (int i = 0; i < points.count; i++) {
		if (i == 0 || i == points.count - 1) {
			b[i] = 0.0;
		}
		else {
			b[i] = 6 / (points.x[i] - points.x[i - 1] + points.x[i + 1] - points.x[i]);
			temp = (points.fx[i + 1] - points.fx[i]) / (points.x[i + 1] - points.x[i]);
			temp -= (points.fx[i] - points.fx[i - 1]) / (points.x[i] - points.x[i - 1]);
			b[i] *= temp;
		}
	}
	return b;
}

double calculateMu(int i, Informations points) {
	double result;

	if (i == points.count - 1) {
		result = 0.0;
	}
	else {
		result = (points.x[i] - points.x[i - 1]);
		result /= (points.x[i] - points.x[i - 1] + points.x[i + 1] - points.x[i]);
	}

	return result;
}

double calculateLambda(int i, Informations points) {
	double lambda;

	if (i == 0) {
		lambda = 0.0;
	}
	else {
		lambda = (points.x[i + 1] - points.x[i]);
		lambda /= (points.x[i] - points.x[i - 1] + points.x[i + 1] - points.x[i]);
	}

	return lambda;
}

double **calculateA(Informations points) {
	double **A = (double **)malloc(sizeof(double) * points.count);
	for (int i = 0; i < points.count; i++) {
		A[i] = (double *)malloc(sizeof(double) * points.count);
	}

	for (int i = 0; i < points.count; i++) {
		for (int j = 0; j < points.count; j++) {
			if (i == j) {
				A[i][j] = 2.0;
			}
			else if (i == j - 1) {
				A[i][j] = calculateLambda(i, points);
			}
			else if (i == j + 1) {
				A[i][j] = calculateMu(i, points);
			}
			else {
				A[i][j] = 0.0;
			}
		}
	}

	return A;
}

void gauss(double **A, double *b, double *x, int N) {
	double m;
	for (int k = 0; k < N; k++) {
		for (int i = k + 1; i < N; i++) {
			m = A[i][k] / A[k][k];
			b[i] -= b[k] * m;
			for (int j = k; j < N; j++) {
				double odj = A[k][j] * m;
				A[i][j] -= odj;
			}
		}
	}
	x[N - 1] = b[N - 1] / A[N - 1][N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = b[i];
		for (int j = i + 1; j < N; j++) {
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}
}

double *calculateM(Informations points) {
	double *b;
	double *x;
	double **A;

	b = calculateB(points);
	A = calculateA(points);
	x = (double *)malloc(sizeof(double) * points.count);

	gauss(A, b, x, points.count);

	for (int i = 0; i < points.count; i++) {
		free(A[i]);
	}
	free(A);
	free(b);
	return x;
}

double cubicSplineInterpolation(double searchX, Informations points) {
	double result;
	int section = 0;
	double b;
	double d;

	double *m = calculateM(points);

	while (searchX > points.x[section + 1]) {
		section++;
	}

	b = (points.fx[section + 1] - points.fx[section]) / (points.x[section + 1] - points.x[section]);
	b -= (2 * m[section] + m[section + 1]) * (points.x[section + 1] - points.x[section]) / 6;

	d = (m[section + 1] - m[section]) / (6 * (points.x[section + 1] - points.x[section]));

	result = points.fx[section];
	result += b * (searchX - points.x[section]);
	result += (m[section] / 2) * (searchX - points.x[section]) * (searchX - points.x[section]);
	result += d * (searchX - points.x[section]) * (searchX - points.x[section]) * (searchX - points.x[section]);

	delete[] m;
	return result;
}