#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "interpolation.h"


int main(int argc, char *argv[]) {
	int skipPoints = atoi(argv[2]);

	Informations originalPoints = readFromCSV(argv[1]);

	if (originalPoints.count == NULL) {
		return 1;
	}

	Informations points = preparePoints(skipPoints, originalPoints.x, originalPoints.fx, originalPoints.count);

	FILE *fp = fopen("resultsLagrange.txt", "w");
	if (fp == NULL) {
		printf("The file could not be created.\n");
		return 1;
	}

	double temp;
	for (int i = 0; i < originalPoints.count; i++) {
		temp = lagrange(originalPoints.x[i], points);
		fprintf(fp, "%.14f\t%.14f\t%.14f\n", originalPoints.x[i], originalPoints.fx[i], temp);
	}
	fclose(fp);

	fp = fopen("resultsSplineInterpolation.txt", "w");
	if (fp == NULL) {
		printf("The file could not be created.\n");
		return 1;
	}

	for (int i = 0; i < originalPoints.count; i++) {
		temp = cubicSplineInterpolation(originalPoints.x[i], points);
		fprintf(fp, "%.14f\t%.14f\t%.14f\n", originalPoints.x[i], originalPoints.fx[i], temp);
	}
	fclose(fp);

	free(points.fx);
	free(points.x);

	free(originalPoints.x);
	free(originalPoints.fx);

	return 0;
}