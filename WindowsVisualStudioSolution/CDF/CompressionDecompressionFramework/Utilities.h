#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <random>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace Utilities {
	double* createMatrixWave(int x, int y, int z, double frequency, double amplitude, double phase, double fx, double fy, double fz);
	void writeWaveToCSV(double* matrix, int x, int y, int z, const std::string& filename);
	double* createMatrixRandom(int x, int y, int z, double minVal, double maxVal);
	double measureTimingOverhead();
	void printMatrix(const double* matrix, int x, int y, int z);
	double calculateMSE(const double* originalData, const double* decompressedData, int size);
	void printComparison(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z);
	size_t calculateOriginalDataBytes(size_t size);
}

#endif // _UTILITIES_H_
