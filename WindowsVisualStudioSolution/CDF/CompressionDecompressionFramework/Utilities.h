#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <random>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>

class Utilities {
public:

    static double* createMatrixWave(int x, int y, int z, double frequency, double amplitude, double phase, double fx, double fy, double fz);
    static void writeWaveToCSV(double* matrix, int x, int y, int z, const std::string& filename);
    static double* createMatrixRandom(int x, int y, int z, double minVal, double maxVal);
    static double measureTimingOverhead();
    static void printMatrix(const double* matrix, int x, int y, int z);
    static double calculateMSE(const double* originalData, const double* decompressedData, int size);
    static void printComparison(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z);
};

#endif // _UTILITIES_H_
