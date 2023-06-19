#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <random>


class Utilities {
public:

    static double* createMatrixWave(int x, int y, int z, double frequency, double amplitude, double phase, double fx, double fy, double fz);
    static void writeWaveToCSV(double* matrix, int x, int y, int z, const std::string& filename);
    //static double* createMatrixWave(int x, int y, int z, double frequency, double amplitude, double phase);
    static double* createMatrixRandom(int x, int y, int z, double minVal, double maxVal);
    static double measureTimingOverhead();
    static void printMatrix(const double* matrix, int x, int y, int z);
    static double calculateMSE(const double* originalData, const double* decompressedData, int size);
};

#endif // _UTILITIES_H_
