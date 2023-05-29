#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <vector>
#include <random>


class Utilities {
public:
    // Generates a 3D matrix filled with a wave function
    static std::vector<double> createMatrixWave(int x, int y, int z, double frequency, double amplitude, double phase);
    // Generates a 3D matrix filled with random values between minVal and maxVal
    static std::vector<double> createMatrixRandom(int x, int y, int z, double minVal, double maxVal);
    static double measureTimingOverhead();
    static void printMatrix(const std::vector<double>& matrix, int x, int y, int z);
};

#endif // _UTILITIES_H_
