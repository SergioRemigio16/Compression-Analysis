#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <random>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>

namespace Utilities {
	/**
	 * This function generates a 3D sinusoidal wave with customizable amplitude, frequency, phase, and propagation rates along each axis.
	 * The generated wave is stored in a one-dimensional array using row-major ordering.
	 * 
	 * @param x The size of the wave along the x dimension.
	 * @param y The size of the wave along the y dimension.
	 * @param z The size of the wave along the z dimension.
	 * @param amplitude The amplitude of the wave, i.e., the maximum displacement of the wave from its equilibrium position.
	 * @param phase The phase of the wave, i.e., the initial angle of a sinusoidal function at its origin.
	 * @param fx The propagation factor along the x dimension. A higher value means the wave propagates faster along this dimension.
	 * @param fy The propagation factor along the y dimension. A higher value means the wave propagates faster along this dimension.
	 * @param fz The propagation factor along the z dimension. A higher value means the wave propagates faster along this dimension.
	 * 
	 * @return A pointer to the dynamically allocated array storing the wave data. Each element represents the wave's value at a point in 3D space.
	 *         The size of the array is x * y * z. The caller of the function is responsible for deallocating this memory.
	 */
	double* createMatrixWave(int x, int y, int z, double amplitude, double phase, double fx, double fy, double fz);
	double* createWave1D(int n, double amplitude, double phase, double fx);
	void writeWaveToCSV(double* matrix, int x, int y, int z, const std::string& filename);
	double* createMatrixRandom(int x, int y, int z, double minVal, double maxVal);
	double measureTimingOverhead();
	void printMatrix(const double* matrix, int x, int y, int z);
	double calculateMSE(const double* originalData, const double* decompressedData, int size);
	/**
	 * Prints comparison of a 3d array
	*/
	void printComparison(const double* originalArray, const double* decompressedArray, int x, int y, int z, int originalMatrixBytes, int compressedSize);
	/**
	 * Prints comparison of a 1d array
	*/
	void printComparison(const double* originalArray, const double* decompressedArray, int n, int originalMatrixBytes, int compressedSize);
	/**
	 * Prints comparison of a 3d array
	*/
	void printError(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z);
	/**
	 * Prints error of a 1d array
	*/
	void printError(const double* originalArray, const double* decompressedArray, int n);
    /// <summary>
	/// Calculate byte size of 
	/// </summary>
	/// <param name="size"></param>
	/// <returns></returns>
	size_t calculateOriginalDataBytes(size_t size);
	double* readBinaryFile(const std::string& fileName, int& size);
}
#endif // _UTILITIES_H_
