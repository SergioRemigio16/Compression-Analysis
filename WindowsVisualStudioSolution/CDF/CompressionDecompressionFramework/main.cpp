/*
A zfp_field is a descriptor for an array of data. It contains information
about the array, including the type of data it contains (e.g., floating point numbers, integers),
the dimensions of the array (how many elements there are in each dimension),
and a pointer to the raw data itself. Essentially, a zfp_field is used to tell the ZFP library
what data you want to compress (or have just decompressed).

A zfp_stream represents a compressed stream and maintains the state for a sequence of
read or write operations. It's associated with a zfp_field for the data being compressed
or decompressed, and a bitstream where the actual compressed data is written or read.
Compression parameters, such as the desired compression rate or precision,
are also set on the zfp_stream.

*/

#include <fstream>
#include <iostream>
#include <chrono>
#include <zfp.h>
#include <vector>
#include <numeric>
#include <cmath>
#include "Utilities.h"
#include "CompressionDecompression.h"

// Define the number of runs for the experiment and the warm-up runs
#define RUNS 10000
#define WARMUP_RUNS 100



// This function runs the compression/decompression experiment on a 3D matrix of provided dimensions
void runExperiment(int x, int y, int z, bool useWave, bool visualizeData) {
	double precision = 0.0;
	// Calculate total size of the 3D matrix
	int size = x * y * z;
	// Create vectors to store the times for compression and decompression, and the MSE values
	std::vector<double> compressTimes(RUNS), decompressTimes(RUNS), mseValues(RUNS);
	// Measure the timing overhead
	double timingOverhead = Utilities::measureTimingOverhead();

	// Run the experiment for the specified number of runs plus the warm-up runs
	for (int i = 0; i < RUNS + WARMUP_RUNS; i++) {
		std::vector<double> originalData(size);

		// Depending on useWave, generate data for each run based on a wave or random distribution
		originalData = useWave ? Utilities::createMatrixWave(x, y, z, 1.0f, 1.0f, 0.0f) :
			Utilities::createMatrixRandom(x, y, z, 0.0f, 1.0f);

		if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
			std::cout << "Original Data:\n";
			Utilities::printMatrix(originalData, x, y, z);
		}

		auto start = std::chrono::high_resolution_clock::now();
		CompressionResult compressionResult = compressData(originalData, x, y, z, precision);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> compressTime = end - start;

		start = std::chrono::high_resolution_clock::now();
		std::vector<double> decompressedData = decompressData(compressionResult);
		end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> decompressTime = end - start;

		// If the run is not a warm-up run, record the times and calculate the MSE
		if (i >= WARMUP_RUNS) {
			compressTimes[i - WARMUP_RUNS] = compressTime.count() - timingOverhead;
			decompressTimes[i - WARMUP_RUNS] = decompressTime.count() - timingOverhead;
			mseValues[i - WARMUP_RUNS] = Utilities::calculateMSE(originalData, decompressedData);
		}

		if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
			std::cout << "Decompressed Data:\n";
			Utilities::printMatrix(decompressedData, x, y, z);
		}
	}

	// Compute the means of the compression times, decompression times, and MSEs
	double meanCompressTime = std::accumulate(compressTimes.begin(), compressTimes.end(), 0.0) / RUNS;
	double meanDecompressTime = std::accumulate(decompressTimes.begin(), decompressTimes.end(), 0.0) / RUNS;
	double meanMSE = std::accumulate(mseValues.begin(), mseValues.end(), 0.0) / RUNS;

	// Print the results
	std::cout << x << "x" << y << "x" << z << "," << meanCompressTime << "," << meanDecompressTime << "," << meanMSE << std::endl;
}

// The main function
int main() {
	// Print the header for the results
	std::cout << "Matrix Size,Mean Compression Time (s),Mean Decompression Time (s),Mean Loss (MSE) - Wave Distribution" << std::endl;
	// Run the experiment for wave distribution
	runExperiment(3, 7, 7, true, false);
	runExperiment(4, 7, 7, true, false);
	runExperiment(5, 7, 7, true, false);
	runExperiment(6, 7, 7, true, false);
	runExperiment(7, 7, 7, true, false);

	// Print the header for the results
	std::cout << "Matrix Size,Mean Compression Time (s),Mean Decompression Time (s),Mean Loss (MSE) - Random Distribution" << std::endl;
	// Run the experiment for random distribution
	runExperiment(3, 7, 7, false, false);
	runExperiment(4, 7, 7, false, false);
	runExperiment(5, 7, 7, false, false);
	runExperiment(6, 7, 7, false, false);
	runExperiment(7, 7, 7, false, false);

	return 0;
}
