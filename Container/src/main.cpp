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



struct CompressionResult {
	int x;
	int y;
	int z;
	/// <summary>
	/// Holds the compressed data
	/// </summary>
	std::vector<unsigned char> buffer;  
	size_t bufsize;
	double precision;
};


CompressionResult compressData(std::vector<double>& originalData, int x, int y, int z, double precision) {
	// Creates a new ZFP field in 3D that will be compressed. The field takes the raw data from 'originalData',
	// specifies that the data type is double, and provides the dimensions of the 3D field (x, y, z).
	zfp_field* field = zfp_field_3d(originalData.data(), zfp_type_double, x, y, z); 
	// Opens a new ZFP stream for compression. The argument is NULL because we don't want to use a pre-existing stream.
	zfp_stream* zfp = zfp_stream_open(NULL);
	// Sets the compression precision of the ZFP stream. In this case, we want lossless compression,
	// so the precision is set to 0.0 (meaning no data will be lost in compression).
	zfp_stream_set_accuracy(zfp, precision);
	// Gets the maximum buffer size required for the compression.
	size_t bufsize = zfp_stream_maximum_size(zfp, field);
	// Create buffer that will hold the compressed data
	std::vector<unsigned char> buffer(bufsize);
	// Creates a bitstream that will be used to hold the compressed data.
	bitstream* stream = stream_open(buffer.data(), bufsize);
	// Associates the bitstream with the ZFP stream, meaning that when we compress our field,
	// the data will be written into this bitstream.
	zfp_stream_set_bit_stream(zfp, stream);
	// Compresses the field with the ZFP stream.
	zfp_compress(zfp, field);

	// Free and close ZFP fields and streams.
	stream_close(stream);
	zfp_stream_close(zfp);
	zfp_field_free(field);

	CompressionResult result;
	result.x = x;
	result.y = y;
	result.z = z;
	result.buffer = std::move(buffer);  // Avoid expensive copy with std::move (transfer ownership)
	result.bufsize = bufsize;
	result.precision = precision;

	return result;
}


void decompressData(std::vector<double>& decompressedData, CompressionResult& result) {
	zfp_stream* zfp = zfp_stream_open(NULL);
	zfp_stream_set_accuracy(zfp, result.precision);
	// Creates a bitstream that will be used to hold the compressed data.
	bitstream* stream = stream_open(result.buffer.data(), result.bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	// Creates a new ZFP field in 3D that will be decompressed.
	zfp_field* dec_field = zfp_field_3d(decompressedData.data(), zfp_type_double, result.x, result.y, result.z);
	// Decompresses the data into the 'dec_field'.
	zfp_decompress(zfp, dec_field);

	// Free and close ZFP fields and streams.
	zfp_field_free(dec_field);
	stream_close(stream);
	zfp_stream_close(zfp);
}


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
		std::vector<double> originalData(size), decompressedData(size);

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
		decompressData(decompressedData, compressionResult);
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

	// Compute the means and standard deviations of the compression times, decompression times, and MSEs
	double meanCompressTime = std::accumulate(compressTimes.begin(), compressTimes.end(), 0.0) / RUNS;
	double meanDecompressTime = std::accumulate(decompressTimes.begin(), decompressTimes.end(), 0.0) / RUNS;
	double meanMSE = std::accumulate(mseValues.begin(), mseValues.end(), 0.0) / RUNS;

	double sq_sum1 = std::inner_product(compressTimes.begin(), compressTimes.end(), compressTimes.begin(), 0.0);
	double stdCompressTime = std::sqrt(sq_sum1 / compressTimes.size() - meanCompressTime * meanCompressTime);

	double sq_sum2 = std::inner_product(decompressTimes.begin(), decompressTimes.end(), decompressTimes.begin(), 0.0);
	double stdDecompressTime = std::sqrt(sq_sum2 / decompressTimes.size() - meanDecompressTime * meanDecompressTime);

	double sq_sum3 = std::inner_product(mseValues.begin(), mseValues.end(), mseValues.begin(), 0.0);
	double stdMSE = std::sqrt(sq_sum3 / mseValues.size() - meanMSE * meanMSE);

	// Print the results
	std::cout << x << "x" << y << "x" << z << "," << meanCompressTime << "," << stdCompressTime << "," << meanDecompressTime << "," << stdDecompressTime << "," << meanMSE << "," << stdMSE << std::endl;
}

// The main function
int main() {
	// Print the header for the results
	std::cout << "Matrix Size,Mean Compression Time (s),STD Compression Time (s),Mean Decompression Time (s),STD Decompression Time (s),Mean Loss (MSE),STD Loss (MSE) - Wave Distribution" << std::endl;
	// Run the experiment for wave distribution
	runExperiment(3, 7, 7, true, true);
	/*
	runExperiment(4, 7, 7, true, false);
	runExperiment(5, 7, 7, true, false);
	runExperiment(6, 7, 7, true, false);
	runExperiment(7, 7, 7, true, false);

	// Print the header for the results
	std::cout << "Matrix Size,Mean Compression Time (s),STD Compression Time (s),Mean Decompression Time (s),STD Decompression Time (s),Mean Loss (MSE),STD Loss (MSE) - Random Distribution" << std::endl;
	// Run the experiment for random distribution
	runExperiment(3, 7, 7, false, false);
	runExperiment(4, 7, 7, false, false);
	runExperiment(5, 7, 7, false, false);
	runExperiment(6, 7, 7, false, false);
	runExperiment(7, 7, 7, false, false);
	*/

	return 0;
}