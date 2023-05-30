#include <iostream>
#include <chrono>
#include <zfp.h>
#include <vector>
#include <numeric>
#include <cmath>
#include "Utilities.h"
#include "CompressionDecompression.h"

#define RUNS 10000
#define WARMUP_RUNS 100

// Run the compression/decompression experiment on a matrix of provided dimensions
void runExperiment(int x, int y, int z, bool useWave, bool visualizeData) {
    int size = x * y * z;
    std::vector<double> compressTimes(RUNS);
    std::vector<double> mseValues(RUNS);
    double timingOverhead = Utilities::measureTimingOverhead();

    for (int i = 0; i < RUNS + WARMUP_RUNS; i++) {
        std::vector<double> originalData;
        std::vector<double> decompressedData(size);

        // Generate data for each run based on a wave or random distribution
        if (useWave)
            originalData = Utilities::createMatrixWave(x, y, z, 1.0f, 1.0f, 0.0f);
        else
            originalData = Utilities::createMatrixRandom(x, y, z, 0.0f, 1.0f);

        if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
            std::cout << "Original Data:\n";
            Utilities::printMatrix(originalData, x, y, z);
        }

        zfp_field* field = zfp_field_3d(originalData.data(), zfp_type_double, x, y, z);
        zfp_stream* zfp = zfp_stream_open(NULL);

        // Using zfp for maximum accuracy
        zfp_stream_set_accuracy(zfp, 0.0);
        size_t bufsize = zfp_stream_maximum_size(zfp, field);

        auto start = std::chrono::high_resolution_clock::now();
        zfp_field* compressedData = compressMatrix(zfp, field, bufsize, decompressedData);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> compressTime = end - start;

	    start = std::chrono::high_resolution_clock::now();
        decompressMatrix(zfp, field, bufsize, compressedData);
        end = std::chrono::high_resolution_clock::now();
        compressTime += end - start;


        if (i >= WARMUP_RUNS) {
            compressTimes[i - WARMUP_RUNS] = compressTime.count() - timingOverhead;
            mseValues[i - WARMUP_RUNS] = Utilities::calculateMSE(originalData, decompressedData);
        }

        if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
            std::cout << "Decompressed Data:\n";
            Utilities::printMatrix(decompressedData, x, y, z);
        }

        zfp_stream_close(zfp);
        zfp_field_free(field);
    }

    double meanCompressTime = std::accumulate(compressTimes.begin(), compressTimes.end(), 0.0) / RUNS;
    double meanMSE = std::accumulate(mseValues.begin(), mseValues.end(), 0.0) / RUNS;

    double sq_sum1 = std::inner_product(compressTimes.begin(), compressTimes.end(), compressTimes.begin(), 0.0);
    double stdCompressTime = std::sqrt(sq_sum1 / compressTimes.size() - meanCompressTime * meanCompressTime);

    double sq_sum2 = std::inner_product(mseValues.begin(), mseValues.end(), mseValues.begin(), 0.0);
    double stdMSE = std::sqrt(sq_sum2 / mseValues.size() - meanMSE * meanMSE);

    std::cout << x << "x" << y << "x" << z << "," << meanCompressTime << "," << stdCompressTime << "," << meanMSE << "," << stdMSE << std::endl;
}

int main() {
    std::cout << "Matrix Size,Mean Compress-Decompress Time (s),STD Compress-Decompress Time (s),Mean Loss (MSE),STD Loss (MSE) - Wave Distribution" << std::endl;
    runExperiment(3, 7, 7, true, false);
    runExperiment(4, 7, 7, true, false);
    runExperiment(5, 7, 7, true, false);
    runExperiment(6, 7, 7, true, false);
    runExperiment(7, 7, 7, true, false);

    std::cout << "Matrix Size,Mean Compress-Decompress Time (s),STD Compress-Decompress Time (s),Mean Loss (MSE),STD Loss (MSE) - Random Distribution" << std::endl;
    runExperiment(3, 7, 7, false, false);
    runExperiment(4, 7, 7, false, false);
    runExperiment(5, 7, 7, false, false);
    runExperiment(6, 7, 7, false, false);
    runExperiment(7, 7, 7, false, false);

    return 0;
}
