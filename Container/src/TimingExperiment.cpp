#include <fstream>
#include <iostream>
#include <chrono>
#include <zfp.h>
#include <vector>
#include <numeric>
#include <cmath>
#include "Utilities.h"
#include "CompressionDecompression.h"

#define RUNS 1000
#define WARMUP_RUNS 100

size_t calculateSize(const std::vector<double>& matrix) {
    return matrix.size() * sizeof(double);
}

size_t calculateSize(const CompressionResult& compressionResult) {
    return sizeof(compressionResult.x) + sizeof(compressionResult.y) + sizeof(compressionResult.z) +
        compressionResult.buffer.size() * sizeof(unsigned char) + sizeof(compressionResult.bufsize) +
        sizeof(compressionResult.rate);
}

void runExperiment(int x, int y, int z, bool useWave, bool visualizeData, int rate) {
    int size = x * y * z;
    std::vector<double> compressTimes(RUNS), decompressTimes(RUNS), mseValues(RUNS);
    std::vector<size_t> originalSizes(RUNS), compressedSizes(RUNS), decompressedSizes(RUNS), totalCompressedSizes(RUNS);
    double timingOverhead = Utilities::measureTimingOverhead();

    for (int i = 0; i < RUNS + WARMUP_RUNS; i++) {
        std::vector<double> originalMatrix(size);
        originalMatrix = useWave ? Utilities::createMatrixWave(x, y, z, 1.0f, 1.0f, 0.0f) :
            Utilities::createMatrixRandom(x, y, z, 0.0f, 1.0f);

        if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
            std::cout << "Original Data:\n";
            Utilities::printMatrix(originalMatrix, x, y, z);
        }

        auto start = std::chrono::high_resolution_clock::now();
        CompressionResult compressionResult = compressMatrixFixedRate(originalMatrix, x, y, z, rate);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> compressTime = end - start;

        start = std::chrono::high_resolution_clock::now();
        std::vector<double> decompressedMatrix = decompressMatrixFixedRate(compressionResult);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> decompressTime = end - start;

        if (i >= WARMUP_RUNS) {
            compressTimes[i - WARMUP_RUNS] = compressTime.count() - timingOverhead;
            decompressTimes[i - WARMUP_RUNS] = decompressTime.count() - timingOverhead;
            mseValues[i - WARMUP_RUNS] = Utilities::calculateMSE(originalMatrix, decompressedMatrix);
            originalSizes[i - WARMUP_RUNS] = calculateSize(originalMatrix);
            compressedSizes[i - WARMUP_RUNS] = compressionResult.bufsize;
            decompressedSizes[i - WARMUP_RUNS] = calculateSize(decompressedMatrix);
            totalCompressedSizes[i - WARMUP_RUNS] = calculateSize(compressionResult);
        }

        if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
            std::cout << "Decompressed Data:\n";
            Utilities::printMatrix(decompressedMatrix, x, y, z);
        }
    }

    double meanCompressTime = std::accumulate(compressTimes.begin(), compressTimes.end(), 0.0) / RUNS;
    double meanDecompressTime = std::accumulate(decompressTimes.begin(), decompressTimes.end(), 0.0) / RUNS;
    double meanMSE = std::accumulate(mseValues.begin(), mseValues.end(), 0.0) / RUNS;
    size_t meanOriginalSize = std::accumulate(originalSizes.begin(), originalSizes.end(), 0.0) / RUNS;
    size_t meanCompressedSize = std::accumulate(compressedSizes.begin(), compressedSizes.end(), 0.0) / RUNS;
    size_t meanDecompressedSize = std::accumulate(decompressedSizes.begin(), decompressedSizes.end(), 0.0) / RUNS;
    size_t meanTotalCompressedSize = std::accumulate(totalCompressedSizes.begin(), totalCompressedSizes.end(), 0.0) / RUNS;

    std::cout << "Rate: " << rate << ", " << x << "x" << y << "x" << z << ", " << meanCompressTime << ", " << meanDecompressTime << ", " << meanMSE
        << ", " << meanOriginalSize << ", " << meanCompressedSize << ", " << meanDecompressedSize << ", " << meanTotalCompressedSize << std::endl;
}

void runExperimentsForRates(int x, int y, int z, bool useWave, bool visualizeData, int minRate, int maxRate) {
    for (int rate = minRate; rate <= maxRate; ++rate) {
        runExperiment(x, y, z, useWave, visualizeData, rate);
    }
}

int runTimingExperiment() {
    int minRate = 1;
    int maxRate = 64;

    std::cout << "Rate, Matrix Size, Mean Compression Time (s), Mean Decompression Time (s), Mean Loss (MSE), Mean Original Size (bytes), Mean Compressed Size (bytes), Mean Decompressed Size (bytes), Mean Total Compressed Size (bytes) - Random Distribution" << std::endl;
    runExperimentsForRates(3, 7, 7, false, false, minRate, maxRate);
    //runExperimentsForRates(4, 7, 7, false, false, minRate, maxRate);
    //runExperimentsForRates(5, 7, 7, false, false, minRate, maxRate);
    //runExperimentsForRates(6, 7, 7, false, false, minRate, maxRate);
    //runExperimentsForRates(7, 7, 7, false, false, minRate, maxRate);

    return 0;
}
