
#include <fstream>
#include <iostream>
#include <chrono>
#include <zfp.h>
#include <numeric>
#include <cmath>
#include <vector>
#include "Utilities.h"
#include "CompressionDecompression.h"
#include "TimingExperiment.h"

#define RUNS 1000
#define WARMUP_RUNS 100

size_t calculateSize(size_t size) {
    return size * sizeof(double);
}

size_t calculateSize(const CompressionResult& compressionResult) {
    return sizeof(compressionResult.x) + sizeof(compressionResult.y) + sizeof(compressionResult.z) +
        compressionResult.buffer.size() * sizeof(unsigned char) + sizeof(compressionResult.bufsize) +
        sizeof(compressionResult.rate);
}

void runExperiment(int x, int y, int z, bool useWave, bool visualizeData, int rate) {
    int size = x * y * z;
    double timingOverhead = Utilities::measureTimingOverhead();

    std::vector<double> compressTimes, decompressTimes, mseValues;
    std::vector<size_t> originalSizes, compressedSizes, totalCompressedSizes;

    // Preallocating memory
    compressTimes.reserve(RUNS);
    decompressTimes.reserve(RUNS);
    mseValues.reserve(RUNS);
    originalSizes.reserve(RUNS);
    compressedSizes.reserve(RUNS);
    totalCompressedSizes.reserve(RUNS);

    for (int i = 0; i < RUNS + WARMUP_RUNS; i++) {
        double* originalMatrix = nullptr;
        if (useWave)
            originalMatrix = Utilities::createMatrixWave(x, y, z, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
        else
            originalMatrix = Utilities::createMatrixRandom(x, y, z, 0.0, 1.0);

        if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
            std::cout << "Original Data:\n";
            Utilities::printMatrix(originalMatrix, x, y, z);
        }

        auto start = std::chrono::high_resolution_clock::now();
        CompressionResult compressionResult = compressMatrixFixedRate(originalMatrix, x, y, z, rate);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> compressTime = end - start;

        start = std::chrono::high_resolution_clock::now();
        double* decompressedMatrix = decompressMatrixFixedRate(compressionResult);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> decompressTime = end - start;

        if (i >= WARMUP_RUNS) {
            compressTimes.push_back(compressTime.count() - timingOverhead);
            decompressTimes.push_back(decompressTime.count() - timingOverhead);
            mseValues.push_back(Utilities::calculateMSE(originalMatrix, decompressedMatrix, size));
            originalSizes.push_back(calculateSize(size));
            compressedSizes.push_back(compressionResult.bufsize);
            totalCompressedSizes.push_back(calculateSize(compressionResult));
        }

        if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
            std::cout << "Decompressed Data:\n";
            Utilities::printMatrix(decompressedMatrix, x, y, z);
        }

        delete[] originalMatrix;
        delete[] decompressedMatrix;
    }

    double meanCompressTime = std::accumulate(compressTimes.begin(), compressTimes.end(), 0.0) / RUNS;
    double meanDecompressTime = std::accumulate(decompressTimes.begin(), decompressTimes.end(), 0.0) / RUNS;
    double meanMSE = std::accumulate(mseValues.begin(), mseValues.end(), 0.0) / RUNS;
    size_t meanOriginalSize = std::accumulate(originalSizes.begin(), originalSizes.end(), 0) / RUNS;
    size_t meanCompressedSize = std::accumulate(compressedSizes.begin(), compressedSizes.end(), 0) / RUNS;
    size_t meanTotalCompressedSize = std::accumulate(totalCompressedSizes.begin(), totalCompressedSizes.end(), 0) / RUNS;

    std::cout << "Rate: " << rate << ", " << x << "x" << y << "x" << z << ", " << meanCompressTime << ", " << meanDecompressTime << ", " << meanMSE
        << ", " << meanOriginalSize << ", " << meanCompressedSize << ", " << meanTotalCompressedSize << std::endl;
}

void runExperimentsForRates(int x, int y, int z, bool useWave, bool visualizeData, int minRate, int maxRate) {
    for (int rate = minRate; rate <= maxRate; ++rate) {
        runExperiment(x, y, z, useWave, visualizeData, rate);
    }
}

int runTimingExperiment() {
    int minRate = 1;
    int maxRate = 64;

    std::cout << "Rate, Matrix Size, Mean Compression Time (s), Mean Decompression Time (s), Mean Loss (MSE), Mean Original Size (bytes), Mean Compressed Size (bytes), Mean Total Compressed Size (bytes) - Random Distribution" << std::endl;
    //runExperimentsForRates(3, 7, 7, true, false, minRate, maxRate);

    int x = 3;
    int y = 3;
    int z = 7;
    double frequency = 1;
    double amplitude = M_PI;
    double phase = 0;

    double* testMatrix = Utilities::createMatrixWave(x, y, z, frequency, amplitude, phase, 1.0, 1.0, 7);
    Utilities::writeWaveToCSV(testMatrix, x, y, z, "wave_data.csv");

    //runExperimentsForRates(4, 7, 7, false, false, minRate, maxRate);
    //runExperimentsForRates(5, 7, 7, false, false, minRate, maxRate);
    //runExperimentsForRates(6, 7, 7, false, false, minRate, maxRate);
    //runExperimentsForRates(7, 7, 7, false, false, minRate, maxRate);

    return 0;
}

