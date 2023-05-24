#include <iostream>
#include <chrono>
#include <zfp.h>
#include <vector>
#include <numeric>
#include <cmath>
#include "Utilities.h"

#define RUNS 10000
#define WARMUP_RUNS 100

// Compress and decompress matrix using ZFP library
void compressAndDecompressMatrix(zfp_stream* zfp, zfp_field* field, size_t bufsize, std::vector<float>& data) {
    void* buffer = malloc(bufsize);
    bitstream* stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);
    zfp_compress(zfp, field);
    zfp_stream_rewind(zfp);

    // link the decompressed field with the data vector
    zfp_field* dec_field = zfp_field_3d(data.data(), zfp_type_float, field->nx, field->ny, field->nz);
    zfp_decompress(zfp, dec_field);
    
    zfp_field_free(dec_field); // free the decompressed field
    stream_close(stream);
    free(buffer);
}



// Calculate Mean Squared Error (MSE) between original and decompressed data
double calculateMSE(const std::vector<float>& originalData, const std::vector<float>& decompressedData) {
    double sum = 0.0;
    for (size_t i = 0; i < originalData.size(); i++) {
        double diff = originalData[i] - decompressedData[i];
        sum += diff * diff;
    }
    return sum / originalData.size();
}

// Measure the overhead of timing process
double measureTimingOverhead() {
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> overhead = end - start;
    return overhead.count();
}

void printMatrix(const std::vector<float>& matrix, int x, int y, int z) {
    for (int i = 0; i < z; i++) {
        std::cout << "z = " << i << ":\n";
        for (int j = 0; j < y; j++) {
            for (int k = 0; k < x; k++) {
                std::cout << matrix[i*x*y + j*x + k] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

// Run the compression/decompression experiment on a matrix of provided dimensions
void runExperiment(int x, int y, int z, bool useWave, bool visualizeData) {
    int size = x * y * z;
    std::vector<double> compressTimes(RUNS);
    std::vector<double> mseValues(RUNS);
    double timingOverhead = measureTimingOverhead();

    for (int i = 0; i < RUNS + WARMUP_RUNS; i++) {
        std::vector<float> originalData;
        std::vector<float> decompressedData(size);

        // Generate data for each run based on a wave or random distribution
        if (useWave)
            originalData = Utilities::createMatrixWave(x, y, z, 1.0f, 1.0f, 0.0f);
        else
            originalData = Utilities::createMatrixRandom(x, y, z, 0.0f, 1.0f);

        if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
            std::cout << "Original Data:\n";
            printMatrix(originalData, x, y, z);
        }

        zfp_field* field = zfp_field_3d(originalData.data(), zfp_type_float, x, y, z);
        zfp_stream* zfp = zfp_stream_open(NULL);

        // Using zfp for maximum accuracy
        zfp_stream_set_accuracy(zfp, 0.0); 
        size_t bufsize = zfp_stream_maximum_size(zfp, field);

        auto start = std::chrono::high_resolution_clock::now();
        compressAndDecompressMatrix(zfp, field, bufsize, decompressedData);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> compressTime = end - start;

        if (i >= WARMUP_RUNS) {
            compressTimes[i - WARMUP_RUNS] = compressTime.count() - timingOverhead;  
            mseValues[i - WARMUP_RUNS] = calculateMSE(originalData, decompressedData);
        }

        if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
            std::cout << "Decompressed Data:\n";
            printMatrix(decompressedData, x, y, z);
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

