#include <iostream>
#include <chrono>
#include <zfp.h>
#include <cstring>
#include "Utilities.h"

#define RUNS 10000

void copyMatrix(const std::vector<float>& src, std::vector<float>& dest) {
    dest = src;
}

void compressAndDecompressMatrix(zfp_stream* zfp, zfp_field* field, size_t bufsize, std::vector<float>& data) {
    void* buffer = malloc(bufsize);

    // Compress the data
    bitstream* stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);
    zfp_compress(zfp, field);

    // Decompress the data
    zfp_stream_rewind(zfp);
    zfp_decompress(zfp, field);

    stream_close(stream);
    free(buffer);
}

void runExperiment(int SIZE, bool useWave) {
    int size = SIZE * SIZE * SIZE;
    std::vector<float> originalData;
    std::vector<float> copyData(size);

    if (useWave)
        originalData = Utilities::createMatrixWave(SIZE, SIZE, SIZE, 1.0f, 1.0f, 0.0f);
    else
        originalData = Utilities::createMatrixRandom(SIZE, SIZE, SIZE, 0.0f, 1.0f);

    double totalCopyTime = 0.0;
    double totalCompressTime = 0.0;

    for (int i = 0; i < RUNS; i++) {
        // Time copying the matrix
        auto start = std::chrono::high_resolution_clock::now();
        copyMatrix(originalData, copyData);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> copyTime = end - start;
        totalCopyTime += copyTime.count();

        // Time compressing and decompressing the matrix
        zfp_field* field = zfp_field_3d(originalData.data(), zfp_type_float, SIZE, SIZE, SIZE);
        zfp_stream* zfp = zfp_stream_open(NULL);
        zfp_stream_set_rate(zfp, 8.0, zfp_type_float, 3, 0); // Most efficient settings
        size_t bufsize = zfp_stream_maximum_size(zfp, field);

        start = std::chrono::high_resolution_clock::now();
        compressAndDecompressMatrix(zfp, field, bufsize, originalData);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> compressTime = end - start;
        totalCompressTime += compressTime.count();

        zfp_stream_close(zfp);
        zfp_field_free(field);
    }

    std::cout << SIZE << "x" << SIZE << "x" << SIZE << "," << totalCopyTime / RUNS << "," << totalCompressTime / RUNS << std::endl;
}

int main() {
    std::cout << "Matrix Size,Copy Time (s),Compress-Decompress Time (s)" << std::endl;
    runExperiment(3, true);
    runExperiment(5, true);
    runExperiment(7, true);

    std::cout << "Matrix Size,Copy Time (s),Compress-Decompress Time (s) - Random Distribution" << std::endl;
    runExperiment(3, false);
    runExperiment(5, false);
    runExperiment(7, false);

    return 0;
}

