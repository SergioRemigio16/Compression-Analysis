#ifndef FFT_ALGORITHMS_H
#define FFT_ALGORITHMS_H

#include <vector>
#include <complex>

// Function to print various types of error between original and decompressed matrices
void printError(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z);

typedef std::vector<std::complex<double>> ComplexVec;
typedef std::pair<double, size_t> MagnitudeIndexPair;

// Structure to hold a complex value and its index
struct FrequencyComponent {
    std::complex<double> value;
    size_t index;
};

// Function to compare magnitudes for sorting
bool compareMagnitude(const MagnitudeIndexPair& a, const MagnitudeIndexPair& b);

// Function to serialize the compressed data to an unsigned char buffer for sending via MPI
unsigned char* serializeCompressedData(const std::vector<FrequencyComponent>& compressedData, size_t& size);

// Function to deserialize the compressed data from an unsigned char buffer after receiving via MPI
std::vector<FrequencyComponent> deserializeCompressedData(const unsigned char* buffer, size_t size);
double* decompressData(const std::vector<FrequencyComponent>& compressedData, size_t x, size_t y, size_t z);
std::vector<FrequencyComponent> compressData(double* originalMatrix, size_t x, size_t y, size_t z, double compressionRatio);

#endif // FFT_ALGORITHMS_H
