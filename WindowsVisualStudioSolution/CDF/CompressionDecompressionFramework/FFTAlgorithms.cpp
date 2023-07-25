#include "FFTAlgorithms.h"

/*
Bitmask Usage Example:

frequency matrix: [1 + 2i, .001 + .001i, .002 + .001i, 3 + 4i]
zero out: [1 + 2i, 0, 0, 3 + 4i]
bitmask: 1001
bytearray: (x), (y), (z), (bitmask),(8 bytes for 0th component),(8 bytes for 3rd component)

decompression:
frequency matrix:  [1 + 2i, 0, 0, 3 + 4i]
notice there were only two components in the bytearray but the bitmask keeps track of when to fill components in.
*/

bool FFTAlgorithms::compareMagnitude(const MagnitudeIndexPair& a, const MagnitudeIndexPair& b) {
    return a.first < b.first;
}

unsigned char* FFTAlgorithms::compressMatrix(double* originalMatrix, int x, int y, int z, double compressionRatio, int& size) {
    // Performing FFT on the original matrix
    ComplexVec complexMatrix(x * y * z);
    fftw_plan p = fftw_plan_dft_r2c_3d(x, y, z, originalMatrix,
        reinterpret_cast<fftw_complex*>(complexMatrix.data()), FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    // Get magnitudes and corresponding indices
    std::vector<MagnitudeIndexPair> magnitudes;
    for (int i = 0; i < x * y * z; ++i) {
        magnitudes.push_back(MagnitudeIndexPair(std::abs(complexMatrix[i]), i));
    }

    // Sort the magnitudes
    std::sort(magnitudes.begin(), magnitudes.end(), compareMagnitude);

    // Determine how many frequency components to keep
    int numToKeep = static_cast<int>(compressionRatio * x * y * z);

    // Prepare the bitmask
    int bitMaskSize = (x * y * z + 7) / 8;
    unsigned char* bitMask = new unsigned char[bitMaskSize];
    memset(bitMask, 0, bitMaskSize);

    // Compute the size of the data and prepare the buffer
    size = bitMaskSize + numToKeep * sizeof(std::complex<double>) + 3 * sizeof(int);
    unsigned char* byteStream = new unsigned char[size];
    unsigned char* dataPtr = byteStream + bitMaskSize + 3 * sizeof(int);

    // Zero out least significant frequency components
    for (int i = 0; i < numToKeep; ++i) {
        int idx = magnitudes[x * y * z - 1 - i].second;
        bitMask[idx / 8] |= (1 << (idx % 8));
        std::memcpy(dataPtr, &complexMatrix[idx], sizeof(std::complex<double>));
        dataPtr += sizeof(std::complex<double>);
    }

    // Copy the bitmask to the byteStream
    std::memcpy(byteStream + 3 * sizeof(int), bitMask, bitMaskSize);

    // Write x, y, z to the start of the bytestream
    std::memcpy(byteStream, &x, sizeof(int));
    std::memcpy(byteStream + sizeof(int), &y, sizeof(int));
    std::memcpy(byteStream + 2 * sizeof(int), &z, sizeof(int));

    delete[] bitMask;

    return byteStream;
}

double* FFTAlgorithms::decompressMatrix(unsigned char* byteStream, int byteStreamSize) {
    // Read x, y, z from the start of the bytestream
    int x, y, z;
    std::memcpy(&x, byteStream, sizeof(int));
    std::memcpy(&y, byteStream + sizeof(int), sizeof(int));
    std::memcpy(&z, byteStream + 2 * sizeof(int), sizeof(int));

    // Calculate the size of the bitmask
    int bitMaskSize = (x * y * z + 7) / 8;

    // Calculate the number of frequency components to keep
    int numToKeep = (byteStreamSize - bitMaskSize - 3 * sizeof(int)) / sizeof(std::complex<double>);

    // Copy the bitmask from the byteStream
    unsigned char* bitMask = new unsigned char[bitMaskSize];
    std::memcpy(bitMask, byteStream + 3 * sizeof(int), bitMaskSize);

    // Prepare to read the data
    unsigned char* dataPtr = byteStream + bitMaskSize + 3 * sizeof(int);

    // Initialize complexMatrix with zeros
    ComplexVec complexMatrix(x * y * z, 0.0);

    // Populate the complexMatrix
    for (int i = 0; i < x * y * z; ++i) {
        if (bitMask[i / 8] & (1 << (i % 8))) {
            std::memcpy(&complexMatrix[i], dataPtr, sizeof(std::complex<double>));
            dataPtr += sizeof(std::complex<double>);
        }
    }

    delete[] bitMask;

    // Perform inverse FFT
    double* decompressedMatrix = new double[x * y * z];
    fftw_plan q = fftw_plan_dft_c2r_3d(x, y, z, reinterpret_cast<fftw_complex*>(complexMatrix.data()),
        decompressedMatrix, FFTW_ESTIMATE);
    fftw_execute(q);
    fftw_destroy_plan(q);

    // FFTW's backwards transform does not normalize the result, so we have to do it manually
    for (int i = 0; i < x * y * z; ++i) {
        decompressedMatrix[i] /= (x * y * z);
    }

    return decompressedMatrix;
}

/*
int main() {
    int x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    double compressionRatio = 0.56; // Compression ratio. 0.9 means keeping 90% of the orignial data

    // Define original matrix
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);

    // Compressing the data
    int byteStreamSize;
    unsigned char* byteStream = FFTAlgorithms::compressData(originalMatrix, x, y, z, compressionRatio, byteStreamSize);

    // Decompressing the received data
    double* decompressedMatrix = FFTAlgorithms::decompressData(byteStream, byteStreamSize);

    // Printing comparison of original and decompressed data
    Utilities::printComparison(originalMatrix, decompressedMatrix, x, y, z);

    // Printing size of original and serialized compressed data
    std::cout << "Size of original matrix (bytes): " << (x * y * z) * sizeof(double) << "\n";
    std::cout << "Size of serialized compressed data (bytes): " << byteStreamSize << "\n";

    // Printing various types of error between original and decompressed data
    Utilities::printError(originalMatrix, decompressedMatrix, x, y, z);

    // Freeing the memory
    delete[] originalMatrix;
    delete[] decompressedMatrix;
    delete[] byteStream;

    return 0;
}

*/

