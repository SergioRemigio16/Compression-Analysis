#include "TimingExperiment.h"
#include "FFTAlgorithms.h"


int main() {
    size_t x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    double compressionRatio = 0.6; // Compression ratio. 0.5 means keeping half of the frequency components.

    // Define original matrix
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);

    // Compressing the data
    std::vector<FFTAlgorithms::FrequencyComponent> compressedData = FFTAlgorithms::compressData(originalMatrix, x, y, z, compressionRatio);

    // Serializing the compressed data for sending via MPI

    size_t serializedSize;
    unsigned char* serializedCompressedData = serializeCompressedData(compressedData, serializedSize);

    // Receiving the serialized data via MPI and deserializing it
    // Note: In this code, we simply reuse the serialized data for demonstration purposes. In a real application,
    // the data would be sent and received via MPI.
    std::vector<FFTAlgorithms::FrequencyComponent> receivedCompressedData = FFTAlgorithms::deserializeCompressedData(serializedCompressedData, serializedSize);

    // Decompressing the received data
    double* decompressedMatrix = decompressData(receivedCompressedData, x, y, z);

    // Printing comparison of original and decompressed data
    Utilities::printComparison(originalMatrix, decompressedMatrix, x, y, z);

    // Printing size of original and serialized compressed data
    std::cout << "Size of original matrix (bytes): " << (x * y * z) * sizeof(double) << "\n";
    std::cout << "Size of serialized compressed data (bytes): " << serializedSize << "\n";

    // Printing various types of error between original and decompressed data
    Utilities::printError(originalMatrix, decompressedMatrix, x, y, z);

    // Freeing the memory
    delete[] serializedCompressedData;
    delete[] originalMatrix;
    delete[] decompressedMatrix;

    return 0;
}
