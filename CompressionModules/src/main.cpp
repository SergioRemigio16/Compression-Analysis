#include "TimingExperiment.h"
#include "Utilities.h"
#include "compression.h"


int main() {
    int x = 30;
    int y = 70;
    int z = 70;
    int n = x * y * z;
    //double* originalMatrix = Utilities::readBinaryFile("../../data/nlsm_uncompressed_dumped_SEND_from_0_to_1.bin", n);
    double* originalMatrix = Utilities::createWave1D(n, 1, 0, 1);
    double threshold = 0.94;

    // Define original matrix
    int originalMatrixBytes = n * sizeof(double);

    // Compress the matrix
    int compressedSize;
    unsigned char* compressedMatrix = FFTAlgorithms::compressMatrix1D(originalMatrix, n, threshold, compressedSize);

    // Decompress the matrix
    double* decompressedMatrix = FFTAlgorithms::decompressMatrix1D(compressedMatrix, compressedSize);

    // Printing comparison of original and decompressed data
    Utilities::printComparison(originalMatrix, decompressedMatrix, n, originalMatrixBytes, compressedSize);
    // Printing various types of error between original and decompressed data
    Utilities::printError(originalMatrix, decompressedMatrix, n);

    // Freeing the memory
    delete[] originalMatrix;
    delete[] compressedMatrix;
    delete[] decompressedMatrix;

    return 0;

    return 0;
}

/*
int main() {
    start();
    return 0;
}
*/


/*
int main() {
    int x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    int n = x * y * z;
    int k = 7; 

    // Define original matrix
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 2.14, 0, 3.14, 3.14, 3.14);
    int originalMatrixBytes = x * y * z * sizeof(double);

    // Compress the matrix
    int compressedSize;
    unsigned char* compressedMatrix = SVDAlgorithms::compressMatrix1dsq(originalMatrix, n, k, compressedSize);

    // Decompress the matrix
    double* decompressedMatrix = SVDAlgorithms::decompressMatrix1dsq(compressedMatrix, compressedSize);

    // Printing comparison of original and decompressed data
    Utilities::printComparison(originalMatrix, decompressedMatrix, x, y, z);

    std::cout << "Original matrix size: " << originalMatrixBytes << " bytes" << std::endl;
    std::cout << "Compressed matrix size: " << compressedSize << " bytes" << std::endl;
    // Printing various types of error between original and decompressed data
    Utilities::printError(originalMatrix, decompressedMatrix, x, y, z);

    // Freeing the memory
    delete[] originalMatrix;
    delete[] compressedMatrix;
    delete[] decompressedMatrix;

    return 0;
}
*/

/*
int main() {
	int x = 3, y = 7, z = 7;
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 2.14, 0, 3.14, 3.14, 3.14);
	// The first degree polynomials for each dimension
	int N = 3;
	int Q = 5;
	int S = 5;
	
	int bufferSize;

	unsigned char* buffer = ChebyshevAlgorithms::compressMatrix(originalMatrix, x, y, z, N, Q, S, bufferSize);
	double* decompressedMatrix = ChebyshevAlgorithms::decompressMatrix(buffer, bufferSize);

	Utilities::printComparison(originalMatrix, decompressedMatrix, x, y, z);
	Utilities::printError(originalMatrix, decompressedMatrix, x, y, z);

	// Printing size of original and serialized compressed data
	std::cout << "Size of original matrix (bytes): " << (x * y * z) * sizeof(double) << "\n";
	std::cout << "Size of serialized compressed data (bytes): " << bufferSize << "\n";

	delete[] decompressedMatrix;
	delete[] originalMatrix;

	return 0;

}
*/
