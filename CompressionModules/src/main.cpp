#include "TimingExperiment.h"
#include "Utilities.h"
#include "compression.h"

int main() {
    start();
    return 0;
}

/*
int main() {
	int x = 3, y = 7, z = 7;
	double* originalMatrix = Utilities::createMatrixWave(x, y, z, 0.1, 0.1, 0.1, 1.0, 1.0, 1);
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

/*

int main() {
    int x = 3, y = 7, z = 30; // Dimensions of matrix. Modify as needed.
    //int k = std::min(y, z) / 2;
    int k = 21;

    // Define original matrix
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);
    int originalMatrixBytes = x * y * z * sizeof(double);

    // Compress the matrix
    int compressedSize;
    unsigned char* compressedMatrix = SVDAlgorithms::compressMatrix(originalMatrix, x, y, z, k, compressedSize);

    // Decompress the matrix
    double* decompressedMatrix = SVDAlgorithms::decompressMatrix(compressedMatrix, compressedSize);

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

/*
int main() {
	//start();
}
*/

