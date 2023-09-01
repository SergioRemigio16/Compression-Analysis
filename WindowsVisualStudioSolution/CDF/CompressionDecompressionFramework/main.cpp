#include "TimingExperiment.h"
#include "FFTAlgorithms.h"
#include "Utilities.h"
#include <boost/math/special_functions/chebyshev.hpp>
#include <utility>
#include <unordered_map>
#include "ChebyshevAlgorithms.h"
#include <Eigen/Core>



int main() {
	int x = 3, y = 7, z = 7;
	double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 1, 1, 1.0, 1.0, 1);
	// The first degree polynomials for each dimension
	int N = 3;
	int Q = 6;
	int S = 7;
	
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

/*
int main() {
	//start();
}
*/


/*
int main() {
	int x = 3, y = 7, z = 7; // Dimensions of matrix.
	double compressionRatio = .8; // Compression ratio. 0.9 means keeping 90% of the orignial data

	// Define original matrix
	double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 1, 1, 1.0, 1.0, 1);

	// Compressing the data
	int byteStreamSize;
	unsigned char* byteStream = FFTAlgorithms::compressMatrix(originalMatrix, x, y, z, compressionRatio, byteStreamSize);

	// Decompressing the received data
	double* decompressedMatrix = FFTAlgorithms::decompressMatrix(byteStream, byteStreamSize);

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

/*
typedef std::vector<std::complex<double>> ComplexVec;

int main() {
	// Initialize your parameters...
	int x = 3;
	int y = 7;
	int z = 7;
	int size;
	double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);
	double compressionRatio = 0.5;

	ComplexVec complexMatrix(x * y * z);
	double* in = new double[x * y * z];
	fftw_complex* out = reinterpret_cast<fftw_complex*>(complexMatrix.data());

	// Create the plans
	fftw_plan p = fftw_plan_dft_r2c_3d(x, y, z, originalMatrix, out, FFTW_ESTIMATE);
	fftw_plan q = fftw_plan_dft_c2r_3d(x, y, z, out, originalMatrix, FFTW_ESTIMATE);

	// Use the methods
	unsigned char* compressed = FFTAlgorithms::compressMatrix(p, originalMatrix, x, y, z, compressionRatio, size);
	double* decompressed = FFTAlgorithms::decompressMatrix(q, compressed, size);

	// Cleanup
	fftw_destroy_plan(p);
	fftw_destroy_plan(q);

	// Delete the dynamically allocated arrays
	delete[] compressed;
	delete[] decompressed;

	return 0;
}
*/



/*
int main() {
	int x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
	double rate = 20; // Compression ratio. 0.9 means keeping 90% of the orignial data

	// Define original matrix
	double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);

	// Compressing the data
	int byteStreamSize;
	unsigned char* byteStream = ZFPAlgorithms::compressMatrix(originalMatrix, x, y, z, rate, byteStreamSize);

	// Decompressing the received data
	double* decompressedMatrix = ZFPAlgorithms::decompressMatrix(byteStream, byteStreamSize);

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
