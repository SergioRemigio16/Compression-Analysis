#include "TimingExperiment.h"
#include "FFTAlgorithms.h"
#include "Utilities.h"
#include <boost/math/special_functions/chebyshev.hpp>
#include <utility>
#include <unordered_map>
#include "ChebyshevAlgorithms.h"
#include <Eigen/Core>


Eigen::MatrixXd kroneckerProduct(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
	Eigen::MatrixXd C(A.rows() * B.rows(), A.cols() * B.cols());
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < A.cols(); ++j) {
			C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
		}
	}
	return C;
}

/**
 * @brief
 *
 * N, Q, and S must be less than or equal to x, y, and z respectively
 *
 * @param originalMatrix A three dimensional array represented in one dimension through row majoring order
 * @param x The number of elements in the x dimension
 * @param y The number of elements in the y dimension
 * @param z The number of elements in the z dimension
 * @param N The number of the coefficients in the x dimension
 * @param Q The number of the coefficients in the y dimension
 * @param S The number of the coefficients in the z dimension
 *
 * @return A pointer to the byte stream containing the compressed data.
 */
double* testing(double* originalMatrix, int x, int y, int z, int N, int Q, int S) {
	if (N > x) {
		N = x;
	}
	if (Q > y) {
		Q = y;
	}
	if (S > z) {
		S = z;
	}
	// Normalize the data inside originalMatrix
	double minVal = *std::min_element(originalMatrix, originalMatrix + x * y * z);
	double maxVal = *std::max_element(originalMatrix, originalMatrix + x * y * z);
	double range = maxVal - minVal;



	// Create normalizedMatrix as Eigen::VectorXd
	Eigen::VectorXd normalizedMatrix(x * y * z);
	for (int i = 0; i < x; ++i) {
		for (int j = 0; j < y; ++j) {
			for (int k = 0; k < z; ++k) {
				normalizedMatrix(i * y * z + j * z + k) = ((originalMatrix[i * y * z + j * z + k] - minVal) / range) * 2.0 - 1.0;
			}
		}
	}

	// Example how to get a chebyshev polynomial 
	double result = ChebyshevAlgorithms::chebyshevT(4, 0.3);

	// Calculate A, B, and C, where A = x x N, B = y x Q, and C = z x S

	Eigen::MatrixXd A(x, N);
	Eigen::MatrixXd B(y, Q);
	Eigen::MatrixXd C(z, S);


	// Populate A using Chebyshev polynomials
	for (int i = 0; i < x; ++i) {
		double normalized_x = static_cast<double>(i) / (x - 1);
		for (int n = 0; n < N; ++n) {
			A(i, n) = ChebyshevAlgorithms::chebyshevT(n, normalized_x);
		}
	}

	// Populate B using Chebyshev polynomials
	for (int j = 0; j < y; ++j) {
		double normalized_y = static_cast<double>(j) / (y - 1);
		for (int q = 0; q < Q; ++q) {
			B(j, q) = ChebyshevAlgorithms::chebyshevT(q, normalized_y);
		}
	}

	// Populate C using Chebyshev polynomials
	for (int k = 0; k < z; ++k) {
		double normalized_z = static_cast<double>(k) / (z - 1);
		for (int s = 0; s < S; ++s) {
			C(k, s) = ChebyshevAlgorithms::chebyshevT(s, normalized_z);
		}
	}


	// Perform SVD on A, B, and C
    Eigen::JacobiSVD<Eigen::MatrixXd> svdA(A, Eigen::ComputeFullU);
    Eigen::JacobiSVD<Eigen::MatrixXd> svdB(B, Eigen::ComputeFullU);
    Eigen::JacobiSVD<Eigen::MatrixXd> svdC(C, Eigen::ComputeFullU);

    // Truncate U matrices
    Eigen::MatrixXd UA_truncated = svdA.matrixU().leftCols(N);
    Eigen::MatrixXd UB_truncated = svdB.matrixU().leftCols(Q);
    Eigen::MatrixXd UC_truncated = svdC.matrixU().leftCols(S);

    // Calculate the Kronecker product with truncated U matrices
    Eigen::MatrixXd UA_UB_kron = kroneckerProduct(UA_truncated, UB_truncated);
    Eigen::MatrixXd full_kron_product_truncated = kroneckerProduct(UA_UB_kron, UC_truncated);

    // Calculate the coefficients using the truncated Kronecker product
    Eigen::VectorXd coefficients_truncated = full_kron_product_truncated.transpose() * normalizedMatrix;

    // Reconstruct the matrix using the coefficients and the truncated Kronecker product
    Eigen::VectorXd reconstructedMatrix = full_kron_product_truncated * coefficients_truncated;

	// Convert Eigen::VectorXd back to double*
	double* reconstructedArray = new double[x * y * z];
	for (int i = 0; i < x * y * z; ++i) {
		reconstructedArray[i] = reconstructedMatrix(i);
	}

	// Unnormalize the data
	for (int i = 0; i < x * y * z; ++i) {
		reconstructedArray[i] = ((reconstructedArray[i] + 1.0) / 2.0) * range + minVal;
	}


	Utilities::printComparison(originalMatrix, reconstructedArray, x, y, z);
	Utilities::printError(originalMatrix, reconstructedArray, x, y, z);

	return reconstructedArray;

}

/*
void decompress_matrix(double* compressedCoeffs, double* decompressedMatrix, int x, int y, int z, int degree) {

}
*/

int main() {
	int x = 21, y = 21, z = 1;
	double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 1, 1, 1.0, 1.0, 1);
	// The first degree polynomials for each dimension
	int N = x - 3;
	int Q = y - 3;
	int S = z;

	double* compressedCoeffs = testing(originalMatrix, x, y, z, N, Q, S);
	//double* decompressedMatrix = new double[x * y * z];
	//decompress_matrix(compressedCoeffs, decompressedMatrix, x, y, z);

	//Utilities::printComparison(originalMatrix, decompressedMatrix, x, y, z);
	//Utilities::printError(originalMatrix, decompressedMatrix, x, y, z);

	delete[] originalMatrix;
	delete[] compressedCoeffs;
	//delete[] decompressedMatrix;

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
