#ifndef _SVDALGORITHMS_H_
#define _SVDALGORITHMS_H_

// Disables Eigen's memory alignment which could lead to extra memory padding.
#define EIGEN_DONT_ALIGN
#include <Eigen/Dense> // Include this header for matrix and vector operations
#include <iostream>

// This namespace encapsulates all methods and types for SVD-based compression.
namespace SVDAlgorithms {
	// This struct stores the result of SVD.
	// U, V are matrices, S is a vector, and x, y, z are the dimensions of the original matrix.
	struct SVDResult {
		Eigen::MatrixXd U; // Left singular vectors
		Eigen::MatrixXd V; // Right singular vectors
		Eigen::VectorXd S; // Singular values
		int x; 
		int y; 
		int z; 
	};

	// This function compresses the original matrix using SVD and returns the compressed version.
	// originalMatrix: pointer to the original matrix
	// x, y, z: dimensions of the original matrix
	// k: the desired number of singular values to retain in the compressed matrix
	//    A smaller k will result in higher compression but less accuracy in the reconstructed data.
	//    k should be less than or equal to min(x*y, z)
	SVDResult compressMatrix(double*& originalMatrix, const int x, const int y, const int z, const int k);

	// This function decompresses the compressed data back into the original form.
	// compressedData: the compressed data from compressMatrix()
	// x, y, z: dimensions of the original matrix
	double* decompressMatrix(const SVDResult& compressedData, const int x, const int y, const int z);

	// This function prints the memory sizes of U, V, S and the total size.
	// U, V: singular vectors
	// S: singular values
	void printByteSize(const SVDResult& compressedData);

	// This function reports the byte size of the compressed data.
	void printByteSizeReport();
}

#endif // _SVDALGORITHMS_H_
