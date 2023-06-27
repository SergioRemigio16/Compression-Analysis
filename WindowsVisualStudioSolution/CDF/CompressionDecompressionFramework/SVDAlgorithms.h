
/**
 * @file SVDAlgorithms.h
 * @brief This header file provides SVD-based data compression and decompression algorithms.
 *
 * No padding. 
 * This file contains methods and types to perform SVD-based compression and decompression,
 * byte size reporting for the compressed data, and calculation of decompressed data bytes.
 */

#ifndef _SVDALGORITHMS_H_
#define _SVDALGORITHMS_H_

// Disables Eigen's memory alignment which could lead to extra memory padding.
#define EIGEN_DONT_ALIGN
#include <Eigen/Dense> // Include this header for matrix and vector operations
#include <iostream>

// This namespace encapsulates all methods and types for SVD-based compression.
namespace SVDAlgorithms {
	
	/**
	 * @struct SVDResult
	 * @brief This structure stores the result of SVD.
	 *
	 * The result includes left singular vectors (U), right singular vectors (V),
	 * singular values (S), and the dimensions of the original matrix (x, y, z).
	 */
	struct SVDResult {
		Eigen::MatrixXd U; // Left singular vectors
		Eigen::MatrixXd V; // Right singular vectors
		Eigen::VectorXd S; // Singular values
		int x; 
		int y; 
		int z; 
	};

	/**
	 * @brief Compresses the original matrix using SVD and returns the compressed version.
	 *
	 * @param originalMatrix Pointer to the original matrix.
	 * @param x X-dimension of the original matrix.
	 * @param y Y-dimension of the original matrix.
	 * @param z Z-dimension of the original matrix.
	 * @param k The desired number of singular values to retain in the compressed matrix.
	 * @return The compressed matrix as SVDResult.
	 */
	SVDResult compressMatrix(double*& originalMatrix, const int x, const int y, const int z, const int k);

	/**
	 * @brief Decompresses the compressed data back into the original form.
	 *
	 * @param compressedData The compressed data from compressMatrix().
	 * @param x X-dimension of the original matrix.
	 * @param y Y-dimension of the original matrix.
	 * @param z Z-dimension of the original matrix.
	 * @return A pointer to the decompressed data.
	 */
	double* decompressMatrix(const SVDResult& compressedData, const int x, const int y, const int z);

	/**
	 * @brief Prints the memory sizes of U, V, S and the total size.
	 *
	 * @param compressedData The compressed data from compressMatrix().
	 */
	void printByteSizeReport(const SVDResult& compressedData);

	/**
	 * @brief Calculates the byte size of the compressed data.
	 *
	 * @param compressedData The compressed data from compressMatrix().
	 */
	size_t calculateCompressedDataBytes(const SVDAlgorithms::SVDResult& compressedData);
}

#endif // _SVDALGORITHMS_H_
