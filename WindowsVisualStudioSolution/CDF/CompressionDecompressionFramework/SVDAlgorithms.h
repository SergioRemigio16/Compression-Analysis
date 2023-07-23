
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

	unsigned char* compressMatrix(double*& originalMatrix, const int x, const int y, const int z, const int k, int& size);

	double* decompressMatrix(unsigned char*& compressedData, int size);

}
#endif // _SVDALGORITHMS_H_
