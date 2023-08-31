/*
  -----------------------------------------------------------------------
  Implementation of data compression algorithm based on the following paper:

  Author: Peter M. Williams
  Title: Image Compression for Neural Networks using Chebyshev Polynomials
  Journal: Artificial Neural Networks
  Publisher: North-Holland
  Year: 1992
  Pages: 1139-1142
  ISBN: 9780444894885
  DOI: https://doi.org/10.1016/B978-0-444-89488-5.50065-8
  Link: https://www.sciencedirect.com/science/article/pii/B9780444894885500658
  
  Abstract: This paper proposes an algorithm for data compression of multi-dimensional images, 
  prior to input into neural network feature detectors, using Chebyshev approximation. 
  Illustrations are given of the fidelity obtainable for geophysical data.
  -----------------------------------------------------------------------
*/

#ifndef _CHEBYSHEVALGORITHMS_H_
#define _CHEBYSHEVALGORITHMS_H_

// Disables Eigen's memory alignment which could lead to extra memory padding.
#define EIGEN_DONT_ALIGN
#include <Eigen/Dense> // Include this header for matrix and vector operations
#include <iostream>
#include <utility>
#include <unordered_map>
#include <boost/math/special_functions/chebyshev.hpp>





// This namespace encapsulates all methods and types for SVD-based compression.
namespace ChebyshevAlgorithms {

	double chebyshevT(int n, double x);

	/*
	unsigned char* compressMatrix(double*& originalMatrix, const int x, const int y, const int z, const int k, int& size);
	double* decompressMatrix(unsigned char*& buffer, int bufferSize);
	*/

}
#endif // _CHEBYSHEVALGORITHMS_H_

