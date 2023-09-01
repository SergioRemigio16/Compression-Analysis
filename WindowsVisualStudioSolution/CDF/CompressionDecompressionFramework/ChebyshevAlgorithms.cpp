#include "ChebyshevAlgorithms.h"

// Type alias for cache key, where first is the degree n and second is the value x
using CacheKey = std::pair<int, double>;

// Custom hash function for CacheKey
struct CacheKeyHash {
	std::size_t operator()(const CacheKey& key) const {
		return std::hash<int>()(key.first) ^ std::hash<double>()(key.second);
	}
};

thread_local std::unordered_map<CacheKey, double, CacheKeyHash> cache;


// Type alias for cache key
struct SVDKey {
	int x;
	int y;
	int z;

	bool operator==(const SVDKey& other) const {
		return x == other.x && y == other.y && z == other.z;
	}
};

// Custom hash function for SVDKey
struct SVDKeyHash {
	std::size_t operator()(const SVDKey& key) const {
		return std::hash<int>()(key.x) ^ std::hash<int>()(key.y) ^ std::hash<int>()(key.z);
	}
};

// SVD struct
struct SVDs {
	Eigen::JacobiSVD<Eigen::MatrixXd> svdA;
	Eigen::JacobiSVD<Eigen::MatrixXd> svdB;
	Eigen::JacobiSVD<Eigen::MatrixXd> svdC;
};

// Thread-local cache for SVD objects
thread_local std::unordered_map<SVDKey, SVDs, SVDKeyHash> svdCache;

Eigen::MatrixXd kroneckerProduct(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
	Eigen::MatrixXd C(A.rows() * B.rows(), A.cols() * B.cols());
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < A.cols(); ++j) {
			C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
		}
	}
	return C;
}

double ChebyshevAlgorithms::chebyshevT(int n, double x) {
	CacheKey key = {n, x};
	
	// Check if value is in cache
	auto iter = cache.find(key);
	if (iter != cache.end()) {
		return iter->second;
	}

	// Calculate Chebyshev polynomial
	double result;
	if (n == 0) {
		result = 1.0;
	} else if (n == 1) {
		result = x;
	} else {
		double T_prev = 1.0;
		double T_curr = x;
		for (int i = 2; i <= n; ++i) {
			double T_next = 2 * x * T_curr - T_prev;
			T_prev = T_curr;
			T_curr = T_next;
		}
		result = T_curr;
	}

	// Store result in cache and return
	cache[key] = result;
	return result;
}


/**
 * @brief
 *
 * N, Q, and S must be less than or equal to x, y, and z respectively. 
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
unsigned char* ChebyshevAlgorithms::compressMatrix(double*& originalMatrix, int x, int y, int z, int N, int Q, int S, int& bufferSize) {
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

	Eigen::JacobiSVD<Eigen::MatrixXd> svdA;
	Eigen::JacobiSVD<Eigen::MatrixXd> svdB;
	Eigen::JacobiSVD<Eigen::MatrixXd> svdC;

	// Define your key based on whatever criteria makes matrices A, B, C identical or close enough
	SVDKey key{ x, y, z };

	// Check if the SVDs for this key are already computed and stored in the cache
	auto it = svdCache.find(key);
	if (it != svdCache.end()) {
		// Cache hit: SVD results are already computed
		SVDs cachedSVDs = it->second;
		svdA = cachedSVDs.svdA;
		svdB = cachedSVDs.svdB;
		svdC = cachedSVDs.svdC;
	} else {
		// Cache miss: Need to perform SVD and store the results in the cache
		svdA = Eigen::JacobiSVD<Eigen::MatrixXd>(A, Eigen::ComputeFullU);
		svdB = Eigen::JacobiSVD<Eigen::MatrixXd>(B, Eigen::ComputeFullU);
		svdC = Eigen::JacobiSVD<Eigen::MatrixXd>(C, Eigen::ComputeFullU);

		SVDs newSVDs{ svdA, svdB, svdC };
		svdCache[key] = newSVDs;
	}

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

	return NULL;

}


double* ChebyshevAlgorithms::decompressMatrix(unsigned char*& buffer, int bufferSize) {
	return NULL;
}



