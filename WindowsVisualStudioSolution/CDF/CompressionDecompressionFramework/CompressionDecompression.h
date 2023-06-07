#ifndef _COMPRESSIONDECOMPRESSION_H_
#define _COMPRESSIONDECOMPRESSION_H_

#include <zfp.h>
#include <vector>
#include <complex>

struct CompressionResult {
	int x;
	int y;
	int z;
	/// <summary>
	/// Holds the compressed data
	/// </summary>
	std::vector<unsigned char> buffer;
	size_t bufsize;
	double precision;
	uint rate;
	CompressionResult() : x(0), y(0), z(0), bufsize(0), rate(0), precision(0.0) {}
};

struct CompressionResultFFT {
	int x, y, z; // Original dimensions
	int compressSize; // Number of strongest frequencies we keep
	std::vector<int> indices; // Indices of the kept frequencies
	std::vector<std::complex<double>> frequencies; // Complex values of the kept frequencies
};


CompressionResult compressMatrixFixedRate(std::vector<double>& originalData, int x, int y, int z, uint rate);
CompressionResult compressMatrixAccuracy(std::vector<double>& originalData, int x, int y, int z, double precision);
std::vector<double> decompressMatrixFixedRate(CompressionResult& result);
std::vector<double> decompressMatrixAccuracy(CompressionResult& result);


CompressionResultFFT compressMatrixFFT(const std::vector<double>& originalMatrix, int x, int y, int z, int compressSize);
std::vector<double> decompressMatrixFFT(const CompressionResultFFT& compressionResult);

#endif // _COMPRESSIONDECOMPRESSION_H_
