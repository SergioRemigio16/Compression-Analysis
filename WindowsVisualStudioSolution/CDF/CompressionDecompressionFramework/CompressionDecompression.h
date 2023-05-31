#ifndef _COMPRESSIONDECOMPRESSION_H_
#define _COMPRESSIONDECOMPRESSION_H_

#include <zfp.h>
#include <vector>

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

CompressionResult compressMatrixFixedRate(std::vector<double>& originalData, int x, int y, int z, uint rate);
CompressionResult compressMatrixAccuracy(std::vector<double>& originalData, int x, int y, int z, double precision);
std::vector<double> decompressMatrixFixedRate(CompressionResult& result);
std::vector<double> decompressMatrixAccuracy(CompressionResult& result);

#endif // _COMPRESSIONDECOMPRESSION_H_
