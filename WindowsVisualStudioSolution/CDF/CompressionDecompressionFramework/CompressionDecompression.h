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
	uint rate;
	CompressionResult() : x(0), y(0), z(0), bufsize(0), rate(0.0) {}
};

CompressionResult compressMatrix(std::vector<double>& originalData, int x, int y, int z, uint rate);
std::vector<double> decompressMatrix(CompressionResult& result);

#endif // _COMPRESSIONDECOMPRESSION_H_
