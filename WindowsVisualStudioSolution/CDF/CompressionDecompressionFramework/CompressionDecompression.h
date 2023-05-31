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
};

CompressionResult compressData(std::vector<double>& originalData, int x, int y, int z, double precision);
std::vector<double> decompressData(CompressionResult& result);

#endif // _COMPRESSIONDECOMPRESSION_H_
