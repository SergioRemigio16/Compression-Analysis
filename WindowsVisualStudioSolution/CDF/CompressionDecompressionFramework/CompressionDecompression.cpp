#include <iostream>
#include <chrono>
#include <zfp.h>
#include <vector>
#include <numeric>
#include <cmath>
#include "Utilities.h"
#include "CompressionDecompression.h"

CompressionResult compressMatrix(std::vector<double>& originalData, int x, int y, int z, uint rate) {
	zfp_field* field = zfp_field_3d(originalData.data(), zfp_type_double, x, y, z);
	zfp_stream* zfp = zfp_stream_open(NULL);

	// Set the compression rate instead of accuracy
	zfp_stream_set_rate(zfp, rate, zfp_type_double, 3, 0);
	size_t bufsize = zfp_stream_maximum_size(zfp, field);
	std::vector<unsigned char> buffer(bufsize);
	bitstream* stream = stream_open(buffer.data(), bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	zfp_compress(zfp, field);
	stream_close(stream);
	zfp_stream_close(zfp);
	zfp_field_free(field);

	CompressionResult result;
	result.x = x;
	result.y = y;
	result.z = z;
	result.buffer = std::move(buffer);
	result.bufsize = bufsize;
	result.rate = rate;

	return result;
}

std::vector<double> decompressMatrix(CompressionResult& result) {
	zfp_stream* zfp = zfp_stream_open(NULL);

	// Set the decompression rate instead of accuracy
	zfp_stream_set_rate(zfp, result.rate, zfp_type_double, 3, 0);
	bitstream* stream = stream_open(result.buffer.data(), result.bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	std::vector<double> decompressedData(result.x * result.y * result.z);
	zfp_field* dec_field = zfp_field_3d(decompressedData.data(), zfp_type_double, result.x, result.y, result.z);
	zfp_decompress(zfp, dec_field);
	zfp_field_free(dec_field);
	stream_close(stream);
	zfp_stream_close(zfp);
	return decompressedData;
}
