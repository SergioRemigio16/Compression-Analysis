
#include <iostream>
#include <chrono>
#include <zfp.h>
#include <vector>
#include <numeric>
#include <cmath>
#include "Utilities.h"
#include "CompressionDecompression.h"

CompressionResult compressData(std::vector<double>& originalData, int x, int y, int z, double precision) {
	// Creates a new ZFP field in 3D that will be compressed. The field takes the raw data from 'originalData',
	// specifies that the data type is double, and provides the dimensions of the 3D field (x, y, z).
	zfp_field* field = zfp_field_3d(originalData.data(), zfp_type_double, x, y, z); 
	// Opens a new ZFP stream for compression. The argument is NULL because we don't want to use a pre-existing stream.
	zfp_stream* zfp = zfp_stream_open(NULL);
	// Sets the compression precision of the ZFP stream. In this case, we want lossless compression,
	// so the precision is set to 0.0 (meaning no data will be lost in compression).
	zfp_stream_set_accuracy(zfp, precision);
	// Gets the maximum buffer size required for the compression.
	size_t bufsize = zfp_stream_maximum_size(zfp, field);
	// Create buffer that will hold the compressed data
	std::vector<unsigned char> buffer(bufsize);
	// Creates a bitstream that will be used to hold the compressed data.
	bitstream* stream = stream_open(buffer.data(), bufsize);
	// Associates the bitstream with the ZFP stream, meaning that when we compress our field,
	// the data will be written into this bitstream.
	zfp_stream_set_bit_stream(zfp, stream);
	// Compresses the field with the ZFP stream.
	zfp_compress(zfp, field);

	// Free and close ZFP fields and streams.
	stream_close(stream);
	zfp_stream_close(zfp);
	zfp_field_free(field);

	CompressionResult result;
	result.x = x;
	result.y = y;
	result.z = z;
	result.buffer = std::move(buffer);  // Avoid expensive copy with std::move (transfer ownership)
	result.bufsize = bufsize;
	result.precision = precision;

	return result;
}

std::vector<double> decompressData(CompressionResult& result) {
	std::vector<double> decompressedData(result.x * result.y * result.z);
	zfp_stream* zfp = zfp_stream_open(NULL);
	zfp_stream_set_accuracy(zfp, result.precision);
	// Creates a bitstream that will be used to hold the compressed data.
	bitstream* stream = stream_open(result.buffer.data(), result.bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	// Creates a new ZFP field in 3D that will be decompressed.
	zfp_field* dec_field = zfp_field_3d(decompressedData.data(), zfp_type_double, result.x, result.y, result.z);
	// Decompresses the data into the 'dec_field'.
	zfp_decompress(zfp, dec_field);

	// Free and close ZFP fields and streams.
	zfp_field_free(dec_field);
	stream_close(stream);
	zfp_stream_close(zfp);
	return decompressedData;
}


