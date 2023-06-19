#include <iostream>
#include <chrono>
#include <complex>
#include <zfp.h>
#include <fftw3.h>
#include <vector>
#include <numeric>
#include <cmath>
#include "Utilities.h"
#include "CompressionDecompression.h"

CompressionResult compressMatrixFixedRate(double*& originalData, int x, int y, int z, uint rate) {
	zfp_field* field = zfp_field_3d(originalData, zfp_type_double, x, y, z);
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

double* decompressMatrixFixedRate(CompressionResult& result) {
	zfp_stream* zfp = zfp_stream_open(NULL);

	// Set the decompression rate instead of accuracy
	zfp_stream_set_rate(zfp, result.rate, zfp_type_double, 3, 0);
	bitstream* stream = stream_open(result.buffer.data(), result.bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	double* decompressedData = new double[result.x * result.y * result.z];
	zfp_field* dec_field = zfp_field_3d(decompressedData, zfp_type_double, result.x, result.y, result.z);
	zfp_decompress(zfp, dec_field);
	zfp_field_free(dec_field);
	stream_close(stream);
	zfp_stream_close(zfp);
	return decompressedData;
}



CompressionResult compressMatrixAccuracy(double*& originalData, int x, int y, int z, double precision) {
	// Creates a new ZFP field in 3D that will be compressed. The field takes the raw data from 'originalData',
	// specifies that the data type is double, and provides the dimensions of the 3D field (x, y, z).
	zfp_field* field = zfp_field_3d(originalData, zfp_type_double, x, y, z);
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
/*
	std::cout << "minbits ";
	std::cout << zfp->minbits << std::endl;
	std::cout << "maxbits ";
	std::cout << zfp->maxbits << std::endl;
	std::cout << "minexp ";
	std::cout << zfp->minexp << std::endl;
	std::cout << "maxprec ";
	std::cout << zfp->maxprec << std::endl;
	std::cout << "stream ";
	std::cout << zfp->stream << std::endl;
	*/


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



double* decompressMatrixAccuracy(CompressionResult& result) {
	zfp_stream* zfp = zfp_stream_open(NULL);
	zfp_stream_set_accuracy(zfp, result.precision);
	// Creates a bitstream that will be used to hold the compressed data.
	bitstream* stream = stream_open(result.buffer.data(), result.bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	// Where the decompressedData will be placed
	double* decompressedData = new double[result.x * result.y * result.z];
	// Creates a new ZFP field in 3D that will be decompressed.
	zfp_field* dec_field = zfp_field_3d(decompressedData, zfp_type_double, result.x, result.y, result.z);
	// Decompresses the data into the 'dec_field'.
	zfp_decompress(zfp, dec_field);

	// Free and close ZFP fields and streams.
	zfp_field_free(dec_field);
	stream_close(stream);
	zfp_stream_close(zfp);
	return decompressedData;
}



// Compress the input data using FFT and keep the COMPRESS_SIZE strongest frequencies.
CompressionResultFFT compressMatrixFFT(const double*& originalMatrix, int x, int y, int z, int compressSize) {
	CompressionResultFFT result;
	result.x = x;
	result.y = y;
	result.z = z;
	result.compressSize = compressSize;

	// Prepare the input data for FFTW.
	int size = x * y * z;
	fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
	for (int i = 0; i < size; ++i) {
		in[i][0] = originalMatrix[i]; // real part
		in[i][1] = 0; // imaginary part
	}

	// Execute FFT using FFTW.
	fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
	fftw_plan p = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	// Select the COMPRESS_SIZE strongest frequencies and store their indices and values.
	for (int i = 0; i < compressSize; ++i) {
		result.indices.push_back(i);
		result.frequencies.push_back(std::complex<double>(out[i][0], out[i][1]));
	}

	// Clean up.
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);

	return result;
}

// Decompress the data using inverse FFT.
double* decompressMatrixFFT(const CompressionResultFFT& compressionResult) {
	int size = compressionResult.x * compressionResult.y * compressionResult.z;

	// Prepare the input data for inverse FFT.
	fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
	for (int i = 0; i < size; ++i) {
		in[i][0] = 0; // real part
		in[i][1] = 0; // imaginary part
	}
	for (int i = 0; i < compressionResult.compressSize; ++i) {
		in[compressionResult.indices[i]][0] = compressionResult.frequencies[i].real();
		in[compressionResult.indices[i]][1] = compressionResult.frequencies[i].imag();
	}

	// Execute inverse FFT.
	fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
	fftw_plan p = fftw_plan_dft_1d(size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	// Copy the decompressed data to the output vector.
	double* decompressedMatrix = new double[size];
	for (int i = 0; i < size; ++i) {
		decompressedMatrix[i] = out[i][0]; // we only need the real part
	}

	// Clean up.
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);

	return decompressedMatrix;
}
