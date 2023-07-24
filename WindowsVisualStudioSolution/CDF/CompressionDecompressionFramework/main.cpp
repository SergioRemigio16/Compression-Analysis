#include "TimingExperiment.h"

#include <zfp.h>
#include <vector>
#include <complex>
#include <iostream>
#include <chrono>
#include <fftw3.h>
#include <numeric>
#include <cmath>

unsigned char* compressData(double* originalData, int x, int y, int z, double rate, int& size) {
	// Initialize a 3D array with original data, using ZFP's special 'field' structure. 
	// The field has the type of double and the dimensions are given by x, y, and z.
	zfp_field* field = zfp_field_3d(originalData, zfp_type_double, x, y, z);
	// Open a new ZFP stream. A ZFP stream is responsible for compressing and decompressing data.
	zfp_stream* zfp = zfp_stream_open(NULL);
	// Set the compression rate for the ZFP stream. The type of the data is double, the dimensionality is 3, and '0' indicates we're not using a user-specified precision.
	zfp_stream_set_rate(zfp, rate, zfp_type_double, 3, 0);

	// Determine the maximum buffer size necessary for this ZFP stream given the input field.
	int bufsize = zfp_stream_maximum_size(zfp, field);
	size = bufsize + 3 * sizeof(int) + sizeof(double); // metadata x,y,z,rate
	// Create a buffer with enough capacity to store the compressed data.
	unsigned char* buffer = new unsigned char[size];
	unsigned char* ptr = buffer;
	memcpy(ptr, &x, sizeof(int));
	ptr += sizeof(int);
	memcpy(ptr, &y, sizeof(int));
	ptr += sizeof(int);
	memcpy(ptr, &z, sizeof(int));
	ptr += sizeof(int);
	memcpy(ptr, &rate, sizeof(double));
	ptr += sizeof(double);
	// Create a bitstream from the buffer to store compressed data.
	bitstream* stream = stream_open(ptr, bufsize);
	// Associate the bitstream with the ZFP stream, so compressed data will go into our buffer.
	zfp_stream_set_bit_stream(zfp, stream);
	// Compress the data. The results will be stored in the buffer we've created.
	zfp_compress(zfp, field);

	// Close the bitstream. All compressed data should now reside in our buffer.
	stream_close(stream);
	// Close the ZFP stream since we're done with compression.
	zfp_stream_close(zfp);
	// Release the memory allocated for the field since we're done with it.
	zfp_field_free(field);

	return buffer;
}


double* decompressData(unsigned char* buffer, int bufferSize) {
	// Deserialize metadata
	int x, y, z;
	memcpy(&x, buffer, sizeof(int));
	buffer += sizeof(int);
	memcpy(&y, buffer, sizeof(int));
	buffer += sizeof(int);
	memcpy(&z, buffer, sizeof(int));
	buffer += sizeof(int);
	double rate;
	memcpy(&rate, buffer, sizeof(double));
	buffer += sizeof(double);

	zfp_stream* zfp = zfp_stream_open(NULL);

	// Set the decompression rate instead of accuracy
	zfp_stream_set_rate(zfp, rate, zfp_type_double, 3, 0);
	bitstream* stream = stream_open(buffer, bufferSize - 3 * sizeof(int) - sizeof(double));
	zfp_stream_set_bit_stream(zfp, stream);
	double* decompressedData = new double[x * y * z];
	zfp_field* dec_field = zfp_field_3d(decompressedData, zfp_type_double, x, y, z);
	zfp_decompress(zfp, dec_field);
	zfp_field_free(dec_field);
	stream_close(stream);
	zfp_stream_close(zfp);
	return decompressedData;
}


int main() {
    int x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    double rate = 8; // Compression ratio. 0.9 means keeping 90% of the orignial data

    // Define original matrix
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);

    // Compressing the data
    int byteStreamSize;
    unsigned char* byteStream = compressData(originalMatrix, x, y, z, rate, byteStreamSize);

    // Decompressing the received data
    double* decompressedMatrix = decompressData(byteStream, byteStreamSize);

    // Printing comparison of original and decompressed data
    Utilities::printComparison(originalMatrix, decompressedMatrix, x, y, z);

    // Printing size of original and serialized compressed data
    std::cout << "Size of original matrix (bytes): " << (x * y * z) * sizeof(double) << "\n";
    std::cout << "Size of serialized compressed data (bytes): " << byteStreamSize << "\n";

    // Printing various types of error between original and decompressed data
    Utilities::printError(originalMatrix, decompressedMatrix, x, y, z);

    // Freeing the memory
    delete[] originalMatrix;
    delete[] decompressedMatrix;
    delete[] byteStream;

    return 0;
}
