#include <iostream>
#include <chrono>
#include <zfp.h>
#include <vector>
#include <numeric>
#include <cmath>
#include "Utilities.h"

// Compress the given matrix using the ZFP library.
// Input: A zfp_stream, a zfp_field, the buffer size, and a vector of data.
// Output: The zfp_field of the compressed matrix.
zfp_field* compressMatrix(zfp_stream* zfp, zfp_field* field, size_t bufsize, std::vector<double>& data) {
    void* buffer = malloc(bufsize);
    bitstream* stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);
    zfp_compress(zfp, field);
    zfp_stream_rewind(zfp);

    // link the decompressed field with the data vector
    zfp_field* dec_field = zfp_field_3d(data.data(), zfp_type_double, field->nx, field->ny, field->nz);

    stream_close(stream);
    free(buffer);

    return dec_field;
}

// Decompress the given matrix using the ZFP library.
// Input: A zfp_stream, a zfp_field, the buffer size, and a zfp_field of compressed data.
// Output: The zfp_field of the decompressed matrix.
zfp_field* decompressMatrix(zfp_stream* zfp, zfp_field* field, size_t bufsize, zfp_field* dec_field) {
    void* buffer = malloc(bufsize);
    bitstream* stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);
    zfp_compress(zfp, field);
    zfp_stream_rewind(zfp);

    // link the decompressed field with the data vector
    zfp_decompress(zfp, dec_field);

    stream_close(stream);
    zfp_field_free(dec_field); // free the decompressed field
    free(buffer);

    return dec_field;
}
