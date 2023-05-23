
#include <iostream>
#include <zfp.h>

int main() {
    // Define dimensions of the 3D array
    uint nx = 10;
    uint ny = 10;
    uint nz = 10;

    // Allocate memory for the 3D array
    float* data = new float[nx * ny * nz];

    // Fill the 3D array with some example data (e.g., a simple pattern)
    for (uint z = 0; z < nz; ++z) {
        for (uint y = 0; y < ny; ++y) {
            for (uint x = 0; x < nx; ++x) {
                data[z * ny * nx + y * nx + x] = x * y * z;
            }
        }
    }

    // Create a zfp field for the 3D array
    zfp_field* field = zfp_field_3d(data, zfp_type_float, nx, ny, nz);

    // Create a zfp compression stream with a specific compression rate
    zfp_stream* zfp = zfp_stream_open(NULL);
    double rate = 8; // Compression rate (bits/value)
    zfp_stream_set_rate(zfp, rate, zfp_type_float, 3, 0);

    // Allocate memory for the compressed data
    size_t bufsize = zfp_stream_maximum_size(zfp, field);
    uchar* buffer = new uchar[bufsize];

    // Create a bit stream and associate it with the compression stream
    bitstream* stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

    // Compress the data
    size_t compressed_size = zfp_compress(zfp, field);
    if (!compressed_size) {
        std::cerr << "Compression failed\n";
        return 1;
    }

    // Cleanup
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);

    std::cout << "Original size: " << nx * ny * nz * sizeof(float) << " bytes" << std::endl;
    std::cout << "Compressed size: " << compressed_size << " bytes" << std::endl;

    delete[] data;
    delete[] buffer;

    return 0;
}
