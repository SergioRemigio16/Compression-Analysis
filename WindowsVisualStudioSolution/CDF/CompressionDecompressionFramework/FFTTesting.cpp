#include <iostream>
#include <fftw3.h>
#include <vector>
#include <fstream>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Constants for data size and compression level
const int DATA_SIZE = 1000;
const int COMPRESS_SIZE = 100; // We'll only keep the strongest 100 frequencies

struct Bin {
    int index;
    double re, im;
};

inline int runFFTTest() {
    // Create input data (for now, just a simple sin wave)
    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * DATA_SIZE);
    for (int i = 0; i < DATA_SIZE; ++i) {
        in[i][0] = sin(2 * M_PI * i / DATA_SIZE); // real part
        in[i][1] = 0; // imaginary part
    }

    // Prepare FFTW
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * DATA_SIZE);
    fftw_plan p = fftw_plan_dft_1d(DATA_SIZE, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Measure compression time
    auto start = std::chrono::high_resolution_clock::now();

    // Execute FFT
    fftw_execute(p);

    // Now the 'out' array contains frequency components. Each component has real and imaginary parts
    // which represent the amplitude and phase of that frequency.

    // Compress the data: select the COMPRESS_SIZE strongest frequency components.
    // For simplicity, we'll just select the first COMPRESS_SIZE components.
    // In a real application, you should select the strongest components based on their amplitude.
    std::vector<Bin> compressed(COMPRESS_SIZE);
    for (int i = 0; i < COMPRESS_SIZE; ++i) {
        compressed[i].index = i;
        compressed[i].re = out[i][0];
        compressed[i].im = out[i][1];
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Compression time: " << duration.count() << " microseconds" << std::endl;

    // Write the compressed data to a file. For simplicity, we use raw binary format.
    std::ofstream compressed_file("compressed.dat", std::ios::binary);
    compressed_file.write(reinterpret_cast<char*>(compressed.data()), sizeof(Bin) * COMPRESS_SIZE);
    compressed_file.close();

    // Get compressed data size
    std::ifstream compressed_file_read("compressed.dat", std::ios::binary | std::ios::ate);
    std::streamsize compressed_size = compressed_file_read.tellg();
    compressed_file_read.close();

    // Measure decompression time
    start = std::chrono::high_resolution_clock::now();

    // At this point, you could send the compressed file to another machine or save it for later use.
    // To decompress and recover the data, you would first read the compressed data:

    std::vector<Bin> read_compressed(COMPRESS_SIZE);
    std::ifstream read_file("compressed.dat", std::ios::binary);
    read_file.read(reinterpret_cast<char*>(read_compressed.data()), sizeof(Bin) * COMPRESS_SIZE);
    read_file.close();

    // Then prepare the data for inverse FFT:

    fftw_complex* recovered_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * DATA_SIZE);
    for (int i = 0; i < DATA_SIZE; ++i) {
        // Initialize all bins to zero
        recovered_in[i][0] = 0;
        recovered_in[i][1] = 0;
    }
    for (const Bin& bin : read_compressed) {
        // Set the bins we have kept
        recovered_in[bin.index][0] = bin.re;
        recovered_in[bin.index][1] = bin.im;
    }

    // Prepare and execute inverse FFT

    fftw_complex* recovered_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * DATA_SIZE);
    fftw_plan p_inverse = fftw_plan_dft_1d(DATA_SIZE, recovered_in, recovered_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p_inverse);

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Decompression time: " << duration.count() << " microseconds" << std::endl;

    // Now, 'recovered_out' should contain data similar to the original input data,
    // except that it has been lossily compressed and then decompressed.

    // Get original data size
    std::streamsize original_size = sizeof(fftw_complex) * DATA_SIZE;
    std::cout << "Original data size: " << original_size << " bytes" << std::endl;
    std::cout << "Compressed data size: " << compressed_size << " bytes" << std::endl;

    // Don't forget to free the memory and destroy the plans
    fftw_destroy_plan(p);
    fftw_destroy_plan(p_inverse);
    fftw_free(in);
    fftw_free(out);
    fftw_free(recovered_in);
    fftw_free(recovered_out);

    return 0;
}