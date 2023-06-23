#include "Utilities.h"

// This macro calculates the index for a 3D point in a 1D array. It works by treating the 1D array as a 3D array,
// with each (i, j, k) triplet corresponding to a unique index in the 1D array.
#define IDX(i, j, k) ((i)*y*z + (j)*z + (k))  // 3D to 1D indexing

/**
 * This function generates a 3D sinusoidal wave with customizable amplitude, frequency, phase, and propagation rates along each axis.
 * The generated wave is stored in a one-dimensional array using row-major ordering.
 * 
 * @param x The size of the wave along the x dimension.
 * @param y The size of the wave along the y dimension.
 * @param z The size of the wave along the z dimension.
 * @param frequency The frequency of the wave, i.e., how many times the wave cycles over a certain distance.
 * @param amplitude The amplitude of the wave, i.e., the maximum displacement of the wave from its equilibrium position.
 * @param phase The phase of the wave, i.e., the initial angle of a sinusoidal function at its origin.
 * @param fx The propagation factor along the x dimension. A higher value means the wave propagates faster along this dimension.
 * @param fy The propagation factor along the y dimension. A higher value means the wave propagates faster along this dimension.
 * @param fz The propagation factor along the z dimension. A higher value means the wave propagates faster along this dimension.
 * 
 * @return A pointer to the dynamically allocated array storing the wave data. Each element represents the wave's value at a point in 3D space.
 *         The size of the array is x * y * z. The caller of the function is responsible for deallocating this memory.
 */
double* Utilities::createMatrixWave(int x, int y, int z, double frequency, double amplitude, double phase, double fx, double fy, double fz) {
    // Allocate a 1D array to hold the wave data. The size of the array is the product of the dimensions x, y, and z.
    double* matrix = new double[x * y * z];

    // Loop over each index in the x dimension.
    for (int i = 0; i < x; i++) {
        // Loop over each index in the y dimension.
        for (int j = 0; j < y; j++) {
            // Loop over each index in the z dimension.
            for (int k = 0; k < z; k++) {
                // Calculate the value of the wave at the point (i, j, k). The value is determined by a sine function, which 
                // is influenced by the frequency, amplitude, and phase input parameters, as well as the propagation factors 
                // fx, fy, and fz for each dimension. The result is a wave that can propagate differently along the x, y, and z dimensions.
                double waveValue = amplitude * sin(frequency * (fx * i + fy * j + fz * k) + phase);

                // Store the wave value in the 1D array. The 3D index (i, j, k) is converted to a 1D index using the IDX macro.
                matrix[IDX(i, j, k)] = waveValue;
            }
        }
    }

    // Return the 1D array containing the wave data.
    return matrix;
}




void Utilities::writeWaveToCSV(double* matrix, int x, int y, int z, const std::string& filename) {
    std::ofstream file;
    file.open(filename);

    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            for (int k = 0; k < z; k++) {
                file << i << ',' << j << ',' << k << ',' << matrix[IDX(i, j, k)] << '\n';
            }
        }
    }

    file.close();
}

double* Utilities::createMatrixRandom(int x, int y, int z, double minVal, double maxVal) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(minVal, maxVal);

    double* matrix = new double[x * y * z];
    for (int i = 0; i < x * y * z; i++) {
        matrix[i] = dis(gen);
    }
    return matrix;
}

double Utilities::measureTimingOverhead() {
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> overhead = end - start;
    return overhead.count();
}

void Utilities::printMatrix(const double* matrix, int x, int y, int z) {
    for (int k = 0; k < z; k++) {
        std::cout << "z = " << k << ":\n";
        for (int i = 0; i < x; i++) {
            for (int j = 0; j < y; j++) {
                std::cout << " " << matrix[IDX(i, j, k)];
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}


double Utilities::calculateMSE(const double* originalData, const double* decompressedData, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        double diff = originalData[i] - decompressedData[i];
        sum += diff * diff;
    }
    return sum / size;
}


void Utilities::printComparison(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z) {
    double epsilon = 1e-9;
    std::cout << std::fixed << std::setprecision(15);  // for consistent number of decimal places
    int width = 20;  // 6 decimal digits, 1 dot, 1 sign and up to 4 whole part digits

    for (int k = 0; k < z; k++) {
        std::cout << "z = " << k << ":\n";
        for (int i = 0; i < x; i++) {
            for (int j = 0; j < y; j++) {
                double original = originalMatrix[IDX(i, j, k)];
                double decompressed = decompressedMatrix[IDX(i, j, k)];

                std::cout << "[";
                std::cout << std::setw(width) << original << ",\n ";
                std::cout << std::setw(width) << decompressed;

                // Check if the numbers differ significantly
                if (std::abs(original - decompressed) > epsilon) {
                    std::cout << "*]\n";
                }
                else {
                    std::cout << "]\n";
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}