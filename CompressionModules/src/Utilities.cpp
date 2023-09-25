#include "Utilities.h"

// This macro calculates the index for a 3D point in a 1D array. It works by treating the 1D array as a 3D array,
// with each (i, j, k) triplet corresponding to a unique index in the 1D array.
#define IDX(i, j, k) ((i)*y*z + (j)*z + (k))  // 3D to 1D indexing

double* Utilities::createMatrixWave(int x, int y, int z, double amplitude, double phase, double fx, double fy, double fz) {
    // Allocate a 1D array to hold the wave data. The size of the array is the product of the dimensions x, y, and z.
    double* matrix = new double[x * y * z];

    // Loop over each index in the x dimension.
    for (int i = 0; i < x; i++) {
        // Loop over each index in the y dimension.
        for (int j = 0; j < y; j++) {
            // Loop over each index in the z dimension.
            for (int k = 0; k < z; k++) {
                // Calculate the value of the wave at the point (i, j, k). The value is determined by a sine function, which 
                // is influenced by the amplitude, and phase input parameters, as well as the propagation factors 
                // fx, fy, and fz for each dimension. The result is a wave that can propagate differently along the x, y, and z dimensions.
                double waveValue = amplitude * sin((fx * i + fy * j + fz * k) + phase);

                // Store the wave value in the 1D array. The 3D index (i, j, k) is converted to a 1D index using the IDX macro.
                matrix[IDX(i, j, k)] = waveValue;
            }
        }
    }

    // Return the 1D array containing the wave data.
    return matrix;
}

double* Utilities::createWave1D(int n, double amplitude, double phase, double fx) {
    double* array = new double[n];
    double next = 0.0001;
    for (int i = 0; i < n; i++) {
        double waveValue = amplitude * sin((fx * next * i) + phase);
        array[i] = waveValue;
    }
    return array;
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

size_t Utilities::calculateOriginalDataBytes(size_t size) {
    return size * sizeof(double);
}


void Utilities::printComparison(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z, int originalMatrixBytes, int compressedSize) {
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

void Utilities::printComparison(const double* originalArray, const double* decompressedArray, int n, int originalMatrixBytes, int compressedSize) {
    double epsilon = 1e-9;
    std::cout << std::fixed << std::setprecision(15);  // for consistent number of decimal places
    int width = 20;  // 6 decimal digits, 1 dot, 1 sign and up to 4 whole part digits

    for (int i = 0; i < n; i++) {
        double original = originalArray[i];
        double decompressed = decompressedArray[i];

        std::cout << "[";
        std::cout << std::setw(width) << original << ", ";
        std::cout << std::setw(width) << decompressed;

        // Check if the numbers differ significantly
        if (std::abs(original - decompressed) > epsilon) {
            std::cout << "*]\n";
        }
        else {
            std::cout << "]\n";
        }
    }
    std::cout << "Original matrix size: " << originalMatrixBytes << " bytes" << std::endl;
    std::cout << "Compressed matrix size: " << compressedSize << " bytes" << std::endl;
    std::cout << std::fixed << std::setprecision(2) << "Compression Ratio: " << (double)originalMatrixBytes/(double)compressedSize << std::endl;
    std::cout << std::fixed << std::setprecision(15);
}


void Utilities::printError(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z) {
	// Computing and printing Mean Squared Error, Mean Absolute Error,
    // Maximum Absolute Error, Root Mean Squared Error and Peak Signal to Noise Ratio.
    double mse = 0;
    for (size_t i = 0; i < x * y * z; ++i) {
        double error = originalMatrix[i] - decompressedMatrix[i];
        mse += error * error;
    }
    mse /= (x * y * z);

    std::cout << "Mean Squared Error: " << mse << "\n";
    // Compute and print Mean Absolute Error
    double mae = 0;
    for (size_t i = 0; i < x * y * z; ++i) {
        double error = std::abs(originalMatrix[i] - decompressedMatrix[i]);
        mae += error;
    }
    mae /= (x * y * z);
    std::cout << "Mean Absolute Error: " << mae << "\n";

    // Compute and print Maximum Absolute Error
    double maxError = 0;
    for (size_t i = 0; i < x * y * z; ++i) {
        double error = std::abs(originalMatrix[i] - decompressedMatrix[i]);
        maxError = std::max(maxError, error);
    }
    std::cout << "Max Absolute Error: " << maxError << "\n";

    // Compute and print Root Mean Squared Error
    double rmse = sqrt(mse);
    std::cout << "Root Mean Squared Error: " << rmse << "\n";

    // Compute and print Peak Signal to Noise Ratio (PSNR)
    double maxOriginalValue = *std::max_element(originalMatrix, originalMatrix + (x * y * z));
    double psnr = 20 * log10(maxOriginalValue / rmse);
    std::cout << "Peak Signal to Noise Ratio (in dB): " << psnr << "\n";
}

void Utilities::printError(const double* originalArray, const double* decompressedArray, int n) {
    // Computing and printing Mean Squared Error
    double mse = 0;
    for (size_t i = 0; i < n; ++i) {
        double error = originalArray[i] - decompressedArray[i];
        mse += error * error;
    }
    mse /= n;

    std::cout << "Mean Squared Error: " << mse << "\n";

    // Compute and print Mean Absolute Error
    double mae = 0;
    for (size_t i = 0; i < n; ++i) {
        double error = std::abs(originalArray[i] - decompressedArray[i]);
        mae += error;
    }
    mae /= n;
    std::cout << "Mean Absolute Error: " << mae << "\n";

    // Compute and print Maximum Absolute Error
    double maxError = 0;
    for (size_t i = 0; i < n; ++i) {
        double error = std::abs(originalArray[i] - decompressedArray[i]);
        maxError = std::max(maxError, error);
    }
    std::cout << "Max Absolute Error: " << maxError << "\n";

    // Compute and print Root Mean Squared Error
    double rmse = std::sqrt(mse);
    std::cout << "Root Mean Squared Error: " << rmse << "\n";

    // Compute and print Peak Signal to Noise Ratio (PSNR)
    double maxOriginalValue = *std::max_element(originalArray, originalArray + n);
    double psnr = 20 * std::log10(maxOriginalValue / rmse);
    std::cout << "Peak Signal to Noise Ratio (in dB): " << psnr << "\n";
}



double* Utilities::readBinaryFile(const std::string& fileName, int& size) {
    std::ifstream inputFile(fileName, std::ios::binary);
    if (!inputFile) {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return nullptr;
    }

    // Read the size of the array from the first 4 bytes
    inputFile.read(reinterpret_cast<char*>(&size), sizeof(int));
    
    if (!inputFile) {
        std::cerr << "Error: Could not read the size of the array from the file." << std::endl;
        return nullptr;
    }

    // Allocate memory for the double array
    double* array = new double[size];
    
    // Read the double values
    inputFile.read(reinterpret_cast<char*>(array), size * sizeof(double));
    
    if (!inputFile) {
        std::cerr << "Error: Could not read the double values from the file." << std::endl;
        delete[] array;
        return nullptr;
    }

    // Close the input file
    inputFile.close();

    return array;
}