#include "TimingExperiment.h"
#include "SVDAlgorithms.h"

#define EIGEN_DONT_ALIGN
#include <Eigen/Dense> 
#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <algorithm>
#include <cstring>  

void printError2(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z) {
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

void serializeInt(std::vector<unsigned char>& byteStream, int value) {
    // Write the integer to the byte stream
    const unsigned char* start = reinterpret_cast<unsigned char*>(&value);
    const unsigned char* end = start + sizeof(int);
    byteStream.insert(byteStream.end(), start, end);
}

int deserializeInt(std::vector<unsigned char>& byteStream) {
    // Read the integer from the byte stream
    int value;
    std::copy(byteStream.begin(), byteStream.begin() + sizeof(int), reinterpret_cast<unsigned char*>(&value));

    // Erase the used data from the byteStream
    byteStream.erase(byteStream.begin(), byteStream.begin() + sizeof(int));

    return value;
}


void serializeMatrix(std::vector<unsigned char>& byteStream, const Eigen::MatrixXd& matrix) {

    // First write the size of the matrix
    serializeInt(byteStream, matrix.rows());
    serializeInt(byteStream, matrix.cols());

    // Then write the data
    const unsigned char* dataStart = reinterpret_cast<const unsigned char*>(matrix.data());
    const unsigned char* dataEnd = dataStart + sizeof(double) * matrix.size();
    byteStream.insert(byteStream.end(), dataStart, dataEnd);
}

Eigen::MatrixXd deserializeMatrix(std::vector<unsigned char>& byteStream) {
    // First read the size of the matrix
    int rows = deserializeInt(byteStream);
    int cols = deserializeInt(byteStream);

    // Then read the data
    Eigen::MatrixXd matrix(rows, cols);
    std::copy(byteStream.begin(), byteStream.begin() + sizeof(double) * matrix.size(), reinterpret_cast<unsigned char*>(matrix.data()));

    // Erase the used data from the byteStream
    byteStream.erase(byteStream.begin(), byteStream.begin() + sizeof(double) * matrix.size());

    return matrix;
}


unsigned char* compressMatrix(double*& originalMatrix, const int x, const int y, const int z, const int k, size_t& size) {
    // Convert to Eigen's VectorXd for easier manipulation
    Eigen::VectorXd eigenMatrix = Eigen::Map<Eigen::VectorXd>(originalMatrix, x * y * z);

    // Reshape your originalMatrix into a 2D matrix
    Eigen::MatrixXd dataMatrix(x * y, z);
    for (int i = 0; i < x * y; i++) {
        for (int j = 0; j < z; j++) {
            dataMatrix(i, j) = eigenMatrix[i * z + j];
        }
    }

    // Perform SVD
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(dataMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Now we will get the first k columns of U and V, and the first k singular values
    Eigen::MatrixXd U = svd.matrixU().leftCols(k);
    Eigen::MatrixXd V = svd.matrixV().leftCols(k);
    Eigen::VectorXd S = svd.singularValues().head(k);

    // Serialize the matrices and metadata into a byte stream
    std::vector<unsigned char> byteStream;
    serializeMatrix(byteStream, U);
    serializeMatrix(byteStream, V);
    serializeMatrix(byteStream, S);
    serializeInt(byteStream, x);
    serializeInt(byteStream, y);
    serializeInt(byteStream, z);

    // Save the size of the byteStream
    size = byteStream.size();

    // Convert the vector to a raw pointer
    unsigned char* rawBytes = new unsigned char[size];
    std::copy(byteStream.begin(), byteStream.end(), rawBytes);

    return rawBytes;
}

double* decompressMatrix(unsigned char*& compressedData, size_t size) {
    // Initialize byte stream with data from the raw pointer
    std::vector<unsigned char> byteStream(compressedData, compressedData + size);

    // Deserialize the matrices and metadata
    Eigen::MatrixXd U = deserializeMatrix(byteStream);
    Eigen::MatrixXd V = deserializeMatrix(byteStream);
    Eigen::MatrixXd S = deserializeMatrix(byteStream);
    int x = deserializeInt(byteStream);
    int y = deserializeInt(byteStream);
    int z = deserializeInt(byteStream);

    // Perform matrix multiplication to decompress the data
    Eigen::MatrixXd decompressedMatrix = U * S.asDiagonal() * V.transpose();

    // Reshape and convert the decompressed matrix to a raw pointer
    double* decompressedData = new double[x * y * z];
    for (int i = 0; i < x * y; i++) {
        for (int j = 0; j < z; j++) {
            decompressedData[i * z + j] = decompressedMatrix(i, j);
        }
    }

    return decompressedData;
}

// A function to calculate the size of a matrix in bytes.
size_t matrixSizeInBytes(int rows, int cols) {
    return sizeof(double) * rows * cols;
}

int main() {
    size_t x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    int k = std::min(y, z) / 2;

    // Define original matrix
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);
    size_t originalMatrixBytes = matrixSizeInBytes(x * y, z);

    std::cout << "Original matrix size: " << originalMatrixBytes << " bytes" << std::endl;

    // Compress the matrix
    size_t compressedSize;
    unsigned char* compressedMatrix = compressMatrix(originalMatrix, x, y, z, k, compressedSize);
    size_t compressedMatrixBytes = matrixSizeInBytes(k, z) * 2 + matrixSizeInBytes(k, k) + 4 * sizeof(int); // U, V, S, and 4 integers (x, y, z, k)

    std::cout << "Compressed matrix size: " << compressedMatrixBytes << " bytes" << std::endl;

    // Decompress the matrix
    double* decompressedMatrix = decompressMatrix(compressedMatrix, compressedSize);

	// Printing comparison of original and decompressed data
    Utilities::printComparison(originalMatrix, decompressedMatrix, x, y, z);

    // Printing various types of error between original and decompressed data
    printError2(originalMatrix, decompressedMatrix, x, y, z);

    // Freeing the memory
    delete[] originalMatrix;
    delete[] compressedMatrix;
    delete[] decompressedMatrix;

    return 0;
}
