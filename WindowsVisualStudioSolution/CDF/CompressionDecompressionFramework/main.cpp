#include "TimingExperiment.h"
#include "SVDAlgorithms.h"

#define EIGEN_DONT_ALIGN
#include <Eigen/Dense> 
#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <algorithm>
#include <cstring>  


unsigned char* compressMatrix(double*& originalMatrix, const int x, const int y, const int z, const int k, int& size) {
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
	int U_rows = U.rows();
    int U_cols = U.cols();
    int V_rows = V.rows();
    int V_cols = V.cols();
    int S_length = S.size();

    int U_size = U.size() * sizeof(double);
    int V_size = V.size() * sizeof(double);
    int S_size = S.size() * sizeof(double);
    int header_size = 8 * sizeof(int); // x, y, z, U_rows, U_cols, V_rows, V_cols, S_size
    size = header_size + U_size + V_size + S_size;

    unsigned char* byteStream = new unsigned char[size];
    unsigned char* ptr = byteStream; 

    // Serialize the metadata into the byte stream
    memcpy(ptr, &x, sizeof(int));
    ptr += sizeof(int);
    memcpy(ptr, &y, sizeof(int));
    ptr += sizeof(int);
    memcpy(ptr, &z, sizeof(int));
    ptr += sizeof(int);
    memcpy(ptr, &U_rows, sizeof(int));
    ptr += sizeof(int);
    memcpy(ptr, &U_cols, sizeof(int));
    ptr += sizeof(int);
    memcpy(ptr, &V_rows, sizeof(int));
    ptr += sizeof(int);
    memcpy(ptr, &V_cols, sizeof(int));
    ptr += sizeof(int);
    memcpy(ptr, &S_length, sizeof(int));
    ptr += sizeof(int);
    // Now serialize the matrix and vector data
    // For U
    // Ret
    for (int i = 0; i < U_rows; ++i)
    {
        for (int j = 0; j < U_cols; ++j)
        {
            double value = U(i, j);
            memcpy(ptr, &value, sizeof(double));
            ptr += sizeof(double);
        }
    }
    // For V
    for (int i = 0; i < V_rows; ++i)
    {
        for (int j = 0; j < V_cols; ++j)
        {
            double value = V(i, j);
            memcpy(ptr, &value, sizeof(double));
            ptr += sizeof(double);
        }
    }

    // For S
    for (int i = 0; i < S_length; ++i)
    {
        double value = S(i);
        memcpy(ptr, &value, sizeof(double));
        ptr += sizeof(double);
    }
    return byteStream;
}

double* decompressMatrix(unsigned char*& compressedData, int size) {
    unsigned char* ptr = compressedData;

    // Deserialize the metadata
    int x, y, z;
    memcpy(&x, ptr, sizeof(int));
    ptr += sizeof(int);
    memcpy(&y, ptr, sizeof(int));
    ptr += sizeof(int);
    memcpy(&z, ptr, sizeof(int));
    ptr += sizeof(int);
    int U_rows, U_cols, V_rows, V_cols, S_length;
    memcpy(&U_rows, ptr, sizeof(int));
    ptr += sizeof(int);
    memcpy(&U_cols, ptr, sizeof(int));
    ptr += sizeof(int);
    memcpy(&V_rows, ptr, sizeof(int));
    ptr += sizeof(int);
    memcpy(&V_cols, ptr, sizeof(int));
    ptr += sizeof(int);
    memcpy(&S_length, ptr, sizeof(int));
    ptr += sizeof(int);

    // Deserialize the matrices 
    Eigen::MatrixXd U(U_rows, U_cols);
    for (int i = 0; i < U_rows; ++i) {
        for (int j = 0; j < U_cols; ++j) {
            double value;
            memcpy(&value, ptr, sizeof(double));
            U(i, j) = value;
            ptr += sizeof(double);
        }
    }
    // For V
    Eigen::MatrixXd V(V_rows, V_cols);
    for (int i = 0; i < V_rows; ++i) {
        for (int j = 0; j < V_cols; ++j) {
            double value;
            memcpy(&value, ptr, sizeof(double));
            V(i, j) = value;
            ptr += sizeof(double);
        }
    }

    // For S
    Eigen::VectorXd S(S_length);
    for (int i = 0; i < S_length; ++i) {
        double value;
        memcpy(&value, ptr, sizeof(double));
        S(i) = value;
        ptr += sizeof(double);
    }

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
int matrixSizeInBytes(int rows, int cols) {
    return sizeof(double) * rows * cols;
}

int main() {
    int x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    int k = std::min(y, z) / 2;

    // Define original matrix
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);
    int originalMatrixBytes = matrixSizeInBytes(x * y, z);

    // Compress the matrix
    int compressedSize;
    unsigned char* compressedMatrix = compressMatrix(originalMatrix, x, y, z, k, compressedSize);

    // Decompress the matrix
    double* decompressedMatrix = decompressMatrix(compressedMatrix, compressedSize);

	// Printing comparison of original and decompressed data
    Utilities::printComparison(originalMatrix, decompressedMatrix, x, y, z);

    std::cout << "Original matrix size: " << originalMatrixBytes << " bytes" << std::endl;
    std::cout << "Compressed matrix size: " << compressedSize << " bytes" << std::endl;
    // Printing various types of error between original and decompressed data
    Utilities::printError(originalMatrix, decompressedMatrix, x, y, z);

    // Freeing the memory
    delete[] originalMatrix;
    delete[] compressedMatrix;
    delete[] decompressedMatrix;

    return 0;
}
