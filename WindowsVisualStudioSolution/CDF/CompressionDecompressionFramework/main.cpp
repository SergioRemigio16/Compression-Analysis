#include "TimingExperiment.h"
#include <Eigen/Dense>
#include "Utilities.h"
#include <iomanip>


#include <iostream>
#define IDX(i, j, k) ((i)*y*z + (j)*z + (k))  // 3D to 1D indexing

using namespace std;
using namespace Eigen;


void printComparison(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z, double epsilon = 1e-9) {
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
                } else {
                    std::cout << "]\n";
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}


int main() {
    // Define the dimensions of your 3D data
    int x = 3; 
    int y = 7;
    int z = 7; 

    // Define your original data with your function
	
    double* originalMatrixData = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);
    // Convert to Eigen's VectorXd for easier manipulation
    VectorXd originalMatrix = Map<VectorXd>(originalMatrixData, x * y * z);

    // Reshape your originalMatrix into a 2D matrix
    MatrixXd dataMatrix(x * y, z);
    for(int i = 0; i < x * y; i++) {
        for(int j = 0; j < z; j++) {
            dataMatrix(i, j) = originalMatrix[i*z + j];
        }
    }

    // Perform SVD
    JacobiSVD<MatrixXd> svd(dataMatrix, ComputeThinU | ComputeThinV);

    // Choose the number of singular values to keep
    // This number k is a parameter you should choose depending on your requirement for compression and data accuracy.
    // A smaller k will result in higher compression but less accuracy in the reconstructed data
    // A larger k will result in less compression but higher accuracy in the reconstructed data
    // k should be less than or equal to min(x*y, z)
    int k = 2 ; // replace with the number of singular values you want to keep

    // Now we will get the first k columns of U and V, and the first k singular values
    MatrixXd U = svd.matrixU().leftCols(k);
    MatrixXd V = svd.matrixV().leftCols(k);
    VectorXd S = svd.singularValues().head(k);

    // Print the size of original and compressed matrix
    cout << "Size of original matrix: " << originalMatrix.size() << endl;
    long long compressedSize = U.cols() * (U.rows() + V.rows() + 1); // U*S*V'
    cout << "Size of compressed matrix: " << compressedSize << endl;

    // Decompress the matrix
    MatrixXd decompressedMatrix = U * S.asDiagonal() * V.transpose();

    // Convert decompressed matrix back to original 3D shape
    double* decompressedMatrixData = new double[x * y * z];
    for (int i = 0; i < x * y; i++) {
        for (int j = 0; j < z; j++) {
            decompressedMatrixData[i * z + j] = decompressedMatrix(i, j);
        }
    }

    // Calculate and print the data loss
    double dataLoss = (dataMatrix - decompressedMatrix).norm() / dataMatrix.norm();
    cout << "Data loss: " << dataLoss << endl;

    // Print decompressed matrix with original 3D shape
    printComparison(originalMatrixData, decompressedMatrixData, x, y, z);

    // Clean up
    delete[] originalMatrixData;
    delete[] decompressedMatrixData;

    return 0;
}



/*
int main() {
	return runTimingExperiment();
}
*/

/*

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fftw3.h>

// Define typedef for convenience
typedef std::vector<std::complex<double>> ComplexVec;

int main() {
    size_t x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    double compressionRatio = 0.1; // Compression ratio. 0.5 means keeping half of the frequency components.

    // Define original matrix
    double* originalMatrix = new double[x * y * z];

    // Initialize the original matrix with some data
    for (size_t i = 0; i < x * y * z; ++i) {
        originalMatrix[i] = i;
    }

    // Print original matrix
    std::cout << "Original Matrix:\n";
    for (size_t i = 0; i < x * y * z; ++i) {
        std::cout << originalMatrix[i] << ' ';
    }
    std::cout << "\n";

    // Perform FFT
    ComplexVec complexMatrix(x * y * z);
    fftw_plan p = fftw_plan_dft_r2c_3d(x, y, z, originalMatrix,
        reinterpret_cast<fftw_complex*>(complexMatrix.data()), FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    // Compress the data by zeroing out the least significant components
    size_t compressedSize = static_cast<size_t>(x * y * z * compressionRatio);
    for (size_t i = compressedSize; i < x * y * z; ++i) {
        complexMatrix[i] = 0;
    }

    // Perform inverse FFT
    double* decompressedMatrix = new double[x * y * z];
    fftw_plan q = fftw_plan_dft_c2r_3d(x, y, z, reinterpret_cast<fftw_complex*>(complexMatrix.data()),
        decompressedMatrix, FFTW_ESTIMATE);
    fftw_execute(q);
    fftw_destroy_plan(q);

    // FFTW's backwards transform does not normalize the result, so we have to do it manually
    for (size_t i = 0; i < x * y * z; ++i) {
        decompressedMatrix[i] /= (x * y * z);
    }

    // Print decompressed matrix
    std::cout << "Decompressed Matrix:\n";
    for (size_t i = 0; i < x * y * z; ++i) {
        std::cout << decompressedMatrix[i] << ' ';
    }
    std::cout << "\n";

    // Print size of original and compressed matrices in bytes
    std::cout << "Size of original matrix (bytes): " << (x * y * z) * sizeof(double) << "\n";
    std::cout << "Size of compressed matrix (bytes): " << compressedSize * sizeof(std::complex<double>) << "\n";

    // Compute and print Mean Squared Error
    double mse = 0;
    for (size_t i = 0; i < x * y * z; ++i) {
        double error = originalMatrix[i] - decompressedMatrix[i];
        mse += error * error;
    }
    mse /= (x * y * z);
    std::cout << "Mean Squared Error: " << mse << "\n";

    // Deallocate memory
    delete[] originalMatrix;
    delete[] decompressedMatrix;

    return 0;
}
*/

