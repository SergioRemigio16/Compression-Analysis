#include "TimingExperiment.h"
#include <Eigen/Dense>
#include "Utilities.h"

using namespace std;
using namespace Eigen;

/*
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
    Utilities::printComparison(originalMatrixData, decompressedMatrixData, x, y, z);

    // Clean up
    delete[] originalMatrixData;
    delete[] decompressedMatrixData;

    return 0;
}
*/

///*
int main() {
	return runTimingExperiment();
}
//*/
