#include "SVDAlgorithms.h"

SVDAlgorithms::SVDResult SVDAlgorithms::compressMatrix(double*& originalMatrix, const int x, const int y, const int z, const int k) {
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
    SVDAlgorithms::SVDResult result;
    result.U = U;
    result.V = V;
    result.S = S;
    result.x = x;
    result.y = y;
    result.z = z;

    return result; // return U, S and V
}

double* SVDAlgorithms::decompressMatrix(const SVDAlgorithms::SVDResult& compressedData, const int x, const int y, const int z) {
    // Decompress the matrix
    Eigen::MatrixXd decompressedMatrix = compressedData.U * compressedData.S.asDiagonal() * compressedData.V.transpose();

    // Convert decompressed matrix back to original 3D shape
    double* decompressedMatrixData = new double[x * y * z];
    for (int i = 0; i < x * y; i++) {
        for (int j = 0; j < z; j++) {
            decompressedMatrixData[i * z + j] = decompressedMatrix(i, j);
        }
    }
    return decompressedMatrixData;
}
void SVDAlgorithms::printByteSizeReport(const SVDAlgorithms::SVDResult& compressedData) {
    // Compute size assuming no padding
    size_t sizeU = compressedData.U.size() * sizeof(double);
    size_t sizeV = compressedData.V.size() * sizeof(double);
    size_t sizeS = compressedData.S.size() * sizeof(double);

    size_t totalSize = sizeU + sizeV + sizeS;

    std::cout << "Size of U: " << sizeU << " bytes" << std::endl;
    std::cout << "Size of V: " << sizeV << " bytes" << std::endl;
    std::cout << "Size of S: " << sizeS << " bytes" << std::endl;
    std::cout << "Total size: " << totalSize << " bytes" << std::endl;
    std::cout << "EIGEN_MAX_STATIC_ALIGN_BYTES: " << EIGEN_MAX_STATIC_ALIGN_BYTES << std::endl;
}

void SVDAlgorithms::calculateDecompressedDataBytes(const SVDAlgorithms::SVDResult& compressedData) {

}
