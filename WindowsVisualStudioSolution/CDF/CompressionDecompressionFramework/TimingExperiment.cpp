#include "TimingExperiment.h"
#include <iomanip>
#include <fstream>

#define RUNS 20
#define WARMUP_RUNS 5

// Parameters for wave matrix creation
struct WaveParams {
    double frequency;
    double amplitude;
    double phase;
    double fx;
    double fy;
    double fz;
};


void runExperiment(int x, int y, int z, bool useWave, bool visualizeData,
    const WaveParams& waveParams, int rate, int k, double compressionRatio, int algorithm, std::ofstream& csv) {

    std::string algorithmName;

    switch (algorithm) {
    case 0: algorithmName = "ZFP"; break;
    case 1: algorithmName = "SVD"; break;
    case 2: algorithmName = "FFT"; break;
    default:
        std::cerr << "Invalid algorithm specified.\n";
        return;
    }


    /*
    std::cout << "Name: " + algorithmName
        << ", Frequency: " << waveParams.frequency
        << ", Amplitude: " << waveParams.amplitude
        << ", Phase: " << waveParams.phase
        << ", Fx: " << waveParams.fx
        << ", Fy: " << waveParams.fy
        << ", Fz: " << waveParams.fz
        << ", Rate: " << rate
        << ", K: " << k
        << ", CompressionRatio: " << std::left << std::setw(4) << compressionRatio;
    */

    int size = x * y * z;
    double timingOverhead = Utilities::measureTimingOverhead();

    std::vector<double> compressTimes, decompressTimes, mseValues;
    std::vector<size_t> originalSizes, compressedSizes;

    // Preallocating memory
    compressTimes.reserve(RUNS);
    decompressTimes.reserve(RUNS);
    mseValues.reserve(RUNS);
    originalSizes.reserve(RUNS);
    compressedSizes.reserve(RUNS);

    for (int i = 0; i < RUNS + WARMUP_RUNS; i++) {
        double* originalMatrix =
            Utilities::createMatrixWave(x, y, z, waveParams.frequency, waveParams.amplitude, waveParams.phase, waveParams.fx, waveParams.fy, waveParams.fz);

        if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
            std::cout << "Original Data:\n";
            Utilities::printMatrix(originalMatrix, x, y, z);
        }

        unsigned char* compressedMatrix = nullptr;
        double* decompressedMatrix = nullptr;
        int compressedSize;
        std::chrono::duration<double> compressTime, decompressTime;


        switch (algorithm) {
        case 0:  // ZFP
        {
            auto start = std::chrono::high_resolution_clock::now();
            compressedMatrix = ZFPAlgorithms::compressMatrix(originalMatrix, x, y, z, rate, compressedSize);
            auto end = std::chrono::high_resolution_clock::now();
            compressTime = end - start;

            start = std::chrono::high_resolution_clock::now();
            decompressedMatrix = ZFPAlgorithms::decompressMatrix(compressedMatrix, compressedSize);
            end = std::chrono::high_resolution_clock::now();
            decompressTime = end - start;
            break;
        }
        case 1:  // SVD
        {
            auto start = std::chrono::high_resolution_clock::now();
            compressedMatrix = SVDAlgorithms::compressMatrix(originalMatrix, x, y, z, k, compressedSize);
            auto end = std::chrono::high_resolution_clock::now();
            compressTime = end - start;

            start = std::chrono::high_resolution_clock::now();
            decompressedMatrix = SVDAlgorithms::decompressMatrix(compressedMatrix, compressedSize);
            end = std::chrono::high_resolution_clock::now();
            decompressTime = end - start;
            break;
        }
        case 2:  // FFT
        {
            auto start = std::chrono::high_resolution_clock::now();
            compressedMatrix = FFTAlgorithms::compressMatrix(originalMatrix, x, y, z, compressionRatio / 100.0, compressedSize);
            auto end = std::chrono::high_resolution_clock::now();
            compressTime = end - start;

            start = std::chrono::high_resolution_clock::now();
            decompressedMatrix = FFTAlgorithms::decompressMatrix(compressedMatrix, compressedSize);
            end = std::chrono::high_resolution_clock::now();
            decompressTime = end - start;
            break;
        }
        }



        if (i >= WARMUP_RUNS) {
            compressTimes.push_back(compressTime.count() - timingOverhead);
            decompressTimes.push_back(decompressTime.count() - timingOverhead);
            mseValues.push_back(Utilities::calculateMSE(originalMatrix, decompressedMatrix, size));
            originalSizes.push_back(size * sizeof(double));
            compressedSizes.push_back(compressedSize);
        }

        if (visualizeData && i == RUNS + WARMUP_RUNS - 1) {
            std::cout << "Decompressed Data:\n";
            Utilities::printMatrix(decompressedMatrix, x, y, z);
        }

        delete[] compressedMatrix;
        delete[] decompressedMatrix;

        delete[] originalMatrix;
    }

    double meanCompressTime = std::accumulate(compressTimes.begin(), compressTimes.end(), 0.0) / RUNS;
    double meanDecompressTime = std::accumulate(decompressTimes.begin(), decompressTimes.end(), 0.0) / RUNS;
    double meanMSE = std::accumulate(mseValues.begin(), mseValues.end(), 0.0) / RUNS;
    size_t meanOriginalSize = std::accumulate(originalSizes.begin(), originalSizes.end(), 0) / RUNS;
    size_t meanCompressedSize = std::accumulate(compressedSizes.begin(), compressedSizes.end(), 0) / RUNS;

    // Writing to CSV
    csv << algorithmName << ',' << waveParams.frequency << ',' << waveParams.amplitude << ',' << waveParams.phase << ',' << waveParams.fx << ','
        << waveParams.fy << ',' << waveParams.fz << ',' << rate << ',' << k << ',' << compressionRatio << ','
        << x << ',' << y << ',' << z << ',' << meanCompressTime << ',' << meanDecompressTime << ',' << meanMSE << ','
        << meanOriginalSize << ',' << meanCompressedSize << std::endl;

    /*
    std::cout << std::left << std::setw(15) << "Matrix size:"
        << std::setw(15) << (std::to_string(x) + "x" + std::to_string(y) + "x" + std::to_string(z) + ",")
        << std::setw(25) << ("meanCompressTime: " + std::to_string(meanCompressTime) + ",")
        << std::setw(30) << ("meanDecompressTime: " + std::to_string(meanDecompressTime) + ",")
        << std::setw(20) << ("meanMSE: " + std::to_string(meanMSE) + ",")
        << std::setw(25) << ("meanOriginalSize: " + std::to_string(meanOriginalSize) + ",")
        << std::setw(25) << ("meanCompressedSize: " + std::to_string(meanCompressedSize))
        << std::endl;
    */
}



void runExperimentsForRatesAndKs(int x, int y, int z, bool useWave, bool visualizeData,
    const WaveParams& waveParams, int minRate, int maxRate, int minK, int maxK,
    double minCompressionRatio, double maxCompressionRatio, double compressionRatioStep,
    int algorithm, std::ofstream& csv) {

    int minDimension = std::min({ x, y, z });
    maxK = std::min(maxK, minDimension);

    for (int rate = minRate; rate <= maxRate; ++rate) {
        for (int k = minK; k <= maxK; ++k) {
            for (double compressionRatio = minCompressionRatio; compressionRatio <= maxCompressionRatio;
                compressionRatio += compressionRatioStep) {

                if ((algorithm == 0 && (k != minK || compressionRatio != minCompressionRatio)) ||   // ZFP
                    (algorithm == 1 && (rate != minRate || compressionRatio != minCompressionRatio)) ||   // SVD
                    (algorithm == 2 && (rate != minRate || k != minK))) {   // FFT
                    continue;
                }

                runExperiment(x, y, z, useWave, visualizeData, waveParams, rate, k, compressionRatio, algorithm, csv);
            }
        }
    }
}

int runTimingExperiment() {
    std::ofstream csv("results.csv");
    csv << "Algorithm,Frequency,Amplitude,Phase,Fx,Fy,Fz,Rate,K,CompressionRatio,X,Y,Z,MeanCompressTime,MeanDecompressTime,MeanMSE,MeanOriginalSize,MeanCompressedSize\n";

    for (int x = 3; x <= 12; ++x) {
        for (double frequency = 1.0; frequency <= 2.0; frequency += 0.1) {
            std::cout << frequency << std::endl;
            for (double amplitude = 0.1; amplitude <= 2.0; amplitude += 0.1) {
                for (double phase = 0.0; phase <= M_PI; phase += M_PI / 2.0) {
                    for (double fx = 0.1; fx <= 1.0; fx += 0.3) {
                        for (double fy = 0.1; fy <= 1.0; fy += 0.3) {
                            for (double fz = 0.1; fz <= 1.0; fz += 0.3) {
                                WaveParams waveParams = { frequency, amplitude, phase, fx, fy, fz };

                                for (int algorithm = 0; algorithm < 3; ++algorithm) {
                                    runExperimentsForRatesAndKs(x, 7, 7, true, false, waveParams, 1, 64, 1, x * 7 * 7, 0.1, 1.0, 0.1, algorithm, csv);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    csv.close();
    return 0;
}
