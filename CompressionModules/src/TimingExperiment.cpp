#include "TimingExperiment.h"
#include <iomanip>
#include <fstream>
#include <thread>
#include <vector>

#define RUNS 5
#define WARMUP_RUNS 1

// Parameters for wave matrix creation
struct WaveParams {
	double amplitude;
	double phase;
	double fx;
	double fy;
	double fz;
};


int runExperiment(int x, int y, int z, bool useWave, const WaveParams& waveParams, 
	int ZFPRate, int SVDK, double FFTRatio, int N, int Q, int S,
	int algorithm, std::ofstream& csv) {

	std::string algorithmName;

	switch (algorithm) {
	case 0: algorithmName = "ZFP"; break;
	case 1: algorithmName = "SVD"; break;
	case 2: algorithmName = "Chebyshev"; break;
	case 3: algorithmName = "FFT"; break;
	default:
		std::cerr << "Invalid algorithm specified.\n";
		return 2;
	}

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

	// This will store the maximum absolute error for each run
	std::vector<double> maxErrors;
	maxErrors.reserve(RUNS);

	// Variables to track the maxError, and the original and decompressed values that generated it
	double maxError = 0;
	double maxErrorOriginalValue = 0, maxErrorDecompressedValue = 0;

	for (int i = 0; i < RUNS + WARMUP_RUNS; i++) {
		double* originalMatrix =
			Utilities::createMatrixWave(x, y, z, waveParams.amplitude, waveParams.phase, waveParams.fx, waveParams.fy, waveParams.fz);

		unsigned char* compressedMatrix = nullptr;
		double* decompressedMatrix = nullptr;
		int compressedSize;
		std::chrono::duration<double> compressTime, decompressTime;


		switch (algorithm) {
		case 0:  // ZFP
		{
			auto start = std::chrono::high_resolution_clock::now();
			compressedMatrix = ZFPAlgorithms::compressMatrix(originalMatrix, x, y, z, ZFPRate, compressedSize);
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
			compressedMatrix = SVDAlgorithms::compressMatrix(originalMatrix, x, y, z, SVDK, compressedSize);
			auto end = std::chrono::high_resolution_clock::now();
			compressTime = end - start;

			start = std::chrono::high_resolution_clock::now();
			decompressedMatrix = SVDAlgorithms::decompressMatrix(compressedMatrix, compressedSize);
			end = std::chrono::high_resolution_clock::now();
			decompressTime = end - start;
			break;
		}
		case 2: // Chebyshev
		{
			auto start = std::chrono::high_resolution_clock::now();
			compressedMatrix = ChebyshevAlgorithms::compressMatrix(originalMatrix, x, y, z, N, Q, S, compressedSize);
			auto end = std::chrono::high_resolution_clock::now();
			compressTime = end - start;

			start = std::chrono::high_resolution_clock::now();
			decompressedMatrix = ChebyshevAlgorithms::decompressMatrix(compressedMatrix, compressedSize);
			end = std::chrono::high_resolution_clock::now();
			decompressTime = end - start;
			break;
		}
		case 3:  // FFT
		{
			auto start = std::chrono::high_resolution_clock::now();
			compressedMatrix = FFTAlgorithms::compressMatrix(originalMatrix, x, y, z, FFTRatio, compressedSize);
			auto end = std::chrono::high_resolution_clock::now();
			compressTime = end - start;

			start = std::chrono::high_resolution_clock::now();
			decompressedMatrix = FFTAlgorithms::decompressMatrix(compressedMatrix, compressedSize);
			end = std::chrono::high_resolution_clock::now();
			decompressTime = end - start;
			break;
		}
		}

		// Send message to stop computing if compression size is worse than original size
		if (compressedSize > size * sizeof(double)) {
			return 1;
		}

		if (i >= WARMUP_RUNS) {
			compressTimes.push_back(compressTime.count() - timingOverhead);
			decompressTimes.push_back(decompressTime.count() - timingOverhead);
			mseValues.push_back(Utilities::calculateMSE(originalMatrix, decompressedMatrix, size));
			originalSizes.push_back(size * sizeof(double));
			compressedSizes.push_back(compressedSize);

			double currentMaxError = 0;
			double currentMaxErrorOriginalValue = 0, currentMaxErrorDecompressedValue = 0;

			for (int j = 0; j < size; ++j) {
				double error = std::abs(originalMatrix[j] - decompressedMatrix[j]);
				if (error > currentMaxError) {
					currentMaxError = error;
					currentMaxErrorOriginalValue = originalMatrix[j];
					currentMaxErrorDecompressedValue = decompressedMatrix[j];
				}
			}

			maxErrors.push_back(currentMaxError);

			// Update global maxError and related values
			if (currentMaxError > maxError) {
				maxError = currentMaxError;
				maxErrorOriginalValue = currentMaxErrorOriginalValue;
				maxErrorDecompressedValue = currentMaxErrorDecompressedValue;
			}
		}

		delete[] compressedMatrix;
		delete[] decompressedMatrix;

		delete[] originalMatrix;
	}

	double meanCompressTime = std::accumulate(compressTimes.begin(), compressTimes.end(), 0.0) / RUNS;
	double meanDecompressTime = std::accumulate(decompressTimes.begin(), decompressTimes.end(), 0.0) / RUNS;
	double meanMSE = std::accumulate(mseValues.begin(), mseValues.end(), 0.0) / RUNS;
	size_t meanCompressedSize = std::accumulate(compressedSizes.begin(), compressedSizes.end(), 0) / RUNS;
	double meanMaxError = std::accumulate(maxErrors.begin(), maxErrors.end(), 0.0) / RUNS;
	int originalSize = size * sizeof(double);


	csv << algorithmName << ','
	<< x << ','
	<< y << ','
	<< z << ','
	<< waveParams.amplitude << ','
	<< waveParams.phase << ','
	<< waveParams.fx << ','
	<< waveParams.fy << ','
	<< waveParams.fz << ',';
	// Conditionally write rate, k, and compressionRatio
	if (algorithm == 0) { // ZFP
		csv << ZFPRate << ",,,,,,";
	} 
	else if (algorithm == 1) { // SVD
		csv << ',' << SVDK << ",,,,,";
	} 
	else if (algorithm == 2) { // Chebyshev
		csv << ",," << N << ',' << Q << ',' << S << ",,";
	} 
	else if (algorithm == 3) { // FFT
		csv << ",,,,," << FFTRatio << ',';
	} 
	else {
		csv << ",,,,,,,"; // If no algorithm matches, make sure to add the commas
	}
	csv << meanCompressTime << ','
	<< meanDecompressTime << ','
	<< meanMSE << ','
	<< originalSize << ','
	<< meanCompressedSize << ','
	<< maxError << ','
	<< maxErrorOriginalValue << ','
	<< maxErrorDecompressedValue
	<< std::endl;

	/*
	std::cout << std::left
		<< std::setw(10) << "Matrix size:"
		<< std::setw(5) << (std::to_string(x) + "x" + std::to_string(y) + "x" + std::to_string(z) + ",")
		<< std::setw(5) << ("Name: " + algorithmName)
		<< std::setw(20) << (", Amplitude: " + std::to_string(waveParams.amplitude))
		<< std::setw(10) << (", Phase: " + std::to_string(waveParams.phase))
		<< std::setw(15) << (", Fx: " + std::to_string(waveParams.fx))
		<< std::setw(15) << (", Fy: " + std::to_string(waveParams.fy))
		<< std::setw(15) << (", Fz: " + std::to_string(waveParams.fz))
		<< std::setw(10) << (", Rate: " + std::to_string(rate))
		<< std::setw(5) << (", K: " + std::to_string(k))
		<< std::setw(30) << (", CompressionRatio: " + std::to_string(compressionRatio))
		<< std::setw(30) << ("meanCompressTime: " + std::to_string(meanCompressTime) + ",")
		<< std::setw(30) << ("meanDecompressTime: " + std::to_string(meanDecompressTime) + ",")
		<< std::setw(20) << ("meanMSE: " + std::to_string(meanMSE) + ",")
		<< std::setw(20) << ("OriginalSize: " + std::to_string(size * sizeof(double)) + ",")
		<< std::setw(30) << ("meanCompressedSize: " + std::to_string(meanCompressedSize))
		<< std::setw(20) << ("meanMaxError: " + std::to_string(meanMaxError) + ",")
		<< std::setw(20) << ("maxError: " + std::to_string(maxError) + ",")
		<< std::setw(20) << ("maxErrorOriginalValue: " + std::to_string(maxErrorOriginalValue) + ",")
		<< std::setw(20) << ("maxErrorDecompressedValue: " + std::to_string(maxErrorDecompressedValue))
		<< std::endl;
		*/

		return 0;
}

void runExperimentsForRatesAndKs(int x, int y, int z, bool useWave, const WaveParams& waveParams, 
	int minRate, int maxRate, int minK, int maxK, 
	int minN, int maxN, int minQ, int maxQ, int minS, int maxS,
	double minCompressionRatio, double maxCompressionRatio, double compressionRatioStep,
	int algorithm, std::ofstream& csv) {
	switch (algorithm) {
		case 0: // ZFP
			for (int rate = minRate; rate <= maxRate; ++rate) {
				int result = runExperiment(x, y, z, useWave, waveParams, rate, 0, 0, 
					0, 0, 0, algorithm, csv);
				if (result == 1) {
					return;
				}
				else if (result == 2)
				{
					// Throwing an exception
					throw std::runtime_error("Invalid algorithm used in experiment.");
				}
			}
			break;
		case 1: // SVD
			for (int k = minK; k <= maxK; ++k) {
				int result = runExperiment(x, y, z, useWave, waveParams, 0, k, 0, 
					0, 0, 0, algorithm, csv);
				if (result == 1) {
					return;
				}
				else if (result == 2)
				{
					// Throwing an exception
					throw std::runtime_error("Invalid algorithm used in experiment.");
				}
			}
			break;
		case 2: // Chebyshev
			for (int N = minN; N <= maxN; N = N + 2) {
				for (int Q = minQ; Q <= maxQ; Q = Q + 2) {
					for (int S = minS; S <= maxS; S = S + 2) {
						int result = runExperiment(x, y, z, useWave, waveParams, 0, 0, 0, N, Q, S, algorithm, csv);
						if (result == 1) {
							return;
						}
						else if (result == 2)
						{
							// Throwing an exception
							throw std::runtime_error("Invalid algorithm used in experiment.");
						}
					}
				}
			}
			break;
		case 3: // FFT
			for (double compressionRatio = minCompressionRatio; compressionRatio <= maxCompressionRatio; compressionRatio += compressionRatioStep) {
				int result = runExperiment(x, y, z, useWave, waveParams, 0, 0, compressionRatio, 
					0, 0, 0, algorithm, csv);
				if (result == 1) {
					return;
				}
				else if (result == 2)
				{
					// Throwing an exception
					throw std::runtime_error("Invalid algorithm used in experiment.");
				}
			}
			break;
		default:
			// Handle error or unknown algorithm
			break;
	}
}


void runTimingExperiment(int x, int y, int z) {
	std::ofstream csv("results_" + std::to_string(x) + ".csv");
	csv.precision(15);

	csv << "Algorithm Name" << ','
		<< "X" << ','
		<< "Y" << ','
		<< "Z" << ','
		<< "Amplitude" << ','
		<< "Phase" << ','
		<< "Fx" << ','
		<< "Fy" << ','
		<< "Fz" << ','
		<< "ZFP Rate" << ','
		<< "SVD K" << ','
		<< "Cheby N" << ','
		<< "Cheby Q" << ','
		<< "Cheby S" << ','
		<< "FFT Ratio" << ','
		<< "Mean Compress Time" << ','
		<< "Mean Decompress Time" << ','
		<< "Mean MSE" << ','
		<< "Original Size" << ','
		<< "Mean Compressed Size" << ','
		<< "Max Error" << ','
		<< "Max Error Original Value" << ','
		<< "Max Error Decompressed Value"
		<< std::endl;

	for (double amplitude = 0.1; amplitude <= 2.0; amplitude += 0.1) {
		std::cout << x << ": " << amplitude << std::endl;
		for (double phase = 0.0; phase <= M_PI; phase += M_PI / 2.0) {
			for (double fx = 0.1; fx <= 1.0; fx += 0.3) {
				for (double fy = 0.1; fy <= 1.0; fy += 0.3) {
					for (double fz = 0.1; fz <= 1.0; fz += 0.3) {
						WaveParams waveParams = {amplitude, phase, fx, fy, fz };

						for (int algorithm = 0; algorithm < 3; ++algorithm) {
							runExperimentsForRatesAndKs(x, y, z, true, waveParams,
								1, 64, //minRate, maxRate
								1, std:: min(x * y, z), //minK, maxK
								1, x, 1, y, 1, z,// minN, maxN, minQ, maxQ, minS, maxS
								0, 0, 0, //FFT 
								algorithm, csv);
						}
					}
				}
			}
		}
	}

	csv.close();
}

void runTimingExperimentFFTOnly() {
	std::ofstream csv("results_FFT.csv");
	csv.precision(15);
	csv << "Algorithm Name" << ','
		<< "X" << ','
		<< "Y" << ','
		<< "Z" << ','
		<< "Amplitude" << ','
		<< "Phase" << ','
		<< "Fx" << ','
		<< "Fy" << ','
		<< "Fz" << ','
		<< "FFT Rate" << ','
		<< "SVD K" << ','
		<< "Cheby N" << ','
		<< "Cheby Q" << ','
		<< "Cheby S" << ','
		<< "FFT Ratio" << ','
		<< "Mean Compress Time" << ','
		<< "Mean Decompress Time" << ','
		<< "Mean MSE" << ','
		<< "Original Size" << ','
		<< "Mean Compressed Size" << ','
		<< "Max Error" << ','
		<< "Max Error Original Value" << ','
		<< "Max Error Decompressed Value"
		<< std::endl;
	for (int x = 3; x <= 3; x++) {
		std::cout << "FFT: " << x << std::endl;
		for (double amplitude = 0.1; amplitude <= 2.0; amplitude += 0.1) {
			for (double phase = 0.0; phase <= M_PI; phase += M_PI / 2.0) {
				for (double fx = 0.1; fx <= 1.0; fx += 0.3) {
					for (double fy = 0.1; fy <= 1.0; fy += 0.3) {
						for (double fz = 0.1; fz <= 1.0; fz += 0.3) {
							WaveParams waveParams = {amplitude, phase, fx, fy, fz };

							for (int algorithm = 3; algorithm < 4; ++algorithm) { // Only the FFT algorithm
								runExperimentsForRatesAndKs(x, 7, 7, true, waveParams, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
									0.1, 1.0, 0.05, // minCompressionRatio, maxCompressionRatio, compressionRatioStep
									algorithm, csv);
							}
						}
					}
				}
			}
		}
	}

	csv.close();
}


void start() {

	std::vector<std::thread> threads;

	int y = 7;
	int z = 7;

	for (int x = 3; x <= 3; ++x) {
		threads.push_back(std::thread(runTimingExperiment, x, y, z));
		//threads.push_back(std::thread(std::bind(runTimingExperiment, x, y, z)));
	}
	threads.push_back(std::thread(runTimingExperimentFFTOnly));


	// Join the threads
	for (auto& th : threads) {
		th.join();
	}
	int i = 0;
}

