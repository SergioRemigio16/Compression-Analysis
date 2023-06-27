#ifndef ZFP_TESTING_H
#define ZFP_TESTING_H

#include <fstream>
#include <iostream>
#include <chrono>
#include <zfp.h>
#include <vector>
#include <numeric>
#include <cmath>
#include "Utilities.h"
#include "ZFPAlgorithms.h"



#define RUNS 1000
#define WARMUP_RUNS 100
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Run experiment for different compression rates
void runExperiment(int x, int y, int z, bool useWave, bool visualizeData, int rate);

// Run experiments for a range of rates
void runExperimentsForRates(int x, int y, int z, bool useWave, bool visualizeData, int minRate, int maxRate);

int runTimingExperiment();

#endif 
